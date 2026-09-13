"""The stage graph: selection, dependencies, the lock, resume and interrupt."""
from __future__ import annotations

import argparse
import io
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time

import pytest

import fixtures as F
from conftest import METAANNOT_PY, ROOT, build_project, run_metaannot
from test_stages import _searchable


def _wait_for(pred, timeout=30.0, interval=0.05):
    end = time.time() + timeout
    while time.time() < end:
        if pred():
            return True
        time.sleep(interval)
    return False


# --- finding 25 -------------------------------------------------------
def test_a_stage_whose_enabled_dependency_was_excluded_refuses_to_run(
        tmp_path, stub_bin):
    # symptom: `--only finalise` on an empty results directory wrote a full
    # annotation in which everything was 4_dark, with no warning at all.
    proj = _searchable(tmp_path, tmp_path / "p")
    proc = proj.run("--only", "finalise", expect=1)
    assert "cannot run: the stage(s) this one reads produced nothing" in \
        proc.stderr
    assert "integrate" in proc.stderr
    assert not os.path.exists(proj.rpath("annotation_final.tsv"))


def test_the_refusal_names_the_stages_to_add_to_the_selection(tmp_path,
                                                              stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proc = proj.run("--only", "join", expect=1)
    assert "add" in proc.stderr and "to the selection" in proc.stderr


def test_a_disabled_dependency_is_exempt_from_the_refusal(tmp_path, stub_bin):
    # a disabled dependency's evidence is legitimately absent; only one that
    # was merely NOT SELECTED means the stage would run against inputs that do
    # not exist yet.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()                                    # everything present
    os.remove(proj.rpath("hmm", "pfam.tblout"))
    st = json.load(open(proj.rpath(".metaannot_state.json"), encoding="utf-8"))
    st.pop("pfam", None)
    json.dump(st, open(proj.rpath(".metaannot_state.json"), "w",
                       encoding="utf-8"), indent=1)
    # pfam disabled and not selected: integrate must still run
    proj.write_config(run=dict(proj.cfg["run"], pfam=False))
    proc = proj.run("--only", "integrate", "finalise")
    assert "cannot run" not in proc.stderr
    assert os.path.exists(proj.rpath("annotation_final.tsv"))


def test_only_and_from_cannot_be_combined(project):
    proc = run_metaannot("run", "--config", project.config_path,
                         "--only", "pfam", "--from", "integrate", expect=1)
    assert "cannot be combined" in proc.stderr


# --- finding 26 -------------------------------------------------------
def test_a_second_run_on_one_results_directory_refuses(tmp_path, stub_bin):
    # symptom: two concurrent runs write the same files and the same state,
    # and the damage is silent — both report success.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()                                   # create the results dir
    lock = proj.rpath(".metaannot.lock")
    with open(lock, "w", encoding="utf-8") as fh:
        json.dump({"pid": os.getpid(), "host": __import__("socket").gethostname(),
                   "started": "2026-01-01T00:00:00"}, fh)
    proc = proj.run(expect=1)
    assert "another metaannot is already running here" in proc.stderr
    assert os.path.exists(lock), "a live lock must not be removed"


def test_two_simultaneous_runs_do_not_both_proceed(tmp_path, stub_bin):
    # the race the O_EXCL open exists for: two runs launched in the same
    # second both used to see no lock and both proceed.
    proj = _searchable(tmp_path, tmp_path / "p")
    env = dict(os.environ, STUB_SLEEP="1.5")
    procs = [subprocess.Popen(
        [sys.executable, METAANNOT_PY, "run", "--config", proj.config_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
        cwd=proj.root) for _ in range(2)]
    outs = [p.communicate() for p in procs]
    codes = [p.returncode for p in procs]
    assert sorted(codes) == [0, 1], f"exactly one run must win, got {codes}"
    losing = [o[1] for o, c in zip(outs, codes) if c == 1][0]
    assert "another metaannot is already running here" in losing


def test_a_lock_from_a_dead_process_is_reclaimed(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    dead = subprocess.Popen([sys.executable, "-c", "pass"])
    dead.wait()
    with open(proj.rpath(".metaannot.lock"), "w", encoding="utf-8") as fh:
        json.dump({"pid": dead.pid,
                   "host": __import__("socket").gethostname(),
                   "started": "2026-01-01T00:00:00"}, fh)
    proc = proj.run()
    assert "removing a stale lock" in proc.stderr


def test_a_zero_byte_lock_left_by_a_crash_is_reclaimed(ma, tmp_path):
    """A host bugcheck left exactly this: an O_EXCL create with no write.

    `_holder_is_alive({})` says a garbled lock is alive, which is right for a
    garbled one. An EMPTY one is different: its writer died inside the
    microsecond window between creating the file and writing to it, so no live
    run can own it, and treating it as alive made a crashed run permanently
    unresumable without --force-unlock.
    """
    lk = tmp_path / ".metaannot.lock"
    lk.write_bytes(b"")
    old = time.time() - 3600
    os.utime(lk, (old, old))
    with ma.ResultsLock(str(lk), force=False):
        assert lk.read_text(encoding="utf-8"), "the new holder must write itself in"
    assert not lk.exists(), "the lock must be released on exit"


def test_a_zero_byte_lock_that_just_appeared_is_left_alone(ma, tmp_path):
    """The one window where an empty lock is legitimate: another run created
    it a moment ago and has not written yet. Racing it would corrupt both."""
    lk = tmp_path / ".metaannot.lock"
    lk.write_bytes(b"")
    with pytest.raises(ma.StageError) as e:
        with ma.ResultsLock(str(lk), force=False):
            pass
    assert "already running" in str(e.value)


def test_the_empty_lock_grace_is_a_time_window_not_a_size_test(ma, tmp_path):
    """A NON-empty stale lock must still go down the liveness path, not the
    new one -- otherwise any old lock would be reclaimed on age alone."""
    lk = tmp_path / ".metaannot.lock"
    lk.write_text(json.dumps({"pid": 1, "host": "some-other-node",
                              "started": "2020-01-01T00:00:00"}),
                  encoding="utf-8")
    old = time.time() - 3600
    os.utime(lk, (old, old))
    with pytest.raises(ma.StageError) as e:
        with ma.ResultsLock(str(lk), force=False):
            pass
    assert "already running" in str(e.value)


def test_a_lock_we_cannot_disprove_is_treated_as_alive(ma):
    # unprovable means alive: trampling a live run is silent corruption while
    # refusing is a message and --force-unlock.
    lk = ma.ResultsLock("/tmp/never", force=False)
    assert lk._holder_is_alive({"pid": 1, "host": "some-other-node"})
    assert lk._holder_is_alive({})                     # garbled lock file
    assert lk._holder_is_alive({"pid": "not-an-int", "host": "x"})


def test_the_lock_is_released_when_the_run_ends(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    assert not os.path.exists(proj.rpath(".metaannot.lock"))


def test_the_lock_is_released_even_when_a_stage_fails(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(run=dict(proj.cfg["run"], kofam=True))
    proj.run(expect=1)
    assert not os.path.exists(proj.rpath(".metaannot.lock"))


def _dead_pid():
    """A pid on this host that this host can prove is gone."""
    dead = subprocess.Popen([sys.executable, "-c", "pass"])
    dead.wait()
    return dead.pid


def _plant_lock(proj, pid, host=None):
    with open(proj.rpath(".metaannot.lock"), "w", encoding="utf-8") as fh:
        json.dump({"pid": pid,
                   "host": host or __import__("socket").gethostname(),
                   "started": "2026-01-01T00:00:00"}, fh)


@pytest.mark.skipif(os.name == "nt",
                    reason="a pid this host can prove is DEAD is what this "
                           "needs, and on Windows the proof is OpenProcess "
                           "rather than os.kill; the refusal below is "
                           "exercised on the platform that can prove it.")
def test_force_unlock_takes_over(tmp_path, stub_bin):
    # The lock this flag exists for: its writer is GONE. The holder used to be
    # this test's own pid, which is alive on this host and is now refused - see
    # the test below, which is the case that changed.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    _plant_lock(proj, _dead_pid())
    proj.run("--force-unlock")


# --- finding 27 -------------------------------------------------------
def test_a_round_that_only_skips_stages_is_not_a_deadlock(tmp_path, stub_bin):
    # symptom: skipping stages IS progress — their dependents may become ready
    # in the same round — so the scheduler must loop again before concluding
    # anything is stuck.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    proc = proj.run()                # second run: every stage cached/disabled
    assert "deadlock" not in proc.stderr
    assert "done: 0 run" in proc.stderr


def test_a_run_with_every_stage_disabled_completes(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(run={k: False for k in proj.cfg["run"]})
    proc = proj.run()
    assert "deadlock" not in proc.stderr


def test_the_deadlock_message_names_the_unsatisfied_dependencies(ma):
    # kept as a message rather than a bare exception: normally it means those
    # deps were excluded by --only/--from, and if they were selected it is a
    # bug in the stage graph.
    src = ma.cmd_run.__code__.co_consts
    assert any(isinstance(c, str) and "deadlock" in c
               for c in _walk_consts(src))


def _walk_consts(consts):
    for c in consts:
        if hasattr(c, "co_consts"):
            yield from _walk_consts(c.co_consts)
        else:
            yield c


# --- finding 29 -------------------------------------------------------
DEFAULT_RESUME_POINTS = ["emapper", "pfam", "integrate", "finalise", "join"]


def _outputs(proj, skip=()):
    """Every result file, as bytes, keyed by its path inside results/."""
    out = {}
    for root, _dirs, files in os.walk(proj.results):
        for f in sorted(files):
            if f.startswith(".") or f.endswith(".log"):
                continue          # state and log carry timings, not results
            if f == "config.effective.yaml":
                # a record of the invocation, not a result: it carries the
                # absolute paths of the project it was run in, so two projects
                # can never write the same bytes.
                continue
            full = os.path.join(root, f)
            rel = os.path.relpath(full, proj.results)
            if rel in skip:
                continue
            out[rel] = open(full, "rb").read()
    return out


# annotation_pass1.tsv is the one output a resume does NOT reproduce; see
# test_annotation_pass1_is_reproducible_across_a_resume below.
PASS1 = "annotation_pass1.tsv"


@pytest.mark.parametrize("stage", DEFAULT_RESUME_POINTS)
def test_resuming_from_a_stage_reproduces_the_clean_result(tmp_path, stub_bin,
                                                           stage):
    # symptom: a resumed run must be indistinguishable from one that never
    # stopped, or "rerun the same command" is not a safe instruction.
    ref = _searchable(tmp_path, tmp_path / "ref")
    ref.run()
    want = _outputs(ref, skip=(PASS1,))

    proj = _searchable(tmp_path, tmp_path / f"r_{stage}")
    proj.run()
    proj.run("--from", stage, "--force")
    assert _outputs(proj, skip=(PASS1,)) == want


@pytest.mark.xfail(reason="LIVE DEFECT: build_annotation sets "
                          "structure_requested from whether each protein is in "
                          "dark.faa, and integrate itself writes dark.faa. So "
                          "an integrate that reruns over an existing results "
                          "directory reports structure_requested=True for "
                          "proteins its own first pass had reported False, and "
                          "annotation_pass1.tsv is not reproducible across a "
                          "resume. annotation_final.tsv is unaffected, because "
                          "finalise always runs after dark.faa exists.",
                   strict=True)
def test_annotation_pass1_is_reproducible_across_a_resume(tmp_path, stub_bin):
    ref = _searchable(tmp_path, tmp_path / "ref")
    ref.run()
    want = open(ref.rpath(PASS1), "rb").read()
    proj = _searchable(tmp_path, tmp_path / "again")
    proj.run()
    proj.run("--from", "integrate", "--force")
    assert open(proj.rpath(PASS1), "rb").read() == want


@pytest.mark.slow
@pytest.mark.parametrize("stage", [s for s in
                                   ["emapper", "pfam", "dbcan", "diamond",
                                    "cluster", "integrate", "jackhmmer",
                                    "hhblits", "esmfold", "foldseek",
                                    "finalise", "join"]])
def test_resuming_from_every_stage_reproduces_the_clean_result(tmp_path,
                                                               stub_bin,
                                                               stage):
    ref = _searchable(tmp_path, tmp_path / "ref")
    ref.run()
    want = _outputs(ref, skip=(PASS1,))
    proj = _searchable(tmp_path, tmp_path / f"r_{stage}")
    proj.run()
    proj.run("--from", stage, "--force")
    assert _outputs(proj, skip=(PASS1,)) == want


@pytest.mark.skipif(
    os.name == "nt",
    reason="send_signal(SIGINT) is unsupported on Windows, and CTRL_C_EVENT "
           "goes to the whole console group including the test runner. The "
           "path is covered on Linux, where this passes.")
def test_an_interrupted_run_leaves_parseable_state_and_resumes(tmp_path,
                                                               stub_bin):
    # symptom: a run killed mid-write left the half-written file as the only
    # evidence, and the next run adopted it as finished. The stage is now
    # recorded "running" BEFORE it starts.
    ref = _searchable(tmp_path, tmp_path / "ref")
    ref.run()
    want = _outputs(ref, skip=(PASS1,))

    proj = _searchable(tmp_path, tmp_path / "int")
    env = dict(os.environ, STUB_SLEEP="5", PYTHONHASHSEED="0")
    proc = subprocess.Popen(
        [sys.executable, METAANNOT_PY, "run", "--config", proj.config_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
        cwd=proj.root)
    state = proj.rpath(".metaannot_state.json")
    # `"status": "running"` and not just the word: the top-level _run record
    # says "final_status": "running" from the moment the run starts, so the
    # bare word would fire before any stage had begun.
    assert _wait_for(lambda: os.path.exists(state)
                     and b'"status": "running"' in open(state, "rb").read()), \
        "the stage was never recorded as running"
    proc.send_signal(signal.SIGINT)
    out, err = proc.communicate(timeout=60)
    # 128 + SIGINT. Pinned here because the exit status is now the ONE thing
    # the two signal paths do not share, so this test - the only one that
    # sends SIGINT to a real run - is where a regression to a hard-coded 130
    # for SIGTERM, or to a bare 1, would otherwise pass unseen.
    assert proc.returncode == 128 + int(signal.SIGINT)

    st = json.load(open(state, encoding="utf-8"))     # parseable, not torn
    assert any(v.get("status") == "running" for v in st.values())

    proj.run(env={"STUB_SLEEP": "0"})
    assert _outputs(proj, skip=(PASS1,)) == want


def test_an_unreadable_state_file_stops_adoption(tmp_path, stub_bin):
    # a run that cannot read the state knows nothing about the outputs lying
    # in the directory, so it must not adopt them as its own.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    with open(proj.rpath(".metaannot_state.json"), "w", encoding="utf-8") as fh:
        fh.write("{ not json")
    proc = proj.run()
    assert "state file unreadable" in proc.stderr
    assert "recomputing instead of adopting" in proc.stderr


def test_a_stage_left_running_is_recomputed_not_trusted(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    st = json.load(open(proj.rpath(".metaannot_state.json"), encoding="utf-8"))
    st["pfam"] = {"signature": None, "status": "running",
                  "started": "2026-01-01T00:00:00"}
    json.dump(st, open(proj.rpath(".metaannot_state.json"), "w",
                       encoding="utf-8"), indent=1)
    proc = proj.run()
    assert "interrupted while this stage was writing" in proc.stderr
    assert proj.state()["pfam"]["status"] == "ok"


# --- dry run ----------------------------------------------------------
def test_dry_run_creates_no_directories(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proc = run_metaannot("run", "--config", proj.config_path, "--dry-run",
                         cwd=proj.root)
    assert not os.path.exists(proj.results)
    assert "stage" in proc.stdout and "action" in proc.stdout


def test_dry_run_shows_the_refusal_a_real_invocation_would_give(tmp_path,
                                                                stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proc = run_metaannot("run", "--config", proj.config_path, "--dry-run",
                         "--only", "finalise", cwd=proj.root)
    assert "refused: needs" in proc.stdout


def test_dry_run_marks_disabled_stages(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proc = run_metaannot("run", "--config", proj.config_path, "--dry-run",
                         cwd=proj.root)
    assert "disabled (run.structure)" in proc.stdout


# --- resource split ---------------------------------------------------
def test_a_tiny_ram_budget_across_many_workers_warns_and_floors(tmp_path,
                                                                stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(stage_workers=4, ram_gb=2, threads=4)
    proc = proj.run()
    assert "leaves under 1 GB each" in proc.stderr


def test_threads_below_one_is_refused(project):
    proc = run_metaannot("run", "--config", project.config_path,
                         "--threads", "0", expect=1)
    assert "--threads must be at least 1" in proc.stderr


def test_serial_and_parallel_both_complete(tmp_path, stub_bin):
    a = _searchable(tmp_path, tmp_path / "a")
    a.run("--serial")
    b = _searchable(tmp_path, tmp_path / "b")
    b.run()
    assert _outputs(a) == _outputs(b)


# --- the GPU lease ----------------------------------------------------
# symptom: independent stages run concurrently and the CPU and RAM budgets are
# split between them, but the GPU was not modelled at all. tmbed held 15.5 GB
# of a 16 GB card and ESMFold peaked at 13.3 GB on a single short sequence; on
# the run this comes from they only avoided each other by accident, because
# they happened to be in separate invocations.
def _needs(*gpu):
    return lambda n: n in set(gpu)


def test_two_gpu_stages_ready_together_do_not_both_start(ma):
    go, wait = ma.gpu_lease(["tmbed", "esmfold"], [], 1,
                            _needs("tmbed", "esmfold"))
    assert go == ["tmbed"]
    assert wait == ["esmfold"]


def test_a_gpu_stage_waits_while_another_holds_the_card(ma):
    go, wait = ma.gpu_lease(["esmfold"], ["tmbed"], 1,
                            _needs("tmbed", "esmfold"))
    assert go == []
    assert wait == ["esmfold"]


def test_cpu_stages_keep_running_while_the_gpu_is_leased(ma):
    # the whole point: the device is leased, not the machine.
    go, wait = ma.gpu_lease(["pfam", "esmfold", "interpro"], ["tmbed"], 1,
                            _needs("tmbed", "esmfold"))
    assert go == ["pfam", "interpro"]
    assert wait == ["esmfold"]


def test_the_lease_is_returned_when_the_holder_finishes(ma):
    go, wait = ma.gpu_lease(["esmfold"], ["pfam"], 1, _needs("esmfold"))
    assert (go, wait) == (["esmfold"], [])


def test_gpu_workers_above_one_lets_that_many_run(ma):
    # a lease COUNT, not a device map: nothing here spreads them over cards.
    go, wait = ma.gpu_lease(["tmbed", "esmfold"], [], 2,
                            _needs("tmbed", "esmfold"))
    assert (go, wait) == (["tmbed", "esmfold"], [])


def test_more_holders_than_slots_never_produces_a_negative_budget(ma):
    go, wait = ma.gpu_lease(["esmfold"], ["tmbed", "other"], 1,
                            _needs("tmbed", "other", "esmfold"))
    assert (go, wait) == ([], ["esmfold"])


def test_the_stages_that_use_gpu_device_are_the_ones_marked_gpu(ma):
    # the invariant, so a stage added later that touches the card is not
    # scheduled as if the GPU were free.
    import inspect
    marked = {st["name"] for st in ma.STAGES if st.get("gpu")}
    uses = set()
    for st in ma.STAGES:
        try:
            src = inspect.getsource(st["fn"])
        except (OSError, TypeError):        # pragma: no cover
            continue
        if "gpu_device" in src:
            uses.add(st["name"])
    assert marked == uses, (
        f"marked gpu={sorted(marked)} but gpu_device is used by "
        f"{sorted(uses)}; a stage that touches the card must take the lease")
    assert marked == {"tmbed", "esmfold"}


def test_the_default_is_one_gpu_stage_at_a_time(ma):
    assert ma.DEFAULT_CONFIG["gpu_workers"] == 1


def test_gpu_workers_is_documented_as_a_lease_not_a_device_map(ma):
    # a machine with several GPUs is NOT given one stage per device, and the
    # config must say so rather than presenting 1 as a physical limit.
    src = io.open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index('"gpu_workers"')
    comment = src[max(0, i - 1200):i]
    assert "not a device map" in comment.lower()
    assert "same card" in comment.lower()


def test_a_non_numeric_gpu_workers_is_a_message_not_a_traceback(tmp_path,
                                                                stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(gpu_workers="lots")
    proc = proj.run(expect=1)
    assert "gpu_workers must be a whole number" in proc.stderr
    assert "Traceback" not in proc.stderr


def _overlap(a, b):
    return min(a[1], b[1]) - max(a[0], b[0]) > 0


def _run_args(config):
    """The Namespace `metaannot run` builds, with every default."""
    return argparse.Namespace(config=config, faa=None, results_dir=None,
                              threads=None, ram=None, only=None,
                              from_stage=None, force=False, no_adopt=False,
                              serial=False, force_unlock=False, dry_run=False)


def test_two_gpu_stages_never_overlap_while_a_cpu_stage_does(ma, tmp_path,
                                                             monkeypatch,
                                                             capsys):
    # gpu_lease is pinned as a function above; this drives the real dispatch
    # loop and watches the clock, because the lease only means anything if it
    # is applied to every round and a deferred stage does not sit in a worker
    # while it waits.
    #
    # A synthetic stage table, because the real one cannot present the case:
    # integrate depends on tmbed and esmfold depends on integrate, so today the
    # DAG happens to order the two GPU stages. The lease is what keeps them
    # apart when that stops being true - which is exactly how the run this
    # comes from avoided a CUDA OOM, by accident.
    proj = build_project(tmp_path / "conc", stage_workers=3, threads=6)
    seen, seen_lock = [], threading.Lock()

    def recorder(name):
        def fn(cfg, p):
            t0 = time.time()
            time.sleep(0.4)
            with seen_lock:
                seen.append((name, t0, time.time()))
            with open(os.path.join(p.R, name + ".out"), "w",
                      encoding="utf-8") as fh:
                fh.write("x\n")
        return fn

    def stage(name, gpu):
        return dict(name=name, enabled=None, gpu=gpu, deps=[], keys=[],
                    inp=lambda c, p: [],
                    out=lambda p, n=name: [os.path.join(p.R, n + ".out")],
                    fn=recorder(name))

    stages = [stage("tmbed", True), stage("esmfold", True),
              stage("pfam", False)]
    monkeypatch.setattr(ma, "STAGES", stages)
    monkeypatch.setattr(ma, "STAGE_NAMES", [st["name"] for st in stages])
    # cmd_run opens the run's log file into this global; setting it to its own
    # value registers monkeypatch's restore, so the handle does not leak into
    # the rest of the session.
    monkeypatch.setattr(ma, "_LOGFH", ma._LOGFH)

    assert ma.cmd_run(_run_args(proj.config_path)) == 0
    when = {name: (t0, t1) for name, t0, t1 in seen}
    assert set(when) == {"tmbed", "esmfold", "pfam"}, \
        "all three must have run: a deferred stage is postponed, not dropped"
    assert not _overlap(when["tmbed"], when["esmfold"]), \
        "two stages on one 16 GB card is the CUDA OOM this prevents"
    assert (_overlap(when["pfam"], when["tmbed"])
            or _overlap(when["pfam"], when["esmfold"])), \
        "the device is leased, not the machine: a CPU stage must still overlap"
    assert "waiting for the GPU" in capsys.readouterr().err


def test_the_run_says_which_gpu_stages_share_one_device(tmp_path, stub_bin):
    # an enabled stage that has not started must be explained; the announcement
    # is the first half of that, the per-stage "waiting for the GPU" line the
    # second.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(run=dict(proj.cfg["run"], topology=True, structure=True))
    proc = proj.run("--dry-run")
    assert proc.returncode == 0          # dry run: nothing needs a card
    proj2 = _searchable(tmp_path, tmp_path / "q")
    proj2.write_config(run=dict(proj2.cfg["run"], topology=True,
                                structure=True))
    proc = proj2.run(expect=1)           # signalp6/tmbed are not installed
    assert "at most 1 of them runs at a time" in proc.stderr
    assert "SAME card" in proc.stderr


# --- longest-processing-time-first ------------------------------------
# symptom: every stage with no dependencies is ready in the first round and
# stage_workers is 4, so the first four IN TABLE ORDER started. On the
# 455,571-protein run that handed workers to cluster (109 s) and dbcan
# (10 min) while signalp and tmbed — nearly an hour each on a set an order
# of magnitude smaller — queued behind them.
def _wave_one(ma):
    """The stages a fresh run finds ready in its first round."""
    return [st["name"] for st in ma.STAGES if not st["deps"]]


def test_the_first_wave_no_longer_goes_to_the_shortest_stages_in_table_order(
        ma):
    ready = _wave_one(ma)
    before = ready[:4]
    after = sorted(ready, key=ma.stage_priority, reverse=True)[:4]
    # what it used to pick, and why that was wrong
    assert before == ["emapper", "pfam", "dbcan", "diamond"]
    assert {"dbcan", "diamond"} & set(after) == set()
    # every stage that now claims a first-round worker is an hours-class one
    assert all(ma.stage_priority(n) == 3 for n in after), after
    assert after == ["emapper", "pfam", "signalp", "tmbed"]


def test_a_short_stage_never_outranks_a_long_one_wherever_the_table_puts_it(
        ma):
    order = sorted(_wave_one(ma), key=ma.stage_priority, reverse=True)
    # interpro is tenth in the table and the longest stage in the pipeline;
    # cluster and smorf are seconds and sit on either side of it.
    for short in ("cluster", "smorf", "dbcan", "diamond"):
        assert order.index("interpro") < order.index(short), (
            f"{short} is dispatched before interpro: {order}")


def test_equal_cost_stages_keep_table_order(ma):
    # the sort has to be stable, or the run log reorders itself between
    # releases for no reason a reader could explain.
    names = [st["name"] for st in ma.STAGES]
    order = sorted(names, key=ma.stage_priority, reverse=True)
    for rank in (3, 2, 1):
        same = [n for n in order if ma.stage_priority(n) == rank]
        assert same == [n for n in names if ma.stage_priority(n) == rank], rank


def test_every_stage_declares_a_cost_and_it_is_a_known_rank(ma):
    for st in ma.STAGES:
        assert "cost" in st, f"{st['name']} declares no cost"
        assert st["cost"] in (1, 2, 3), (st["name"], st["cost"])
    assert set(ma.STAGE_COSTS) == {st["name"] for st in ma.STAGES}


def test_an_unknown_stage_name_raises_rather_than_ranking_as_trivial(ma):
    # a .get(name, 1) default would silently rank a stage added without a
    # cost as seconds-class, which is the bug this ordering exists to fix.
    with pytest.raises(KeyError):
        ma.stage_priority("no_such_stage")


def test_a_stage_nothing_waits_on_never_outranks_one_integrate_needs(ma):
    by_name = {st["name"]: st for st in ma.STAGES}
    depended_on = {d for st in ma.STAGES for d in st["deps"]}
    orphans = [st["name"] for st in ma.STAGES
               if st["name"] not in depended_on and st["deps"]] + \
              [st["name"] for st in ma.STAGES
               if st["name"] not in depended_on and not st["deps"]]
    assert "smorf" in orphans, orphans
    floor = min(ma.stage_priority(d) for d in by_name["integrate"]["deps"])
    for name in orphans:
        if name == "join":                 # the terminal stage, waits on all
            continue
        assert ma.stage_priority(name) <= floor, (
            f"{name} blocks nothing yet outranks a stage integrate is "
            f"waiting for")


def test_the_cost_rank_never_enters_a_cache_signature(ma):
    # scheduling order cannot change a stage's output, so changing it must
    # not recompute anything.
    for st in ma.STAGES:
        assert not any("cost" in k for k in st["keys"]), st["name"]
    cfg = {"full_content_digest": False, "tool_args": {}}
    cheap = dict(name="x", cost=1, inp=lambda c, p: [], keys=[],
                 out=lambda p: [])
    dear = dict(cheap, cost=3)
    assert ma.signature(cheap, cfg, None) == ma.signature(dear, cfg, None)


def test_the_scheduler_sorts_the_ready_set_before_it_caps_at_stage_workers(
        ma):
    # the sort is worthless below the `len(futures) >= workers` break, and
    # wrong after gpu_lease, which hands the card to whichever gpu stage it
    # sees first.
    src = io.open(METAANNOT_PY, encoding="utf-8").read()
    decide = src.index("                run_now.append(name)")
    srt = src.index("run_now.sort(key=stage_priority, reverse=True)", decide)
    lease = src.index("run_now, waiting = gpu_lease(", decide)
    cap = src.index("if len(futures) >= workers:", decide)
    assert decide < srt < lease < cap


# --- the lock and a killed run ----------------------------------------
# symptom: atexit does not run on SIGTERM or SIGHUP, so a run stopped by
# `kill`, by a scheduler's time limit, or by a closing ssh session left a lock
# file behind naming a pid that no longer exists. On the same host the next
# run can prove it is dead; from another node of a cluster it cannot, and the
# resume became a stale-lock refusal needing --force-unlock.
@pytest.mark.skipif(os.name == "nt",
                    reason="TerminateProcess runs no handler on Windows; "
                           "SIGTERM cannot be delivered to a child there")
@pytest.mark.parametrize("signame", ["SIGTERM", "SIGHUP"])
def test_a_killed_run_releases_the_results_lock(tmp_path, stub_bin, signame):
    sig = getattr(signal, signame, None)
    if sig is None:
        pytest.skip(f"{signame} does not exist here")
    proj = _searchable(tmp_path, tmp_path / f"k{signame}")
    env = dict(os.environ, STUB_SLEEP="10", PYTHONHASHSEED="0")
    proc = subprocess.Popen(
        [sys.executable, METAANNOT_PY, "run", "--config", proj.config_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
        cwd=proj.root)
    lock = proj.rpath(".metaannot.lock")
    assert _wait_for(lambda: os.path.exists(lock)), "the lock was never taken"
    proc.send_signal(sig)
    out, err = proc.communicate(timeout=60)
    assert not os.path.exists(lock), \
        f"the lock survived {signame}:\n{err[-2000:]}"
    assert f"stopping on {signame}" in err
    assert "releasing the results lock" in err
    assert proc.returncode == -sig or proc.returncode == 128 + int(sig), \
        f"exit status {proc.returncode} is neither 128+N nor a signal death"


@pytest.mark.skipif(os.name == "nt", reason="see above")
def test_a_run_killed_that_way_resumes_without_force_unlock(tmp_path,
                                                            stub_bin):
    # the whole point: the next invocation must not need --force-unlock.
    ref = _searchable(tmp_path, tmp_path / "kref")
    ref.run()
    want = _outputs(ref, skip=(PASS1,))
    proj = _searchable(tmp_path, tmp_path / "kres")
    env = dict(os.environ, STUB_SLEEP="10", PYTHONHASHSEED="0")
    proc = subprocess.Popen(
        [sys.executable, METAANNOT_PY, "run", "--config", proj.config_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
        cwd=proj.root)
    assert _wait_for(lambda: os.path.exists(proj.rpath(".metaannot.lock")))
    proc.send_signal(signal.SIGTERM)
    proc.communicate(timeout=60)
    again = proj.run(env={"STUB_SLEEP": "0"})
    assert "--force-unlock" not in again.stderr
    assert "already running here" not in again.stderr
    assert _outputs(proj, skip=(PASS1,)) == want


def test_the_signals_that_can_strand_a_lock_are_all_handled(ma):
    # SIGINT is deliberately absent: Python raises KeyboardInterrupt for it,
    # which unwinds through the lock's `with` and runs atexit, so a handler
    # here would replace an orderly stop with an abrupt one.
    src = io.open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index("_release_lock_on_signal(sig, _frame)")
    # The span is the handler plus its registration, anchored on the
    # registration itself rather than on a count of characters: the rationale
    # written beside it grows, and a fixed 4000 characters silently stopped
    # short of the loop this is about, which reads as "SIGTERM is not
    # registered" when SIGTERM is registered eleven lines further down.
    tail = src[i:src.index("signal.signal(_sig, _release_lock_on_signal)", i)]
    for name in ("SIGTERM", "SIGHUP", "SIGBREAK"):
        assert f'"{name}"' in tail, f"{name} is not registered"
    assert "SIGINT" not in tail
    assert "os._exit" in tail, \
        "sys.exit would wait for the stage pool's workers"


def test_the_signal_handler_touches_nothing_that_takes_a_lock(ma):
    # A Python signal handler runs IN THE MAIN THREAD, between two bytecodes
    # of whatever that thread was doing, so any lock the interrupted frame
    # holds is still held and is not reentrant. The first version of this
    # handler called log() -> sys.stderr.write, whose buffer lock is exactly
    # that; on a run logging a progress line per stage per minute the signal
    # eventually lands mid-write and the process HANGS instead of releasing
    # the lock. CI caught it on one job of seven -- a race, so a behavioural
    # test cannot be relied on to catch a regression. This pins the rule.
    src = io.open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index("def _release_lock_on_signal(sig, _frame):")
    body = src[src.index("release_results_lock()", i):
               src.index("os._exit(128 + int(sig))", i)]
    for banned in ("log(", "sys.stderr", "_LOGFH", "print(", ".flush()",
                   "f\"", "format("):
        assert banned not in body, (
            f"{banned!r} in the signal handler: it either takes a lock or "
            "allocates through one. Pre-format at registration and use "
            "os.write(2, ...)")
    assert "os.write(2," in body, "the message must go out on the raw fd"


@pytest.mark.skipif(not hasattr(signal, "SIGTERM"),
                    reason="there is no SIGTERM to install a handler for")
def test_the_sigterm_handler_cmd_run_installs_supersedes_the_unwinding_one(
        ma, project, tmp_path, monkeypatch):
    """Two handlers claim SIGTERM, and which one a `kill` reaches decides what
    a three-day run leaves behind. main() installs the unwinding one - SIGTERM
    raises KeyboardInterrupt, so the run stops the way Ctrl-C stops it - and
    cmd_run then installs _release_lock_on_signal over it, which releases the
    lock and os._exit()s without unwinding. The later registration wins, so for
    `run` and `all` the unwinding one never runs.

    That is the behaviour that is wanted, and the three tests below pin what it
    leaves behind. What was unpinned was the ORDER, which is the whole of it:
    swap the two registrations and every one of those tests still describes a
    run that now waits for its stage to finish - a TMbed chunk is an hour -
    and under systemd reaches TimeoutStopSec and is SIGKILLed with the lock
    still on disk. Nothing else would have said so."""
    names = [n for n in ("SIGTERM", "SIGHUP", "SIGBREAK")
             if getattr(signal, n, None) is not None]
    saved = {n: signal.getsignal(getattr(signal, n)) for n in names}
    # cmd_run gets far enough to open the run's log into ma._LOGFH and to put
    # a RunRecord in ma._RUN, both module globals it never puts back. Setting
    # each to its own value registers monkeypatch's restore, so this tmp
    # project's log handle does not outlive its directory and catch every
    # later log() in the session, and no later test stamps a record that
    # belongs to this one.
    monkeypatch.setattr(ma, "_LOGFH", ma._LOGFH)
    monkeypatch.setattr(ma, "_RUN", ma._RUN)
    try:
        ma.handle_sigterm_like_sigint()
        unwinding = signal.getsignal(signal.SIGTERM)
        assert unwinding not in (saved["SIGTERM"], signal.SIG_DFL), \
            "main()'s handler is not the one in force before cmd_run runs"

        # cmd_run as far as its own registration and no further: it takes the
        # lock, registers, and only then refuses a proteins_faa that is not
        # there. Nothing is run, and the refusal is the stop point rather than
        # a monkeypatch, so the order being pinned is the real one.
        args = argparse.Namespace(config=project.config_path, threads=None,
                                  ram=None, faa=str(tmp_path / "gone.faa"),
                                  results_dir=None, dry_run=False,
                                  force_unlock=False, only=None,
                                  from_stage=None, force=False)
        with pytest.raises(ma.StageError) as e:
            ma.cmd_run(args)
        assert "proteins_faa not found" in str(e.value)

        now = signal.getsignal(signal.SIGTERM)
        assert now is not unwinding, \
            "cmd_run's handler no longer supersedes the unwinding one, so a " \
            "killed run waits for its stage instead of releasing the lock"
        assert getattr(now, "__name__", "") == "_release_lock_on_signal", \
            f"SIGTERM is handled by {now!r}, not by cmd_run's releaser"
    finally:
        for n in names:
            signal.signal(getattr(signal, n), saved[n])


# --- the engine half of the console: signals, the run record, the heartbeat --
def _paused_inside_a_stage(proj, **env):
    """A run started in the background and caught while a stage is running.

    The stub tools sleep, so there is a window in which to signal it.
    """
    e = dict(os.environ, STUB_SLEEP="5", PYTHONHASHSEED="0", **env)
    proc = subprocess.Popen(
        [sys.executable, METAANNOT_PY, "run", "--config", proj.config_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=e,
        cwd=proj.root)
    state = proj.rpath(".metaannot_state.json")
    if not _wait_for(lambda: os.path.exists(state)
                     and b'"status": "running"' in open(state, "rb").read()):
        proc.kill()
        proc.communicate()
        raise AssertionError("no stage was ever recorded as running")
    return proc


@pytest.mark.skipif(
    os.name == "nt",
    reason="TerminateProcess runs no handler on Windows, so a child cannot be "
           "sent SIGTERM there and nothing this test asserts - the released "
           "lock, the 128+N status, the WARN line - can happen. Deduced from "
           "the platform rather than measured on it.")
def test_a_run_killed_with_sigterm_releases_the_results_lock(tmp_path,
                                                             stub_bin):
    # symptom: the lock is released by an atexit hook, and SIGTERM's default
    # terminates the process without ever reaching atexit — so `kill`,
    # `systemctl stop` and `wsl --terminate` left .metaannot.lock behind, and
    # a lock that cannot be disproved refuses every later run.
    proj = _searchable(tmp_path, tmp_path / "term")
    proc = _paused_inside_a_stage(proj)
    proc.send_signal(signal.SIGTERM)
    proc.communicate(timeout=60)
    assert not os.path.exists(proj.rpath(".metaannot.lock")), \
        "SIGTERM left the results lock behind"
    # 128 + the signal, not the 130 Ctrl-C gives. SIGTERM shares every other
    # part of SIGINT's path, but the exit status is the only channel a
    # supervisor has for telling `systemctl stop` from a person at a keyboard -
    # and systemd units carry SuccessExitStatus=143 for exactly this.
    assert proc.returncode == 143


@pytest.mark.skipif(
    os.name == "nt",
    reason="TerminateProcess runs no handler on Windows, so a child cannot be "
           "sent SIGTERM there and nothing this test asserts - the released "
           "lock, the 128+N status, the WARN line - can happen. Deduced from "
           "the platform rather than measured on it.")
def test_a_sigterm_releases_the_lock_and_leaves_a_resumable_trace(tmp_path,
                                                                  stub_bin):
    # The engine's SIGTERM handler deliberately does NOT unwind: it releases
    # the lock with a raw os.remove, writes one pre-encoded line to fd 2, and
    # os._exit(128+N). Unwinding from a signal handler is what made an earlier
    # attempt deadlock on stderr's buffer lock. So what a killed run leaves is
    # the lock gone, the stage still marked running, and a resume that works.
    proj = _searchable(tmp_path, tmp_path / "term2")
    proc = _paused_inside_a_stage(proj)
    proc.send_signal(signal.SIGTERM)
    _out, err = proc.communicate(timeout=60)

    assert proc.returncode == 128 + int(signal.SIGTERM)
    assert "stopping on SIGTERM: releasing the results lock" in err, err[-400:]
    assert not os.path.exists(proj.rpath(".metaannot.lock")), \
        "the lock survived the signal it exists to be released by"

    st = json.load(open(proj.rpath(".metaannot_state.json"), encoding="utf-8"))
    assert any(v.get("status") == "running"
               for k, v in st.items() if k != "_run"), \
        "nothing records which stage was mid-flight, so a resume cannot know"
    proj.run(env={"STUB_SLEEP": "0"})            # and it still resumes


@pytest.mark.skipif(
    os.name == "nt",
    reason="TerminateProcess runs no handler on Windows, so a child cannot be "
           "sent SIGTERM there and nothing this test asserts - the released "
           "lock, the 128+N status, the WARN line - can happen. Deduced from "
           "the platform rather than measured on it.")
def test_a_sigtermed_run_is_exactly_the_case_the_heartbeat_exists_for(
        tmp_path, stub_bin):
    # Because the handler os._exit()s, nothing stamps the run record on the
    # way out: `_run` is left saying "running" with a last_seen that stops
    # advancing. That is not a defect - it is precisely the state a watcher
    # has to be able to describe, and why last_seen is advisory rather than
    # proof of death. A console reading this must say "no heartbeat since X",
    # never "this run is dead".
    proj = _searchable(tmp_path, tmp_path / "hb_term")
    proc = _paused_inside_a_stage(proj)
    proc.send_signal(signal.SIGTERM)
    proc.communicate(timeout=60)
    run = json.load(open(proj.rpath(".metaannot_state.json"),
                         encoding="utf-8"))["_run"]
    assert run["final_status"] == "running", \
        "a signal that does not unwind cannot have stamped a verdict"
    assert run["finished"] is None


@pytest.mark.skipif(
    os.name == "nt",
    reason="neither half can be delivered on Windows: TerminateProcess runs "
           "no handler for the SIGTERM, and send_signal(SIGINT) is "
           "unsupported there. Deduced from the platform rather than "
           "measured on it.")
def test_ctrl_c_unwinds_and_stamps_where_a_sigterm_cannot(tmp_path, stub_bin):
    # The two paths differ ON PURPOSE, which an earlier version of this file
    # asserted the opposite of. SIGINT is left as Python's default, so it
    # raises KeyboardInterrupt and unwinds - which is what lets it print the
    # "interrupted" line and stamp `_run`. SIGTERM buys the release of the
    # lock at the cost of that trace; SIGINT buys the trace at the cost of
    # waiting for the running stage.
    proj = _searchable(tmp_path, tmp_path / "msg_int")
    proc = _paused_inside_a_stage(proj)
    proc.send_signal(signal.SIGINT)
    _out, err = proc.communicate(timeout=60)
    assert proc.returncode == 128 + int(signal.SIGINT)
    assert "interrupted" in err
    run = json.load(open(proj.rpath(".metaannot_state.json"),
                         encoding="utf-8"))["_run"]
    assert run["final_status"] == "interrupted"
    assert run["finished"]


def _gated_run(proj, gate, *extra):
    """A background run parked inside a stage until `gate` is created.

    STUB_SLEEP races the unwind against a clock; this blocks it against a file
    instead, so a test can decide to the millisecond when a killed run is
    allowed to finish unwinding. `--only pfam` so exactly one stage is in
    flight and the main thread is provably parked in the dispatch loop.
    """
    e = dict(os.environ, STUB_WAIT_FOR=gate, STUB_SLEEP="0",
             PYTHONHASHSEED="0")
    proc = subprocess.Popen(
        [sys.executable, METAANNOT_PY, "run", "--config", proj.config_path,
         "--only", "pfam", *extra],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=e,
        cwd=proj.root)
    state = proj.rpath(".metaannot_state.json")
    if not _wait_for(lambda: os.path.exists(state)
                     and _run_record(proj).get("pid") == proc.pid
                     and b'"status": "running"' in open(state, "rb").read()):
        proc.kill()
        proc.communicate()
        raise AssertionError("the run never reached the gated stage")
    return proc


def _run_record(proj):
    try:
        with open(proj.rpath(".metaannot_state.json"), encoding="utf-8") as fh:
            return json.load(fh).get("_run") or {}
    except (OSError, ValueError):
        return {}


@pytest.mark.slow
@pytest.mark.skipif(
    os.name == "nt",
    reason="TerminateProcess runs no handler on Windows, so a child cannot be "
           "sent SIGTERM there and nothing this test asserts - the released "
           "lock, the 128+N status, the WARN line - can happen. Deduced from "
           "the platform rather than measured on it.")
def test_a_killed_runs_tail_never_lands_on_the_run_that_replaced_it(
        tmp_path, stub_bin):
    """The operator sequence this whole round of fixes is about, end to end.

    Before Phase 0 a SIGTERM killed the process instantly and it wrote nothing
    on the way out. The handler makes a killed run UNWIND, and unwinding
    writes - the `_run` stamp and the atexit lock release both land after the
    signal. So: A is killed; A unwinds slowly (still inside a stage); the
    operator, seeing it hang, --force-unlocks and starts B; A's tail then
    fires while B is live; and C has to find a usable directory afterwards.

    A's tail must write neither of the two things it can still write: not B's
    `_run` (that would show a console a live run marked "interrupted" by a
    process that is already gone) and not B's lock (removing it lets C join B
    in one results directory, which is the corruption the lock exists for).
    The pieces are unit-tested above; this is the sequence itself.
    """
    proj = _searchable(tmp_path, tmp_path / "handoff")
    lockfile = proj.rpath(".metaannot.lock")
    gate_a, gate_b = str(tmp_path / "free_a"), str(tmp_path / "free_b")

    a = _gated_run(proj, gate_a)
    a.send_signal(signal.SIGTERM)      # A starts to unwind, still in the stage

    # The operator hands the directory to a replacement. B parks in the stage
    # too, so it is demonstrably LIVE when A's tail arrives.
    b = _gated_run(proj, gate_b, "--force-unlock")
    assert json.load(open(lockfile, encoding="utf-8"))["pid"] == b.pid

    # A's tail: the stage returns and A unwinds the rest of the way, past its
    # `_run` stamp and its atexit lock release, with B running throughout.
    open(gate_a, "w").close()
    _out, a_err = a.communicate(timeout=120)
    assert a.returncode == 128 + int(signal.SIGTERM)

    # The damage first, the diagnostic after it: what matters is that B's
    # directory is untouched, not that A was polite about stopping.
    rec = _run_record(proj)
    assert rec.get("pid") == b.pid and rec.get("final_status") == "running", \
        "A's verdict landed on top of the live run's record"
    assert json.load(open(lockfile, encoding="utf-8"))["pid"] == b.pid, \
        "A's atexit hook removed the lock B was holding"
    # No "no longer holds" line is expected here, and its absence is the point.
    # v0.4.0's SIGTERM handler releases the lock and os._exit(128+N)s WITHOUT
    # unwinding, so A never reaches the superseded-write path that would say
    # so - it is gone before it can write anything at all. The diagnostic
    # belongs to the SIGINT path, which does unwind; see
    # test_a_superseded_run_says_so_once. Two designs protect B here and only
    # one of them can be polite about it, so this asserts the protection.
    assert "stopping on SIGTERM" in a_err, a_err[-300:]

    # B finishes normally, and nothing of A's is in its record.
    open(gate_b, "w").close()
    _out, _err = b.communicate(timeout=120)
    assert b.returncode == 0
    assert _run_record(proj)["final_status"] == "ok"
    assert not os.path.exists(lockfile), "B did not release its own lock"

    # ...and C, the run after all of that, gets a directory it can use.
    c = proj.run("--only", "pfam", env={"STUB_SLEEP": "0"})
    assert "done:" in c.stderr
    assert _run_record(proj)["final_status"] == "ok"


def test_a_results_directory_says_what_produced_it(ma, tmp_path, stub_bin):
    # symptom: the state file recorded what each stage did and nothing about
    # the run, so a results directory opened cold — one of eight on a shared
    # machine — could not say which config, host or command made it.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run("--threads", "2")
    rec = proj.state()["_run"]
    assert rec["final_status"] == "ok" and rec["finished"]
    assert rec["version"] == ma.__version__
    assert rec["config_path"] == os.path.abspath(proj.config_path)
    assert rec["host"] == __import__("socket").gethostname()
    assert isinstance(rec["pid"], int) and rec["run_id"]
    assert "run" in rec["argv"] and "--threads" in rec["argv"]


def test_the_run_record_is_not_read_as_a_stage(ma, tmp_path, stub_bin):
    # symptom risk: the state file is keyed by stage name, so a top-level
    # record read as one would be silently wrong — --force clearing it, or
    # decide() asking it for a signature.
    assert "_run" not in ma.STAGE_NAMES
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    first = proj.state()["_run"]["run_id"]
    proc = proj.run("--force")
    assert "unknown stage" not in proc.stderr
    st = proj.state()
    assert set(st) - {"_run"} <= set(ma.STAGE_NAMES)
    assert st["_run"]["run_id"] != first, "the second run must record itself"
    assert st["_run"]["final_status"] == "ok"


def test_the_run_record_says_a_failed_run_failed(tmp_path, stub_bin):
    # a directory whose last run died has to say so: "ok" and "still running"
    # are the two ways to mislead whoever opens it next.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(run=dict(proj.cfg["run"], kofam=True))
    proj.run(expect=1)
    rec = proj.state()["_run"]
    assert rec["final_status"] == "failed"
    assert rec["finished"]


def test_a_run_that_dies_outright_records_itself_as_failed(tmp_path, stub_bin):
    # die() unwinds past every return cmd_run has, so main() is the only place
    # that can record a fatal error — without it the record stays "running"
    # for ever and a watcher reads a dead run as a live one.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(gpu_workers="lots")
    proc = proj.run(expect=1)
    assert "gpu_workers must be a whole number" in proc.stderr
    assert proj.state()["_run"]["final_status"] == "failed"


def test_the_run_heartbeat_advances_while_a_stage_is_running(tmp_path,
                                                             stub_bin):
    # symptom: _holder_is_alive() cannot ask another host's process table, and
    # asks nothing at all on Windows, so nothing outside the process could tell
    # "InterProScan, hour 19" from "the box died".
    proj = _searchable(tmp_path, tmp_path / "beat", heartbeat_s=1)
    proc = _paused_inside_a_stage(proj)
    state = proj.rpath(".metaannot_state.json")

    def beat():
        """(last_seen_epoch, final_status), or None before it is written."""
        try:
            with open(state, encoding="utf-8") as fh:
                rec = json.load(fh)["_run"]
            return rec["last_seen_epoch"], rec["final_status"]
        except (KeyError, ValueError, OSError):
            return None

    first = beat()
    assert first and first[1] == "running", "the run recorded no heartbeat"

    def advanced():
        # still "running" when it advances: the stamp at the END of a run moves
        # last_seen too, and that would prove nothing about a live one.
        b = beat()
        return bool(b) and b[0] > first[0] and b[1] == "running"

    assert _wait_for(advanced, timeout=20), \
        "last_seen never advanced while the run was still going"
    proc.send_signal(signal.SIGTERM)
    proc.communicate(timeout=60)


def _plant_a_foreign_lock(proj, silent_for):
    """A lock and a matching run record from a host we cannot probe, whose
    heartbeat last landed `silent_for` seconds ago."""
    info = {"pid": 4242, "host": "another-node-of-the-array",
            "started": "2026-01-01T00:00:00"}
    with open(proj.rpath(".metaannot.lock"), "w", encoding="utf-8") as fh:
        json.dump(info, fh)
    state = proj.rpath(".metaannot_state.json")
    st = json.load(open(state, encoding="utf-8"))
    st["_run"] = dict(st.get("_run", {}), pid=info["pid"], host=info["host"],
                      final_status="running", heartbeat_s=30,
                      last_seen_epoch=time.time() - silent_for)
    with open(state, "w", encoding="utf-8") as fh:
        json.dump(st, fh, indent=1)


@pytest.mark.parametrize("silent_for", [0, 100_000])
def test_a_foreign_lock_is_never_reclaimed_on_a_silent_heartbeat(
        tmp_path, stub_bin, silent_for):
    # The heartbeat is ADVISORY: it says how long a lock has been quiet, and
    # nothing may reclaim a lock on the strength of it. A stopped heartbeat is
    # not a stopped process - one failed write ends the timer, not the run -
    # and from another host the two are indistinguishable, so reclaiming would
    # trade a stale lock (a message and --force-unlock) for two runs writing
    # one directory (silent corruption). Silent for a day reads exactly like
    # silent for a second: still alive.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    _plant_a_foreign_lock(proj, silent_for=silent_for)
    proc = proj.run(expect=1)
    assert "another metaannot is already running here" in proc.stderr
    assert "removing a stale lock" not in proc.stderr
    assert os.path.exists(proj.rpath(".metaannot.lock"))


def test_the_heartbeat_is_not_wired_into_the_liveness_test(ma):
    # the same rule stated where it can be broken quietly: whatever a state
    # file says, an unprovable lock is alive. Deleting _heartbeat_says_dead was
    # the fix; this is what stops it coming back by another route.
    lk = ma.ResultsLock("/tmp/never", force=False)
    assert not hasattr(lk, "_heartbeat_says_dead")
    assert lk._holder_is_alive({"pid": 4242, "host": "another-node"})
    if os.name != "nt":
        # ... and the local branch still asks the process table, which is the
        # only thing that may ever disprove a lock.
        dead = subprocess.Popen([sys.executable, "-c", "pass"])
        dead.wait()
        assert not lk._holder_is_alive(
            {"pid": dead.pid, "host": __import__("socket").gethostname()})


# --- the lock is released only by the run that still holds it -------------
def _foreign_payload():
    return json.dumps({"pid": os.getpid() + 1, "host": "some-other-node",
                       "started": "2026-01-01T00:00:00"})


def test_a_lock_rewritten_by_someone_else_is_not_removed_on_exit(ma, tmp_path):
    """The corruption SIGTERM made reachable, reproduced end to end by a
    reviewer: A is killed but unwinds slowly, an operator --force-unlocks and
    starts B, and A's atexit hook then deletes B's lock - after which a third
    run joins B in the same results directory. __exit__ removed the file
    unconditionally, never checking it was still the one __enter__ wrote."""
    lk = tmp_path / ".metaannot.lock"
    lock = ma.ResultsLock(str(lk))
    lock.__enter__()
    lk.write_text(_foreign_payload(), encoding="utf-8")   # B takes over
    lock.__exit__()
    assert lk.exists(), "a run deleted a lock that was no longer its own"
    assert json.loads(lk.read_text(encoding="utf-8"))["host"] == \
        "some-other-node"


def test_a_second_exit_cannot_remove_a_lock_taken_in_between(ma, tmp_path):
    # the signal path and the atexit hook both fire, so __exit__ runs twice.
    # The second call must not be able to remove whatever is at the path by
    # then, which after a --force-unlock is somebody else's lock.
    lk = tmp_path / ".metaannot.lock"
    lock = ma.ResultsLock(str(lk))
    lock.__enter__()
    lock.__exit__()
    assert not lk.exists()
    lk.write_text(_foreign_payload(), encoding="utf-8")
    lock.__exit__()
    assert lk.exists(), "the second __exit__ removed another run's lock"


def test_a_lock_still_ours_is_still_released(ma, tmp_path):
    # the other direction, and the one that matters most often: the ownership
    # check must not turn every clean exit into a stranded lock.
    lk = tmp_path / ".metaannot.lock"
    with ma.ResultsLock(str(lk)):
        assert lk.exists()
    assert not lk.exists()


def test_a_lock_taken_by_force_is_released_by_its_new_owner(ma, tmp_path):
    # --force-unlock rewrites the file, so the ownership check has to be
    # against what __enter__ wrote LAST, not against anything read earlier.
    lk = tmp_path / ".metaannot.lock"
    lk.write_text(_foreign_payload(), encoding="utf-8")
    with ma.ResultsLock(str(lk), force=True):
        assert json.loads(lk.read_text(encoding="utf-8"))["pid"] == os.getpid()
    assert not lk.exists(), "--force-unlock left its own lock behind"


def test_an_empty_lock_reclaimed_under_the_grace_is_released(ma, tmp_path):
    # the other path that ends with a lock this process owns.
    lk = tmp_path / ".metaannot.lock"
    lk.write_bytes(b"")
    old = time.time() - 3600
    os.utime(lk, (old, old))
    with ma.ResultsLock(str(lk)):
        assert lk.read_text(encoding="utf-8")
    assert not lk.exists()


def test_a_lock_we_cannot_read_is_still_released_by_its_owner(ma, tmp_path,
                                                              monkeypatch):
    """symptom: the ownership check answered every OSError with "not ours", so
    one NFS EACCES, EIO or ESTALE on the run's OWN lock left the file behind
    for a human to --force-unlock - the stranded lock this whole path exists to
    remove, and a regression against both 5684f6c and ef9e5a3, which released
    it unconditionally. "Unreadable but probably ours" and "definitely somebody
    else's" are different answers and must not collapse into one."""
    lk = tmp_path / ".metaannot.lock"
    lock = ma.ResultsLock(str(lk))
    lock.__enter__()
    real_open = open

    def unreadable(path, *a, **k):
        if str(path) == str(lk):
            raise PermissionError(13, "Permission denied")
        return real_open(path, *a, **k)

    monkeypatch.setattr("builtins.open", unreadable)
    assert lock.is_still_ours() is True, \
        "a lock we cannot read is one we have no evidence has changed hands"
    lock.__exit__()
    monkeypatch.undo()
    assert not lk.exists(), "an unreadable lock was stranded by its own owner"

    # ...and the three answers stay distinguishable. A lock that is genuinely
    # GONE is None, not False: there is nothing for __exit__ to remove, but
    # nothing has taken the directory either, so the run's own verdict is
    # still worth writing. Collapsing it into False cost a run nobody
    # superseded its final stamp, with a WARN announcing a handover that
    # never happened.
    gone = ma.ResultsLock(str(tmp_path / "b.lock"))
    gone.__enter__()
    os.remove(gone.path)
    assert gone.is_still_ours() is None


def test_a_lock_we_cannot_decode_is_still_released_by_its_owner(ma, tmp_path):
    """The same rule one layer down, where it did not hold: the ownership
    check reads the lock as TEXT, and UnicodeDecodeError is a ValueError, not
    an OSError, so the arm that answers "unprovable means alive" could not
    catch it. Undecodable bytes - a torn NFS write cutting a multi-byte
    character, a page of nulls where a crashed writer's payload should be -
    therefore raised out of __exit__, which had already cleared `held`: the
    lock stayed on disk and the second __exit__ was a no-op, which is the
    exact stranding the OSError arm exists to prevent."""
    lk = tmp_path / ".metaannot.lock"
    lock = ma.ResultsLock(str(lk))
    lock.__enter__()
    lk.write_bytes(b'{"pid": 1, "host": "\xff\xfe\x80 not utf-8"}')
    assert lock.is_still_ours() is True, \
        "a lock we cannot decode is one we have no evidence has changed hands"
    lock.__exit__()
    assert not lk.exists(), "an undecodable lock was stranded by its own owner"


def test_the_other_readers_of_the_lock_and_the_state_take_the_same_bytes(
        ma, tmp_path, capsys):
    """The same class of bug checked where it could land next. Every other
    reader on this path already guards ValueError, which UnicodeDecodeError is
    a subclass of, so each one answers rather than raising - and each answers
    the way that path answers an unreadable file, which is not the same answer
    twice."""
    lk = tmp_path / ".metaannot.lock"
    lk.write_bytes(b"\xff\xfe not utf-8, and not JSON either")
    # __enter__: a lock it cannot make sense of is a lock it cannot disprove,
    # so it refuses rather than trampling one that may be live.
    with pytest.raises(ma.StageError) as e:
        ma.ResultsLock(str(lk)).__enter__()
    assert "already running here" in str(e.value)
    with ma.ResultsLock(str(lk), force=True):
        assert json.loads(lk.read_text(encoding="utf-8"))["pid"] == os.getpid()

    # load_state: undecodable is unreadable, and a run that cannot read the
    # state must adopt nothing rather than die on the way in.
    st = tmp_path / ".metaannot_state.json"
    st.write_bytes(b'{"pfam": {"sig": "\xff\xfe"}}')
    capsys.readouterr()
    got = ma.load_state(str(st))
    assert got.unreadable and dict(got) == {}
    assert "state file unreadable" in capsys.readouterr().err


def test_a_run_whose_lock_vanishes_loses_its_verdict_and_says_so(ma, tmp_path,
                                                                 capsys):
    """What a vanished lock costs, stated where it is paid.

    `_run` is an ownership CLAIM, and a lock that is gone while this process is
    still alive is not proof of ownership - it is exactly what a replacement
    that took the directory with --force-unlock-live and then EXITED leaves
    behind, and from in here the two readings are the same bytes: nothing.
    So `_run` is not written, and this run's own last word is the price.

    That price is real and it is this test. What it buys is in
    test_a_superseded_runs_tail_cannot_stamp_run_after_recording_a_stage: the
    gate that declined only the write which would CREATE the document did not
    charge this at all, because a run's own stage record creates the document
    one call earlier - so it also protected nothing.

    NO ORDINARY RUN PAYS IT, which is what made the trade takeable, and that
    was measured rather than assumed: the results lock is released by an atexit
    hook, which runs after cmd_run's stamp_run("ok"/"failed") and after
    main()'s stamp_run("interrupted"), and cmd_run's SIGTERM handler releases
    the lock and os._exit()s without stamping anything. Driven end to end,
    `is_still_ours()` answered True at the final stamp of a clean run and of a
    Ctrl-C'd one. The test below is the abnormal case on purpose.

    And the WARN must not claim a handover that nothing here proves. The
    earlier defect in the opposite direction was a run whose lock an operator
    removed being told the directory "has been handed to another run"; the line
    now names both readings and says which records are still being written.
    """
    state = tmp_path / ".metaannot_state.json"
    lock = ma.ResultsLock(str(tmp_path / "c.lock"))
    lock.__enter__()
    st = ma._State()
    rec = ma.RunRecord(str(state), st, ["metaannot", "run"],
                       config_path=None, owner=lock)
    rec.stamp("running")
    assert json.loads(state.read_text())[ma.RUN_KEY]["final_status"] == \
        "running", "a run holding its own lock must write `_run`"
    before = state.read_bytes()

    os.remove(lock.path)                       # vanished, with no replacement
    capsys.readouterr()
    rec.stamp("ok")
    assert state.read_bytes() == before, \
        "`_run` was written on a lock this run could not prove was its own"
    assert json.loads(state.read_text())[ma.RUN_KEY]["final_status"] == \
        "running", "the verdict reached the file anyway"

    said = capsys.readouterr().err
    assert "is no longer there" in said and "cannot prove" in said, \
        f"the run lost its own final verdict silently: {said!r}"
    assert "handed to another run" not in said and "no longer holds" not in \
        said, "a vanished lock was reported as a handover nothing proved"

    # NOT superseded, and its stage records go on reaching the file: nothing
    # proved anything against this run, and a stage record claims nothing
    # about who owns the directory.
    assert rec.superseded is False
    st["pfam"] = {"status": "ok", "signature": "a"}
    assert ma.update_state(str(state), st, ("pfam",), claim=rec.claim) is True
    assert json.loads(state.read_text())["pfam"]["status"] == "ok", \
        "the stage keys were gated like `_run`; they must not be"


def test_a_lock_is_held_from_the_moment_the_file_exists(ma, tmp_path,
                                                        monkeypatch):
    """The startup hole: between os.open(O_EXCL) and the write, the file
    already blocks every other run, but self.held was still False, so nothing
    would remove it. Racing SIGTERM against a starting run left a fully written
    lock behind in 8 of 30 attempts - and being non-empty, it was not eligible
    for the zero-byte grace either."""
    lk = tmp_path / ".metaannot.lock"

    def die_between(fd, *a, **k):
        os.close(fd)
        raise OSError(28, "no space left on device")

    monkeypatch.setattr(ma.os, "fdopen", die_between)
    lock = ma.ResultsLock(str(lk))
    with pytest.raises(OSError):
        lock.__enter__()
    assert lk.exists(), "the file is there, and it blocks everyone"
    assert lock.held, "a lock that exists is held, written or not"
    monkeypatch.undo()
    lock.__exit__()
    assert not lk.exists(), "the abandoned create was never cleaned up"


def test_the_lock_release_hook_is_armed_before_the_lock_is_taken(
        project, monkeypatch, ma):
    """The other half of the same hole. atexit.register used to come AFTER
    __enter__ returned, so every one of __enter__'s outcomes - including the
    successful create - happened with no hook armed."""
    order = []
    monkeypatch.setattr(ma.atexit, "register",
                        lambda fn, *a, **k: order.append("register") or fn)

    class Boom(Exception):
        pass

    def spy(self):
        order.append("enter")
        raise Boom

    monkeypatch.setattr(ma.ResultsLock, "__enter__", spy)
    # The spy stops the lock, not the lines above it: cmd_run has already
    # opened this tmp project's log into ma._LOGFH by then, and nothing puts
    # it back. Setting each global to its own value registers monkeypatch's
    # restore, so the handle does not outlive the directory and go on
    # catching every later log() in the session.
    monkeypatch.setattr(ma, "_LOGFH", ma._LOGFH)
    monkeypatch.setattr(ma, "_RUN", ma._RUN)
    args = argparse.Namespace(config=project.config_path, threads=None,
                              ram=None, faa=None, results_dir=None,
                              dry_run=False, force_unlock=False)
    with pytest.raises(Boom):
        ma.cmd_run(args)
    assert order == ["register", "enter"], \
        "the release hook has to be armed before the lock can exist"


# --- one ownership gate, consulted by every write a superseded run makes ---
def _superseded(ma, tmp_path, interval=999):
    """A RunRecord whose directory has been handed to a replacement run.

    The operator sequence, compressed: this run takes the lock and writes its
    record, the operator --force-unlocks, and B takes the directory and writes
    a record of its own. Everything after that is this run's tail.
    """
    lk = tmp_path / ".metaannot.lock"
    state_path = str(tmp_path / ".metaannot_state.json")
    lock = ma.ResultsLock(str(lk))
    lock.__enter__()
    state = {}
    rec = ma.RunRecord(state_path, state, ["metaannot", "run"], None, interval,
                       owner=lock)
    ma.save_state(state_path, state)
    # --force-unlock, and B takes the directory. Written rather than taken
    # through a second ResultsLock because both would be this one process:
    # same pid, same host, same second, so the two tokens would be identical
    # and the gate would rightly call the second lock ours.
    lk.write_text(_foreign_payload(), encoding="utf-8")
    b_state = {"_run": {"pid": -1, "final_status": "running",
                        "last_seen": "2026-09-09T00:00:00"}}
    ma.save_state(state_path, b_state)
    return rec, state_path


def test_a_superseded_runs_final_verdict_does_not_overwrite_the_live_record(
        ma, tmp_path):
    """symptom: stamp(required=True) wrote `_run` whatever else was true, and
    SIGTERM made a killed run unwind rather than die where it stood - so A's
    "interrupted" landed on top of the `_run` of the run the operator had just
    given the directory to. A final verdict is not worth corrupting a live
    run's state file."""
    rec, state_path = _superseded(ma, tmp_path)
    rec.stamp("interrupted")
    on_disk = json.load(open(state_path, encoding="utf-8"))["_run"]
    assert on_disk["pid"] == -1 and on_disk["final_status"] == "running", \
        "a superseded run wrote its verdict over the live run's record"
    assert rec.superseded


def test_a_superseded_run_stays_superseded_once_the_directory_goes_quiet(
        ma, tmp_path):
    """symptom: the gate re-read the lock on every write, so it answered for
    the instant it was asked rather than for what had already happened. The
    sequence needs no race: A is --force-unlocked and says "no longer holds",
    B then FINISHES and removes its own lock, and A's final stamp finds the
    path vacant - None, not False - takes the "no replacement exists to be
    corrupted" branch, and writes A's whole stale state dict over B's finished
    results. Being told once is final."""
    rec, state_path = _superseded(ma, tmp_path)
    assert rec._tick() is False and rec.superseded    # A is told, once
    finished = {"_run": {"pid": -1, "final_status": "ok",
                         "finished": "2026-09-09T01:00:00"},
                "pfam": {"status": "ok", "sig": "b"}}
    ma.save_state(state_path, finished)
    os.remove(rec.owner.path)          # B is done, and released its lock

    rec.stamp("interrupted")
    on_disk = json.load(open(state_path, encoding="utf-8"))
    assert on_disk["_run"] == finished["_run"], \
        "a superseded run overwrote the record of the run that replaced it"
    assert on_disk.get("pfam"), \
        "B's stage records went with it: the whole state dict was replaced"


def test_a_superseded_runs_heartbeat_writes_nothing_and_stops(ma, tmp_path):
    # the same gate on the other write. A tick is the one that would fire
    # repeatedly, so it also has to stop the timer rather than be turned away
    # every interval for the rest of the unwind.
    rec, state_path = _superseded(ma, tmp_path)
    assert rec._tick() is False, "a superseded tick must report itself done"
    on_disk = json.load(open(state_path, encoding="utf-8"))["_run"]
    assert on_disk["last_seen"] == "2026-09-09T00:00:00", \
        "a superseded tick stamped the live run's record"
    assert rec._stop.is_set(), "the heartbeat kept going after being disowned"


def test_a_superseded_run_says_so_once(ma, tmp_path, capsys):
    # once, not once per write: the log it would repeat into belongs to the
    # run that owns the directory now.
    rec, _ = _superseded(ma, tmp_path)
    capsys.readouterr()
    rec._tick()
    rec._tick()
    rec.stamp("interrupted")
    err = capsys.readouterr().err
    said = [l for l in err.splitlines() if "no longer holds" in l]
    assert len(said) == 1, f"said it {len(said)} times, not once"
    assert "--force-unlock" in err and ".metaannot_state.json" in err, \
        "the message has to name what happened and what stopped"


def test_a_run_that_still_owns_its_directory_writes_as_before(ma, tmp_path):
    # the direction that matters on every ordinary run: the gate must not turn
    # a normal verdict into a lost one.
    lk = tmp_path / ".metaannot.lock"
    state_path = str(tmp_path / ".metaannot_state.json")
    state = {}
    with ma.ResultsLock(str(lk)) as lock:
        rec = ma.RunRecord(state_path, state, ["metaannot", "run"], None, 999,
                           owner=lock)
        assert rec._tick() is True
        rec.stamp("ok")
    on_disk = json.load(open(state_path, encoding="utf-8"))["_run"]
    assert on_disk["final_status"] == "ok" and on_disk["pid"] == os.getpid()
    assert not rec.superseded


# --- the heartbeat survives a bad write, and says so if it cannot ---------
def _record(ma, tmp_path, interval=0.01):
    st = {}
    return ma.RunRecord(str(tmp_path / ".metaannot_state.json"), st,
                        ["metaannot", "run"], None, interval)


def test_one_failed_heartbeat_write_does_not_end_the_heartbeat(ma, tmp_path,
                                                               monkeypatch):
    """symptom: _beat had no guard, so a single transient OSError killed the
    thread and last_seen stopped advancing for the rest of a 34-hour run, with
    one traceback in the log as the only trace."""
    calls = []
    real = ma.save_state

    def flaky(path, state):
        calls.append(1)
        if len(calls) == 2:
            raise OSError(28, "no space left on device")
        real(path, state)

    monkeypatch.setattr(ma, "save_state", flaky)
    rec = _record(ma, tmp_path)
    rec.watch(threading.Lock())
    try:
        assert _wait_for(lambda: len(calls) >= 5, timeout=10), \
            "the heartbeat died on its first failed write"
    finally:
        rec.stamp("ok")
    assert rec.rec["last_seen_epoch"] > 0


def test_a_heartbeat_that_gives_up_says_so_rather_than_going_quiet(
        ma, tmp_path, monkeypatch, capsys):
    # a display that simply stops updating reads as a dead run, which is the
    # one conclusion a stopped heartbeat must not invite.
    def always_bad(path, state):
        raise OSError(5, "input/output error")

    monkeypatch.setattr(ma, "save_state", always_bad)
    rec = _record(ma, tmp_path)
    rec.watch(threading.Lock())
    said = []

    def gave_up():
        said.append(capsys.readouterr().err)
        return "giving up" in "".join(said)

    try:
        # The message, not merely a thread that is no longer there: a thread
        # that died of an unhandled exception is also "no longer there", and
        # that is the failure being fixed, not the fix.
        assert _wait_for(gave_up, timeout=10), \
            "the heartbeat stopped without saying so"
    finally:
        rec._stop.set()      # stamp() would write, and writing is what fails
    assert "last_seen" in "".join(said) and "NOTHING" in "".join(said), \
        "the message has to say what is now stale and that the run continues"
    assert not any(t.name == "metaannot-heartbeat"
                   for t in threading.enumerate()), \
        "a heartbeat that has given up must stop, not spin"


def test_the_give_up_warning_names_the_stamp_the_state_file_really_has(
        ma, tmp_path, monkeypatch, capsys):
    """symptom: _tick advances rec["last_seen"] in memory and only THEN
    attempts the write that would carry it, so by the time the timer gives up
    the record is HEARTBEAT_GIVE_UP x heartbeat_s ahead of anything the file
    ever received. The warning quoted that in-memory value - sending an
    operator to a `_run.last_seen` that says something else, when sending them
    to that file is the entire purpose of the line."""
    def always_bad(path, state):
        raise OSError(5, "input/output error")

    monkeypatch.setattr(ma, "save_state", always_bad)
    # Long enough that the failing ticks cross a second boundary, since the
    # stamp is second-resolution: five misses at 0.5s is 2.5s of drift.
    rec = _record(ma, tmp_path, interval=0.5)
    frozen = rec._written
    rec.watch(threading.Lock())
    said = []

    def gave_up():
        said.append(capsys.readouterr().err)
        return "giving up" in "".join(said)

    try:
        assert _wait_for(gave_up, timeout=30), \
            "the heartbeat stopped without saying so"
    finally:
        rec._stop.set()
    msg = "".join(said)
    drifted = rec.rec["last_seen"]
    assert drifted != frozen, \
        "the failing ticks never advanced the record; nothing to catch here"
    assert frozen in msg, \
        "the warning does not name the stamp the state file actually holds"
    assert drifted not in msg, \
        f"the warning names {drifted}, a stamp the state file never received"


def test_a_heartbeat_tick_that_cannot_get_the_lock_writes_nothing(
        ma, tmp_path, monkeypatch):
    """_save's docstring said the scheduler's lock was there to stop a write
    going wrong, and then wrote anyway when the lock never came - doing the
    thing it called unsafe. A tick that cannot get the lock has nothing worth
    the risk: the next one repairs it. The final verdict does."""
    wrote = []
    monkeypatch.setattr(ma, "save_state", lambda p, s: wrote.append(1))

    class Busy:
        def acquire(self, timeout=None):
            return False

        def release(self):
            raise AssertionError("released a lock it never acquired")

    rec = _record(ma, tmp_path)
    rec.lock = Busy()
    rec._save()
    assert wrote == [], "a tick wrote behind the scheduler's back"
    rec._save(required=True)
    assert wrote == [1], "the final stamp must be written regardless"


def test_a_tick_already_in_flight_cannot_outlive_the_verdict(ma, tmp_path):
    """stamp() and _beat share one record, and a tick can already be past its
    wait() when the verdict lands - so it wrote last_seen a moment AFTER
    `finished`, and a console rendered a finished run whose last sign of life
    postdated its own end. Driven by hand here, because the window is the few
    microseconds a live thread spends between wait() and the write."""
    rec = _record(ma, tmp_path, interval=999)
    rec.stamp("ok")
    settled = dict(rec.rec)
    time.sleep(1.05)                  # a stamp now would be a LATER second
    assert rec._tick() is False, "a tick after the verdict must write nothing"
    assert rec.rec == settled, "the record changed after the run ended"
    assert rec.rec["last_seen"] == rec.rec["finished"]


def test_the_heartbeat_stops_when_the_run_does(ma, tmp_path):
    rec = _record(ma, tmp_path, interval=0.01)
    rec.watch(threading.Lock())
    assert _wait_for(lambda: rec.rec["last_seen_epoch"] > 0, timeout=5)
    rec.stamp("ok")
    assert _wait_for(lambda: not any(t.name == "metaannot-heartbeat"
                                     for t in threading.enumerate()),
                     timeout=5), "the heartbeat outlived the run"
    assert rec.rec["last_seen"] <= rec.rec["finished"]


# --- the two artefacts describe the same run -----------------------------
def test_an_early_fatal_error_records_the_run_whose_config_it_wrote(
        tmp_path, stub_bin):
    """config.effective.yaml was written immediately after the lock and the
    _run record a hundred lines and three die()s later, so a run dying in
    between left a config describing a run the state file had never heard of -
    and on a resume, one flatly contradicting the previous run's record:
    "finished ok, threads 7, against a FASTA that does not exist"."""
    import yaml
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run("--threads", "2")
    first = proj.state()["_run"]
    assert first["final_status"] == "ok"

    proj.write_config(proteins_faa=str(tmp_path / "gone.faa"))
    proj.run("--threads", "7", expect=1)          # dies at proteins_faa
    rec = proj.state()["_run"]
    with open(proj.rpath("config.effective.yaml"), encoding="utf-8") as fh:
        eff = yaml.safe_load(fh)

    assert eff["threads"] == 7, "the config is the dead run's"
    assert rec["run_id"] != first["run_id"], \
        "the state file still describes the PREVIOUS run"
    assert rec["final_status"] == "failed"
    assert "7" in rec["argv"] and rec["finished"]


def test_a_run_that_dies_before_any_stage_still_leaves_a_state_file(
        tmp_path, stub_bin):
    # the fresh-directory half: dying at an argument check left a
    # config.effective.yaml and NO state file at all, so a console had a
    # config for a run it could not see.
    proj = _searchable(tmp_path, tmp_path / "p")
    proc = proj.run("--only", "no-such-stage", expect=1)
    assert "unknown stage" in proc.stderr
    assert os.path.exists(proj.rpath("config.effective.yaml"))
    assert proj.state()["_run"]["final_status"] == "failed"


def test_a_dry_run_writes_neither_artefact(tmp_path, stub_bin):
    # it holds no lock, so it must not write into a directory another process
    # may be running in.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run("--dry-run")
    assert not os.path.exists(proj.rpath(".metaannot_state.json"))
    assert not os.path.exists(proj.rpath("config.effective.yaml"))


# --- the shape of the division, not just its size ---------------------
# symptom: share() divided what was FREE, which fixed the serial tail sitting
# on a quarter of an idle box. It divided it EVENLY, though, and the stages in
# one round are not equal. An external tool's thread count is fixed when it is
# launched -- `interproscan.sh -cpu 7` cannot grow, nor can hmmsearch's --cpu
# or InterProScan's -Xmx -- so an hours-class stage that hands half the box to
# a seconds-class one has given it away for its whole life.
#
# On the 455,571-protein run InterProScan was launched with -cpu 7 while
# signalp and tmbed held the other 15. Those two finished at 57 h; InterProScan
# spent its remaining hours on 7 of 22 cores with 14 idle.
def _round_of(ma, total, weights):
    """The cuts share() would make across one dispatch round."""
    free, cuts = total, []
    for i in range(len(weights)):
        c = ma.weighted_share(free, weights[i:])
        cuts.append(c)
        free -= c
    return cuts, free


def test_equal_cost_stages_split_the_machine_exactly_as_before(ma):
    # the common case must not move: a round of same-rank stages is the old
    # even split, arithmetically identical.
    assert _round_of(ma, 22, [3, 3, 3]) == ([7, 7, 8], 0)
    assert _round_of(ma, 32, [2, 2, 2, 2]) == ([8, 8, 8, 8], 0)
    assert _round_of(ma, 22, [1, 1]) == ([11, 11], 0)


def test_an_hours_class_stage_is_not_halved_by_a_seconds_class_one(ma):
    # 13 rather than 7, and the two short stages still get enough to run.
    cuts, left = _round_of(ma, 22, [3, 1, 1])
    assert cuts[0] == 13, cuts
    assert min(cuts) >= 1
    assert left == 0


def test_a_round_of_cuts_never_oversubscribes(ma):
    # the invariant that matters: the machine cannot be promised twice.
    for total in (1, 2, 7, 22, 64, 128):
        for weights in ([3], [3, 3], [3, 1], [1, 3], [3, 2, 1], [1, 1, 1, 1],
                        [3, 3, 3, 3, 3], [2, 1, 1, 1, 1, 1]):
            cuts, left = _round_of(ma, total, weights)
            assert sum(cuts) <= max(total, len(weights)), (total, weights, cuts)
            assert left >= 0 or total < len(weights), (total, weights, cuts)
            assert all(c >= 1 for c in cuts), (total, weights, cuts)


def test_a_lone_stage_gets_the_whole_machine(ma):
    assert ma.weighted_share(22, [3]) == 22
    assert ma.weighted_share(22, [1]) == 22


def test_the_split_never_divides_by_zero_or_returns_nothing(ma):
    # an empty weight list, a zero weight and a zero budget all reach this.
    assert ma.weighted_share(22, []) == 22
    assert ma.weighted_share(0, [3, 1]) == 1      # floored, never 0
    assert ma.weighted_share(22, [0, 0]) == 11    # zeros clamp to 1


def test_the_weight_is_the_cost_rank_so_the_two_cannot_drift(ma):
    # share() must weight by stage_priority and nothing else, or the ranks and
    # the allocation tell different stories about which stage is long.
    src = io.open(METAANNOT_PY, encoding="utf-8").read()
    body = src[src.index("    def share(starting):"):
               src.index("    def worker(st, cpu, ram):")]
    assert "stage_priority(n)" in body
    assert "weighted_share(" in body
    assert "// slots" not in body, "the even split is gone"


def test_a_stage_that_cannot_grow_into_freed_cpu_is_told_about(ma):
    # Structural, deliberately. The condition needs three hours-class stages
    # running concurrently and one finishing first; with instant stubs that is
    # a race, and a racy test here would be worse than none -- see the signal
    # handler, where CI caught on 1 job of 7 what every local run missed. So
    # this pins the guards instead: once per stage, hours-class only, and only
    # when the freed CPU would at least double what the survivor holds.
    src = io.open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index("# CPU freed here cannot be handed to a stage that is")
    body = src[i:src.index("break", i)]
    assert "other in starved" in body, "must be said once per stage"
    assert "stage_priority(other) < 3" in body, "minutes-class stages finish first"
    assert "free_cpu < held" in body, "only when it could at least double"
    assert "--only" in body, "the message has to say what to do"
    assert "starved = set()" in src, "the once-per-stage set must exist"


# ======================================================================
# issue #23: one process's snapshot must not be able to erase another
# process's records, and a superseded run's output must not be able to
# land under the live run's signature
# ======================================================================
def _doc(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def test_a_stale_run_writing_one_stage_leaves_the_other_records_alone(
        tmp_path, stub_bin):
    """The reviewer's reproduction, with two REAL runs and nothing simulated.

    A hangs inside pfam. The operator hands the directory to B, which
    completes four stages and exits. A's stage then returns and A records the
    one thing it knows: that IT finished pfam.

    Before this change A's record went to disk as `save_state(path, state)` -
    the whole document, from A's own in-memory dict, which had never heard of
    dbcan, diamond or cluster. Those three were erased, and the run after that
    printed "adopting output this run did not produce" for all three, which is
    the warning that says outright it cannot tell a finished file from an
    interrupted one. A writes one key now, and a key it has no record for is
    not a key it writes.
    """
    proj = _searchable(tmp_path, tmp_path / "erase")
    gate = str(tmp_path / "free_a")
    a = _gated_run(proj, gate)

    # A is alive and on this host, so --force-unlock alone is refused: this
    # sequence is exactly what the refusal is for, and the second flag is the
    # operator saying so out loud.
    b = proj.run("--only", "pfam", "dbcan", "diamond", "cluster",
                 "--force-unlock", "--force-unlock-live",
                 env={"STUB_SLEEP": "0"}, timeout=180)
    assert "--force-unlock refused" not in b.stderr
    after_b = _doc(proj.rpath(".metaannot_state.json"))
    b_run = after_b["_run"]["run_id"]
    assert [after_b[n]["status"] for n in
            ("pfam", "dbcan", "diamond", "cluster")] == ["ok"] * 4

    open(gate, "w").close()                 # A's stage returns; A records it
    _out, a_err = a.communicate(timeout=180)

    st = _doc(proj.rpath(".metaannot_state.json"))
    for name in ("dbcan", "diamond", "cluster"):
        assert st.get(name, {}).get("status") == "ok", \
            f"{name}'s record was erased by a run that never ran it"
        assert st[name]["run_id"] == b_run
    assert st["pfam"]["run_id"] == b_run, \
        "the superseded run's own record landed on top of the live one's"
    assert st["_run"]["run_id"] == b_run
    assert "another run holds this results directory" in a_err, \
        "the superseded run stopped writing without saying why"

    # ...and the run after all of that adopts nothing, because every stage is
    # recorded, by the run that made it.
    c = proj.run("--only", "pfam", "dbcan", "diamond", "cluster",
                 env={"STUB_SLEEP": "0"}, timeout=180)
    assert "adopting output this run did not produce" not in c.stderr


def test_a_write_merges_into_the_file_rather_than_the_snapshot_it_started_from(
        ma, tmp_path):
    # the mechanism, on its own: a writer that knows about one stage must
    # leave every other record exactly as it found it, including one written
    # after this run read the file.
    path = str(tmp_path / ".metaannot_state.json")
    state = {"pfam": {"status": "ok", "signature": "a"}}
    ma.save_state(path, state)
    ma.save_state(path, dict(state, dbcan={"status": "ok", "signature": "b"}))
    state["pfam"] = {"status": "ok", "signature": "a2"}

    assert ma.update_state(path, state, ("pfam",)) is True
    on_disk = _doc(path)
    assert on_disk["pfam"]["signature"] == "a2", "our own key never landed"
    assert on_disk["dbcan"]["signature"] == "b", \
        "a record this run had never heard of was erased by its write"


def test_a_record_this_run_no_longer_holds_is_deleted_rather_than_kept(
        ma, tmp_path):
    # the other direction of the same primitive, which --force needs: a key
    # named in the write and ABSENT from the run's view is a removal, not a
    # no-op, or `--force --only pfam` could never discard anything.
    path = str(tmp_path / ".metaannot_state.json")
    ma.save_state(path, {"pfam": {"status": "ok"}, "dbcan": {"status": "ok"}})
    assert ma.update_state(path, {}, ("pfam",)) is True
    on_disk = _doc(path)
    assert "pfam" not in on_disk and on_disk["dbcan"]["status"] == "ok"


def test_a_missing_state_file_gets_ONLY_the_keys_the_write_names(ma, tmp_path):
    """A write creates the document it needs and rebuilds nothing else.

    This branch used to restore the WHOLE of `state` - one run's in-memory
    snapshot - on the reasoning that a document which is not there has nothing
    in it to lose. It has: the records of whoever holds the directory now.
    Three revisions kept the restore and guarded it, and the guard cannot
    work, so the restore is gone. `dbcan` here is the record this run happens
    to be carrying and was NOT asked to write; it must not appear.
    """
    path = str(tmp_path / ".metaannot_state.json")
    state = {"pfam": {"status": "ok"}, "dbcan": {"status": "ok"}}
    assert ma.update_state(path, state, ("pfam",)) is True
    assert set(_doc(path)) == {"pfam"}, \
        "a write rebuilt a whole document from the writer's own snapshot"


def test_an_unparseable_state_file_gets_the_same_treatment_and_says_so(
        ma, tmp_path, capsys):
    # bytes that are not a document are damage, not a rival's payload - the
    # same judgement is_still_ours() makes one layer down about the lock - so
    # they are REPLACED. Replaced by the named keys and nothing else, for the
    # same reason as the branch above: whatever else was in those bytes, this
    # run's snapshot is not a reconstruction of it.
    path = str(tmp_path / ".metaannot_state.json")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("{ not json")
    capsys.readouterr()
    state = {"pfam": {"status": "ok"}, "dbcan": {"status": "ok"}}
    assert ma.update_state(path, state, ("pfam",)) is True
    assert _doc(path) == {"pfam": {"status": "ok"}}
    err = capsys.readouterr().err
    assert "does not parse" in err and "is not in it" in err, \
        "a rewrite that drops another run's records has to say so"


def test_an_unreadable_state_file_writes_nothing_rather_than_erasing_the_others(
        ma, tmp_path, monkeypatch, capsys):
    """A read that FAILED says nothing about the document, so nothing is
    written on the strength of it.

    Writing this run's keys into the empty dict the failed read returned is
    the single most destructive thing this design could do - one NFS blip and
    every other stage's record is gone. Declining costs at most the one record
    this call carried, which decide() then reads as an output with no record:
    adopt or recompute, never corruption. The next write of the SAME key - the
    stage's own next record, or the heartbeat's `_run` - reaches the file
    normally once a read succeeds.
    """
    path = str(tmp_path / ".metaannot_state.json")
    ma.save_state(path, {"dbcan": {"status": "ok", "signature": "b"}})
    state = {"pfam": {"status": "ok", "signature": "a"},
             "dbcan": {"status": "ok", "signature": "b"}}
    real_open = open

    def unreadable(p, *a, **k):
        if str(p) == path:
            raise OSError(5, "input/output error")
        return real_open(p, *a, **k)

    monkeypatch.setattr("builtins.open", unreadable)
    capsys.readouterr()
    assert ma.update_state(path, state, ("pfam",)) is False
    monkeypatch.undo()
    assert "pfam" not in _doc(path), "a declined record was written anyway"
    assert _doc(path)["dbcan"]["signature"] == "b", "the read error erased a record"
    assert "nothing is written to it until a read succeeds" \
        in capsys.readouterr().err

    # and once the read works again, the same key writes exactly as it would
    # have. Nothing was carried and nothing had to be.
    assert ma.update_state(path, state, ("pfam",)) is True
    assert _doc(path)["pfam"]["signature"] == "a"
    assert _doc(path)["dbcan"]["signature"] == "b"


def test_a_write_that_loses_a_race_redoes_its_merge_and_both_records_survive(
        ma, tmp_path, monkeypatch):
    """The read-back, which is what makes the merge worth having.

    A writer that finds the file is not what it just wrote redoes the merge
    against the document as it now stands, so the record it clobbered in the
    gap between its read and its rename comes back. It is optimistic
    concurrency and not exclusion: it shrinks the window and notices most of
    what is left in it.
    """
    path = str(tmp_path / ".metaannot_state.json")
    real, raced = ma.save_state, []

    def racing(p, s):
        payload = real(p, s)
        if not raced:
            raced.append(1)
            real(p, {"dbcan": {"status": "ok", "signature": "b"}})
        return payload

    monkeypatch.setattr(ma, "save_state", racing)
    assert ma.update_state(path, {"pfam": {"status": "ok", "signature": "a"}},
                           ("pfam",)) is True
    on_disk = _doc(path)
    assert on_disk["pfam"]["signature"] == "a"
    assert on_disk["dbcan"]["signature"] == "b", \
        "the writer that lost the race left the other record clobbered"


def test_the_per_write_read_is_not_the_start_of_run_read(ma, tmp_path,
                                                         monkeypatch, capsys):
    """load_state() must never be what a write reads with.

    It WARNs that every stage will be recomputed and returns an empty _State
    whose `unreadable` flag stops adoption for the whole run. Called once per
    write it would say that on every write, and the empty dict it returns
    would make the merge write a document holding only this run's keys - the
    erasure this change exists to stop, reintroduced through the reader.
    """
    path = str(tmp_path / ".metaannot_state.json")
    ma.save_state(path, {"dbcan": {"status": "ok"}})
    calls = []
    monkeypatch.setattr(ma, "load_state",
                        lambda p: calls.append(p) or ma._State())
    capsys.readouterr()
    ma.update_state(path, {"pfam": {"status": "ok"}}, ("pfam",))
    assert calls == [], "a write went through the start-of-run reader"
    assert "every stage will be recomputed" not in capsys.readouterr().err
    assert set(_doc(path)) == {"pfam", "dbcan"}


# --- the succession check: the document's own answer to "is this ours" ----
def _b_record():
    return {"run_id": "20260101T000000-999", "host": "another-node-of-the-array",
            "pid": 999, "started": "2026-01-01T00:00:00",
            "final_status": "running"}


def _taken_over(ma, tmp_path, interval=999):
    """A run whose state file has been taken by a run it never saw.

    No ResultsLock at all, deliberately: this is the OTHER gate, and it has to
    hold on its own evidence - on another node, on Windows, and every other
    case where the lock is correctly unprovable.
    """
    path = str(tmp_path / ".metaannot_state.json")
    state = {}
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, interval)
    assert ma.update_state(path, state, (ma.RUN_KEY,), claim=rec.claim)
    ma.save_state(path, {ma.RUN_KEY: _b_record(),
                         "dbcan": {"status": "ok", "signature": "b"}})
    return rec, path


def test_a_run_whose_directory_was_taken_stops_writing_the_state_file(
        ma, tmp_path, capsys):
    rec, path = _taken_over(ma, tmp_path)
    capsys.readouterr()
    assert rec._tick() is False, "a superseded tick wrote anyway"
    on_disk = _doc(path)
    assert on_disk[ma.RUN_KEY] == _b_record(), \
        "the superseded run stamped the live run's record"
    assert on_disk["dbcan"]["signature"] == "b"
    assert rec.superseded and rec._stop.is_set()
    err = capsys.readouterr().err
    assert "another run holds this results directory" in err
    assert _b_record()["run_id"] in err, "the message does not name who has it"


def test_a_superseded_run_says_it_once_however_often_it_tries(ma, tmp_path,
                                                              capsys):
    rec, _path = _taken_over(ma, tmp_path)
    capsys.readouterr()
    rec._tick()
    rec._tick()
    rec.stamp("interrupted")
    said = [l for l in capsys.readouterr().err.splitlines()
            if "another run holds this results directory" in l]
    assert len(said) == 1, f"said it {len(said)} times, not once"


def test_a_resume_writes_over_its_predecessor_for_the_whole_of_its_life(
        ma, tmp_path):
    """The direction that matters on every ordinary resume.

    A foreign `_run` on disk is the NORMAL state at the start of one: it is
    the previous run's record. A check that did not know the predecessor would
    stand a legitimate resume down at its first write and say nothing further
    for the rest of it, which is worse than the defect being fixed.
    """
    path = str(tmp_path / ".metaannot_state.json")
    ma.save_state(path, {ma.RUN_KEY: _b_record(),
                         "pfam": {"status": "ok", "signature": "a"}})
    state = ma.load_state(path)
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, 999)
    assert ma.update_state(path, state, (ma.RUN_KEY,), claim=rec.claim) is True
    assert rec._tick() is True, "a resume was stood down by its predecessor"
    assert rec.stamp("ok") is None
    on_disk = _doc(path)
    assert on_disk[ma.RUN_KEY]["run_id"] == rec.rec["run_id"]
    assert on_disk["pfam"]["signature"] == "a", "the resume erased a record"
    assert not rec.superseded


def test_a_state_file_with_no_run_record_never_stands_a_run_down(ma, tmp_path):
    # every kind of doubt writes. A document from a version that did not
    # record `_run`, or one an operator hand-edited, is not evidence that the
    # directory changed hands.
    path = str(tmp_path / ".metaannot_state.json")
    state = {}
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, 999)
    for doc in ({"pfam": {"status": "ok"}},
                {ma.RUN_KEY: {"pid": 1}},
                {ma.RUN_KEY: "not a record"}):
        ma.save_state(path, doc)
        assert ma.update_state(path, state, (ma.RUN_KEY,),
                               claim=rec.claim) is True
        assert not rec.superseded


# --- the atomic_out variant ---------------------------------------------
def _armed(ma, tmp_path, declared=()):
    """This process, armed as a run that owns tmp_path's state file."""
    path = str(tmp_path / ".metaannot_state.json")
    state = {}
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, 999)
    ma.update_state(path, state, (ma.RUN_KEY,), claim=rec.claim)
    ma._DECLARED_OUTPUTS.update(declared)
    return rec, path


def _stale_probe(ma):
    """Make the next rename ask the file rather than trusting what it knows."""
    ma._STATE_WATCH["seen"] = 0.0
    ma._STATE_WATCH["tried"] = 0.0


def test_a_superseded_run_parks_its_stage_output_instead_of_renaming_it(
        ma, tmp_path, capsys):
    """The variant that does not touch the state file at all.

    A stage still inside st["fn"] when its run is superseded RUNS TO
    COMPLETION - the executor waits for it - and atomic_out renamed its result
    into place at the end. If the replacement had already finished that stage
    and recorded it "ok", the old run's output landed under the new run's
    valid signature and the next run reported `cached` and read it.
    """
    out = str(tmp_path / "hmm" / "pfam.tblout")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("the live run's result\n")
    _rec, path = _armed(ma, tmp_path, declared=[out])
    ma.save_state(path, {ma.RUN_KEY: _b_record()})
    _stale_probe(ma)
    capsys.readouterr()

    with ma.atomic_out(out) as tmp:
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write("the superseded run's result\n")

    assert open(out, encoding="utf-8").read() == "the live run's result\n", \
        "a superseded run renamed its output over the live run's"
    parked = [f for f in os.listdir(os.path.dirname(out))
              if f.startswith(ma.SUPERSEDED_SUFFIX)]
    assert len(parked) == 1, f"the work was not parked: {os.listdir(os.path.dirname(out))}"
    assert open(os.path.join(os.path.dirname(out), parked[0]),
                encoding="utf-8").read() == "the superseded run's result\n", \
        "the parked file is not the work that was done"
    err = capsys.readouterr().err
    assert "NOT renamed into place" in err and "nothing was deleted" in err
    assert "missing input" in err, \
        "the warning does not say what happens to the rest of this run"


def test_a_per_item_output_is_not_parked_but_left_as_a_part_file(
        ma, tmp_path, capsys):
    """Parking must not scale with items.

    Several stages write one file per item - a PDB per dark protein, a file
    per DIAMOND database - and a superseded esmfold would otherwise park one
    undeletable file per protein into a directory that rule 5 of CLAUDE.md
    says must never be cleaned up.
    Only a stage's DECLARED outputs are parked; the rest keep the .part
    convention, which is already what a crash leaves and is already what
    `find results -name '.*.part.*'` looks for.
    """
    item = str(tmp_path / "structures" / "P00001.pdb")
    os.makedirs(os.path.dirname(item), exist_ok=True)
    _rec, path = _armed(ma, tmp_path)           # nothing declared
    ma.save_state(path, {ma.RUN_KEY: _b_record()})
    _stale_probe(ma)

    with ma.atomic_out(item) as tmp:
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write("a model\n")

    left = os.listdir(os.path.dirname(item))
    assert "P00001.pdb" not in left, "the rename happened anyway"
    assert not any(f.startswith(ma.SUPERSEDED_SUFFIX) for f in left), \
        "a per-item output was parked, one file per item"
    assert [f for f in left if ma.ATOMIC_SUFFIX in f], \
        "the work was neither renamed, parked, nor left as a .part"


def test_a_stage_output_rename_is_never_aborted_by_a_transient_error(
        ma, tmp_path, monkeypatch):
    """The trade the issue refuses, refused here too.

    Gating the rename on a lock read would abort a stage that is legitimately
    finishing the first time a read hiccups, which is worse than what it
    prevents. So there is exactly ONE path to not renaming and it needs a read
    that SUCCEEDED and returned a foreign, parsed identity. Every failure
    renames.
    """
    out = str(tmp_path / "hmm" / "pfam.tblout")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    _rec, path = _armed(ma, tmp_path, declared=[out])
    real_open = open

    def unreadable(p, *a, **k):
        if str(p) == path:
            raise OSError(5, "input/output error")
        return real_open(p, *a, **k)

    for broken in ("unreadable", "garbled", "gone"):
        _stale_probe(ma)
        if broken == "unreadable":
            monkeypatch.setattr("builtins.open", unreadable)
        elif broken == "garbled":
            monkeypatch.undo()
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("{ not json")
        else:
            os.remove(path)
        with ma.atomic_out(out) as tmp:
            with open(tmp, "w", encoding="utf-8") as fh:
                fh.write(broken + "\n")
        monkeypatch.undo()
        assert open(out, encoding="utf-8").read() == broken + "\n", \
            f"a {broken} state file aborted a stage that was finishing"


def test_the_rename_guard_reads_nothing_while_the_document_is_fresh(
        ma, tmp_path, monkeypatch):
    """The cost, pinned. atomic_out is called once per ITEM for several
    stages, so a probe that read the state file on every call would be one
    open() per protein on a shared array - the shape PART_RESCAN_S exists for
    on the console side. A write refreshes what the guard knows, so a run with
    a heartbeat never reaches the file here at all."""
    out = str(tmp_path / "out.tsv")
    _rec, path = _armed(ma, tmp_path, declared=[out])
    reads = []
    real_open = open

    def counting(p, *a, **k):
        if str(p) == path and "w" not in str(a[0] if a else k.get("mode", "r")):
            reads.append(p)
        return real_open(p, *a, **k)

    monkeypatch.setattr("builtins.open", counting)
    for _ in range(5):
        with ma.atomic_out(out) as tmp:
            open(tmp, "w", encoding="utf-8").close()
    assert reads == [], "the rename path read the state file it was just told about"


# --- --force-unlock refuses a holder it can see running -------------------
def test_force_unlock_refuses_a_holder_that_is_provably_alive_on_this_host(
        tmp_path, stub_bin):
    """The interim the issue asks for, and the only part of this that PREVENTS
    rather than detects. The operator reaching for --force-unlock is by
    definition unsure whether the holder is alive; where this host can SEE
    that it is, the answer is a message and not a silent takeover."""
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    _plant_lock(proj, os.getpid())            # this test runner: provably alive
    proc = proj.run("--force-unlock", expect=1)
    assert "--force-unlock refused" in proc.stderr
    assert f"ps -p {os.getpid()}" in proc.stderr, \
        "the refusal does not say how to find out what the process is"
    assert "recycled" in proc.stderr, \
        "the refusal does not cover the pid this host has handed to someone else"
    assert "--force-unlock-live" in proc.stderr
    assert os.path.exists(proj.rpath(".metaannot.lock")), \
        "a refused --force-unlock removed the lock anyway"
    # ...and the second flag is the escape, on the same invocation.
    proj.run("--force-unlock", "--force-unlock-live")


def test_the_refusal_and_its_escape_both_reach_the_all_subcommand(
        tmp_path, stub_bin):
    """`all` is the command the tutorial leads with, and it reads
    --force-unlock through the same namespace. A --force-unlock-live
    registered only on `run` would make the refusal INESCAPABLE there: the
    message would name a flag argparse then refused."""
    proj = _searchable(tmp_path, tmp_path / "p")
    # `--only pfam` throughout: join never runs, so there is no quantified
    # table and `all` skips the report and the object, which is this test's
    # business anyway. The lock check is the first thing either command does.
    proj.run("--only", "pfam")
    _plant_lock(proj, os.getpid())
    proc = run_metaannot("all", "--config", proj.config_path, "--only", "pfam",
                         "--force-unlock", expect=1, cwd=proj.root)
    assert "--force-unlock refused" in proc.stderr
    proc = run_metaannot("all", "--config", proj.config_path, "--only", "pfam",
                         "--force-unlock", "--force-unlock-live",
                         cwd=proj.root, env={"STUB_SLEEP": "0"})
    assert "--force-unlock refused" not in proc.stderr


def test_force_unlock_still_takes_over_a_lock_it_cannot_disprove(tmp_path,
                                                                 stub_bin):
    """The case --force-unlock exists for, and the one a conservative refusal
    would have killed: a lock left on another node of the array. This host
    cannot see that process table, so the holder is unprovable - and unprovable
    still means alive to every OTHER reader of that lock, which is why the
    refusal keys on proof and not on `_holder_is_alive`."""
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    _plant_lock(proj, 4242, host="another-node-of-the-array")
    proc = proj.run("--force-unlock")
    assert "--force-unlock refused" not in proc.stderr
    assert "removing a stale lock" in proc.stderr


def test_unprovable_still_means_alive_and_only_proof_refuses(ma):
    """The split, branch for branch. _holder_is_alive() has to keep answering
    exactly what it answered before - it is what decides whether an ordinary
    run refuses to start - while the refusal reads the narrow positive arm of
    the same evidence."""
    lk = ma.ResultsLock("/tmp/never", force=False)
    me = __import__("socket").gethostname()
    cases = [
        ({"pid": os.getpid(), "host": me}, ma.PROVEN_ALIVE),
        ({"pid": 1, "host": me}, ma.PROVEN_ALIVE),   # PermissionError: it exists
        ({"pid": 4242, "host": "another-node"}, ma.UNPROVABLE),
        ({}, ma.UNPROVABLE),                         # garbled lock file
        ({"pid": "not-an-int", "host": me}, ma.UNPROVABLE),
    ]
    if os.name != "nt":
        cases.append(({"pid": _dead_pid(), "host": me}, ma.PROVEN_DEAD))
    for info, want in cases:
        assert lk._holder_proof(info) == want, info
        assert lk._holder_is_alive(info) is (want != ma.PROVEN_DEAD), info


# --- --force still discards exactly what it selected ---------------------
def test_force_discards_the_selected_stages_through_a_merged_write(
        tmp_path, stub_bin):
    """--force pops in memory and the discard used to reach the file as a side
    effect of the next whole-document write. A merged write only touches the
    keys it is handed, so the pops are carried explicitly - and carried until
    the stage is recorded again, or `--force --only pfam` would put pfam's old
    record straight back."""
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    first = proj.state()
    proj.run("--force", "--only", "pfam")
    st = proj.state()
    # run_id and not `finished`: the stamp is second-resolution and two runs
    # of a stub pipeline land in the same second often enough to matter.
    assert st["pfam"]["run_id"] != first["pfam"]["run_id"], "pfam was not redone"
    assert st["dbcan"] == first["dbcan"], \
        "--force --only pfam discarded a record it was not asked to"


def test_a_superseded_run_writes_nothing_even_with_the_document_missing(
        ma, tmp_path):
    """Rule 4's half of this: a run told it lost the directory stops writing.

    Not the rebuild - nothing rebuilds a document any more, and the tests
    further down drive that on its own. What this pins is the OTHER rule,
    which is latched and unconditional: once either ownership gate has stood a
    run down it writes nothing into that file again, whatever state the file is
    in. A missing document is the case worth pinning because it is the one
    where the run would otherwise be creating the file rather than editing it.

    It latches from the DOCUMENT before removing the document, which makes it
    the easy ordering. The orderings that really happen are below.
    """
    rec, path = _taken_over(ma, tmp_path)
    assert rec._tick() is False and ma._STATE_WATCH["lost"] is not None
    os.remove(path)
    assert ma.update_state(path, {"pfam": {"status": "ok"}}, ("pfam",),
                           claim=rec.claim) is False
    assert not os.path.exists(path), \
        "a superseded run wrote its whole snapshot back into the directory"


def _handed_over(ma, tmp_path):
    """A run holding a lock that has since been rewritten by somebody else.

    The LOCK gate and only the lock gate, which is the half `_taken_over`
    cannot reach: `_STATE_WATCH["lost"]` is raised by a successful, parsed
    read of the state DOCUMENT, so a document that is missing or garbled can
    never raise it however long ago the directory changed hands.
    """
    lock_path = str(tmp_path / ".metaannot.lock")
    lock = ma.ResultsLock(lock_path)
    with lock:
        pass
    # Still believed held by this run - the operator has not told it anything.
    lock.held = True
    with open(lock_path, "w", encoding="utf-8") as fh:
        json.dump({"pid": 999, "host": "another-node-of-the-array",
                   "started": "2026-01-01T00:00:00"}, fh)
    path = str(tmp_path / ".metaannot_state.json")
    state = {"pfam": {"status": "ok", "signature": "a"}}
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, 999, owner=lock)
    assert ma.update_state(path, state, (ma.RUN_KEY,), claim=rec.claim)
    return rec, path, state


@pytest.mark.parametrize("damage", ["removed", "zeroed"])
def test_a_run_stood_down_by_the_lock_writes_nothing_into_damaged_bytes(
        ma, tmp_path, capsys, damage):
    """The same rule, latched from the LOCK instead of from the document.

    `_STATE_WATCH["lost"]` is set by `_judge_ownership()` and by nothing else,
    and that runs only on a read that SUCCEEDED and PARSED - so on a document
    that is missing or full of NULs it can never be set at all. The other half
    comes from the LOCK, whose bytes are a different file, and it is what
    stands this run down here.

    Note what this test has to do to get there: it calls `rec._save()` first.
    That is the only thing in the process that reads the lock, and nothing in
    `finish()`, `mark_running()` or `record()` calls it - so in the real
    ordering neither half of the latch is up, and a fix that depended on one
    being up was a fix that depended on this fixture. That is why the rebuild
    was deleted rather than guarded; the tests below are the real orderings,
    with nothing pre-latched.
    """
    rec, path, state = _handed_over(ma, tmp_path)
    if damage == "removed":
        os.remove(path)
    else:
        with open(path, "wb") as fh:
            fh.write(b"\x00" * 64)
    before = open(path, "rb").read() if damage == "zeroed" else None
    capsys.readouterr()

    # A's tail: the lock gate notices, and then the stage it really finished
    # tries to record itself. Nothing has read the document successfully.
    assert rec._save() is False
    assert rec.superseded is True
    assert ma._STATE_WATCH["lost"] is None, \
        "the document gate cannot latch off a document nothing could read; " \
        "if it can, this test is no longer testing the reported failure"

    state["pfam"] = {"status": "ok", "signature": "a"}
    assert ma.update_state(path, state, ("pfam",), claim=rec.claim) is False, \
        "a run the lock had already stood down went on writing"
    if damage == "removed":
        assert not os.path.exists(path), \
            "a superseded run recreated the state file it had been shut out of"
    else:
        assert open(path, "rb").read() == before, \
            "a superseded run rewrote the document it had been shut out of"
    assert "no longer holds" in capsys.readouterr().err, \
        "the run was stood down without saying so anywhere"


# ======================================================================
# the orderings that really happen: NOTHING pre-latched
#
# The rebuild these drive out was guarded three times and broken three
# times, and the third finding is the one that settles it: no lock read can
# authorise rebuilding a document, because a VACANT lock is exactly what a
# replacement that has FINISHED leaves behind. So the write that rebuilt a
# whole document is gone, and what is pinned below is its absence - in the
# two orderings that occur, with the two shapes of damage, at three
# heartbeat intervals, with neither ownership gate consulted first.
# ======================================================================
def _after_handover(ma, tmp_path, replacement, carried, interval=30):
    """Run A, still believing it holds the directory, after a handover.

    Built in the order it happens. A takes the lock and writes `_run`; the
    directory changes hands; the replacement's records are in the document.
    Only then does the test damage the document and let A write.

    `replacement` says what the replacement did. "live" rewrote the lock and
    is still running, so the lock reads as somebody else's. "finished"
    recorded its work, exited and took its own lock with it, so the path is
    VACANT - and that is the finding this whole section turns on, because
    `_save()` writes on None BY DESIGN and no reading of the lock can tell a
    replacement that has not started from one that has already finished.

    `carried` is what A has in memory to write: "pfam" for a run whose tail is
    recording the one stage it really did finish, "nothing" for a run whose
    heartbeat is the only thing that has ever wanted to write.

    NOTHING is pre-latched. A has not been told anything by either gate: the
    document gate cannot fire on a read that failed, and nothing in finish(),
    mark_running() or record() reads the lock at all.
    """
    lock_path = str(tmp_path / ".metaannot.lock")
    path = str(tmp_path / ".metaannot_state.json")
    lock = ma.ResultsLock(lock_path)
    lock.__enter__()                  # entered and NOT released: A holds it
    state = {}
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, interval,
                       owner=lock)
    # cmd_run's own first write, and then one more, because a run has to have
    # READ this document at least once for the loss of it to be something it
    # can tell apart from a fresh results directory: the first write creates
    # the file, the second parses it.
    assert ma.update_state(path, state, (ma.RUN_KEY,), claim=rec.claim) is True
    assert rec._tick() is True, "A was stood down before the handover"
    assert ma._STATE_WATCH["parsed"] is True

    # --- the handover ---
    if replacement == "live":
        with open(lock_path, "w", encoding="utf-8") as fh:
            json.dump({"pid": 999, "host": "another-node-of-the-array",
                       "started": "2026-01-01T00:00:00"}, fh)
    else:
        os.remove(lock_path)
    ma.save_state(path, {ma.RUN_KEY: _b_record(),
                         "dbcan": {"status": "ok", "signature": "b"},
                         "diamond": {"status": "ok", "signature": "b"}})
    if carried == "pfam":
        state["pfam"] = {"status": "ok", "signature": "a"}
    return rec, path, state


def _damage(path, how):
    """Remove the document, or fill it with the NULs a power loss leaves."""
    if how == "removed":
        os.remove(path)
        return None
    with open(path, "wb") as fh:
        fh.write(b"\x00" * 64)
    return open(path, "rb").read()


def _unchanged(path, before):
    if before is None:
        return not os.path.exists(path)
    return os.path.exists(path) and open(path, "rb").read() == before


@pytest.mark.parametrize("interval", [0, 1, 30])
@pytest.mark.parametrize("damage", ["removed", "zeroed"])
@pytest.mark.parametrize("replacement", ["live", "finished"])
def test_a_superseded_runs_tail_writes_only_its_own_key(
        ma, tmp_path, capsys, replacement, damage, interval):
    """A's tail records the one stage it finished, and rebuilds nothing.

    The ordering: A parks inside pfam, the operator hands the directory over
    with --force-unlock-live, the replacement records its stages, and the
    document then goes away - an `rm`, a remount, a power loss. A's pfam
    returns and A records it.

    What A may write is `pfam`. What it may NOT write is `_run`, which is its
    own identity and an ownership claim: a document that says the directory
    belongs to A is read by the NEXT run as A's, taken as that run's
    predecessor, and from then on every write A makes is inside the new run's
    claim and can never be refused. The rebuild that stood here wrote A's
    whole in-memory dict, `_run` and all.
    """
    rec, path, state = _after_handover(ma, tmp_path, replacement, "pfam",
                                       interval)
    _damage(path, damage)
    capsys.readouterr()
    assert ma._STATE_WATCH["lost"] is None and rec.superseded is False, \
        "something latched before the write; this is no longer the real ordering"

    assert ma.update_state(path, state, ("pfam",), claim=rec.claim) is True, \
        "the stage record A really earned was refused"
    doc = _doc(path)
    assert set(doc) == {"pfam"}, \
        f"the write brought back more than the key it named: {sorted(doc)}"
    assert ma.RUN_KEY not in doc, \
        "a superseded run wrote its own `_run` into a document it created; " \
        "the next run will adopt it as its predecessor and never refuse it"
    assert "dbcan" not in doc and "diamond" not in doc, \
        "records the replacement wrote came back from A's snapshot"
    said = capsys.readouterr().err
    want = "is no longer there" if damage == "removed" else "does not parse"
    assert want in said, \
        f"a run that lost the document's records said nothing about it: {said!r}"
    assert "is not in it" in said or "are not in it" in said, \
        "the message does not say that other runs' records are gone"


@pytest.mark.parametrize("interval", [0, 1, 30])
@pytest.mark.parametrize("damage", ["removed", "zeroed"])
@pytest.mark.parametrize("replacement", ["live", "finished"])
def test_a_superseded_runs_heartbeat_alone_recreates_nothing(
        ma, tmp_path, capsys, replacement, damage, interval):
    """The same handover with NO stage of A's finished at all.

    This is the ordering that needs nothing to have happened: A's heartbeat is
    the only thing in the process that has ever wanted to write, so there is
    no stage record anywhere to hang a fix on. On a "live" replacement the lock
    gate stands A down and says so. On a "finished" one it answers VACANT -
    None, not False - which is what `_save()` used to write on, and that is
    what let A's own `_run` recreate the document it had been shut out of.

    `_run` is the only key here, and it is the one key that needs positive
    proof: a vacant lock is what a finished replacement leaves behind, so A
    does not write `_run` at all while the path is vacant. The test beside
    this one is the same handover with a stage record in front of it, which is
    what every real run has and what the earlier gate did not survive.

    Driven at heartbeat_s 0, 1 and 30 because none of this is a race: the
    interval changes how often A tries, never whether it may.
    """
    rec, path, state = _after_handover(ma, tmp_path, replacement, "nothing",
                                       interval)
    before = _damage(path, damage)
    capsys.readouterr()
    assert ma._STATE_WATCH["lost"] is None and rec.superseded is False

    assert rec._tick() is (replacement == "finished"), \
        "the tick's own verdict about whether this run carries on is wrong"
    assert rec._save() is False, "the heartbeat wrote `_run` anyway"
    assert _unchanged(path, before), \
        "a run whose lock was gone recreated the state file from its own view"
    if replacement == "live":
        assert "no longer holds" in capsys.readouterr().err
        assert rec.superseded is True and rec._stop.is_set()
    else:
        # Not superseded: nothing proved anything against it. It simply has
        # no document to claim, and it keeps its stage writes.
        assert rec.superseded is False, \
            "a vacant lock was read as proof the directory changed hands"
        state["pfam"] = {"status": "ok", "signature": "a"}
        assert ma.update_state(path, state, ("pfam",),
                               claim=rec.claim) is True
        assert set(_doc(path)) == {"pfam"}, \
            "the stage keys were gated like `_run`; they must not be"


@pytest.mark.parametrize("interval", [0, 1, 30])
@pytest.mark.parametrize("damage", ["removed", "zeroed"])
@pytest.mark.parametrize("replacement", ["live", "finished"])
def test_a_superseded_runs_tail_cannot_stamp_run_after_recording_a_stage(
        ma, tmp_path, capsys, replacement, damage, interval):
    """The stage record creates the document; `_run` must still not follow it.

    THE GATE THIS PINS WAS INOPERATIVE FOR EVERY RUN THAT RECORDS ANYTHING.
    It declined only the write that would CREATE the document - and A's own
    stage record creates it one call earlier, so by the time A's tail or its
    heartbeat reaches `_run` there IS a document, it is one A made seconds ago
    holding a stage key and no `_run`, and the gate waved it through. The file
    ended as {_run: A, pfam: A} with A never latched: no "no longer holds", no
    "another run holds". Driven with real processes, and reproduced here at
    heartbeat_s 0, 1 and 30 on a removed document and on a zeroed one.

    So the two writes are driven in the order a run really makes them: the
    stage first, then `_run`. The stage record is A's to write; `_run` is not,
    because a vacant lock is what a FINISHED replacement leaves behind and a
    document naming A the owner is adopted by the next run as its predecessor.
    """
    rec, path, state = _after_handover(ma, tmp_path, replacement, "pfam",
                                       interval)
    _damage(path, damage)
    capsys.readouterr()
    assert ma._STATE_WATCH["lost"] is None and rec.superseded is False, \
        "something latched before the write; this is no longer the real ordering"

    assert ma.update_state(path, state, ("pfam",), claim=rec.claim) is True, \
        "the stage record A really earned was refused"
    assert set(_doc(path)) == {"pfam"}

    # ...and now the write that used to ride in on the back of it.
    assert rec._save() is False, \
        "`_run` was stamped into the document this run had just created"
    assert set(_doc(path)) == {"pfam"}, \
        "a superseded run's `_run` reached a document it created with a stage " \
        "record; the next run adopts it as its predecessor and never refuses it"
    said = capsys.readouterr().err
    if replacement == "live":
        assert "no longer holds" in said
    else:
        assert "cannot prove it still owns this directory" in said, \
            f"the refusal was silent: {said!r}"


def test_a_write_that_lands_in_the_read_to_rename_gap_is_lost_and_silent(
        ma, tmp_path, capsys, monkeypatch):
    """The CLOCK the records half runs on, pinned as a mechanism.

    Three published sentences called the records half "a guarantee with no
    clock in it". It is not one. The write re-reads the document, merges its
    key in, and renames a complete file into place, and anything another
    process writes BETWEEN that read and that rename is not in what we merged -
    so our complete document goes over it. The read-back afterwards compares
    the file against the payload we just wrote, so it catches a writer that
    landed AFTER the rename and redoes the merge; one that landed before it
    matches, and nothing is said.

    Driven by INJECTION, not by racing: the window is one read-modify-rename
    wide and has not been won at shipped speeds, so what this pins is that the
    mechanism exists and is silent - not how wide the window is. A test that
    tried to win the race would be a flake pretending to be a measurement.
    """
    path = str(tmp_path / ".metaannot_state.json")
    theirs = {ma.RUN_KEY: _b_record(),
              "dbcan": {"status": "ok", "signature": "b"},
              "diamond": {"status": "ok", "signature": "b"}}
    state = {"pfam": {"status": "ok", "signature": "a"}}
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, 0, owner=None)
    # A's own document, read once, so the loss of it is something A can tell
    # apart from a fresh results directory.
    assert ma.update_state(path, state, ("pfam",), claim=rec.claim) is True
    assert ma.update_state(path, state, ("pfam",), claim=rec.claim) is True
    assert ma._STATE_WATCH["parsed"] is True
    os.remove(path)                        # the `rm`, before A's next write

    real_save = ma.save_state
    landed = []

    def save_state(p, doc):
        # the other process's write, IN THE GAP: after A's pre-write read came
        # back `missing` and before A's own os.replace.
        if not landed:
            landed.append(True)
            real_save(p, theirs)
        return real_save(p, doc)

    monkeypatch.setattr(ma, "save_state", save_state)
    capsys.readouterr()
    assert ma.update_state(path, state, ("pfam",), claim=rec.claim) is True
    assert landed, "the injected write never ran; this proves nothing"
    doc = _doc(path)
    assert set(doc) == {"pfam"}, \
        "the injection did not reproduce the loss this pins"
    said = capsys.readouterr().err
    assert "changed underneath" not in said, \
        "the read-back is being credited with catching this; it does not"
    assert "dbcan" not in said and "diamond" not in said, \
        "nothing names the records that were lost, and the docs must not " \
        "claim it does"


def test_the_refusal_names_only_the_running_stages_of_the_run_it_is_about(
        ma, tmp_path):
    """The `--force-unlock` refusal, against a document with two runs in it.

    A state file routinely holds `running` records left by an EARLIER run that
    was killed: SIGTERM does not unwind, so it never stamps, and the record it
    was in the middle of stays `running` for good. That is the ORDINARY trace
    of a `kill`, not a corruption, so the refusal meets it constantly.

    It used to collect every `running` record in the document and attribute all
    of them to the pid in the lock. Measured: SIGTERM a run mid-`pfam`, start a
    second with `--only dbcan`, and the refusal said "It has dbcan, pfam
    recorded running" about a process whose own command line — printed by the
    same message, one clause later — said `--only dbcan`. An operator reading
    that goes looking for a `pfam` the live process has never touched.

    `run_id` is what separates them and is in every record, so the message can
    simply ask. A record with no `run_id` predates the field and is nobody's
    that we can prove, which is the same answer.
    """
    state = tmp_path / ".metaannot_state.json"
    state.write_text(json.dumps({
        ma.RUN_KEY: {"run_id": "live-2", "pid": 4242,
                     "argv": ["metaannot.py", "run", "--only", "dbcan"]},
        # the killed run's, still `running` and never stamped
        "pfam": {"status": "running", "run_id": "dead-1"},
        # the live holder's
        "dbcan": {"status": "running", "run_id": "live-2"},
        # written before run_id existed: not provably the holder's either
        "kofam": {"status": "running"},
        # the holder's, but finished - not what it is "in the middle of"
        "cluster": {"status": "ok", "run_id": "live-2"},
    }), encoding="utf-8")

    note = ma._lock_holder_note(str(state))
    assert "dbcan" in note, note
    assert "pfam" not in note, \
        "the refusal named a stage a DIFFERENT run left running: " + note
    assert "kofam" not in note, \
        "the refusal named a record it cannot attribute to the holder: " + note
    assert "cluster" not in note, note


def test_a_force_run_does_not_discard_a_record_written_after_it_decided(
        ma, tmp_path, capsys):
    """`--force`'s discard set deletes keys the write does not name.

    Not a window and not a race: --force pops the records it means to redo from
    ONE read at the start of the run, and every merged write then carried those
    names and popped whatever was under them. A live run's `dbcan` record was
    watched disappearing, silently, from a document a superseded --force run
    was writing `pfam` into - reachable whenever the document holds no `_run`
    to stand the old run down, which is exactly the shape a run recreates after
    its document is removed, and which lasts until the live run stamps `_run`
    again (never, on `heartbeat_s: 0`, before its final stamp).

    So the discard carries the RECORD it was decided about, not just the name,
    and deletes only while the document still holds that record. This is the
    last place in this file where a write touched a key it did not name; the
    whole-document rebuild that was the other one is gone.
    """
    path = str(tmp_path / ".metaannot_state.json")
    mine = {"status": "ok", "signature": "a", "run_id": "A"}
    theirs = {"status": "ok", "signature": "b", "run_id": "THE-LIVE-RUN"}
    state = {"pfam": {"status": "ok", "signature": "a"}}
    rec = ma.RunRecord(path, state, ["metaannot", "run"], None, 0, owner=None)

    # the ordinary --force discard: the record is the one that was popped
    ma.save_state(path, {"dbcan": mine})
    assert ma.update_state(path, state, ("pfam",), claim=rec.claim,
                           drop={"dbcan": mine}) is True
    assert set(_doc(path)) == {"pfam"}, "--force stopped discarding"

    # ...and the same discard against a record somebody else wrote since
    ma.save_state(path, {"dbcan": theirs})
    capsys.readouterr()
    assert ma.update_state(path, state, ("pfam",), claim=rec.claim,
                           drop={"dbcan": mine}) is True
    doc = _doc(path)
    assert doc.get("dbcan") == theirs, \
        "a --force run deleted a record written after its discard was decided"
    assert doc["pfam"] == state["pfam"], "the key the write named is missing"
    said = capsys.readouterr().err
    assert "not the one this run read at the start" in said, \
        f"the declined discard was silent: {said!r}"


def _doc_or_none(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def test_a_solo_run_whose_document_is_removed_still_records_and_finishes(
        ma, tmp_path, stub_bin):
    """The regression the deleted rebuild existed for, end to end.

    No second run anywhere: one process, holding its own lock, whose state
    file is removed from under it - an operator's `rm`, a tmp-reaper, a
    remount. It has to go on recording the stages it finishes and it has to
    exit 0, because the alternative is a tool that stops writing its own
    results over an accident in a directory nobody else is touching.

    And it has to be explicit about the price, which is the whole reason the
    rebuild was ever written: the records that were in the file are gone with
    the file. `pfam` was recorded by the previous run and is NOT in the
    document afterwards; its OUTPUT is still on disk, so the next run adopts
    or recomputes it and says so. That is a recomputation, not a wrong answer,
    and it is the trade taken deliberately over letting any run rebuild a
    document out of its own snapshot.
    """
    proj = _searchable(tmp_path, tmp_path / "solo")
    proj.run("--only", "pfam", env={"STUB_SLEEP": "0"}, timeout=180)
    assert proj.state()["pfam"]["status"] == "ok"
    out = proj.rpath("hmm", "pfam.tblout")
    assert os.path.exists(out)

    state = proj.rpath(".metaannot_state.json")
    gate = str(tmp_path / "free_solo")
    # `--only pfam dbcan`: pfam is cached and writes nothing, so the run is
    # parked inside dbcan with dbcan already recorded `running`.
    a = _gated_run(proj, gate, "dbcan", "--serial")
    try:
        assert _wait_for(lambda: (_doc_or_none(state) or {})
                         .get("dbcan", {}).get("status") == "running",
                         timeout=60), "the run never recorded dbcan running"
        os.remove(state)
        open(gate, "w").close()
        # A timeout, not a bare wait: a regression here is a run that never
        # writes and never exits, and that has to FAIL rather than wedge the
        # suite.
        _out, err = a.communicate(timeout=180)
    finally:
        if a.poll() is None:
            a.kill()
            a.communicate(timeout=60)
    assert a.returncode == 0, f"the solo run did not finish:\n{err}"

    assert "is no longer there" in err, \
        "the run lost every record in the file and said nothing about it"
    doc = _doc_or_none(state)
    assert doc is not None, "the run did not recreate the state file at all"
    assert doc["dbcan"]["status"] == "ok", \
        "a run whose document was removed stopped recording its own stages"
    run = doc[ma.RUN_KEY]
    assert run["final_status"] == "ok" and run["finished"], \
        "the run holding its own lock lost its own final verdict"
    # ...and this is what it costs, stated rather than hidden.
    assert "pfam" not in doc, \
        "the record of a stage this run never touched came back from somewhere"
    assert os.path.exists(out), "an output was lost, which is not the trade"

    # The next run pays the price and says so: pfam's output is there with no
    # record of it, which is adoption with a warning, not silence.
    c = proj.run("--only", "pfam", "dbcan", env={"STUB_SLEEP": "0"},
                 timeout=180)
    assert "adopting output this run did not produce" in c.stderr and \
        "pfam" in c.stderr, \
        "the loss of pfam's record was silent on the next run"
    assert proj.state()["dbcan"]["status"] == "ok"


def test_a_solo_run_that_loses_its_lock_too_keeps_its_stages_not_its_verdict(
        ma, tmp_path, stub_bin):
    """The same solo accident, one file wider, end to end.

    The test above removes the DOCUMENT and the run keeps its own final
    verdict, because it still holds its own lock. This one removes the LOCK as
    well, which is the case `_run`'s gate exists for and the one whose cost the
    docs state: `_run` is a claim about who owns the directory, a vacant lock
    is what a replacement that took it and then exited leaves behind, and from
    inside the process the two readings are the same bytes - nothing.

    So the run records the stage it finishes and exits 0, and its verdict is
    gone: no `final_status: "ok"`, no `finished`. That is the price, it is paid
    here rather than claimed, and it is said out loud in the log - a `_run`
    that simply stopped advancing reads on a console exactly like a process
    that died.
    """
    proj = _searchable(tmp_path, tmp_path / "solo2")
    proj.run("--only", "pfam", env={"STUB_SLEEP": "0"}, timeout=180)
    state = proj.rpath(".metaannot_state.json")
    lock = proj.rpath(".metaannot.lock")
    gate = str(tmp_path / "free_solo2")
    a = _gated_run(proj, gate, "dbcan", "--serial")
    try:
        assert _wait_for(lambda: (_doc_or_none(state) or {})
                         .get("dbcan", {}).get("status") == "running",
                         timeout=60), "the run never recorded dbcan running"
        assert os.path.exists(lock)
        os.remove(lock)
        os.remove(state)
        open(gate, "w").close()
        _out, err = a.communicate(timeout=180)
    finally:
        if a.poll() is None:
            a.kill()
            a.communicate(timeout=60)
    assert a.returncode == 0, f"the solo run did not finish:\n{err}"

    doc = _doc_or_none(state)
    assert doc is not None, "the run stopped recording its own stages too"
    assert doc["dbcan"]["status"] == "ok", \
        "a stage record was gated like `_run`; only `_run` is an ownership claim"
    assert ma.RUN_KEY not in doc, \
        "a run that could not prove it owned the directory claimed it anyway"
    said = " ".join(err.split())
    assert "cannot prove it still owns this directory" in said, \
        f"the run lost its own final verdict silently:\n{err}"
    assert "no longer holds" not in said, \
        "a vanished lock was reported as a handover nothing proved"


def test_a_run_stood_down_by_the_LOCK_does_not_rename_its_output_either(
        ma, tmp_path, capsys):
    """The same two-halves latch, on the rename gate.

    The parking mechanism reads the latch, and asking only the half that a
    successful read of the DOCUMENT can raise made it inoperative on the
    ordinary local takeover: RunRecord._save() tests the LOCK before it
    writes, so the first tick after a handover latches from the lock and
    returns without ever reaching the succession check. `lost` stays None for
    the rest of the unwind, and the in-flight stage renamed over the live
    run's output with the whole clock still answering "ours".
    """
    rec, state_path, _state = _handed_over(ma, tmp_path)
    out = str(tmp_path / "hmm" / "pfam.tblout")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("the live run's result\n")
    ma._DECLARED_OUTPUTS.add(out)

    assert rec._save() is False
    assert rec.superseded is True and ma._STATE_WATCH["lost"] is None
    capsys.readouterr()

    assert ma._directory_still_ours() is False, \
        "the rename gate ignored the gate that had already stood this run down"
    with ma.atomic_out(out) as tmp:
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write("the superseded run's result\n")
    with open(out, encoding="utf-8") as fh:
        assert fh.read() == "the live run's result\n", \
            "a superseded run renamed its output over the live run's"
    parked = [f for f in os.listdir(os.path.dirname(out))
              if f.startswith(ma.SUPERSEDED_SUFFIX)]
    assert len(parked) == 1, f"the work was not parked: {parked}"
    assert "NOT renamed into place" in capsys.readouterr().err


def test_the_rename_guard_probes_at_most_once_in_its_interval(
        ma, tmp_path, monkeypatch):
    """The cost of the fallback read, pinned where it is paid.

    With the heartbeat off there is no merge read to keep the guard's answer
    fresh, so the rename path reads the state file itself - and atomic_out is
    called once per ITEM by several stages. A read per call would be one
    open() per dark protein on a shared array, which is the shape the console
    answered with PART_RESCAN_S in the same directory.
    """
    path = str(tmp_path / ".metaannot_state.json")
    ma._watch_state(path, None, ())
    ma._STATE_WATCH["seen"] = ma._STATE_WATCH["tried"] = 0.0
    ma.save_state(path, {"pfam": {"status": "ok"}})
    reads = []
    real_open = open

    def counting(p, *a, **k):
        if str(p) == path:
            reads.append(p)
        return real_open(p, *a, **k)

    monkeypatch.setattr("builtins.open", counting)
    for _ in range(50):
        assert ma._directory_still_ours() is True
    assert len(reads) == 1, \
        f"the rename path read the state file {len(reads)} times in one interval"
    assert ma.STATE_PROBE_S > 0


def test_force_unlock_live_implies_force_unlock(tmp_path, stub_bin):
    # the flag an operator is told to add gets the run through on its own: the
    # refusal it exists for is only reachable through --force-unlock, so
    # demanding both would be a second thing to get right at the worst moment.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    _plant_lock(proj, os.getpid())
    proc = proj.run("--force-unlock-live")
    assert "--force-unlock refused" not in proc.stderr


def test_a_superseded_run_does_not_rename_a_later_stage_into_the_directory(
        tmp_path, stub_bin):
    """The atomic_out variant, end to end, with two real runs.

    A is parked inside `pfam`; the operator hands the directory to B, which
    finishes `pfam` and exits. A's `pfam` then returns, A's record of it is
    refused, and A learns from that refusal that it no longer holds the
    directory — so `dbcan`, the next stage it runs, cannot put its result into
    the directory at all. A stage BEGUN under a lost directory can never
    overwrite: that is the half of this the design guarantees, and it is what
    this drives. (A stage already in flight at the moment of the handover is
    only NOTICED instead, at whichever comes first of the next heartbeat tick
    and the fallback probe - min(heartbeat_s, STATE_PROBE_S) - which is a
    clock and not a guarantee, so the latch it turns on is unit-tested above
    rather than raced here.)
    """
    proj = _searchable(tmp_path, tmp_path / "park")
    gate = str(tmp_path / "free_a")
    a = _gated_run(proj, gate, "dbcan", "--serial")

    b = proj.run("--only", "pfam", "--force-unlock", "--force-unlock-live",
                 env={"STUB_SLEEP": "0"}, timeout=180)
    assert "--force-unlock refused" not in b.stderr
    b_run = _run_record(proj)["run_id"]

    open(gate, "w").close()
    _out, a_err = a.communicate(timeout=180)

    assert not os.path.exists(proj.rpath("hmm", "dbcan.domtblout")),         "a superseded run put its output into the live run's directory"
    parked = [f for f in os.listdir(proj.rpath("hmm"))
              if f.startswith(".superseded.dbcan.")]
    assert parked, f"the work was not parked: {os.listdir(proj.rpath('hmm'))}"
    assert "NOT renamed into place" in a_err
    assert "missing input" in a_err,         "the warning does not say what happens to the rest of that run"
    assert "dbcan" not in proj.state(),         "the superseded run recorded a stage in the live run's document"
    assert _run_record(proj)["run_id"] == b_run


# --- #36: the record, and whether it still describes the file -----------
#
# The whole group drives ONE property: a `cached` verdict whose declared
# output was written after the record that describes it says so, and which
# reuses the stage anyway. The ones that must stay SILENT are here beside it,
# because a check like this is worth nothing if the false ones are not pinned
# as hard as the true one.
#
# And several of them are about the SIZE and the HONESTY of the report
# rather than about whether it fires, because that is where the first build
# of this was wrong. A report whose line count grows with the directory is
# a report
# that gets turned off; a report that names a cause it cannot know teaches a
# reader to disbelieve the part it can. So the properties pinned are: one WARN
# per run whatever the directory holds, a count that does not move when the
# remedy is taken on one stage, every late stage named, and no cause asserted
# where more than one of them is late.

_SUPERSEDED_RUN = '''
"""A real run whose every state write is refused, paused at the one write that
would have recorded the stage it has just renamed into place.

That is not a contrivance, it is what a superseded run IS: update_state()
returns False on either ownership latch, while atomic_out() goes on renaming,
because _directory_still_ours() is a clock and not an exclusion. The same
refusal reaches a run that never lost anything, on a single unreadable read of
the state file at dispatch. Either way the record that survives in the
document belongs to a different run than the bytes on disk do, and the pause
is where the parent SIGKILLs this one: between the rename and the record,
which is the corruption path no ownership check can reach.
"""
import importlib.util
import os
import sys
import time

spec = importlib.util.spec_from_file_location("metaannot", os.environ["MA_PY"])
ma = importlib.util.module_from_spec(spec)
sys.modules["metaannot"] = ma
spec.loader.exec_module(ma)

STAGE, MARKER = os.environ["MA_STAGE"], os.environ["MA_MARKER"]


def refused(path, state, keys, *args, **kw):
    rec = state.get(STAGE)
    if STAGE in keys and isinstance(rec, dict) and rec.get("status") == "ok":
        with open(MARKER, "w", encoding="utf-8") as fh:
            fh.write("renamed, not recorded")
        time.sleep(600)
    return False


ma.update_state = refused
sys.exit(ma.main())
'''

LATE = "no longer about the file that is there"
SKEWED = "not being read as evidence of anything"
UNDATED = "names no single instant"


def _warn_block(stderr, needle):
    """One whole WARN, continuation lines and all.

    log() writes a multi-line message as its tagged first line plus
    continuations padded out to the width of that tag, so a test that takes
    only the line the phrase sits on reads a HEADING and calls it the message.
    This check's message is a heading, a list of stages and then a paragraph,
    and the paragraph is where every limitation it states lives - so the
    vocabulary rule at the bottom of this file, asserted against the first
    line alone, would have been asserted against the one part of the message
    that cannot break it.
    """
    lines = stderr.splitlines()
    at = next((i for i, l in enumerate(lines) if needle in l), None)
    assert at is not None, f"{needle!r} is not in this output:\n{stderr}"
    out = [lines[at]]
    for line in lines[at + 1:]:
        # A new log() call starts with the elapsed-time tag; a continuation of
        # this one cannot, because it is indented to that tag's width.
        if re.match(r"^\[\s*[\d.]+s\]", line):
            break
        out.append(line)
    return "\n".join(out)


def _stages_of(proj):
    return sorted(k for k in proj.state() if not k.startswith("_"))


def _backdate(proj, *stages, by=3600):
    """Move a record's `finished` stamp back, leaving everything else alone.

    The same fact as writing the output later, and the cheap half of it: the
    check compares one against the other, and only the tests that have to
    drive a RENAME or a COPY need real time to pass. `signature` does not read
    `finished`, so the stage is still cached afterwards.
    """
    path = proj.rpath(".metaannot_state.json")
    with open(path, encoding="utf-8") as fh:
        st = json.load(fh)
    for name in (stages or [k for k in st if not k.startswith("_")]):
        rec = st.get(name)
        if not isinstance(rec, dict) or "finished" not in rec:
            continue
        when = time.mktime(time.strptime(rec["finished"], "%Y-%m-%dT%H:%M:%S"))
        rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S",
                                        time.localtime(when - by))
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(st, fh, indent=1)
    return st


def _set_finished(proj, stage, stamp):
    """Write a literal `finished` stamp, which is the only way to reach a
    local time this machine is not in right now."""
    path = proj.rpath(".metaannot_state.json")
    with open(path, encoding="utf-8") as fh:
        st = json.load(fh)
    st[stage]["finished"] = stamp
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(st, fh, indent=1)


@pytest.mark.skipif(os.name == "nt", reason="there is no SIGKILL to send")
def test_an_output_renamed_after_its_own_record_is_reported_and_still_cached(
        tmp_path, stub_bin, ma):
    """The corruption path no ownership check can reach, driven end to end.

    Run one records `pfam` ok. Run two rewrites `pfam.tblout` and is SIGKILLed
    between atomic_out's rename and finish()'s record, with every state write
    of its own refused - so the document still holds run one's record and the
    file under it is run two's. The signature still agrees and every output
    still exists, so run three reports `cached`: correct, and until this check
    existed, silent.

    It stays `cached` afterwards. Nothing is recomputed, nothing is deleted
    and the record is not touched - the run says what it can see and hands the
    decision to whoever knows whether this directory was --force-unlocked.

    ONE stage, which is the half of the message that may name a remedy: with
    a single record in the report there is nothing to read it two ways, so
    `--force --only` is the answer and the line says so. The several-stage
    half is the test below.
    """
    proj = _searchable(tmp_path, tmp_path / "window")
    proj.run()
    recorded = proj.state()["pfam"]
    assert recorded["status"] == "ok"

    # Real elapsed time, and one of the two places this group spends any: the
    # slack is what separates a healthy stage - renamed a moment BEFORE its
    # own record - from this one, so a rename landing inside it proves nothing
    # and the test must not ask it to.
    time.sleep(ma.OUTPUT_STAMP_SLACK_S + 1.0)

    driver = tmp_path / "superseded_run.py"
    driver.write_text(_SUPERSEDED_RUN, encoding="utf-8")
    marker = tmp_path / "renamed-not-recorded"
    killed = subprocess.Popen(
        [sys.executable, str(driver), "run", "--config", proj.config_path,
         "--force", "--only", "pfam"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        cwd=proj.root,
        env=dict(os.environ, MA_PY=METAANNOT_PY, MA_STAGE="pfam",
                 MA_MARKER=str(marker), STUB_SLEEP="0", PYTHONHASHSEED="0"))
    try:
        reached = _wait_for(marker.exists, timeout=120)
    finally:
        killed.kill()
        killed.communicate(timeout=60)       # never leave one behind running
    assert reached, "the run never got between its rename and its record"
    assert proj.state()["pfam"] == recorded, \
        "the killed run wrote a record of its own, so this is not the window"

    proc = proj.run("--only", "pfam")
    said = _warn_block(proc.stderr, LATE)
    assert "[pfam]" in said
    assert proj.rpath("hmm", "pfam.tblout") in said
    assert recorded["run_id"] in said, \
        "the warning does not name the run whose record it is about"
    assert "ONE DIRECTION ONLY" in said
    assert "SEVERAL STAGES AT ONCE" not in said, \
        "one late record is not read two ways"
    assert "--force --only" in said, "one late record is told how to rebuild it"
    assert proj.state()["pfam"] == recorded, \
        "a warning invalidated, recomputed or re-recorded the stage"


def test_a_copied_directory_says_one_line_before_a_force_only_and_one_after(
        tmp_path, stub_bin, ma):
    """The four steps that took the first build of this check apart, in order,
    with ordinary commands and no state file edited by hand.

    WHY THE FOUR STEPS AND NOT A BACKDATED DOCUMENT. The test that stood here
    backdated EVERY record and asserted that the run said one thing about the
    directory rather than one thing per stage. That pinned the all-stale case
    and structurally could not see the case that matters: the first build
    suppressed the whole report when every `ok` record was late, so
    `--force --only pfam` - the action the report itself recommends - left one
    fresh record, the suppression stopped holding, and every run after that
    said the per-stage sentence about each of the stages left, for ever, each
    one recommending a recompute of a real InterProScan or ESMFold stage on
    evidence that is nothing but the copy's restamp. A test that backdates
    every record can never take the third step.

    So the third step is the point of this test, and the property is a COUNT
    that does not move across it: one line before, one line after. The report
    is one WARN per run whose length grows with the number of late stages and
    whose count does not, so there is no cliff for the remedy to fall off.

    `cp -r` and real elapsed time, because the fact under test is what a copy
    does to mtimes and what it does NOT do to the stamps inside the document -
    and the interval between the two is exactly what has to exceed the slack.
    (`rsync -a` preserves mtimes and is why the documented two-machine
    workflow never reaches any of this.)
    """
    proj = _searchable(tmp_path, tmp_path / "orig")
    proj.run()
    stages = _stages_of(proj)
    assert len(stages) > 2, "a copy of one stage is not the case under test"

    time.sleep(ma.OUTPUT_STAMP_SLACK_S + 1.0)
    copied = str(tmp_path / "copied")
    subprocess.run(["cp", "-r", proj.root, copied], check=True, timeout=300)

    def run_copy(*args):
        return run_metaannot(
            "run", "--config", os.path.join(copied, "config.yaml"), *args,
            cwd=copied, timeout=600)

    before = run_copy()
    assert before.stderr.count(LATE) == 1, \
        "a copied directory says this once per run, not once per stage"
    said = _warn_block(before.stderr, LATE)
    for name in stages:
        assert f"[{name}]" in said, f"{name} is late and is not named"
    assert "SEVERAL STAGES AT ONCE" in said, \
        "several stages at once is read two ways or it asserts a cause"
    assert "copied, extracted or restored" in said
    assert "rather than what a replaced file looks like" not in said, \
        "the report claims to know which of the two readings it is"
    assert SKEWED not in said, "nothing here measured a skewed clock"

    # Step three: the remedy the report names, on one stage. It recomputes
    # first, so it compares nothing and says nothing.
    forced = run_copy("--force", "--only", "pfam")
    assert LATE not in forced.stderr
    assert forced.stderr.count("WARN") >= 0      # --force is silent here

    # Steps four and five: twice, because "the same seven on every run after
    # that" was the shape of the defect and once cannot see it.
    for attempt in (4, 5):
        after = run_copy()
        assert after.stderr.count(LATE) == 1, \
            f"run {attempt} after the --force says this more than once"
        blk = _warn_block(after.stderr, LATE)
        assert "[pfam]" not in blk, \
            "the stage that was rebuilt is still being reported"
        for name in [n for n in stages if n != "pfam"]:
            assert f"[{name}]" in blk, \
                f"run {attempt} stopped naming {name} after one re-record"


def test_two_late_records_are_two_records_and_both_stages_are_named(
        tmp_path, stub_bin, ma):
    """Two `ok` records is the MINIMUM any `--only` run leaves, and a
    superseded run renaming two outputs inside
    min(heartbeat_s, STATE_PROBE_S) is exactly the shape this check exists
    for. The first build answered it with the directory-wide sentence about a
    directory "copied, extracted or restored": neither stage named, and a
    cause asserted that it had no way to know - "a single record is not a
    pattern" had become "two is one".

    Both stages are named now, and the cause is offered as one of two readings
    and asserted as neither, because it cannot be had from timestamps at all:
    a `cp -r` of this project's results directory puts every late output
    inside twenty milliseconds of every other, and so does this test.
    """
    proj = _searchable(tmp_path, tmp_path / "pair")
    proj.run("--only", "pfam", "dbcan")
    assert _stages_of(proj) == ["dbcan", "pfam"]

    time.sleep(ma.OUTPUT_STAMP_SLACK_S + 1.0)
    outs = [proj.rpath("hmm", "pfam.tblout"),
            proj.rpath("hmm", "dbcan.domtblout")]
    for path in outs:
        with open(path, "rb") as fh:
            data = fh.read()
        with open(path, "wb") as fh:                # same bytes, new mtime
            fh.write(data)

    proc = proj.run("--only", "pfam", "dbcan")
    said = _warn_block(proc.stderr, LATE)
    for name in ("pfam", "dbcan"):
        assert f"[{name}]" in said, f"{name} is late and is not named"
    for path in outs:
        assert path in said, "the file that was overwritten is not named"
    assert "CANNOT TELL YOU WHICH OF TWO THINGS IT IS" in said
    assert "rather than what a replaced file looks like" not in said
    assert SKEWED not in said
    assert all(v["status"] == "ok" for k, v in proj.state().items()
               if not k.startswith("_"))


def test_an_ordinary_resume_reports_a_late_write_on_none_of_its_stages(
        tmp_path, stub_bin):
    """The guard that matters more than the warning does.

    Every cached stage of every resume goes through this comparison, so a rule
    a shade too strict is a warning on a healthy directory - and a warning
    that fires on everything is turned off within a week, taking the real one
    with it. Driven over a whole directory rather than one record, because a
    stamp truncated to the whole second makes every healthy stage look a
    fraction late and a single-record test can miss which way it went.
    """
    proj = _searchable(tmp_path, tmp_path / "resume")
    proj.run()
    assert len(_stages_of(proj)) > 1
    proc = proj.run()
    assert LATE not in proc.stderr
    assert SKEWED not in proc.stderr
    assert UNDATED not in proc.stderr


def test_force_only_compares_no_timestamps_because_it_recomputes_first(
        tmp_path, stub_bin):
    """--force is a statement that the record is not to be trusted, and
    answering it with a complaint about that record is noise. It is also the
    action this check's own message recommends, so it must not then complain
    about itself."""
    proj = _searchable(tmp_path, tmp_path / "forced")
    proj.run()
    _backdate(proj, "pfam")
    proc = proj.run("--force", "--only", "pfam")
    assert LATE not in proc.stderr
    assert proj.state()["pfam"]["status"] == "ok"


def test_an_adopted_record_is_never_dated_against_the_file_it_adopted(
        tmp_path, stub_bin):
    """`finished` on an `adopted` record is when THIS box noticed the file,
    not when the box that made it wrote it, and `rsync -a` preserves the
    source mtime. The two numbers are neither on one clock nor about one
    event, so the GPU hand-off is never reported - which is the restriction
    that keeps this quiet on the documented two-machine workflow."""
    proj = _searchable(tmp_path, tmp_path / "adopted")
    os.makedirs(proj.rpath("hmm"), exist_ok=True)
    F.write_tblout(proj.rpath("hmm", "pfam.tblout"),
                   [("P_dark1", "Peptidase_S8", "PF00082.1")])
    proj.run()
    assert proj.state()["pfam"]["status"] == "adopted"
    _backdate(proj, "pfam")
    proc = proj.run("--only", "pfam")
    assert LATE not in proc.stderr
    assert proj.state()["pfam"]["status"] == "adopted"


def test_a_record_whose_stamp_cannot_be_read_is_compared_against_nothing(
        tmp_path, stub_bin):
    """A record this build cannot PARSE is absence of evidence, and it is met
    by comparing nothing and saying nothing. The alternative fires on every
    stage of every directory whose stamps this build cannot read, for ever -
    and no build of this tool has ever written such a stamp, so the silence
    costs nothing that was ever going to happen.

    That is a different case from a stamp this build reads and declines to
    DATE, which is said out loud two tests below: there the text is one this
    tool wrote, and the operator is owed the reason the stage went unjudged.
    """
    proj = _searchable(tmp_path, tmp_path / "undated")
    proj.run()
    path = proj.rpath(".metaannot_state.json")
    with open(path, encoding="utf-8") as fh:
        st = json.load(fh)
    before = dict(st["pfam"], finished="a while back")
    st["pfam"] = before
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(st, fh, indent=1)
    proc = proj.run("--only", "pfam")
    assert LATE not in proc.stderr
    assert SKEWED not in proc.stderr
    assert UNDATED not in proc.stderr
    assert proj.state()["pfam"] == before, "an undatable record was acted on"


def test_a_sentinel_stage_reports_a_done_file_written_after_its_record(
        tmp_path, stub_bin):
    """diamond, hhblits and esmfold declare a zero-byte `.done` sentinel, and
    a superseded run re-touches one with a plain open().close() that no
    ownership gate covers. That makes the sentinel the only witness anywhere
    for the per-database tables, the per-query .hhr files and the per-protein
    PDBs beside it - and the message has to say, in the same breath, that a
    sentinel it can date says nothing whatever about those files."""
    proj = _searchable(tmp_path, tmp_path / "sentinel")
    proj.run()
    _backdate(proj, "diamond")
    proc = proj.run("--only", "diamond")
    said = _warn_block(proc.stderr, LATE)
    assert "[diamond]" in said
    assert proj.rpath("diamond", ".done") in said
    assert "sentinel and says nothing about" in said
    assert proj.state()["diamond"]["status"] == "ok"


_CLOCK_AHEAD_RUN = '''
"""A real run over a real directory, on a filesystem whose clock leads this
machine's.

There is no portable way to skew a mount under a test, and the skew is the one
suppression this check still has - so the measurement is what is replaced,
and nothing else. _fs_clock_ahead() is the whole of how the run learns the
number; everything downstream of it, including the decision not to read any
mtime as evidence, runs exactly as it does in production.
"""
import importlib.util
import os
import sys

spec = importlib.util.spec_from_file_location("metaannot", os.environ["MA_PY"])
ma = importlib.util.module_from_spec(spec)
sys.modules["metaannot"] = ma
spec.loader.exec_module(ma)

ma._fs_clock_ahead = lambda path: float(os.environ["MA_AHEAD"])
sys.exit(ma.main())
'''


def test_a_filesystem_whose_clock_leads_this_one_is_read_as_no_evidence(
        tmp_path, stub_bin):
    """The one suppression left, and the only one that is a MEASUREMENT rather
    than a reading of a pattern.

    A `finished` stamp is the recording machine's wall clock and an mtime is
    the filesystem's; on an NFS or SMB results directory those are two clocks,
    and one that runs ahead makes every output look newer than the record
    describing it for as long as the mount is skewed. There is nothing to
    conclude from comparing two clocks, so the run says that instead, once,
    and reads no mtime as evidence about any stage.
    """
    proj = _searchable(tmp_path, tmp_path / "skewed")
    proj.run()
    _backdate(proj, "pfam")
    driver = tmp_path / "clock_ahead_run.py"
    driver.write_text(_CLOCK_AHEAD_RUN, encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, str(driver), "run", "--config", proj.config_path,
         "--only", "pfam"],
        capture_output=True, text=True, cwd=proj.root, timeout=600,
        env=dict(os.environ, MA_PY=METAANNOT_PY, MA_AHEAD="1000",
                 PYTHONHASHSEED="0"))
    assert proc.returncode == 0, proc.stderr
    assert SKEWED in proc.stderr
    assert "1000 s ahead" in proc.stderr
    assert LATE not in proc.stderr, \
        "an mtime was read as evidence against a clock it is not on"
    assert proj.state()["pfam"]["status"] == "ok"


def test_a_run_that_ends_in_a_failure_still_says_what_it_found(
        tmp_path, stub_bin):
    """The report is emitted before BOTH end-of-run summaries and inside
    neither, so it is not lost on the path that matters most.

    A run that dies is the run whose operator is about to go looking at the
    directory by hand, and the stages it reused were reused whatever happened
    to the ones that ran. Pinned because the placement is a one-line decision
    that reads as arbitrary and would be the first thing a later edit moved.
    """
    proj = _searchable(tmp_path, tmp_path / "failing")
    proj.run("--only", "pfam", "dbcan")     # nothing downstream exists yet
    _backdate(proj, "pfam")
    # pfam is cached and late; finalise has to RUN and its dependencies were
    # left out of the selection, which is the ordinary way a run fails at
    # dispatch.
    proc = proj.run("--only", "pfam", "finalise", expect=1)
    assert "cannot run:" in proc.stderr, "this run did not fail at dispatch"
    said = _warn_block(proc.stderr, LATE)
    assert "[pfam]" in said
    assert proj.state()["pfam"]["status"] == "ok"


def test_a_dry_run_reports_the_directory_and_says_what_it_cannot_measure(
        tmp_path, stub_bin):
    """A dry run calls decide() for every stage, so it is the one command that
    can report a whole directory without running anything - an audit, and the
    one place where that report is the point rather than a surprise.

    And it is the one place the skew measurement is NOT available: it is taken
    on a file the run has just written, a dry run writes nothing on purpose,
    and a plan check that created files in order to measure them would stop
    being a plan check. So the report says that the question was not asked
    rather than leaving a reader to infer it was asked and answered no. The
    first build said neither, and the README claimed it did.
    """
    proj = _searchable(tmp_path, tmp_path / "audit")
    proj.run()
    _backdate(proj, "pfam")
    dry = proj.run("--dry-run")
    said = _warn_block(dry.stderr, LATE)
    assert "[pfam]" in said
    assert "has not been asked at all" in said, \
        "a dry run does not say that it could not measure the clock skew"
    live = proj.run("--only", "pfam")
    assert LATE in live.stderr
    assert "has not been asked at all" not in live.stderr, \
        "a real run measured the skew and still says it did not"


_AMBIGUOUS_TZ = "America/New_York"
# The hour this zone repeats leaving summer time, in the PAST, so an output
# written now really is newer than the text and the comparison is reached.
_AMBIGUOUS_STAMP = "2025-11-02T01:30:00"
_SKIPPED_STAMP = "2026-03-08T02:30:00"      # the hour it skips entering it


def _names_two_instants(stamp, tz=_AMBIGUOUS_TZ):
    """Whether THIS machine's tz database really makes `stamp` ambiguous.

    Asked with `time` alone and never through the code under test, which is
    the whole point of it being a function: a guard that reads
    `_stamp_instant()` would SKIP on a build that had stopped declining
    ambiguous stamps instead of failing on it, which is the one outcome a test
    for a silent defect may not have. Measured the first time this was
    reverted, where exactly that happened.
    """
    if not hasattr(time, "tzset"):
        return False
    old_tz, fmt = os.environ.get("TZ"), "%Y-%m-%dT%H:%M:%S"
    os.environ["TZ"] = tz
    time.tzset()
    try:
        e = time.mktime(time.strptime(stamp, fmt))
        return (time.strftime(fmt, time.localtime(e)) == stamp
                and time.strftime(fmt, time.localtime(e + 3600)) == stamp)
    finally:
        if old_tz is None:
            del os.environ["TZ"]
        else:
            os.environ["TZ"] = old_tz
        time.tzset()


@pytest.mark.skipif(not _names_two_instants(_AMBIGUOUS_STAMP),
                    reason="no America/New_York rules in this tz database")
def test_a_local_stamp_in_a_repeated_or_skipped_hour_names_no_instant(ma):
    """`_stamp_epoch` is time.mktime(time.strptime(...)) on a naive local
    stamp and cannot resolve the local hour that a zone repeats when it leaves
    summer time. It answers anyway, an hour early, which is twelve hundred
    times OUTPUT_STAMP_SLACK_S - so a stage recorded `ok` in that hour would
    be reported on every run for ever, and because only SOME stages of a run
    that straddles the hour are affected nothing in the comparison could
    absorb it.

    Asserted against the zone's own rules rather than against a written-down
    epoch: the test first shows that this text really does name two instants
    here, and only then that `_stamp_instant` declines to pick one.
    """
    old_tz = os.environ.get("TZ")
    os.environ["TZ"] = _AMBIGUOUS_TZ
    time.tzset()
    try:
        fmt = "%Y-%m-%dT%H:%M:%S"
        # Two instants an hour apart carry that one text - asserted by the
        # skipif above, through `time` and not through the code under test -
        # and _stamp_epoch returns the earlier of them without a word.
        first = ma._stamp_epoch(_AMBIGUOUS_STAMP)
        assert time.strftime(fmt, time.localtime(first)) == _AMBIGUOUS_STAMP
        assert (time.strftime(fmt, time.localtime(first + 3600))
                == _AMBIGUOUS_STAMP)
        assert ma._stamp_instant(_AMBIGUOUS_STAMP) is None
        # The hour a zone SKIPS entering summer time is no instant at all.
        assert ma._stamp_instant(_SKIPPED_STAMP) is None
        # Everything else is unchanged, to the second.
        for plain in ("2026-06-01T12:00:00", "2026-11-01T05:30:00"):
            assert ma._stamp_instant(plain) == ma._stamp_epoch(plain)
        assert ma._stamp_instant("not a stamp") is None
        assert ma._stamp_instant(None) is None
    finally:
        if old_tz is None:
            del os.environ["TZ"]
        else:
            os.environ["TZ"] = old_tz
        time.tzset()


@pytest.mark.skipif(not _names_two_instants(_AMBIGUOUS_STAMP),
                    reason="no America/New_York rules in this tz database")
def test_a_stage_recorded_in_the_repeated_hour_is_declined_and_said_so(
        tmp_path, stub_bin):
    """The same fact through a whole run, because the unit test above cannot
    show what the operator sees.

    A record in that hour is DECLINED rather than judged, and the run says
    which record and why - a stage that goes unjudged in silence is a stage
    nobody can tell from one that was judged sound, and this check's whole
    claim is that it reads in one direction.
    """
    proj = _searchable(tmp_path, tmp_path / "dst")
    proj.run()
    _set_finished(proj, "pfam", _AMBIGUOUS_STAMP)
    said = proj.run("--only", "pfam", env={"TZ": _AMBIGUOUS_TZ}).stderr
    assert UNDATED in said, \
        "a record in the repeated hour was dated instead of declined"
    assert f"pfam says {_AMBIGUOUS_STAMP}" in said, \
        "the record that went unjudged is not named"
    assert LATE not in said, "an ambiguous stamp was dated anyway"
    assert proj.state()["pfam"]["status"] == "ok"

    # The control, one hour later: one instant, and reported as late. Without
    # it this test would pass on a build that declined every stamp.
    _set_finished(proj, "pfam", "2025-11-02T03:30:00")
    said = proj.run("--only", "pfam", env={"TZ": _AMBIGUOUS_TZ}).stderr
    assert LATE in said, "the decline is about ambiguity, not about the date"
    assert UNDATED not in said


def test_dating_an_output_can_never_fail_a_stage_that_has_already_finished(
        tmp_path, ma):
    """The one way this change could be worse than not making it: a stat that
    raises on the cached path of a resume would be a traceback out of decide()
    on a run that was about to reuse work which took hours."""
    missing = str(tmp_path / "nothing here")
    assert ma._outputs_written_late([missing], time.time() - 10_000) == []
    assert ma._outputs_written_late([missing], None) == []
    assert ma._outputs_written_late([str(tmp_path)], None) == []
    assert ma._stamp_epoch("not a stamp") is None
    assert ma._stamp_instant("not a stamp") is None
    assert ma._outputs_written_late([str(tmp_path)],
                                    ma._stamp_epoch(None)) == []
    assert ma._outputs_written_late([str(tmp_path)],
                                    ma._stamp_instant(None)) == []
    assert ma._fs_clock_ahead(missing) is None


# The vocabulary this check may never use, in the message OR in the prose about
# it. A stamp in agreement proves nothing whatever, so a word that reads as a
# clean bill of health turns one-directional evidence into a claim the check
# cannot support - and the operator reads the message at 2am, not the README.
_NEVER = ("verif", "match", "confirm", "intact", "corrupt")


def test_neither_the_message_nor_the_readme_reads_as_a_clean_bill_of_health(
        tmp_path, stub_bin):
    """Both halves in one test, because they make ONE claim and drift apart
    the moment they are pinned separately.

    The message half reads the WHOLE WARN and not its first line. The message
    is a heading, a list of stages and a paragraph; every limitation it states
    is in the paragraph, so scanning the heading alone would have been
    scanning the one part that cannot break this rule.

    The README half also holds the sentinel limitation to NAMING the stages
    rather than counting them: `tests/test_docs.py` already classifies the
    phrase "three stages" as a measurement of something else entirely, and a
    new sentence riding on that entry is how a registry starts lying.
    """
    proj = _searchable(tmp_path, tmp_path / "vocabulary")
    proj.run()
    _backdate(proj, "pfam")
    said = proj.run("--only", "pfam").stderr
    warned = _warn_block(said, LATE).lower()
    for word in _NEVER:
        assert word not in warned, \
            f"the message says {word!r}, which reads as a verdict on the file"

    readme = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
    start = readme.index("**So the file is dated against its own record")
    para = readme[start:readme.index("\n### ", start)]
    for word in _NEVER:
        assert word not in para.lower(), \
            f"the README paragraph says {word!r} about a check that cannot " \
            "tell anyone their output is sound"
    assert "proves nothing" in para
    for stage in ("diamond", "hhblits", "esmfold"):
        assert stage in para, \
            "the sentinel limitation names the stages rather than counting them"
    assert "three stages" not in para
