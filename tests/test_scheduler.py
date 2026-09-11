"""The stage graph: selection, dependencies, the lock, resume and interrupt."""
from __future__ import annotations

import argparse
import io
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time

import pytest

import fixtures as F
from conftest import METAANNOT_PY, build_project, run_metaannot
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


def test_force_unlock_takes_over(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    with open(proj.rpath(".metaannot.lock"), "w", encoding="utf-8") as fh:
        json.dump({"pid": os.getpid(),
                   "host": __import__("socket").gethostname(),
                   "started": "2026-01-01T00:00:00"}, fh)
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


def test_a_run_whose_lock_vanishes_still_records_how_it_ended(ma, tmp_path):
    # symptom: the ownership gate answered False for a missing lock, and
    # RunRecord refused to write on anything but True - so an operator's `rm`,
    # a tmp-reaper or a remount left a live, unsuperseded run stranded at
    # final_status "running" with finished null, and the log claimed the
    # directory had been handed to another run.
    state = tmp_path / ".metaannot_state.json"
    lock = ma.ResultsLock(str(tmp_path / "c.lock"))
    lock.__enter__()
    rec = ma.RunRecord(str(state), ma._State(), ["metaannot", "run"],
                       config_path=None, owner=lock)
    rec.stamp("running")
    os.remove(lock.path)                       # vanished, with no replacement
    assert rec.stamp("ok") is not False, "the run lost its own final verdict"
    got = json.loads(state.read_text())[ma.RUN_KEY]
    assert got["final_status"] == "ok"
    assert got["finished"], "finished was never stamped"


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
