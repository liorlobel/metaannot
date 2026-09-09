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
    assert _wait_for(lambda: os.path.exists(state)
                     and b"running" in open(state, "rb").read()), \
        "the stage was never recorded as running"
    proc.send_signal(signal.SIGINT)
    out, err = proc.communicate(timeout=60)

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
