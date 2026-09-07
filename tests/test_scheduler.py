"""The stage graph: selection, dependencies, the lock, resume and interrupt."""
from __future__ import annotations

import json
import os
import shutil
import signal
import subprocess
import sys
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
