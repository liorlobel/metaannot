"""The console itself, driven directly.

`console/console.py` is a plain module: the readers are pure functions, the
model turns a state file into rows, and the server is one class. So these drive
it in process, and reach for a real UNIX socket only where the socket is the
thing under test.

The ugly cases are the point. A watcher pointed at a three-day job meets a state
file mid-rewrite, a log that rotated, a heartbeat four hours old and a results
directory that has nothing in it yet, and none of those may produce a traceback,
a blank page, or a claim that a run is dead.
"""
import copy
import http.client
import json
import os
import pathlib
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time

import pytest

from conftest import METAANNOT_PY


# ----------------------------------------------------------------------
# scaffolding
# ----------------------------------------------------------------------

@pytest.fixture(scope="session")
def contract(console):
    """One describe --json for the whole file: it is a fork, not a fixture to
    take lightly, and the console itself only runs it once."""
    return console.Contract(console.run_describe(sys.executable, METAANNOT_PY))


@pytest.fixture
def results(tmp_path, contract):
    """An empty results directory, and helpers to put things in it."""
    root = tmp_path / "results"
    root.mkdir()

    class R:
        path = str(root)

        def state(self, obj):
            self.write(contract.state_name, json.dumps(obj))

        def raw_state(self, text):
            self.write(contract.state_name, text)

        def log(self, text):
            return self.write(contract.log_name, text)

        def lock(self, obj):
            self.write(contract.lock_name, json.dumps(obj))

        def write(self, name, text):
            path = os.path.join(self.path, name)
            mode = "wb" if isinstance(text, bytes) else "w"
            with open(path, mode) as fh:
                fh.write(text)
            return path

    return R()


def stamp(epoch):
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(epoch))


def running_state(contract, age=0, heartbeat=30, started_ago=7200, extra=None):
    now = time.time()
    state = {contract.run_key: {
        "run_id": "20260101T000000-1", "version": "0.3.0",
        "config_path": "/nowhere/config.yaml", "argv": ["metaannot.py"],
        "host": "lab-fedora", "pid": 4242, "started": stamp(now - started_ago),
        "last_seen": stamp(now - age), "last_seen_epoch": now - age,
        "heartbeat_s": heartbeat, "finished": None, "final_status": "running"}}
    state.update(extra or {})
    return state


def view_of(console, results, contract, **kw):
    project = console.Project(0, results.path, contract)
    return console.project_view(project, sys.executable, METAANNOT_PY, **kw)


@pytest.fixture
def short_sock():
    """AF_UNIX allows 104 bytes of path on macOS; pytest's tmp_path spends
    most of that before the filename."""
    d = tempfile.mkdtemp(prefix="mac", dir="/tmp")
    try:
        yield os.path.join(d, "s.sock")
    finally:
        shutil.rmtree(d, ignore_errors=True)


class UnixConn(http.client.HTTPConnection):
    def __init__(self, path):
        http.client.HTTPConnection.__init__(self, "localhost")
        self._path = path

    def connect(self):
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(30)
        self.sock.connect(self._path)


# ----------------------------------------------------------------------
# the state file, in every state a running engine can leave it
# ----------------------------------------------------------------------

def test_a_missing_state_file_renders_a_page(console, results, contract):
    v = view_of(console, results, contract)
    assert v["state_status"] == "missing"
    assert len(v["rows"]) == 21
    html = console.render_project(v, contract)
    # The accurate sentence LEADS the pinned block now; it used to be last, in
    # muted small text, under "There is no run record in the state file.
    # Either nothing has run here, or it was written by a metaannot older than
    # the one that added the record" — about a file that is not there, and
    # with a second branch that cannot happen.
    assert "There is no state file in this directory" in html
    assert "written by a metaannot older" not in html
    assert html.index("There is no state file in this directory") < \
        html.index("<thead>")
    assert "Traceback" not in html


def test_an_empty_state_file_renders_a_page(console, results, contract):
    results.raw_state("")
    v = view_of(console, results, contract)
    assert v["state_status"] == "empty"
    assert console.render_project(v, contract)


def test_a_truncated_state_file_renders_a_page(console, results, contract):
    results.raw_state('{"emapper": {"status": "o')
    v = view_of(console, results, contract)
    assert v["state_status"] == "bad"
    assert "not valid JSON" in v["state_detail"]
    assert console.render_project(v, contract)


def test_a_state_file_that_is_not_an_object_is_refused_not_rendered(
        console, results, contract):
    results.raw_state("[1, 2, 3]")
    v = view_of(console, results, contract)
    assert v["state_status"] == "bad"
    assert "top level is list" in v["state_detail"]


def test_an_unreadable_state_file_falls_back_to_the_last_good_read(
        console, results, contract):
    """The table is never blanked over a transient read, and the page says
    which moment it is showing."""
    results.state(running_state(contract, extra={
        "emapper": {"signature": "a", "status": "ok", "seconds": 5.0,
                    "finished": stamp(time.time())}}))
    project = console.Project(0, results.path, contract)
    first = console.project_view(project, sys.executable, METAANNOT_PY)
    assert first["state_cached"] is False
    assert [r["state"] for r in first["rows"]][0] == "ok"

    results.raw_state("{ this is not json")
    second = console.project_view(project, sys.executable, METAANNOT_PY)
    assert second["state_cached"] is True
    assert [r["state"] for r in second["rows"]][0] == "ok"    # not blanked
    assert second["state_good_at"] == pytest.approx(first["taken"], abs=5)
    html = console.render_project(second, contract)
    assert "showing the last good read" in html


def test_a_state_file_replaced_under_us_is_read_exactly_once_per_view(
        console, results, contract, monkeypatch):
    """os.replace() swaps a whole new file in. Reading it twice inside one
    request is how a header and the table under it end up describing two
    different moments, so the view takes one snapshot and derives everything
    from that."""
    results.state(running_state(contract))
    calls = []
    real = console.read_json_stable

    def counting(path, limit, *a, **kw):
        if path.endswith(contract.state_name):
            calls.append(path)
        return real(path, limit, *a, **kw)

    monkeypatch.setattr(console, "read_json_stable", counting)
    view_of(console, results, contract)
    assert calls == [os.path.join(results.path, contract.state_name)]


def test_a_state_file_that_appears_between_attempts_is_picked_up(
        console, results, contract, monkeypatch):
    """The retry covers a file that momentarily reads as gone - but only for a
    project whose state file has been read before. A directory where nothing
    has ever run must not cost 100 ms of sleep on every poll."""
    seen = {"n": 0}
    real = console.read_json_once
    target = os.path.join(results.path, contract.state_name)

    def flaky(path, limit):
        if path == target:
            seen["n"] += 1
            if seen["n"] < 3:
                return "missing", None, None
        return real(path, limit)

    results.state(running_state(contract))
    monkeypatch.setattr(console, "read_json_once", flaky)
    status, obj, _ = console.read_json_stable(target, console.STATE_MAX,
                                              retry_missing=True)
    assert status == "ok" and seen["n"] == 3
    assert contract.run_key in obj

    seen["n"] = 0
    status, _, _ = console.read_json_stable(target, console.STATE_MAX)
    assert status == "missing" and seen["n"] == 1     # not retried


def test_a_state_record_with_an_unknown_status_becomes_a_question_mark(
        console, results, contract):
    results.state({"emapper": {"status": "levitating"}})
    v = view_of(console, results, contract)
    row = next(r for r in v["rows"] if r["name"] == "emapper")
    assert row["state"] == "bad" and row["label"] == "?"
    assert "unreadable" in row["detail"]


def test_a_state_file_over_the_cap_is_refused_rather_than_read_whole(
        console, results, contract):
    results.raw_state(" " * (console.STATE_MAX + 10))
    v = view_of(console, results, contract)
    assert v["state_status"] == "bad"
    assert "larger than" in v["state_detail"]


# ----------------------------------------------------------------------
# the log tail, by byte offset
# ----------------------------------------------------------------------

def test_a_missing_log_says_so_and_does_not_raise(console, results, contract):
    tail = console.tail_log(os.path.join(results.path, contract.log_name))
    assert tail["ok"] is False and tail["note"] == "no log file"
    assert tail["text"] == ""


def test_a_zero_byte_log_reads_as_empty(console, results, contract):
    path = results.log("")
    tail = console.tail_log(path)
    assert tail["ok"] and tail["size"] == 0 and tail["text"] == ""
    assert console.parse_log(tail["text"]) == []


def test_a_log_shorter_than_the_window_is_read_whole(console, results):
    path = results.log("[    0.0s] INFO  one\n[    1.0s] WARN  two\n")
    tail = console.tail_log(path)
    assert tail["text"].startswith("[    0.0s] INFO  one")
    assert tail["off"] == tail["size"]
    assert [l["level"] for l in console.parse_log(tail["text"])] == ["INFO",
                                                                    "WARN"]


def test_a_log_longer_than_the_window_is_never_read_whole(console, results):
    """A three-day log is hundreds of megabytes. The read is bounded by the
    window and starts on a line boundary, so no half-line is ever shown."""
    body = "".join("[%7.1fs] INFO  line %d\n" % (i, i) for i in range(60000))
    path = results.log(body)
    assert os.path.getsize(path) > console.LOG_WINDOW
    tail = console.tail_log(path)
    assert len(tail["text"]) <= console.LOG_WINDOW
    assert tail["text"].startswith("[")            # no partial first line
    assert tail["text"].endswith("line 59999\n")
    assert tail["off"] == tail["size"]


def test_a_log_that_grew_between_two_reads_delivers_only_the_delta(
        console, results):
    path = results.log("[    0.0s] INFO  first\n")
    first = console.tail_log(path)
    with open(path, "a") as fh:
        fh.write("[    1.0s] INFO  second\n")
    second = console.tail_log(path, first["off"], first["ident"])
    assert second["reset"] is False and second["note"] is None
    assert second["text"] == "[    1.0s] INFO  second\n"
    assert second["off"] == second["size"]


def test_a_truncated_log_resets_and_says_so(console, results):
    path = results.log("[    0.0s] INFO  " + "x" * 5000 + "\n")
    first = console.tail_log(path)
    results.log("[    0.0s] INFO  fresh\n")           # same name, shorter
    second = console.tail_log(path, first["off"], first["ident"])
    assert second["reset"] is True
    assert "shorter" in second["note"] or "replaced" in second["note"]
    assert "fresh" in second["text"]


def test_a_rotated_log_resets_and_says_so(console, results, contract):
    path = results.log("[    0.0s] INFO  old\n")
    first = console.tail_log(path)
    os.rename(path, path + ".1")                      # a new inode takes over
    results.log("[    0.0s] INFO  brand new\n")
    second = console.tail_log(path, first["off"], first["ident"])
    assert second["reset"] is True
    assert "replaced" in second["note"]
    assert second["ident"] != first["ident"]


def test_a_log_that_outran_the_window_says_the_middle_is_missing(
        console, results):
    path = results.log("[    0.0s] INFO  start\n")
    first = console.tail_log(path)
    with open(path, "a") as fh:
        fh.write("[    1.0s] INFO  " + "y" * (console.LOG_WINDOW + 1000) + "\n")
    second = console.tail_log(path, first["off"], first["ident"])
    assert second["reset"] is True
    assert "not shown" in second["note"]
    assert len(second["text"]) <= console.LOG_WINDOW


def test_a_latin1_byte_in_the_log_costs_one_character_not_the_page(
        console, results, contract):
    """opener()'s docstring records a 0xa0 killing `integrate` after three
    hours of InterProScan. A console that raised UnicodeDecodeError on that
    very log would be a bitter joke."""
    results.write(contract.log_name,
                  b"[    0.0s] WARN  caf\xa0 description\n")
    tail = console.tail_log(os.path.join(results.path, contract.log_name))
    assert "caf" in tail["text"]
    assert console.parse_log(tail["text"])[0]["level"] == "WARN"


def test_the_log_pane_counts_only_what_is_in_view(console, results, contract):
    body = "".join("[%7.1fs] INFO  line %d\n" % (i, i) for i in range(400))
    results.log(body + "[  999.0s] FATAL stage 'cluster' failed: boom\n")
    v = view_of(console, results, contract)
    assert len(v["log"]["lines"]) == console.LOG_LINES
    assert v["log"]["fatal"] == 1
    html = console.render_project(v, contract)
    assert "not for the whole" in html      # the honest footnote


# ----------------------------------------------------------------------
# the heartbeat, which never returns a verdict
# ----------------------------------------------------------------------

# Verdicts the console must never reach. Its own sentence "that is not
# evidence the run is dead" is the opposite of a verdict, so the list is of
# claims rather than of the word.
FORBIDDEN = ("has died", "the run died", "presumed dead", "reclaim",
             "the process is gone", "no longer running", "appears to be dead",
             "probably dead")


def test_a_stale_heartbeat_never_says_the_run_is_dead(console, results,
                                                      contract):
    results.state(running_state(contract, age=4 * 3600 + 720))
    results.log("[    0.0s] INFO  working\n")
    v = view_of(console, results, contract)
    assert v["heartbeat"]["band"] == "long"
    text = console.render_project(v, contract).lower()
    for phrase in FORBIDDEN:
        assert phrase not in text, phrase
    assert "not evidence the run is dead" in text
    assert "unprovable means alive" in text


def test_the_staleness_bands_come_from_the_records_own_heartbeat_s(
        console, results, contract):
    for age, band in ((30, "fresh"), (200, "late"), (5000, "long")):
        results.state(running_state(contract, age=age))
        v = view_of(console, results, contract, log=False)
        assert v["heartbeat"]["band"] == band, (age, band)


def test_a_finished_run_shows_no_staleness_alarm_however_old(
        console, results, contract):
    """The easy bug in this design: a run that finished on Tuesday has a log
    that last moved on Tuesday, and that is not late."""
    now = time.time()
    state = running_state(contract)
    state[contract.run_key].update(final_status="ok", finished=stamp(now - 4 * 86400),
                                   last_seen_epoch=now - 4 * 86400)
    results.state(state)
    results.log("[    0.0s] INFO  done\n")
    os.utime(os.path.join(results.path, contract.log_name),
             (now - 4 * 86400, now - 4 * 86400))
    v = view_of(console, results, contract)
    assert v["heartbeat"]["band"] == "none"
    assert v["heartbeat"]["lines"] == []
    assert "Finished ok" in v["heartbeat"]["verdict"]
    assert "late" not in console.render_project(v, contract)


def test_a_heartbeat_from_the_future_is_named_as_clock_skew(
        console, results, contract):
    results.state(running_state(contract, age=-3600))
    v = view_of(console, results, contract, log=False)
    assert v["heartbeat"]["band"] == "skew"
    assert "future" in v["heartbeat"]["verdict"]
    assert "ago" not in v["heartbeat"]["verdict"]   # no negative age


def test_a_run_record_with_no_epoch_says_so_rather_than_guessing(
        console, results, contract):
    state = running_state(contract)
    del state[contract.run_key]["last_seen_epoch"]
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert v["heartbeat"]["band"] == "unknown"
    assert "no machine-readable heartbeat" in v["heartbeat"]["verdict"]


def test_a_state_file_with_no_run_record_still_renders(console, results,
                                                       contract):
    """Not hypothetical: the server README records a run reporting a version
    that matches no release."""
    results.state({"emapper": {"signature": "a", "status": "ok",
                               "seconds": 1.0, "finished": stamp(time.time())}})
    v = view_of(console, results, contract)
    assert v["run"] is None
    assert "no run record" in v["heartbeat"]["verdict"]
    assert "carries no run record" in console.render_project(v, contract)


def test_a_fatal_line_under_a_running_record_is_called_out(console, results,
                                                           contract):
    results.state(running_state(contract, age=9000))
    results.log("[    1.0s] FATAL stage 'cluster' failed: boom\n")
    v = view_of(console, results, contract)
    assert any("killed-mid-write" in line for line in v["heartbeat"]["lines"])


# ----------------------------------------------------------------------
# the lock, and the decision the console refuses to take
# ----------------------------------------------------------------------

def test_no_lock_file_is_stated_plainly(console, results, contract):
    v = view_of(console, results, contract, log=False)
    assert v["lock"]["held"] is False
    assert "Nothing holds this directory" in v["lock"]["text"]


def test_a_lock_hands_over_the_ps_line_and_refuses_the_decision(
        console, results, contract):
    results.lock({"pid": 21877, "host": "lab-fedora",
                  "started": "2026-09-07T22:10:41"})
    results.state(running_state(contract, age=20000))
    v = view_of(console, results, contract)
    assert v["lock"]["ps"] == "ps -p 21877 -o pid,etime,stat,args"
    html = console.render_project(v, contract)
    assert "ps -p 21877 -o pid,etime,stat,args" in html
    assert "lab-fedora" in html            # which box to run it on
    assert "offers no --force-unlock" in html


def test_a_garbled_lock_file_does_not_break_the_page(console, results,
                                                     contract):
    results.write(contract.lock_name, "{not json")
    v = view_of(console, results, contract, log=False)
    assert v["lock"]["held"] is None
    assert "unreadable" in v["lock"]["text"]


# ----------------------------------------------------------------------
# the stage table
# ----------------------------------------------------------------------

def test_all_twenty_one_stages_are_always_shown(console, results, contract):
    results.state(running_state(contract))
    v = view_of(console, results, contract, log=False)
    assert [r["name"] for r in v["rows"]] == contract.stage_names


def test_an_ok_record_shows_seconds_and_invents_no_start_time(
        console, results, contract):
    now = time.time()
    results.state(running_state(contract, extra={
        "emapper": {"signature": "a", "status": "ok", "seconds": 13260.0,
                    "finished": stamp(now)}}))
    row = next(r for r in view_of(console, results, contract, log=False)["rows"]
               if r["name"] == "emapper")
    assert row["took"] == "3h41m"          # the engine's own shape
    assert "started" not in row["detail"]


def test_an_adopted_record_explains_the_dash_rather_than_leaving_a_gap(
        console, results, contract):
    results.state(running_state(contract, extra={
        "pfam": {"signature": "a", "status": "adopted",
                 "finished": stamp(time.time())}}))
    row = next(r for r in view_of(console, results, contract, log=False)["rows"]
               if r["name"] == "pfam")
    assert row["took"] == "—"
    assert "reuse records no duration" in row["detail"]


def test_a_record_from_an_earlier_run_is_labelled_as_one(console, results,
                                                         contract):
    """The state file is cumulative, so a stage disabled this time keeps last
    week's green record. Showing that as this run's success is a lie by
    omission."""
    now = time.time()
    results.state(running_state(contract, started_ago=3600, extra={
        "smorf": {"signature": "a", "status": "ok", "seconds": 44.0,
                  "finished": stamp(now - 5 * 86400)}}))
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "smorf")
    assert row["carried"] is not None
    assert "from an earlier run" in console.render_project(v, contract)


def test_a_disabled_dependency_does_not_block_a_stage(console, results,
                                                      contract, monkeypatch):
    """The engine's own unmet_deps() says a disabled dependency is fine — its
    evidence is legitimately absent — and finish() puts a skipped stage in
    `done`. A console that reported it as blocking would be wrong about half
    the table on a default config."""
    off = {name: False for name in
           {st["enabled"] for st in contract.stages if st["enabled"]}}
    off["eggnog"] = True
    monkeypatch.setattr(console.Project, "config_run",
                        lambda self, py, sc, st: (off, "test"))
    results.state(running_state(contract, extra={
        "emapper": {"signature": "a", "status": "ok", "seconds": 1.0,
                    "finished": stamp(time.time())}}))
    v = view_of(console, results, contract, log=False)
    rows = {r["name"]: r for r in v["rows"]}
    assert rows["pfam"]["state"] == "off"
    assert rows["integrate"]["state"] == "next"     # nothing left to wait for
    assert "do not block it" not in rows["integrate"]["detail"]
    # and the stage that did run is not called off, whatever the config says
    assert rows["emapper"]["state"] == "ok"


def test_an_unknown_config_never_claims_a_stage_is_off(console, results,
                                                       contract):
    """With no effective config and no readable config beside it, the console
    reports waiting or not-yet-reached, never OFF."""
    results.state(running_state(contract))
    v = view_of(console, results, contract, log=False)
    assert v["cfg_run_known"] is False
    assert not any(r["state"] == "off" for r in v["rows"])
    assert "cannot tell" in console.render_project(v, contract)


def test_a_running_stage_shows_the_part_file_it_is_writing(console, results,
                                                           contract):
    """The engine writes work in progress as `.<stem>.<pid>.<tid>.part<ext>` in
    the output's own directory, and that file growing is a more direct sign of
    work than any heartbeat."""
    stage = next(s for s in contract.stages if s["outputs"])
    out = stage["outputs"][0]
    os.makedirs(os.path.join(results.path, os.path.dirname(out)),
                exist_ok=True)
    stem = os.path.splitext(os.path.basename(out))[0]
    part = os.path.join(results.path, os.path.dirname(out),
                        ".%s.4242.140.part.tsv" % stem)
    with open(part, "wb") as fh:
        fh.write(b"z" * 4096)
    results.state(running_state(contract, extra={
        stage["name"]: {"signature": None, "status": "running",
                        "started": stamp(time.time() - 3600)}}))
    v = view_of(console, results, contract)
    row = next(r for r in v["rows"] if r["name"] == stage["name"])
    assert row["part"] is not None and row["part"]["size"] == 4096
    assert row["took"] == "1h00m"
    html = console.render_project(v, contract)
    assert "more direct sign of work than any heartbeat" in html


def test_a_failed_stage_carries_its_error_into_the_row(console, results,
                                                       contract):
    results.state(running_state(contract, extra={
        "cluster": {"signature": None, "status": "failed",
                    "error": "mmseqs exited 1", "finished": stamp(time.time())}}))
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "cluster")
    assert row["state"] == "failed" and "mmseqs exited 1" in row["detail"]
    assert v["bucket"] == "running"        # a live run stays in the live group
    assert "mmseqs exited 1" in console.render_project(v, contract)


# ----------------------------------------------------------------------
# more than one project, which is the case the doc names
# ----------------------------------------------------------------------

def test_the_index_sorts_and_summarises_several_projects(console, tmp_path,
                                                         contract):
    now = time.time()
    def make(name, state, log_age):
        d = tmp_path / name / "results"
        d.mkdir(parents=True)
        with open(str(d / contract.state_name), "w") as fh:
            json.dump(state, fh)
        with open(str(d / contract.log_name), "w") as fh:
            fh.write("[    0.0s] INFO  x\n")
        os.utime(str(d / contract.log_name), (now - log_age, now - log_age))
        return str(d)

    live = running_state(contract, age=10)
    done = running_state(contract)
    done[contract.run_key].update(final_status="ok", finished=stamp(now - 100))
    dead = running_state(contract)
    dead[contract.run_key].update(final_status="failed",
                                  finished=stamp(now - 50))
    paths = [make("finished", done, 3600), make("live", live, 5),
             make("broken", dead, 60)]
    empty = tmp_path / "nothing" / "results"
    empty.mkdir(parents=True)
    with open(str(empty / contract.log_name), "w") as fh:
        fh.write("")
    paths.append(str(empty))

    projects = [console.Project(i, p, contract) for i, p in enumerate(paths)]
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    order = [i["bucket"] for i in view["projects"]]
    assert order == ["running", "failed", "done", "none"]
    assert "1 running" in view["summary"] and "1 failed" in view["summary"]
    assert "1 never started" in view["summary"]
    html = console.render_index(view, contract)
    for name in ("live", "broken", "finished", "nothing"):
        assert name in html
    # the parent directory names the project, because eight rows reading
    # "results" is the failure
    assert ">results<" not in html


def test_discovery_finds_results_directories_without_descending_into_them(
        console, tmp_path, contract):
    for name in ("ds01", "ds02"):
        d = tmp_path / name / "results"
        (d / "eggnog").mkdir(parents=True)
        with open(str(d / contract.state_name), "w") as fh:
            fh.write("{}")
        # a decoy that would be found if the walk descended into a results dir
        with open(str(d / "eggnog" / contract.log_name), "w") as fh:
            fh.write("")
    (tmp_path / ".hidden" / "results").mkdir(parents=True)
    with open(str(tmp_path / ".hidden" / "results" / contract.log_name),
              "w") as fh:
        fh.write("")
    found = console.discover([str(tmp_path)], contract)
    assert sorted(found) == [str(tmp_path / "ds01" / "results"),
                             str(tmp_path / "ds02" / "results")]


# ----------------------------------------------------------------------
# the socket: 0700, one console, and nothing deleted by accident
# ----------------------------------------------------------------------

def bare_console(console, contract):
    return console.Console(contract, [], sys.executable, METAANNOT_PY)


def test_the_socket_is_0700_even_under_a_wide_umask(console, contract,
                                                    short_sock):
    """bind() creates the inode with 0777 & ~umask, so under umask 0 a socket
    that was not created under 0o077 would come out 0777. The mode is never
    briefly wider, because nothing widens it and narrows it back — there is no
    chmod anywhere in the console (see test_console_contract)."""
    old = os.umask(0)
    try:
        fd = console.take_socket_name(short_sock)
        srv = console.ConsoleServer(short_sock, console.Handler,
                                    bare_console(console, contract))
        try:
            assert os.stat(short_sock).st_mode & 0o7777 == 0o700
            assert os.umask(0) == 0        # and it put the umask back
        finally:
            srv.server_close()
            os.close(fd)
    finally:
        os.umask(old)


def test_a_second_console_on_the_same_socket_is_refused(console, short_sock):
    fd = console.take_socket_name(short_sock)
    try:
        with pytest.raises(console.Refuse) as e:
            console.take_socket_name(short_sock)
        assert "already holds" in str(e.value)
    finally:
        os.close(fd)


def test_a_stale_socket_is_probed_and_removed(console, contract, short_sock):
    """server_close() leaves the inode behind, and so does a SIGKILL. The
    recovery is a connect probe, not a blind unlink."""
    fd = console.take_socket_name(short_sock)
    srv = console.ConsoleServer(short_sock, console.Handler,
                                bare_console(console, contract))
    srv.server_close()
    os.close(fd)
    assert os.path.exists(short_sock)      # still there: not self-healing
    fd = console.take_socket_name(short_sock)   # probe says nothing listens
    try:
        assert not os.path.exists(short_sock)
    finally:
        os.close(fd)


def test_a_live_socket_is_never_unlinked(console, contract, short_sock):
    fd = console.take_socket_name(short_sock)
    srv = console.ConsoleServer(short_sock, console.Handler,
                                bare_console(console, contract))
    thread = threading.Thread(target=srv.serve_forever,
                              kwargs={"poll_interval": 0.05}, daemon=True)
    thread.start()
    try:
        os.close(fd)                       # the flock goes; the listener stays
        with pytest.raises(console.Refuse) as e:
            console.take_socket_name(short_sock)
        assert "already serving" in str(e.value)
        assert os.path.exists(short_sock)
    finally:
        srv.shutdown()
        srv.server_close()


def test_a_path_that_is_not_a_socket_is_never_deleted(console, short_sock):
    """`--socket ~/thesis.docx` must not delete the thesis."""
    import pathlib
    thesis = pathlib.Path(os.path.dirname(short_sock)) / "thesis.docx"
    thesis.write_text("chapter one")
    with pytest.raises(console.Refuse) as e:
        console.take_socket_name(str(thesis))
    assert "refusing to remove it" in str(e.value)
    assert thesis.read_text() == "chapter one"


def test_an_over_long_socket_path_is_named_not_left_to_bind(console, tmp_path):
    long = str(tmp_path / ("x" * 120))
    with pytest.raises(console.Refuse) as e:
        console.take_socket_name(long)
    assert "too long for AF_UNIX" in str(e.value)
    assert "--socket" in str(e.value)


def test_runtime_dir_skips_a_directory_that_is_not_private(console, tmp_path,
                                                           capsys):
    wide = tmp_path / "wide"
    wide.mkdir(mode=0o755)
    mine = tmp_path / "mine"
    found = console.runtime_dir(candidates=[str(wide), str(mine)])
    assert found == str(mine)
    assert os.stat(found).st_mode & 0o077 == 0
    assert "not a private directory" in capsys.readouterr().err


# ----------------------------------------------------------------------
# the server: GET and HEAD, and nothing that acts
# ----------------------------------------------------------------------

@pytest.fixture
def served(console, contract, results, short_sock):
    """A running console over a real UNIX socket."""
    results.state(running_state(contract, age=5))
    results.log("[    0.0s] INFO  hello\n")
    con = console.Console(contract,
                          [console.Project(0, results.path, contract)],
                          sys.executable, METAANNOT_PY)
    stop, box = threading.Event(), {}
    thread = threading.Thread(
        target=console.serve, args=(con, short_sock),
        kwargs={"stop": stop, "ready": lambda s: box.setdefault("srv", s)},
        daemon=True)
    thread.start()
    for _ in range(200):
        if "srv" in box:
            break
        time.sleep(0.02)
    assert "srv" in box
    try:
        yield short_sock
    finally:
        stop.set()
        thread.join(10)


def get(sock, route, method="GET"):
    conn = UnixConn(sock)
    try:
        conn.request(method, route)
        resp = conn.getresponse()
        return resp.status, dict(resp.getheaders()), resp.read()
    finally:
        conn.close()


def test_the_page_comes_back_over_the_socket(served):
    status, headers, body = get(served, "/")
    assert status == 200
    assert b"metaannot console" in body and b"read-only" in body
    assert int(headers["Content-Length"]) == len(body)
    assert headers["X-Content-Type-Options"] == "nosniff"
    assert "default-src 'self'" in headers["Content-Security-Policy"]
    assert not any(h.lower().startswith("access-control") for h in headers)


def test_keep_alive_serves_two_requests_on_one_connection(served):
    conn = UnixConn(served)
    try:
        for route in ("/", "/p/0"):
            conn.request("GET", route)
            resp = conn.getresponse()
            body = resp.read()
            assert resp.status == 200 and len(body) == int(
                resp.getheader("Content-Length"))
    finally:
        conn.close()


def test_head_returns_the_headers_and_no_body(served):
    status, headers, body = get(served, "/", method="HEAD")
    assert status == 200 and body == b""
    assert int(headers["Content-Length"]) > 0


def test_there_is_no_endpoint_that_writes(served):
    """No do_POST exists, so the stdlib answers 501 on its own — which is also
    why CSRF has nothing to attack."""
    for method in ("POST", "PUT", "DELETE", "PATCH"):
        status, _, _ = get(served, "/", method=method)
        assert status == 501, method


def test_no_endpoint_accepts_a_filesystem_path(served):
    """Projects are numbered. A path parameter would make the console a
    whole-filesystem read oracle for anyone who reaches the socket."""
    for route in ("/p/../../etc/passwd", "/api/project?p=/etc/passwd",
                  "/api/log?p=/etc/passwd", "/p/999", "/api/project?p=-1"):
        status, _, body = get(served, route)
        assert status == 404, route
        assert b"passwd" not in body


def test_the_json_endpoints_answer_with_json(served, contract):
    status, headers, body = get(served, "/api/contract")
    assert status == 200 and headers["Content-Type"].startswith(
        "application/json")
    doc = json.loads(body)
    assert doc["state_name"] == contract.state_name
    assert doc["stage_names"] == contract.stage_names

    status, _, body = get(served, "/api/projects")
    assert json.loads(body)["projects"][0]["index"] == 0

    status, _, body = get(served, "/api/log?p=0")
    tail = json.loads(body)
    assert tail["ok"] and tail["lines"][0]["text"].endswith("hello")


def test_the_log_endpoint_appends_from_the_offset_it_was_given(served,
                                                               results):
    _, _, body = get(served, "/api/log?p=0")
    first = json.loads(body)
    with open(os.path.join(results.path, "metaannot.log"), "a") as fh:
        fh.write("[    1.0s] WARN  and then this\n")
    _, _, body = get(served, "/api/log?p=0&off=%d&id=%s"
                     % (first["off"], first["ident"]))
    second = json.loads(body)
    assert second["reset"] is False
    assert [l["text"] for l in second["lines"]] == [
        "[    1.0s] WARN  and then this"]


def test_a_broken_route_returns_a_page_not_a_blank_tab(served, monkeypatch,
                                                       console):
    """On AF_UNIX `client_address` is '', and the inherited address_string()
    would raise IndexError from inside send_response — after the socket closed,
    delivering the browser zero bytes. Anything that raises must still produce
    a status line."""
    monkeypatch.setattr(console, "render_index",
                        lambda *a, **k: 1 / 0)
    status, _, body = get(served, "/")
    assert status == 500
    assert b"failed to render" in body
    assert b"Nothing was written" in body


def test_the_connection_cap_refuses_rather_than_thread_bombs(served, console):
    """This box is running twenty-one stages under a CPU and RAM budget."""
    open_conns = []
    try:
        for _ in range(console.MAX_CONNS):
            conn = UnixConn(served)
            conn.request("GET", "/api/contract")
            conn.getresponse().read()
            open_conns.append(conn)
        extra = UnixConn(served)
        try:
            extra.request("GET", "/api/contract")
            assert extra.getresponse().status == 503
        except (BrokenPipeError, ConnectionResetError,
                http.client.RemoteDisconnected):
            # Also a refusal, and the faster one: the canned 503 goes out and
            # the connection closes without a byte of the request being parsed,
            # which can beat the client's own write.
            pass
        finally:
            extra.close()
    finally:
        for conn in open_conns:
            conn.close()


def test_the_canned_503_carries_its_own_correct_content_length(console):
    """It is written as raw bytes with no framework behind it, and a wrong
    Content-Length there is a browser tab that hangs rather than an error."""
    head, _, body = console.BUSY_RESPONSE.partition(b"\r\n\r\n")
    length = int(dict(
        line.split(b": ", 1) for line in head.split(b"\r\n")[1:]
    )[b"Content-Length"])
    assert length == len(body)
    assert b"Connection: close" in head


def test_the_socket_is_removed_on_shutdown(console, contract, short_sock):
    con = console.Console(contract, [], sys.executable, METAANNOT_PY)
    stop, box = threading.Event(), {}
    thread = threading.Thread(
        target=console.serve, args=(con, short_sock),
        kwargs={"stop": stop, "ready": lambda s: box.setdefault("srv", s)},
        daemon=True)
    thread.start()
    for _ in range(200):
        if "srv" in box:
            break
        time.sleep(0.02)
    assert os.path.exists(short_sock)
    stop.set()
    thread.join(10)
    assert not thread.is_alive()
    assert not os.path.exists(short_sock)


# ----------------------------------------------------------------------
# startup
# ----------------------------------------------------------------------

def test_a_missing_engine_refuses_rather_than_guesses(console, tmp_path):
    args = console.build_parser().parse_args(
        ["--metaannot", str(tmp_path / "nope.py")])
    with pytest.raises(console.Refuse) as e:
        console.build(args)
    assert "no such engine" in str(e.value)


def test_an_engine_that_cannot_describe_refuses_with_the_command(console,
                                                                 tmp_path):
    broken = tmp_path / "broken.py"
    broken.write_text("import sys\nsys.stderr.write('kaboom\\n')\n"
                      "sys.exit(3)\n")
    args = console.build_parser().parse_args(["--metaannot", str(broken)])
    with pytest.raises(console.Refuse) as e:
        console.build(args)
    assert "exited 3" in str(e.value) and "describe --json" in str(e.value)
    assert "kaboom" in str(e.value)


def test_the_banner_names_the_forward_the_operator_has_to_type(console,
                                                               contract):
    import io
    con = console.Console(contract, [], sys.executable, METAANNOT_PY)
    out = io.StringIO()
    console.banner(con, "/run/user/1000/metaannot.sock", out)
    text = out.getvalue()
    assert "ssh -N -L 8080:/run/user/1000/metaannot.sock" in text
    assert "$UID" not in text              # resolved, not left to expand
    assert "mode 0700, uid %d" % os.getuid() in text
    assert "read-only" in text
    assert "no token" in text


# ----------------------------------------------------------------------
# under a live run: a state file being replaced while the page is served
# ----------------------------------------------------------------------

def test_the_page_survives_a_state_file_being_replaced_underneath_it(
        console, results, contract, tmp_path):
    """The engine writes state as temp-file-then-os.replace, concurrently with
    every stage that finishes. This renders the page in a tight loop while a
    writer does exactly that, several hundred times, and nothing may raise, and
    no render may come back empty."""
    path = os.path.join(results.path, contract.state_name)
    project = console.Project(0, results.path, contract)
    results.state(running_state(contract))
    stop = threading.Event()
    trouble = []

    def writer():
        n = 0
        while not stop.is_set():
            n += 1
            state = running_state(contract, age=n % 100, extra={
                "emapper": {"signature": "a", "status": "ok", "seconds": 1.0,
                            "finished": stamp(time.time())}})
            tmp = "%s.%d.tmp" % (path, n)
            try:
                with open(tmp, "w") as fh:
                    json.dump(state, fh, indent=1, sort_keys=True)
                os.replace(tmp, path)
            except OSError:                 # pragma: no cover - the teardown
                return

    thread = threading.Thread(target=writer, daemon=True)
    thread.start()
    try:
        for _ in range(300):
            try:
                view = console.project_view(project, sys.executable,
                                            METAANNOT_PY, log=False)
                html = console.render_project(view, contract)
            except Exception as e:          # noqa: BLE001 - that is the point
                trouble.append(repr(e))
                break
            assert len(view["rows"]) == 21
            assert "metaannot console" in html
    finally:
        stop.set()
        thread.join(5)
    assert not trouble, trouble


def test_the_log_endpoint_survives_a_log_being_appended_underneath_it(
        console, served, results):
    """Every byte the writer wrote must arrive exactly once, in order, across
    however many polls it takes - that is what the offset is for."""
    path = os.path.join(results.path, "metaannot.log")
    stop, written = threading.Event(), []

    def writer():
        i = 0
        while not stop.is_set() and i < 2000:
            with open(path, "a") as fh:
                fh.write("[%7.1fs] INFO  line %d\n" % (i, i))
            written.append(i)
            i += 1

    thread = threading.Thread(target=writer, daemon=True)
    thread.start()
    _, _, body = get(served, "/api/log?p=0")
    first = json.loads(body)
    seen, off, ident = [], first["off"], first["ident"]
    for _ in range(60):
        _, _, body = get(served, "/api/log?p=0&off=%d&id=%s" % (off, ident))
        delta = json.loads(body)
        assert delta["reset"] is False, delta["note"]
        seen += [l["text"] for l in delta["lines"]]
        off, ident = delta["off"], delta["ident"]
        if len(seen) > 200:
            break
    stop.set()
    thread.join(10)
    numbers = [int(t.rsplit(" ", 1)[1]) for t in seen]
    assert numbers == list(range(numbers[0], numbers[0] + len(numbers)))


# ----------------------------------------------------------------------
# several projects, served at once, which is the case the doc names
# ----------------------------------------------------------------------

def test_several_projects_are_served_from_one_socket(console, contract,
                                                     tmp_path, short_sock):
    now = time.time()
    kinds = {"live": running_state(contract, age=8),
             "broken": running_state(contract),
             "finished": running_state(contract)}
    kinds["broken"][contract.run_key].update(final_status="failed",
                                             finished=stamp(now - 60))
    kinds["finished"][contract.run_key].update(final_status="ok",
                                               finished=stamp(now - 600))
    paths = []
    for name, state in kinds.items():
        d = tmp_path / name / "results"
        d.mkdir(parents=True)
        with open(str(d / contract.state_name), "w") as fh:
            json.dump(state, fh)
        with open(str(d / contract.log_name), "w") as fh:
            fh.write("[    0.0s] INFO  %s\n" % name)
        paths.append(str(d))

    con = console.Console(contract,
                          [console.Project(i, p, contract)
                           for i, p in enumerate(paths)],
                          sys.executable, METAANNOT_PY)
    stop, box = threading.Event(), {}
    thread = threading.Thread(
        target=console.serve, args=(con, short_sock),
        kwargs={"stop": stop, "ready": lambda s: box.setdefault("srv", s)},
        daemon=True)
    thread.start()
    try:
        for _ in range(200):
            if "srv" in box:
                break
            time.sleep(0.02)
        status, _, body = get(short_sock, "/")
        assert status == 200
        page = body.decode()
        for name in kinds:
            assert name in page
            assert 'href="/p/' in page
        assert "1 running" in page and "1 failed" in page
        # each project page is its own, and names its own log line
        for i, name in enumerate(kinds):
            _, _, one = get(short_sock, "/p/%d" % i)
            assert name.encode() in one
            assert b"hello" not in one
        _, _, body = get(short_sock, "/api/projects")
        assert len(json.loads(body)["projects"]) == 3
    finally:
        stop.set()
        thread.join(10)


# ----------------------------------------------------------------------
# the page under its own Content-Security-Policy
# ----------------------------------------------------------------------

def test_no_page_emits_an_inline_style_attribute(console, served, contract):
    """The console serves `Content-Security-Policy: default-src 'self'`, which
    blocks inline styles — silently. Column widths written as style attributes
    were dropped by the browser and every table came out with five equal
    columns, with nothing in the HTML or the tests to show for it. Every style
    decision therefore lives in app.css.
    """
    for route in ("/", "/?f=1", "/p/0", "/p/0?f=1"):
        _, _, body = get(served, route)
        assert b"style=" not in body, route
        assert b"onclick=" not in body and b"onload=" not in body, route
    _, _, js = get(served, "/app.js")
    assert b".style." not in js, ("CSSOM writes are a second place style is "
                                 "decided; use classes and `hidden`")
    _, headers, _ = get(served, "/")
    assert "'unsafe-inline'" not in headers["Content-Security-Policy"]


def test_every_column_class_the_tables_use_has_a_width_rule(console, served):
    """The classes that replaced the blocked style attributes. A class with no
    rule behind it is the same silent nothing the inline styles were."""
    import re as _re
    _, _, css = get(served, "/app.css")
    css = css.decode()
    for route in ("/", "/p/0"):
        _, _, body = get(served, route)
        heads = _re.findall(r'<th class="([^"]+)"', body.decode())
        assert heads, route
        for cls in heads:
            rule = _re.search(r'\.%s\s*\{([^}]*)\}' % _re.escape(cls), css)
            assert rule and "width" in rule.group(1), (route, cls)
    # and the three utility classes that were inline styles too
    for cls in ("note", "failmark", "path"):
        assert _re.search(r'\.%s\s*[,{]' % cls, css), cls


def test_a_directory_written_by_another_engine_version_says_so(
        console, results, contract):
    """The console asks whatever --metaannot points at for the stage graph, and
    that need not be the engine that made the directory it is looking at."""
    state = running_state(contract)
    state[contract.run_key]["version"] = "0.1.0-ancient"
    results.state(state)
    v = view_of(console, results, contract, log=False)
    html = console.render_project(v, contract)
    assert "written by metaannot 0.1.0-ancient" in html
    assert "the stage list below comes from metaannot %s" % contract.version \
        in html
    # and the matching case says nothing at all
    state[contract.run_key]["version"] = contract.version
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert "written by metaannot" not in console.render_project(v, contract)


def test_a_refused_connection_does_not_hand_out_a_slot(console, contract,
                                                       short_sock):
    """The refusal path must not return a permit it never took: one slot per
    refusal would raise the very cap the refusal is enforcing, and a client
    that retries in a loop would walk it up indefinitely."""
    fd = console.take_socket_name(short_sock)
    srv = console.ConsoleServer(short_sock, console.Handler,
                                console.Console(contract, [], sys.executable,
                                                METAANNOT_PY))
    try:
        for _ in range(console.MAX_CONNS):
            assert srv._slots.acquire(blocking=False)
        for _ in range(5):
            a, b = socket.socketpair()
            try:
                srv.process_request(a, "")      # refused: the cap is reached
                assert b.recv(4096).startswith(b"HTTP/1.1 503")
            finally:
                b.close()
        assert not srv._slots.acquire(blocking=False), \
            "five refusals handed out a slot each"
    finally:
        srv.server_close()
        os.close(fd)


def test_a_socket_that_cannot_be_bound_refuses_with_the_way_out(console,
                                                                contract,
                                                                short_sock):
    """An NFS home that will not hold an AF_UNIX socket is the case this is
    for, and there the operator needs the flag, not an OSError."""
    bad = os.path.join(short_sock, "nested", "s.sock")   # a path under a file
    open(short_sock, "w").close()
    con = console.Console(contract, [], sys.executable, METAANNOT_PY)
    with pytest.raises(console.Refuse) as e:
        console.serve(con, bad)
    assert "cannot create" in str(e.value) or "could not bind" in str(e.value)
    assert "--socket" in str(e.value)


# ======================================================================
# M1 review fixes. Every test below fails against the console as it was
# reviewed; the three security ones are the reviewers' own reproductions.
# ======================================================================

# ----------------------------------------------------------------------
# 1. a failed stage is the first thing in the pinned block
# ----------------------------------------------------------------------

def failed_mid_run(contract, ago=7200, error="CUDA out of memory"):
    """The reviewer's measured case: esmfold died two hours ago, the run
    carries on, and the heartbeat is fresh."""
    now = time.time()
    state = running_state(contract, age=21)
    state["esmfold"] = {"signature": None, "status": "failed",
                        "error": error, "finished": stamp(now - ago)}
    state["emapper"] = {"signature": "a", "status": "ok", "seconds": 12.0,
                        "finished": stamp(now - 3 * ago)}
    return state


def test_a_failed_stage_leads_the_pinned_vitals_block(console, results,
                                                      contract):
    """The single thing an operator opens this page to learn. It used to be
    343 px below the fold behind fifteen green OK rows, while the block that
    is pinned to the top talked about run ids and a fresh heartbeat."""
    results.state(failed_mid_run(contract))
    results.log("[    1.0s] INFO  carrying on\n")
    v = view_of(console, results, contract)
    assert [f["name"] for f in v["failures"]] == ["esmfold"]

    html = console.render_project(v, contract)
    # stage, time and the engine's own error string, in the pinned block
    vitals = html.split('class="panel vitals', 1)[1]
    assert "esmfold" in vitals
    assert "CUDA out of memory" in vitals
    # and before the heartbeat sentence, not after it
    assert html.index("CUDA out of memory") < html.index("Heartbeat stamped")
    # and before the run id, the pid and the host: those are identifiers,
    # this is the answer
    assert html.index("CUDA out of memory") < html.index("run 20260101T000000")
    # and before the stage table
    assert html.index("CUDA out of memory") < html.index("<thead>")


def test_a_fresh_heartbeat_never_reads_as_reassurance_over_a_dead_stage(
        console, results, contract):
    results.state(failed_mid_run(contract))
    v = view_of(console, results, contract, log=False)
    assert v["heartbeat"]["band"] == "fresh"
    html = console.render_project(v, contract)
    assert "says nothing about the failure above" in html
    # the block is banded red, not green, whatever the heartbeat says
    assert 'class="panel vitals sticky band-fail"' in html
    assert "band-fresh" not in html


def test_the_failure_block_names_the_time_and_carries_no_verdict(
        console, results, contract):
    results.state(failed_mid_run(contract, ago=7200))
    v = view_of(console, results, contract, log=False)
    f = v["failures"][0]
    assert f["error"] == "CUDA out of memory"
    assert 7000 < f["age"] < 7400
    html = console.render_project(v, contract)
    assert "2h00m ago" in html or "1h59m ago" in html
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase


# ----------------------------------------------------------------------
# 2. nothing from a file reaches a command string unvalidated
# ----------------------------------------------------------------------

def test_a_lock_pid_that_is_not_a_number_produces_no_command_to_paste(
        console, results, contract):
    """The reviewer's reproduction, verbatim. HTML escaping made this render
    cleanly, which is what invited the paste."""
    results.lock({"pid": "21877; rm -rf /", "host": "lab-fedora",
                  "started": "2026-09-07T22:10:41"})
    results.state(running_state(contract, age=20000))
    v = view_of(console, results, contract)
    assert v["lock"]["ps"] is None
    assert v["lock"]["pid_ok"] is False
    html = console.render_project(v, contract)
    # The served page rendered `ps -p 21877; rm -rf / -o pid,etime,stat,args`
    # inside a <code> block under "To settle it yourself, on that host:".
    assert "ps -p" not in html
    assert "rm -rf / -o pid" not in html
    assert "<code>" not in html
    assert "To settle it yourself" not in html
    # The raw field is still SHOWN, as data, so the operator can see what is
    # actually in the lock file - just never as something to paste.
    assert "Held by pid 21877; rm -rf / on lab-fedora" in html
    assert "is not a number" in html


@pytest.mark.parametrize("pid", [
    "21877; rm -rf /", "$(id)", "1 2", "-1", "0", "", None, True,
    ["21877"], {"pid": 1}, 1.5, "0x21", " 21877\n; ls",
])
def test_no_lock_pid_shape_but_a_positive_integer_becomes_a_command(
        console, pid):
    assert console.valid_pid(pid) is None


@pytest.mark.parametrize("pid,want", [(21877, 21877), ("21877", 21877),
                                      (" 21877 ", 21877), (21877.0, 21877)])
def test_a_real_pid_still_hands_over_the_ps_line(console, pid, want):
    assert console.valid_pid(pid) == want


# ----------------------------------------------------------------------
# 3. every file this console opens is one it does not control
# ----------------------------------------------------------------------

def test_a_fifo_in_place_of_a_state_file_does_not_block_a_read(console,
                                                               tmp_path,
                                                               contract):
    d = tmp_path / "fifo" / "results"
    d.mkdir(parents=True)
    os.mkfifo(str(d / contract.state_name))
    project = console.Project(0, str(d), contract)
    done = []

    def read():
        done.append(console.project_view(project, sys.executable,
                                         METAANNOT_PY, log=False))

    t = threading.Thread(target=read, daemon=True)
    t.start()
    t.join(20)
    assert done, "reading a FIFO state file never returned"
    v = done[0]
    assert v["state_status"] == "bad"
    assert "not a regular file" in (v["state_detail"] or "")


def test_a_fifo_state_file_does_not_wedge_every_route_including_app_css(
        console, tmp_path, contract, short_sock):
    """The reviewer's reproduction. os.open(path, O_RDONLY) on a FIFO blocks
    forever in a serving thread; the index reads every project's state file,
    so exactly MAX_CONNS index loads exhausted the semaphore and every route
    after that, /app.css included, was 503 until the process was restarted."""
    d = tmp_path / "wedge" / "results"
    d.mkdir(parents=True)
    os.mkfifo(str(d / contract.state_name))
    con = console.Console(contract, [console.Project(0, str(d), contract)],
                          sys.executable, METAANNOT_PY)
    stop, box = threading.Event(), {}
    thread = threading.Thread(
        target=console.serve, args=(con, short_sock),
        kwargs={"stop": stop, "ready": lambda s: box.setdefault("srv", s)},
        daemon=True)
    thread.start()
    try:
        for _ in range(200):
            if "srv" in box:
                break
            time.sleep(0.02)
        assert "srv" in box

        def quick(route, timeout=5):
            conn = UnixConn(short_sock)
            try:
                conn.connect()
                conn.sock.settimeout(timeout)
                conn.request("GET", route)
                r = conn.getresponse()
                r.read()
                return r.status
            finally:
                conn.close()

        # One index load must come back at all. A wedge shows up right here.
        assert quick("/") == 200
        assert quick("/p/0") == 200
        for _ in range(console.MAX_CONNS + 4):
            assert quick("/") == 200
        assert quick("/app.css") == 200
    finally:
        stop.set()
        thread.join(10)


def test_a_fifo_in_place_of_a_log_is_not_read_either(console, tmp_path,
                                                     contract):
    d = tmp_path / "logfifo" / "results"
    d.mkdir(parents=True)
    os.mkfifo(str(d / contract.log_name))
    out = []
    t = threading.Thread(target=lambda: out.append(
        console.tail_log(os.path.join(str(d), contract.log_name))), daemon=True)
    t.start()
    t.join(20)
    assert out, "reading a FIFO log never returned"
    assert out[0]["ok"] is False
    assert "not a regular file" in out[0]["note"]


def test_a_directory_in_place_of_a_state_file_is_refused_not_read(
        console, results, contract):
    os.mkdir(os.path.join(results.path, contract.state_name))
    v = view_of(console, results, contract, log=False)
    assert v["state_status"] == "bad"
    assert "not a regular file" in (v["state_detail"] or "")


# ----------------------------------------------------------------------
# symlinks: the console's own files, at a path it was handed
# ----------------------------------------------------------------------

def test_the_socket_lock_is_never_opened_through_a_symlink(console,
                                                          short_sock):
    """The reviewer planted `x.sock.lock -> .../PLANTED_BY_ATTACKER` and the
    call created that file, owned by the console's uid. O_CREAT|O_RDWR with no
    O_NOFOLLOW follows a symlink both ways: onto an existing file, and into
    creating one at a target that does not exist yet.

    short_sock, not tmp_path: on macOS a pytest tmp_path is already past the
    104-byte AF_UNIX limit, so take_socket_name refuses on the length before
    it ever opens the lock and the test would pass without testing anything.
    """
    sock = short_sock
    target = pathlib.Path(os.path.dirname(sock)) / "PLANTED_BY_ATTACKER"
    os.symlink(str(target), sock + ".lock")
    with pytest.raises(console.Refuse) as e:
        console.take_socket_name(sock)
    assert not target.exists(), "a file was created at the symlink's target"
    assert "lock" in str(e.value)


def test_the_socket_lock_refuses_an_existing_file_reached_by_symlink(
        console, short_sock):
    sock = short_sock
    victim = pathlib.Path(os.path.dirname(sock)) / "someone_elses_file"
    victim.write_text("mine", encoding="utf-8")
    os.symlink(str(victim), sock + ".lock")
    with pytest.raises(console.Refuse):
        console.take_socket_name(sock)
    assert victim.read_text(encoding="utf-8") == "mine"


def test_runtime_dir_does_not_accept_a_symlink_to_a_private_directory(
        console, tmp_path):
    """The build report says the runtime directory must lstat as a 0700
    directory this uid owns. os.stat resolved the link first, so the checks
    described the target while the socket would be bound at the link."""
    real = tmp_path / "real"
    real.mkdir(mode=0o700)
    link = tmp_path / "link"
    os.symlink(str(real), str(link))
    with pytest.raises(console.Refuse):
        console.runtime_dir(candidates=[str(link)])


def test_an_explicit_socket_in_a_world_writable_directory_is_refused(
        console, tmp_path):
    """The bind-failure message tells the operator to pass
    `--socket /tmp/metaannot-console-<uid>/metaannot.sock`, and an explicit
    --socket used to skip every check runtime_dir applies to its own pick."""
    d = tmp_path / "loose"
    d.mkdir(mode=0o777)
    os.chmod(str(d), 0o777)
    args = console.build_parser().parse_args(
        ["--socket", str(d / "s.sock"), "--metaannot", METAANNOT_PY])
    with pytest.raises(console.Refuse) as e:
        console.socket_path(args)
    assert "private directory" in str(e.value)
    # NOT `mkdir -m 700 -p <the parent they passed>`. For `--socket
    # /tmp/pwn.sock` that read "Try: mkdir -m 700 -p /tmp" — advice the
    # operator must not take, and a silent no-op on a directory that already
    # exists, so a diligent one runs it, sees nothing happen and retries into
    # the same refusal. The advice names a directory of the console's own.
    assert "mkdir -m 700" not in str(e.value)
    assert str(d) not in str(e.value).split("Leave --socket off")[1]
    assert "~/.cache/metaannot-console" in str(e.value)


def test_an_explicit_socket_in_a_private_directory_is_accepted(console,
                                                               tmp_path):
    d = tmp_path / "tight"
    d.mkdir(mode=0o700)
    os.chmod(str(d), 0o700)
    args = console.build_parser().parse_args(
        ["--socket", str(d / "s.sock"), "--metaannot", METAANNOT_PY])
    assert console.socket_path(args) == str(d / "s.sock")


# ----------------------------------------------------------------------
# a corrupt directory is not an empty one, on either page
# ----------------------------------------------------------------------

def three_kinds(tmp_path, contract):
    """One never-started directory, one empty state file, one corrupt."""
    paths = []
    for name, body in (("never", None), ("blank", ""),
                       ("corrupt", "{ this is not json")):
        d = tmp_path / name / "results"
        d.mkdir(parents=True)
        with open(str(d / contract.log_name), "w") as fh:
            fh.write("[    0.0s] INFO  x\n")
        if body is not None:
            with open(str(d / contract.state_name), "w") as fh:
                fh.write(body)
        paths.append(str(d))
    return paths


def test_the_index_tells_a_corrupt_directory_from_an_empty_one(console,
                                                               tmp_path,
                                                               contract):
    """The front page said NO RECORD and counted all three under
    "never started" — the same badge, strip and sentence for a directory
    nothing has run in and one whose state file will not parse."""
    paths = three_kinds(tmp_path, contract)
    projects = [console.Project(i, p, contract) for i, p in enumerate(paths)]
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    by_name = {i["name"]: i for i in view["projects"]}
    assert by_name["never"]["bucket"] == "none"
    assert by_name["blank"]["bucket"] == "blank"
    assert by_name["corrupt"]["bucket"] == "unreadable"
    assert "1 never started" in view["summary"]
    assert "1 with an unreadable state file" in view["summary"]
    assert "1 with an empty state file" in view["summary"]

    html = console.render_index(view, contract)
    assert "UNREADABLE" in html and "EMPTY STATE" in html
    assert "will not parse" in html
    assert "not an empty directory" in html


def test_a_corrupt_state_file_leads_the_project_page(console, results,
                                                     contract):
    """The pinned block led with "There is no run record in the state file.
    Either nothing has run here, or it was written by a metaannot older than
    the one that added the record." — both explanations wrong — and the true
    one was the fourth line down in muted small text."""
    results.raw_state('{"emapper": {"status": "ok", "fin')
    v = view_of(console, results, contract, log=False)
    assert v["state_trouble"] is not None
    html = console.render_project(v, contract)
    assert "Either nothing has run here" not in html
    assert "Unterminated string" in html
    assert "not an empty directory" in html
    # first, in the pinned block, above the heartbeat sentence
    assert html.index("Unterminated string") < html.index("no run record")
    assert 'band-trouble' in html


def test_an_empty_state_file_is_not_reported_as_a_fresh_directory(
        console, results, contract):
    results.raw_state("")
    v = view_of(console, results, contract, log=False)
    html = console.render_project(v, contract)
    assert "a write that did not complete" in html
    assert "There is no state file here yet" not in html


def test_an_unreadable_results_directory_is_not_reported_as_an_empty_one(
        console, tmp_path, contract):
    """chmod 000: the state and lock readers said "Permission denied" while
    the vitals line said "There is no log file in this directory", because
    stat_or_none folded EACCES into None."""
    if os.getuid() == 0:
        pytest.skip("root reads everything")
    d = tmp_path / "locked" / "results"
    d.mkdir(parents=True)
    with open(str(d / contract.log_name), "w") as fh:
        fh.write("[    0.0s] INFO  x\n")
    os.chmod(str(d), 0o000)
    try:
        project = console.Project(0, str(d), contract)
        v = console.project_view(project, sys.executable, METAANNOT_PY,
                                 log=False)
        assert v["log_why"] and v["log_why"] != "missing"
        html = console.render_project(v, contract)
        assert "There is no log file in this directory" not in html
        assert "cannot be read" in html
    finally:
        os.chmod(str(d), 0o700)


# ----------------------------------------------------------------------
# a run record with no final_status
# ----------------------------------------------------------------------

def test_a_run_record_with_no_final_status_is_never_called_finished(
        console, results, contract):
    """`if status != "running"` turned an ABSENT key into "Finished ? at ?."
    — the strongest possible claim in the one paragraph whose doctrine is
    that unprovable means alive — and switched the staleness apparatus off
    while a stage row counted up in RUN."""
    now = time.time()
    state = running_state(contract, age=120)
    del state[contract.run_key]["final_status"]
    state["interpro"] = {"signature": None, "status": "running",
                         "started": stamp(now - 7000)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    hb = v["heartbeat"]
    assert "Finished" not in (hb["verdict"] or "")
    assert hb["running"] is True
    assert hb["band"] in ("fresh", "late", "long")
    assert hb["beat_age"] is not None
    assert any("no final_status" in line for line in hb["lines"])
    html = console.render_project(v, contract)
    assert "Finished ? at ?" not in html
    assert "unprovable means alive" in html


@pytest.mark.parametrize("value", [None, "", 0])
def test_an_unusable_final_status_reads_as_still_going(console, value):
    run = {"final_status": value} if value is not None else {}
    assert console.run_is_over(run) is False


def test_a_real_final_status_still_reads_as_over(console):
    assert console.run_is_over({"final_status": "ok"}) is True
    assert console.run_is_over({"final_status": "running"}) is False


# ----------------------------------------------------------------------
# a finished run is not waiting on anything
# ----------------------------------------------------------------------

def test_a_finished_run_has_no_present_tense_waiting_rows(console, results,
                                                          contract):
    """heartbeat_view suppresses the whole staleness apparatus once the run is
    over; stage_rows had no equivalent notion, so the vitals said "Finished ok
    at … · ran 2d 04h" while the table under it said "waiting on pfam,
    signalp, dbcan"."""
    now = time.time()
    state = running_state(contract)
    state[contract.run_key].update(final_status="ok",
                                   finished=stamp(now - 3600))
    state["emapper"] = {"signature": "a", "status": "ok", "seconds": 3.0,
                        "finished": stamp(now - 3700)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    for r in v["rows"]:
        assert r["state"] != "wait", r["name"]
        assert r["label"] != "WAIT", r["name"]
        assert r["label"] != "NEXT", r["name"]
    body = console.render_project(v, contract).split("<tbody>", 1)[1]
    assert "waiting on" not in body
    assert "no record yet" not in body
    assert "This run ended" in body


def test_a_live_run_still_says_what_a_stage_is_waiting_on(console, results,
                                                          contract):
    results.state(running_state(contract, age=10))
    v = view_of(console, results, contract, log=False)
    assert any(r["label"] == "WAIT" for r in v["rows"])
    assert "waiting on" in console.render_project(v, contract)


# ----------------------------------------------------------------------
# /api/log discards nothing in silence
# ----------------------------------------------------------------------

def test_the_log_endpoint_says_when_a_burst_outran_the_pane(served, console,
                                                            contract):
    """500 lines appended between two polls returned 200 lines with note=None
    and reset=False; the offset advanced past all 500, so 300 lines were gone
    from the live pane with nothing said."""
    status, _, body = get(served, "/api/log?p=0")
    first = json.loads(body)
    assert first["ok"] and first["note"] is None

    path = os.path.join(json.loads(get(served, "/api/project?p=0")[2])["path"],
                        contract.log_name)
    with open(path, "a") as fh:
        for i in range(500):
            fh.write("[    1.0s] INFO  new line %d\n" % i)

    status, _, body = get(
        served, "/api/log?p=0&off=%d&id=%s" % (first["off"], first["ident"]))
    d = json.loads(body)
    assert d["reset"] is False
    assert len(d["lines"]) == console.LOG_LINES
    assert d["dropped"] == 300
    assert d["note"], "300 lines vanished and nothing was said"
    assert "500 lines arrived" in d["note"]
    assert "300 oldest" in d["note"]


def test_the_log_endpoint_is_silent_when_nothing_was_dropped(served):
    status, _, body = get(served, "/api/log?p=0")
    d = json.loads(body)
    assert d["dropped"] == 0 and d["note"] is None


def test_log_lines_reports_what_it_left_out(console):
    text = "".join("[    1.0s] INFO  line %d\n" % i for i in range(500))
    lines, dropped = console.log_lines(text)
    assert len(lines) == console.LOG_LINES
    assert dropped == 300
    assert lines[0]["text"].endswith("line 300")
    assert "300 oldest" in console.dropped_note(dropped, reset=False)
    assert "300 older ones" in console.dropped_note(dropped, reset=True)


def test_a_log_window_with_no_line_break_is_not_reported_as_empty(console,
                                                                  results):
    """300 KB of a single line rendered "— no log lines —" with note=None,
    which is exactly what an EMPTY log looks like."""
    out = console.tail_log(results.log("x" * 300000))
    assert out["ok"] is True
    assert out["text"] == ""
    assert out["note"] and "no line break" in out["note"]


def test_an_empty_log_still_says_nothing_about_line_breaks(console, results):
    out = console.tail_log(results.log(""))
    assert out["note"] is None


# ----------------------------------------------------------------------
# the console does not contradict the engine about what will run
# ----------------------------------------------------------------------

def effective(results, contract, run_block):
    """Write a config.effective.yaml with a literal `run:` block."""
    body = ["proteins_faa: /nowhere/p.faa", "quant_table: /nowhere/q.tsv",
            "run:"]
    for key, text in run_block.items():
        body.append("  %s:%s" % (key, (" " + text) if text else ""))
    results.write(contract.effective_name, "\n".join(body) + "\n")


@pytest.mark.parametrize("literal,shown", [
    ("0", "0"), ("", "empty"), ('""', '""'), ("[]", "[]"),
    ("false", "false"), ("no", "false"), ("off", "false"),
])
def test_a_falsy_run_flag_is_reported_off_not_next(console, results, contract,
                                                   literal, shown):
    """metaannot.py skips on `not cfg["run"].get(st["enabled"], False)` —
    plain truthiness — while this file used `is False`, so `run.pfam: 0`,
    `: null`, `: ""` and `: []` each made the engine skip a stage the console
    announced as NEXT. `run:\\n  pfam:` is the realistic form."""
    results.state(running_state(contract, age=10))
    effective(results, contract, {"pfam": literal})
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "pfam")
    assert row["state"] == "off", row["detail"]
    assert row["label"] == "OFF"
    assert shown in row["detail"]
    assert "run.pfam" in row["detail"]


def test_a_true_run_flag_is_still_not_off(console, results, contract):
    results.state(running_state(contract, age=10))
    effective(results, contract, {"pfam": "true"})
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "pfam")
    assert row["state"] != "off"


def test_a_stage_with_a_record_is_never_called_off(console, results, contract):
    now = time.time()
    state = running_state(contract, age=10)
    state["pfam"] = {"signature": "a", "status": "ok", "seconds": 2.0,
                     "finished": stamp(now - 10)}
    results.state(state)
    effective(results, contract, {"pfam": "0"})
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "pfam")
    assert row["state"] == "ok"


# ----------------------------------------------------------------------
# a stage record that is not a dict
# ----------------------------------------------------------------------

@pytest.mark.parametrize("value", ["a string", 42, ["a", "list"], True])
def test_a_stage_record_that_is_not_an_object_is_a_bad_row(console, results,
                                                           contract, value):
    """It was silently downgraded to "not reached", which reads as nothing
    being wrong. The `bad` row type exists for exactly this."""
    state = running_state(contract, age=10)
    state["dbcan"] = value
    results.state(state)
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "dbcan")
    assert row["state"] == "bad"
    assert row["label"] == "?"
    assert "not an object" in row["detail"]
    assert "no dependency is blocking it" not in row["detail"]


def test_a_stage_the_contract_does_not_know_is_named_not_dropped(console,
                                                                 results,
                                                                 contract):
    """A directory written by a newer engine. These records were dropped with
    no mention whenever _run.version happened to match."""
    state = running_state(contract, age=10)
    state["quantum_fold"] = {"status": "ok", "seconds": 1.0}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert v["unknown_stages"] == ["quantum_fold"]
    html = console.render_project(v, contract)
    assert "quantum_fold" in html
    assert "no stage for" in html


# ----------------------------------------------------------------------
# the part file is evidence of work only while it is growing
# ----------------------------------------------------------------------

def part_file_project(results, contract, part_age, beat_age):
    """A running interpro whose part file is `part_age` old."""
    now = time.time()
    state = running_state(contract, age=beat_age)
    state["interpro"] = {"signature": None, "status": "running",
                         "started": stamp(now - 3600)}
    for dep in ("emapper",):
        state[dep] = {"signature": "a", "status": "ok", "seconds": 1.0,
                      "finished": stamp(now - 3600)}
    results.state(state)
    out = None
    for st in contract.stages:
        if st["name"] == "interpro" and st["outputs"]:
            out = st["outputs"][0]
            break
    assert out
    full = os.path.join(results.path, out)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    stem = os.path.splitext(os.path.basename(full))[0]
    part = os.path.join(os.path.dirname(full),
                        ".%s.21877.140.part.tsv" % stem)
    with open(part, "w") as fh:
        fh.write("x" * 431923)
    os.utime(part, (now - part_age, now - part_age))
    return part


def test_a_growing_part_file_is_still_offered_as_the_best_evidence(
        console, results, contract):
    part_file_project(results, contract, part_age=5, beat_age=140)
    v = view_of(console, results, contract, log=False)
    html = console.render_project(v, contract)
    assert "more direct sign of work than any heartbeat" in html


def test_a_stalled_part_file_is_not_offered_as_reassurance(console, results,
                                                           contract):
    """Observed: "No heartbeat for 2m18s … so that is late", then two lines
    later "interpro is writing … last grew 26m43s ago. That is a more direct
    sign of work than any heartbeat." A 26-minute-old file offered as stronger
    evidence of work than a 2-minute-old heartbeat inverts the meaning."""
    part_file_project(results, contract, part_age=1603, beat_age=138)
    v = view_of(console, results, contract, log=False)
    html = console.render_project(v, contract)
    assert "more direct sign of work than any heartbeat" not in html
    assert "STALER than every other clock" in html
    assert "stopped producing bytes" in html


def test_the_part_cache_does_not_rescan_a_directory_that_has_not_changed(
        console, results, contract, monkeypatch):
    """newest_part_file does a full scandir per running stage per poll, with
    no early exit; for esmfold that directory is structures/, one file per
    dark protein."""
    part = part_file_project(results, contract, part_age=5, beat_age=10)
    out = os.path.join(os.path.dirname(part),
                       os.path.basename(part)[1:].replace(
                           ".21877.140.part", ""))
    scans = []
    real = console.newest_part_file
    monkeypatch.setattr(console, "newest_part_file",
                        lambda p: (scans.append(p), real(p))[1])
    cache = console.PartCache()
    first = cache.newest(out)
    assert first and first["size"] == 431923
    for _ in range(20):
        again = cache.newest(out)
    assert len(scans) == 1, "the directory was re-scanned %d times" % len(scans)
    # and it still reports the file growing
    with open(part, "a") as fh:
        fh.write("more")
    assert cache.newest(out)["size"] == 431927
    assert len(scans) == 1


# ----------------------------------------------------------------------
# --root is scanned again, and the page says which
# ----------------------------------------------------------------------

def test_a_dataset_launched_after_startup_appears_without_a_restart(
        console, tmp_path, contract):
    """discover() ran inside build(), before the server existed, and the
    footer said "Snapshot <fresh time>" beside a list fixed at startup: you
    kick off dataset five, refresh the tab, and it is silently missing."""
    first = tmp_path / "ds01" / "results"
    first.mkdir(parents=True)
    with open(str(first / contract.state_name), "w") as fh:
        fh.write("{}")
    args = console.build_parser().parse_args(
        ["--root", str(tmp_path), "--python", sys.executable,
         "--metaannot", METAANNOT_PY])
    con = console.build(args)
    con.rescan_s = 0
    assert [p.name for p in con.projects] == ["ds01"]

    later = tmp_path / "Brand_new_run" / "results"
    later.mkdir(parents=True)
    with open(str(later / contract.state_name), "w") as fh:
        fh.write("{}")
    view = con.index()
    assert sorted(i["name"] for i in view["projects"]) == ["Brand_new_run",
                                                           "ds01"]
    # and the one that was already there keeps its number, so a bookmark holds
    assert con.project("0").name == "ds01"
    assert con.project("1").name == "Brand_new_run"
    assert "re-scanned" in console.render_index(view, contract)


def test_a_console_with_no_root_says_the_list_cannot_change(console, results,
                                                            contract):
    con = console.Console(contract,
                          [console.Project(0, results.path, contract)],
                          sys.executable, METAANNOT_PY)
    html = console.render_index(con.index(), contract)
    assert "does not change while the console runs" in html
    assert "re-scanned" not in html


def test_the_banner_says_whether_the_watched_list_can_change(console, results,
                                                             contract):
    import io
    con = console.Console(contract,
                          [console.Project(0, results.path, contract)],
                          sys.executable, METAANNOT_PY)
    out = io.StringIO()
    console.banner(con, "/tmp/x.sock", out=out)
    assert "this list is fixed for the life of the process" in out.getvalue()

    con2 = console.Console(contract, [], sys.executable, METAANNOT_PY,
                           roots=["/data/runs"])
    out2 = io.StringIO()
    console.banner(con2, "/tmp/x.sock", out=out2)
    assert "keeps its number" in out2.getvalue()


def test_discovery_follows_a_symlinked_dataset_directory(console, tmp_path,
                                                          contract):
    """`<root>/DatasetA -> /mnt/nvme/DatasetA` is an ordinary way to keep one
    dataset on another mount; is_dir(follow_symlinks=False) excluded it and
    the operator got an index that was simply short, with no message."""
    elsewhere = tmp_path / "other_mount" / "DatasetA" / "results"
    elsewhere.mkdir(parents=True)
    with open(str(elsewhere / contract.state_name), "w") as fh:
        fh.write("{}")
    root = tmp_path / "runs"
    (root / "DatasetB" / "results").mkdir(parents=True)
    with open(str(root / "DatasetB" / "results" / contract.state_name),
              "w") as fh:
        fh.write("{}")
    os.symlink(str(tmp_path / "other_mount" / "DatasetA"),
               str(root / "DatasetA"))

    found = console.discover([str(root)], contract)
    assert sorted(found) == [str(root / "DatasetA" / "results"),
                             str(root / "DatasetB" / "results")]


def test_discovery_follows_a_symlink_straight_at_a_results_directory(
        console, tmp_path, contract):
    elsewhere = tmp_path / "other_mount" / "DatasetA" / "results"
    elsewhere.mkdir(parents=True)
    with open(str(elsewhere / contract.state_name), "w") as fh:
        fh.write("{}")
    root = tmp_path / "runs"
    root.mkdir()
    os.symlink(str(elsewhere), str(root / "DatasetC"))
    assert console.discover([str(root)], contract) == \
        [str(root / "DatasetC")]


def test_discovery_still_refuses_to_loop_through_a_symlink(console, tmp_path,
                                                           contract):
    root = tmp_path / "loop"
    (root / "a").mkdir(parents=True)
    os.symlink(str(root), str(root / "a" / "back"))
    assert console.discover([str(root)], contract) == []


# ----------------------------------------------------------------------
# the polling cadence actually adapts
# ----------------------------------------------------------------------

def test_the_running_flag_lives_inside_the_fragment_the_poller_replaces(
        console, results, contract):
    """data-running was an attribute of #live, and the fragment endpoints
    return that element's CONTENTS, so box.getAttribute('data-running') after
    innerHTML = t always returned the page-load value. A tab opened before a
    run started polled every 60 s through the whole run."""
    results.state(running_state(contract, age=10))
    v = console.decorate(view_of(console, results, contract, log=False),
                         contract)
    frag = console.render_project_body(v)
    assert 'id="livestate"' in frag
    assert 'data-running="1"' in frag
    assert "data-sig=" in frag

    html = console.render_project(v, contract)
    live = html.split('<div id="live"', 1)[1].split(">", 1)[0]
    assert "data-running" not in live, "still on an element never replaced"
    assert "#livestate" in console.PAGE_JS


def test_the_index_fragment_carries_the_flag_too(console, results, contract):
    results.state(running_state(contract, age=10))
    con = console.Console(contract,
                          [console.Project(0, results.path, contract)],
                          sys.executable, METAANNOT_PY)
    frag = console.render_index_body(con.index())
    assert 'id="livestate"' in frag and 'data-running="1"' in frag


def test_the_change_signature_ignores_the_clock(console, results, contract):
    """`t !== box.innerHTML` was true on every poll because the fragment
    embeds relative ages, so `quiet` never incremented and the 15 s back-off
    branch was unreachable."""
    results.state(running_state(contract, age=10))
    results.log("[    0.0s] INFO  x\n")
    project = console.Project(0, results.path, contract)
    a = console.project_view(project, sys.executable, METAANNOT_PY, log=False)
    b = console.project_view(project, sys.executable, METAANNOT_PY, log=False)
    assert a["taken"] != b["taken"]
    assert console.project_sig(a) == console.project_sig(b)

    results.log("[    1.0s] INFO  y\n[    2.0s] INFO  z\n")
    c = console.project_view(project, sys.executable, METAANNOT_PY, log=False)
    assert console.project_sig(c) != console.project_sig(a)


def test_the_change_signature_carries_no_content_from_a_watched_file(
        console, results, contract):
    results.lock({"pid": 21877, "host": "lab-fedora",
                  "started": "2026-09-07T22:10:41"})
    results.state(running_state(contract, age=10))
    v = console.decorate(view_of(console, results, contract, log=False),
                         contract)
    sig = console.project_sig(v)
    assert re.fullmatch(r"[0-9a-f]{16}", sig)
    assert "lab-fedora" not in sig and "21877" not in sig


# ----------------------------------------------------------------------
# --interval, checked
# ----------------------------------------------------------------------

@pytest.mark.parametrize("bad", ["0", "-1", "0.1", "nan", "inf", "banana"])
def test_a_useless_interval_is_refused_at_the_command_line(console, bad):
    """`--interval 0` made period() return 0 and the page run
    setTimeout(pull, 0) — a busy fetch loop against the NFS mount the comments
    in this file worry about."""
    with pytest.raises(SystemExit):
        console.build_parser().parse_args(["--interval", bad])


def test_a_usable_interval_is_accepted(console):
    args = console.build_parser().parse_args(["--interval", "0.5"])
    assert args.interval == 0.5
    assert console.build_parser().parse_args([]).interval == 3.0


# ----------------------------------------------------------------------
# the SIGHUP an operator asked for is left alone
# ----------------------------------------------------------------------

SIGHUP_PROBE = '''
import importlib.util, signal, sys, threading
spec = importlib.util.spec_from_file_location("c", sys.argv[1])
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
signal.signal(signal.SIGHUP, signal.SIG_IGN)          # what nohup does
con = m.Console(m.Contract(m.run_describe(sys.executable, sys.argv[2])), [],
                sys.executable, sys.argv[2])
stop = threading.Event()


def ready(srv):
    print("HUP=%s" % ("IGN" if signal.getsignal(signal.SIGHUP)
                      is signal.SIG_IGN else "REPLACED"))
    print("TERM=%s" % ("IGN" if signal.getsignal(signal.SIGTERM)
                       is signal.SIG_IGN else "INSTALLED"))
    sys.stdout.flush()
    stop.set()


m.serve(con, sys.argv[3], ready=ready, stop=stop)
'''


def test_a_sighup_already_ignored_is_left_ignored(console, short_sock,
                                                  tmp_path):
    """signal.signal(SIGHUP, …) unconditionally replaced the SIG_IGN nohup
    installs, so `nohup python3 console/console.py … &` died at logout — the
    one command an operator reaches for on a box without tmux."""
    probe = tmp_path / "probe.py"
    probe.write_text(SIGHUP_PROBE, encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, str(probe), console.__file__, METAANNOT_PY,
         short_sock], capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr
    assert "HUP=IGN" in proc.stdout, proc.stdout
    # and the ones nobody ignored are still installed
    assert "TERM=INSTALLED" in proc.stdout, proc.stdout


# ----------------------------------------------------------------------
# the banner never precedes a refusal
# ----------------------------------------------------------------------

def test_a_refusal_is_never_printed_under_a_success_banner(console, contract,
                                                           short_sock):
    """banner() was printed before take_socket_name(), so a refusal such as
    "thesis.docx exists and is not a socket" appeared AFTER a complete success
    banner ending "Ctrl-C to stop."."""
    d = os.path.dirname(short_sock)
    os.chmod(d, 0o700)
    thesis = pathlib.Path(d) / "t.docx"
    thesis.write_text("chapter one", encoding="utf-8")
    said = []
    with pytest.raises(console.Refuse) as e:
        console.serve(bare_console(console, contract), str(thesis),
                      announce=lambda: said.append(1))
    assert "is not a socket" in str(e.value)
    assert said == [], "the banner was printed before the socket was taken"
    assert thesis.read_text(encoding="utf-8") == "chapter one"


def test_main_prints_the_refusal_and_no_banner(console, short_sock, capsys,
                                               monkeypatch):
    d = os.path.dirname(short_sock)
    os.chmod(d, 0o700)
    thesis = pathlib.Path(d) / "t.docx"
    thesis.write_text("chapter one", encoding="utf-8")
    said = []
    monkeypatch.setattr(console, "banner", lambda *a, **k: said.append(1))
    code = console.main(["--socket", str(thesis), "--metaannot", METAANNOT_PY,
                         "--python", sys.executable])
    assert code == 2
    assert said == [], "the banner went out before the socket was taken"
    assert "is not a socket" in capsys.readouterr().err


# ----------------------------------------------------------------------
# the small ones the reviewers batched
# ----------------------------------------------------------------------

def test_a_shared_path_prefix_is_trimmed_from_the_left_not_the_right(console):
    """With /data/runs/<dataset>/results every row showed the identical shared
    prefix and CSS ellipsis cut off exactly the discriminating segment."""
    assert console.elide_path("/data/runs/Cohort_B/results") == \
        "…/Cohort_B/results"
    assert console.elide_path("/short") == "/short"
    assert "text-overflow" not in console.PAGE_CSS.split(".path {", 1)[1] \
        .split("}", 1)[0]


def test_the_index_carries_a_visible_legend_for_the_strip(console, results,
                                                          contract):
    results.state(running_state(contract, age=10))
    con = console.Console(contract,
                          [console.Project(0, results.path, contract)],
                          sys.executable, METAANNOT_PY)
    html = console.render_index(con.index(), contract)
    assert 'class="legend"' in html
    for word in ("finished", "running", "failed", "not reached",
                 "off in the config"):
        assert word in html
    assert ".legend" in console.PAGE_CSS


def test_newest_activity_never_names_a_directory_nobody_ran_in(console,
                                                               tmp_path,
                                                               contract):
    """"newest activity Nm ago (X)" could name a project whose badge is
    NO RECORD — a directory someone merely touched."""
    now = time.time()
    live = tmp_path / "live" / "results"
    live.mkdir(parents=True)
    with open(str(live / contract.state_name), "w") as fh:
        json.dump(running_state(contract, age=5), fh)
    with open(str(live / contract.log_name), "w") as fh:
        fh.write("[    0.0s] INFO  x\n")
    os.utime(str(live / contract.log_name), (now - 30, now - 30))
    os.utime(str(live / contract.state_name), (now - 30, now - 30))

    touched = tmp_path / "touched" / "results"
    touched.mkdir(parents=True)
    with open(str(touched / contract.log_name), "w") as fh:
        fh.write("")                     # touch, and nothing else

    projects = [console.Project(i, str(p), contract)
                for i, p in enumerate([live, touched])]
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    assert "(live)" in view["summary"]
    assert "(touched)" not in view["summary"]


def test_a_config_path_that_is_not_a_regular_file_is_never_handed_to_the_engine(
        console, results, contract, tmp_path):
    """`_run.config_path`, read out of the watched directory's own state file,
    becomes an argv element of `describe --json --config <that path>`."""
    fifo = tmp_path / "wedge.yaml"
    os.mkfifo(str(fifo))
    state = running_state(contract, age=10)
    state[contract.run_key]["config_path"] = str(fifo)
    results.state(state)
    forks = []
    project = console.Project(0, results.path, contract)
    real = console.run_describe

    def counting(*a, **kw):
        forks.append(kw.get("config"))
        return real(*a, **kw)

    console.run_describe = counting
    try:
        run, note = project.config_run(sys.executable, METAANNOT_PY, state)
    finally:
        console.run_describe = real
    assert run is None
    assert "no readable config" in note
    assert forks == [], "the engine was pointed at %r" % forks


def test_the_engines_own_reason_survives_into_the_page(console, results,
                                                       contract):
    """`str(e).splitlines()[0]` is "the engine exited 1:" — a message carrying
    no diagnostic content, for a distinction that governs 9 of 21 rows."""
    results.state(running_state(contract, age=10))
    results.write(contract.effective_name, "run:\n  pfam: [unclosed\n")
    v = view_of(console, results, contract, log=False)
    assert v["cfg_run_known"] is False
    assert "the engine exited" in v["cfg_note"]
    # the sentence the operator can act on, not just "the engine exited 1:"
    assert "not valid YAML" in v["cfg_note"] or "indented" in v["cfg_note"], \
        v["cfg_note"]
    assert "not valid YAML" in console.render_project(v, contract)


def test_one_line_keeps_the_engines_stderr_and_drops_the_command_echo(console):
    e = console.Refuse("the engine exited 1:\n"
                       "  /usr/bin/python3 metaannot.py describe --json\n"
                       "config.yaml: not valid YAML at line 3, column 1\n"
                       "  expected ',' or ']', but got '<stream end>'")
    got = console.one_line(e)
    assert "the engine exited 1:" in got
    assert "not valid YAML at line 3" in got
    assert "expected ',' or ']'" in got
    assert "metaannot.py describe" not in got


def test_a_state_file_that_is_still_corrupt_is_not_retried_again(console,
                                                                 results,
                                                                 contract,
                                                                 monkeypatch):
    """read_json_stable short-circuited only on ok and on missing-when-not-
    retrying, so eight directories with a corrupt state file put most of a
    second into every index poll, for a file that will still be corrupt."""
    results.raw_state("{ this is not json")
    project = console.Project(0, results.path, contract)
    tries = []
    real = console.read_json_once
    monkeypatch.setattr(console, "read_json_once",
                        lambda p, l: (tries.append(p), real(p, l))[1])
    project.read_state()
    assert len(tries) == console.STATE_RETRIES     # a torn write deserves this
    tries.clear()
    for _ in range(5):
        project.read_state()
    # Not read AT ALL, not merely read once: the short-circuit shortened the
    # retry loop and left the read, so an 8 MB state file already known corrupt
    # was still read whole on every poll of every affected project.
    assert tries == [], "an 8 MB corrupt state file is still read every poll"

    # and a file that is rewritten is given the full budget again
    tries.clear()
    results.raw_state("{ still not json, but new")
    project.read_state()
    assert len(tries) == console.STATE_RETRIES


def test_one_project_forks_describe_once_however_many_tabs_open_at_once(
        console, results, contract):
    """config_run released the lock before calling run_describe, so N
    concurrent first-requests forked N identical `describe --config`."""
    results.state(running_state(contract, age=10))
    results.write(contract.effective_name,
                  "proteins_faa: /nowhere/p.faa\nquant_table: /nowhere/q.tsv\n"
                  "run:\n  pfam: false\n")
    project = console.Project(0, results.path, contract)
    forks = []
    real = console.run_describe

    def slow(*a, **kw):
        forks.append(1)
        time.sleep(0.3)
        return real(*a, **kw)

    console.run_describe = slow
    try:
        state = project.read_state()["state"]
        threads = [threading.Thread(
            target=lambda: project.config_run(sys.executable, METAANNOT_PY,
                                              state)) for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(60)
    finally:
        console.run_describe = real
    assert len(forks) == 1, "%d forks for one project" % len(forks)


# ======================================================================
# round 3: the read-only guarantee, the guard on untrusted numbers, the
# fork that wedged the pool, the pinned block, and the stalled run the
# index could not show
# ======================================================================

# ----------------------------------------------------------------------
# 1. the socket never lives in, and never creates, a watched directory
# ----------------------------------------------------------------------

@pytest.fixture
def loose_tmp():
    """A short-pathed scratch directory, and its contents afterwards."""
    d = tempfile.mkdtemp(prefix="mcs", dir="/tmp")
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_a_socket_inside_a_watched_directory_is_refused(console, loose_tmp,
                                                        contract):
    """The reviewer\'s reproduction, verbatim.

    `console.py /tmp/mcres/results --socket /tmp/mcres/results/console.sock`
    printed "metaannot console 0.1.0 - read-only. Nothing is written to any
    results directory." while recording os.mkdir of the results directory, an
    O_CREAT|O_RDWR open of console.sock.lock inside it and a bind of the socket
    beside that. The directory\'s mtime, ctime, nlink and size all moved and
    the lock file survived a clean SIGTERM.
    """
    res = os.path.join(loose_tmp, "results")
    os.mkdir(res, 0o700)
    with open(os.path.join(res, contract.state_name), "w") as fh:
        fh.write("{}")
    before = sorted(os.listdir(res)), os.stat(res).st_mtime_ns
    args = console.build_parser().parse_args(
        [res, "--socket", os.path.join(res, "console.sock"),
         "--metaannot", METAANNOT_PY])
    with pytest.raises(console.Refuse) as e:
        console.socket_path(args)
    assert "is inside" in str(e.value) and res in str(e.value)
    assert "console.sock.lock" in str(e.value)
    assert (sorted(os.listdir(res)), os.stat(res).st_mtime_ns) == before


def test_a_socket_inside_a_scanned_root_is_refused(console, loose_tmp):
    """--root is the same guarantee one level up: the console must not put its
    socket in a tree it is walking for results directories."""
    args = console.build_parser().parse_args(
        ["--root", loose_tmp, "--socket", os.path.join(loose_tmp, "c.sock"),
         "--metaannot", METAANNOT_PY])
    with pytest.raises(console.Refuse) as e:
        console.socket_path(args)
    assert "scanning for results directories" in str(e.value)


def test_a_socket_beside_a_watched_directory_is_still_allowed(console,
                                                              loose_tmp,
                                                              contract):
    """The refusal is about being INSIDE one, not about being near one: a 0700
    directory next to the watched tree is a perfectly good place to bind."""
    res = os.path.join(loose_tmp, "results")
    os.mkdir(res, 0o700)
    with open(os.path.join(res, contract.state_name), "w") as fh:
        fh.write("{}")
    home = os.path.join(loose_tmp, "sockets")
    os.mkdir(home, 0o700)
    args = console.build_parser().parse_args(
        [res, "--socket", os.path.join(home, "c.sock"),
         "--metaannot", METAANNOT_PY])
    assert console.socket_path(args) == os.path.join(home, "c.sock")


def test_an_explicit_socket_never_creates_its_own_parent(console, loose_tmp):
    """The more reachable half of the same defect, and the worse one.

    runtime_dir\'s os.makedirs made the parent, so `--socket
    <root>/DatasetZ/results/c.sock` produced a results directory at 0700 with
    c.sock and c.sock.lock in it - a results directory the CONSOLE made, which
    a run started into that path afterwards then finds. No unusual mode on
    anything and no unusual umask required.
    """
    target = os.path.join(loose_tmp, "DatasetZ", "results")
    args = console.build_parser().parse_args(
        ["--socket", os.path.join(target, "c.sock"), "--metaannot",
         METAANNOT_PY])
    with pytest.raises(console.Refuse) as e:
        console.socket_path(args)
    assert "does not exist" in str(e.value)
    assert not os.path.exists(os.path.join(loose_tmp, "DatasetZ"))
    assert os.listdir(loose_tmp) == []


def test_nothing_is_created_for_an_explicit_socket_even_under_umask_zero(
        console, loose_tmp):
    """os.makedirs applies its mode to the LEAF only, so under umask 0
    `--socket /tmp/a/b/m.sock` left /tmp/a at 0777 with only b at 0700 - and
    /tmp/a is not sticky, so another account can rename b away and substitute
    its own directory, after which `ssh -L` reaches their socket."""
    old = os.umask(0)
    try:
        args = console.build_parser().parse_args(
            ["--socket", os.path.join(loose_tmp, "a", "b", "m.sock"),
             "--metaannot", METAANNOT_PY])
        with pytest.raises(console.Refuse):
            console.socket_path(args)
    finally:
        os.umask(old)
    assert os.listdir(loose_tmp) == []


def test_an_exposed_ancestor_is_refused_not_only_the_immediate_parent(
        console, loose_tmp):
    """Only the socket\'s own directory was ever lstat-checked. A world-
    writable directory ANYWHERE above it is the same hole: the attacker
    replaces the 0700 directory rather than the socket."""
    wide = os.path.join(loose_tmp, "wide")
    tight = os.path.join(wide, "tight")
    os.mkdir(wide, 0o777)
    os.chmod(wide, 0o777)                     # umask does not get a vote
    os.mkdir(tight, 0o700)
    args = console.build_parser().parse_args(
        ["--socket", os.path.join(tight, "m.sock"), "--metaannot",
         METAANNOT_PY])
    with pytest.raises(console.Refuse) as e:
        console.socket_path(args)
    assert "sits under" in str(e.value) and wide in str(e.value)
    assert "not sticky" in str(e.value)


def test_a_group_writable_ancestor_is_reported_and_not_refused(console,
                                                               loose_tmp,
                                                               capsys):
    """A home directory at 0775 is the default on any Linux with umask 002 and
    per-user groups, where that group has exactly one member. Refusing there
    would send a diligent operator round the same loop the old "Try: mkdir -m
    700 -p /tmp" message sent them round, so it is reported and not judged -
    this console cannot count the group\'s members."""
    wide = os.path.join(loose_tmp, "shared")
    tight = os.path.join(wide, "tight")
    os.mkdir(wide, 0o770)
    os.chmod(wide, 0o770)
    os.mkdir(tight, 0o700)
    args = console.build_parser().parse_args(
        ["--socket", os.path.join(tight, "m.sock"), "--metaannot",
         METAANNOT_PY])
    assert console.socket_path(args) == os.path.join(tight, "m.sock")
    err = capsys.readouterr().err
    assert "is 0700, but" in err and wide in err
    assert "cannot tell how many members" in err


def test_a_sticky_ancestor_is_accepted_because_sticky_is_the_protection(
        console, loose_tmp):
    """/tmp is 1777 and every short_sock in this file lives under it. The rule
    is not "no other account may write here"; it is "no other account may
    replace what is here", which is what the sticky bit means."""
    tight = os.path.join(loose_tmp, "tight")
    os.mkdir(tight, 0o700)
    args = console.build_parser().parse_args(
        ["--socket", os.path.join(tight, "m.sock"), "--metaannot",
         METAANNOT_PY])
    assert console.socket_path(args) == os.path.join(tight, "m.sock")


def test_the_refusal_never_tells_the_operator_to_widen_anything(console,
                                                                loose_tmp,
                                                                contract):
    """"Try: mkdir -m 700 -p /tmp" was advice the operator must not take, and
    a silent no-op on a directory that exists, so a diligent one runs it, sees
    nothing happen, and retries into the same refusal. Pointed at a 0755
    results directory it read "mkdir -m 700 -p <results>" - the console
    explaining how to make the violation it had just refused."""
    res = os.path.join(loose_tmp, "results")
    os.mkdir(res, 0o755)
    with open(os.path.join(res, contract.state_name), "w") as fh:
        fh.write("{}")
    for argv in ([res, "--socket", os.path.join(res, "c.sock")],
                 ["--socket", "/tmp/pwn.sock"]):
        args = console.build_parser().parse_args(
            argv + ["--metaannot", METAANNOT_PY])
        with pytest.raises(console.Refuse) as e:
            console.socket_path(args)
        msg = str(e.value)
        assert "mkdir -m 700" not in msg
        # the only mkdir it offers is of a directory of the console\'s own
        for line in msg.splitlines():
            if "mkdir" in line:
                assert "metaannot-console" in line, line
        assert "chmod 700 ~/.cache/metaannot-console" in msg


def test_even_the_default_socket_path_is_refused_inside_a_watched_directory(
        console, loose_tmp, contract, monkeypatch):
    """$XDG_RUNTIME_DIR is somebody\'s environment, and a results directory at
    0700 satisfies every test runtime_dir applies. The guarantee is not
    "unless you asked for it"."""
    res = os.path.join(loose_tmp, "results")
    os.mkdir(res, 0o700)
    with open(os.path.join(res, contract.state_name), "w") as fh:
        fh.write("{}")
    monkeypatch.setenv("XDG_RUNTIME_DIR", res)
    args = console.build_parser().parse_args(
        [res, "--metaannot", METAANNOT_PY])
    assert not args.socket
    with pytest.raises(console.Refuse) as e:
        console.socket_path(args)
    assert "is inside" in str(e.value)
    assert sorted(os.listdir(res)) == [contract.state_name]


def test_a_socket_owned_by_another_account_is_never_unlinked(console,
                                                             short_sock,
                                                             monkeypatch):
    """take_socket_name unlinked a stale socket on S_ISSOCK alone, never on who
    owns it."""
    srv = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    srv.bind(short_sock)
    srv.close()                                 # a socket inode, nothing on it
    monkeypatch.setattr(console.os, "getuid", lambda: 999999)
    with pytest.raises(console.Refuse) as e:
        console.take_socket_name(short_sock)
    assert "not by this account" in str(e.value)
    assert os.path.exists(short_sock)


# ----------------------------------------------------------------------
# 2. one guard, for every number that comes out of a file
# ----------------------------------------------------------------------

@pytest.mark.parametrize("val", [
    float("inf"), float("-inf"), float("nan"), 1e18, -1e18,
    pytest.param(10 ** 400, id="an int too large to be a float"),
    "12", None, True, False, [1], {"a": 1},
])
def test_no_unusable_number_survives_the_guard(console, val):
    assert console.epoch_of(val) is None


@pytest.mark.parametrize("val", [0, 1.5, -3, 1e9])
def test_a_usable_number_comes_back_unchanged(console, val):
    assert console.finite(val) == float(val)


@pytest.mark.parametrize("stamp_text", [
    "1899-01-01T00:00:00", "0001-01-01T00:00:00", "9999-12-31T23:59:59",
    "not a stamp", "", None, 42,
])
def test_a_stamp_outside_the_range_comes_back_none_not_an_exception(
        console, stamp_text):
    assert console.stamp_epoch(stamp_text) is None


@pytest.mark.parametrize("bad", [float("inf"), float("nan"), "?", None])
def test_the_formatters_answer_a_question_mark_rather_than_raising(console,
                                                                   bad):
    assert console.dur(bad) == "?"
    assert console.clock(bad) == "?"
    assert console.human_bytes(bad) == "?"


def test_an_epoch_out_of_range_is_a_question_mark_not_an_oserror(console):
    """`last_seen_epoch: 1e18` was OSError(22) out of time.localtime. 1e18 is
    a fine number of BYTES, which is why the bound belongs to the reader of
    epochs and not to the guard."""
    assert console.clock(1e18) == "?"
    assert console.dur(1e18) == "?"          # over the hundred-year ceiling
    assert console.human_bytes(1e18) != "?"


def state_with(contract, run_extra=None, stage=None):
    run = {"run_id": "r1", "version": "0.3.0", "host": "lab-fedora", "pid": 7,
           "started": "2026-01-01T00:00:00", "heartbeat_s": 30,
           "last_seen_epoch": time.time(), "final_status": "running"}
    run.update(run_extra or {})
    state = {contract.run_key: run}
    state.update(stage or {})
    return state


@pytest.mark.parametrize("label,text", [
    ("a run stamp before the platform range",
     '{"%(run)s": {"final_status": "running", "started": '
     '"1899-01-01T00:00:00", "heartbeat_s": 30, "last_seen_epoch": %(now)f}}'),
    ("an epoch too large for time_t",
     '{"%(run)s": {"final_status": "running", "started": '
     '"2026-01-01T00:00:00", "heartbeat_s": 30, "last_seen_epoch": 1e18}}'),
    ("Infinity, which json accepts",
     '{"%(run)s": {"final_status": "running", "started": '
     '"2026-01-01T00:00:00", "heartbeat_s": 30, "last_seen_epoch": Infinity}}'),
    ("-Infinity",
     '{"%(run)s": {"final_status": "running", "started": '
     '"2026-01-01T00:00:00", "heartbeat_s": 30, "last_seen_epoch": -Infinity}}'),
    ("NaN",
     '{"%(run)s": {"final_status": "running", "started": '
     '"2026-01-01T00:00:00", "heartbeat_s": 30, "last_seen_epoch": NaN}}'),
    ("an infinite heartbeat interval",
     '{"%(run)s": {"final_status": "running", "started": '
     '"2026-01-01T00:00:00", "heartbeat_s": Infinity, '
     '"last_seen_epoch": %(now)f}}'),
    ("a stage that took Infinity seconds",
     '{"%(run)s": {"final_status": "running", "started": '
     '"2026-01-01T00:00:00", "heartbeat_s": 30, "last_seen_epoch": %(now)f}, '
     '"pfam": {"status": "ok", "seconds": Infinity, "finished": '
     '"2026-01-01T01:00:00"}}'),
    ("a stage that finished in year 1",
     '{"%(run)s": {"final_status": "running", "started": '
     '"2026-01-01T00:00:00", "heartbeat_s": 30, "last_seen_epoch": %(now)f}, '
     '"pfam": {"status": "failed", "error": "x", "finished": '
     '"0001-01-01T00:00:00"}}'),
])
def test_an_exotic_number_never_500s_a_page(console, results, contract, label,
                                            text):
    """Every one of these was confirmed live over the socket. Python\'s json
    ACCEPTS NaN and Infinity, so the "not valid JSON" gate never saw them, and
    one odd file 500\'d / and /api/projects for EVERY watched project."""
    results.raw_state(text % {"run": contract.run_key, "now": time.time()})
    results.log("[    0.0s] INFO  x\n")
    project = console.Project(0, results.path, contract)
    v = console.project_view(project, sys.executable, METAANNOT_PY, log=False)
    assert console.render_project(v, contract)
    iv = console.index_view([project], sys.executable, METAANNOT_PY)
    assert console.render_index(iv, contract)
    assert "Traceback" not in console.render_index(iv, contract)


def test_a_deeply_nested_state_file_is_refused_rather_than_500ing(console,
                                                                  results,
                                                                  contract):
    """json.loads raises RecursionError, which is not a ValueError, so a 40 KB
    file of `{"a":` + `[`*20000 + `]`*20000 + `}` 500\'d the index,
    /api/projects and its own project page for every dataset."""
    results.raw_state('{"a":' + "[" * 20000 + "]" * 20000 + "}")
    project = console.Project(0, results.path, contract)
    snap = project.read_state()
    assert snap["status"] == "bad"
    assert "nested" in snap["detail"]
    v = console.project_view(project, sys.executable, METAANNOT_PY, log=False)
    assert console.render_project(v, contract)
    iv = console.index_view([project], sys.executable, METAANNOT_PY)
    assert iv["projects"][0]["bucket"] == "unreadable"


def test_the_index_survives_one_bad_directory_among_good_ones(console,
                                                              tmp_path,
                                                              contract):
    """The point of the finding: the page 500s for EVERY project, not just the
    one with the odd file."""
    paths = []
    for name, body in (("good", None), ("odd", '{"%s": {"final_status": '
                                               '"running", "started": '
                                               '"1899-01-01T00:00:00"}}'
                                               % contract.run_key)):
        d = tmp_path / name / "results"
        d.mkdir(parents=True)
        with open(str(d / contract.state_name), "w") as fh:
            fh.write(body if body is not None
                     else json.dumps(state_with(contract)))
        paths.append(str(d))
    projects = [console.Project(i, p, contract) for i, p in enumerate(paths)]
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    html = console.render_index(view, contract)
    assert "good" in html and "odd" in html


# ----------------------------------------------------------------------
# 3. describe --config never blocks the request pool again
# ----------------------------------------------------------------------

def test_a_slow_describe_config_blocks_one_thread_briefly_and_no_others(
        console, results, contract, monkeypatch):
    """A 2.4 MB config.effective.yaml made one /p/3 take 61.6 s while holding a
    MAX_CONNS permit; 33 concurrent requests then 503\'d /app.css and every
    healthy project. S_ISREG does not close that - a regular file can be as
    slow as it likes - so the fork runs off the request thread and exactly one
    request waits CONFIG_WAIT for it."""
    results.state(running_state(contract, age=10))
    results.write(contract.effective_name,
                  "proteins_faa: /nowhere/p.faa\nquant_table: /nowhere/q.tsv\n"
                  "run:\n  pfam: false\n")
    monkeypatch.setattr(console, "CONFIG_WAIT", 0.2)
    real = console.run_describe

    def slow(*a, **kw):
        if kw.get("config"):
            time.sleep(2.0)
        return real(*a, **kw)

    monkeypatch.setattr(console, "run_describe", slow)
    project = console.Project(0, results.path, contract)
    state = project.read_state()["state"]

    started = time.time()
    run, note = project.config_run(sys.executable, METAANNOT_PY, state)
    first = time.time() - started
    assert run is None and "still being read" in note
    assert first < 1.5, "one request thread waited %.1fs" % first

    times = []
    for _ in range(5):
        at = time.time()
        project.config_run(sys.executable, METAANNOT_PY, state)
        times.append(time.time() - at)
    assert max(times) < 0.2, times      # nobody else waits at all

    for _ in range(100):
        run, note = project.config_run(sys.executable, METAANNOT_PY, state)
        if run is not None:
            break
        time.sleep(0.1)
    assert run == {"pfam": False} or run.get("pfam") is False
    assert "read from" in note


def test_an_enormous_config_is_never_handed_to_the_engine(console, results,
                                                          contract,
                                                          monkeypatch):
    """A courtesy on top of the real fix: no config this engine writes is a
    megabyte, and parsing one costs the machine the run is on."""
    results.state(running_state(contract, age=10))
    results.write(contract.effective_name,
                  "run:\n  pfam: true\n" + "# padding\n" * 200000)
    forks = []
    monkeypatch.setattr(console, "run_describe",
                        lambda *a, **kw: forks.append(kw.get("config")))
    project = console.Project(0, results.path, contract)
    run, note = project.config_run(sys.executable, METAANNOT_PY,
                                   project.read_state()["state"])
    assert run is None and forks == []
    assert "MB" in note and "not a config this engine wrote" in note


# ----------------------------------------------------------------------
# 4. the pinned block, which must not push the page off the screen
# ----------------------------------------------------------------------

def many_failures(contract, n, chars=417, ago=7200):
    """The reviewer\'s measured shape: one wrong `db.*` root fails pfam, dbcan,
    ncbifam, kofam and jackhmmer in the same run, and metaannot caps a stage
    error at 500 characters."""
    now = time.time()
    state = running_state(contract, age=21, started_ago=ago + 3600)
    for i, name in enumerate(contract.stage_names[:n]):
        state[name] = {"signature": None, "status": "failed",
                       "error": "CUDA out of memory. " + "x" * chars,
                       "finished": stamp(now - ago - i)}
    return state


def test_the_pinned_block_is_bounded_however_many_stages_fail(console,
                                                              results,
                                                              contract):
    """Measured at 1440x860: five failures made the vitals block 799 px and
    eight made it 1123 px, with the heartbeat line 245 px below the fold and
    the whole stage table under that. The fix for "the failure is below the
    fold" must not push everything else below the fold."""
    results.state(many_failures(contract, 8))
    v = view_of(console, results, contract, log=False)
    assert len(v["failures"]) == 8
    html = console.render_project(v, contract)
    block = html.split('<div class="dead">', 1)[1].split("</div>", 1)[0]
    assert len(block) < 2000, len(block)
    assert block.count("failed —") == 3
    assert "and 5 more, in the table below" in html
    # the error is truncated in the block and whole in the row below it
    assert "x" * 417 not in block
    assert "more characters)" in block
    # and the pinned block itself stays short - by what goes in it, not by a
    # cap: a `max-height` on a `position:sticky` block put its own tail in a
    # scroller the page could not reach. See
    # test_the_pinned_block_is_never_a_nested_scroll_container.
    assert len(visible_text(pinned_block(html)[1])) < PIN_TEXT_MAX


def test_a_four_megabyte_error_string_does_not_become_the_page(console,
                                                               results,
                                                               contract):
    """render_failures printed the error verbatim and render_row printed it
    again as `detail`: a 4 MB state file produced an 8.03 MB fragment, which
    the page re-fetches every three seconds."""
    state = running_state(contract, age=10)
    state["pfam"] = {"signature": None, "status": "failed",
                     "error": "z" * (4 * 1024 * 1024),
                     "finished": stamp(time.time() - 60)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    frag = console.render_project_body(console.decorate(v, contract))
    assert len(frag) < 100 * 1024, len(frag)
    assert "z" * 1000 not in frag


def test_a_failure_that_predates_the_run_is_not_announced_as_this_runs(
        console, results, contract):
    """Verified live on a reused results directory: a run started an hour ago,
    and the top of the page reading "1 stage has failed in this run - esmfold
    failed - 2026-09-01T23:41:21 (8d 00h ago)" with the run\'s own start time
    on the very next line. The state file is cumulative and adopted reuse is a
    headline feature, so a reused directory is metaannot\'s normal mode."""
    now = time.time()
    state = running_state(contract, age=20, started_ago=3600)
    state["esmfold"] = {"signature": None, "status": "failed",
                        "error": "CUDA out of memory",
                        "finished": stamp(now - 8 * 86400)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert v["failures"] == []
    html = console.render_project(v, contract)
    assert "failed in this run" not in html
    # still visible, still red, still explained - as an earlier run\'s
    assert "from an earlier run" in html
    assert "CUDA out of memory" in html
    # and the index says the same thing
    iv = console.index_view([console.Project(0, results.path, contract)],
                            sys.executable, METAANNOT_PY)
    assert iv["projects"][0]["failed"] == []
    assert iv["projects"][0]["carried_failed"] == ["esmfold"]
    assert "failed in an earlier run" in console.render_index(iv, contract)


def test_a_failure_inside_this_run_is_still_announced(console, results,
                                                      contract):
    """The other half: the filter is on `carried`, not on age."""
    state = running_state(contract, age=20, started_ago=86400)
    state["esmfold"] = {"signature": None, "status": "failed",
                        "error": "CUDA out of memory",
                        "finished": stamp(time.time() - 7200)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert [f["name"] for f in v["failures"]] == ["esmfold"]
    assert "1 stage has failed in this run" in console.render_project(v,
                                                                     contract)


def test_the_failure_block_is_past_tense_over_a_run_that_has_ended(console,
                                                                   results,
                                                                   contract):
    """final_status failed, finished 26 h ago, and the page still said "The
    engine carries on with every stage that does not depend on a casualty".
    Nothing is carrying on; the run is over.

    What this pins is the SUPPRESSION, so it names the live sentence that
    actually ships - whatever that sentence currently is. Asserting the retired
    wording is absent stopped being a test of anything the moment the wording
    was retired: no branch of render_failures can emit it, so the assertion
    could not fail whichever branch the page took. The reintroduction guard
    lives with the block that replaced it.
    """
    now = time.time()
    state = running_state(contract, age=26 * 3600, started_ago=30 * 3600)
    state[contract.run_key].update(final_status="failed",
                                   finished=stamp(now - 26 * 3600))
    state["esmfold"] = {"signature": None, "status": "failed",
                        "error": "CUDA out of memory",
                        "finished": stamp(now - 27 * 3600)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    html = console.render_project(v, contract)
    assert "1 stage failed in this run" in html
    assert "has failed in this run" not in html
    assert "This run has ended" in html
    # the live half of the same ternary, which is the thing being suppressed
    assert "Nothing new starts after the first failure" not in html
    assert "still look busy" not in html
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase


def test_the_outcome_of_a_finished_run_leads_its_identifiers(console, results,
                                                             contract):
    """"Finished failed at 2026-09-08T21:54:42 - ran 4h00m" is the headline;
    `run 20260101T000000-1 - metaannot 0.3.0 - pid 4242` is a set of
    identifiers. They were the other way round."""
    now = time.time()
    state = running_state(contract, age=100, started_ago=4 * 3600 + 100)
    state[contract.run_key].update(final_status="failed",
                                   finished=stamp(now - 100))
    results.state(state)
    html = console.render_project(view_of(console, results, contract,
                                          log=False), contract)
    assert html.index("Finished failed at") < html.index("run 20260101T000000")


# ----------------------------------------------------------------------
# 5. the index could not show a stalled run - the P1 wound, one page over
# ----------------------------------------------------------------------

def two_runs(tmp_path, contract, ages):
    paths = []
    for name, age in ages:
        d = tmp_path / name / "results"
        d.mkdir(parents=True)
        with open(str(d / contract.state_name), "w") as fh:
            json.dump(running_state(contract, age=age), fh)
        with open(str(d / contract.log_name), "w") as fh:
            fh.write("[    0.0s] INFO  x\n")
        paths.append(str(d))
    return [console_project(contract, i, p) for i, p in enumerate(paths)]


def console_project(contract, index, path):
    import metaannot_console
    return metaannot_console.Project(index, path, contract)


def test_the_index_shows_a_stalled_run_without_calling_it_dead(console,
                                                               tmp_path,
                                                               contract):
    """Cohort_D, heartbeat 5h18m dead, was badged RUNNING with an EMPTY
    annotation cell, byte-identical markup to Pilot_A at 13 s, and the summary
    - which summarise()\'s docstring calls "the tmux replacement, in one
    sentence" - said "10 running"."""
    projects = two_runs(tmp_path, contract,
                        [("Pilot_A", 13), ("Cohort_D", 5 * 3600 + 18 * 60)])
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    by = {i["name"]: i for i in view["projects"]}
    assert by["Cohort_D"]["beat_band"] == "long"
    assert by["Pilot_A"]["beat_band"] == "fresh"
    # still RUNNING: this console does not decide that a run is over
    assert by["Cohort_D"]["bucket"] == "running"
    html = console.render_index(view, contract)
    assert "no heartbeat for 5h18m" in html
    assert "2 running" in view["summary"]
    assert "Cohort_D has not stamped a heartbeat for 5h18m" in view["summary"]
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase


def test_the_stalest_running_project_sorts_first_not_last(console, tmp_path,
                                                          contract):
    """The sort key was (bucket, -moved), so the one project in the group that
    had silently stopped sorted to the BOTTOM of it."""
    projects = two_runs(tmp_path, contract,
                        [("fresh_one", 5), ("stalled", 6 * 3600),
                         ("also_fresh", 9)])
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    assert [i["name"] for i in view["projects"]][0] == "stalled"


def test_a_healthy_index_says_nothing_about_heartbeats(console, tmp_path,
                                                       contract):
    projects = two_runs(tmp_path, contract, [("a", 5), ("b", 9)])
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    assert "heartbeat" not in view["summary"]
    assert "no heartbeat for" not in console.render_index(view, contract)


def test_the_summary_of_many_hurt_projects_is_still_one_sentence(console):
    assert console.name_list(["a", "b"]) == "a, b"
    assert console.name_list(["a", "b", "c", "d", "e"]) == "a, b, c and 2 more"


def test_the_poller_is_pinned_at_the_floor_while_anything_is_wrong(console,
                                                                   results,
                                                                   contract):
    """The documented back-off could not fire on a healthy run - the engine
    rewrites the state file for every 30 s heartbeat, so data-sig flips every
    30 s - and could fire only once the heartbeat STOPPED, which is the moment
    you are watching for it to come back."""
    results.state(running_state(contract, age=5))
    ok = console.render_project(view_of(console, results, contract, log=False),
                                contract)
    assert 'data-alert="0"' in ok
    results.state(running_state(contract, age=4 * 3600))
    late = console.render_project(view_of(console, results, contract,
                                          log=False), contract)
    assert 'data-alert="1"' in late
    assert "if (alarm) return floor;" in console.PAGE_JS


# ----------------------------------------------------------------------
# 6. the small ones, each with the sentence that named it
# ----------------------------------------------------------------------

def test_a_log_line_cannot_forge_a_second_line(console):
    """str.splitlines() splits on \\r, \\x0b, \\x0c, \\x1c-\\x1e, \\x85, U+2028 and
    U+2029, so one tool-emitted line containing U+2028 rendered as two, the
    second parsed as a well-formed engine line the engine never wrote."""
    forged = ("[    1.0s] INFO   emapper | eggnog said \u2028"
              "[    2.0s] INFO  all stages finished ok\n")
    lines = console.parse_log(forged)
    assert len(lines) == 1
    assert lines[0]["text"].count("all stages finished ok") == 1
    assert console.log_lines(forged)[0] == lines


def test_control_characters_never_reach_the_served_page(console, results,
                                                        contract):
    """NUL bytes and ANSI escape sequences from the log were written raw into
    the served HTML body - verified present in the response bytes. Escaping is
    correct for markup and does nothing about what is underneath it."""
    results.state(running_state(contract, age=5))
    results.log("[    1.0s] FATAL  \x1b[31mred\x1b[0m and a \x00 byte\n")
    html = console.render_project(view_of(console, results, contract), contract)
    assert "\x1b" not in html and "\x00" not in html
    assert "\ufffd" in html                       # shown, as a replacement
    assert "red" in html


def test_a_lock_pid_is_printed_as_the_number_the_command_uses(console,
                                                              results,
                                                              contract):
    """valid_pid accepts "0021877" and builds `ps -p 21877`, while the prose
    one clause earlier read "Held by pid 0021877" - two different numbers on
    one line, in the block whose whole purpose is handing over an exact
    command."""
    results.lock({"pid": "0021877", "host": "lab-fedora",
                  "started": "2026-09-07T22:10:41"})
    results.state(running_state(contract, age=20))
    html = console.render_project(view_of(console, results, contract),
                                  contract)
    assert "Held by pid 21877 on lab-fedora" in html
    assert "ps -p 21877" in html
    assert "0021877" not in html


def test_a_corrupt_state_file_does_not_produce_twenty_one_positive_rows(
        console, results, contract):
    """The page led with "Nothing below is known about this directory" and then
    made 21 positive claims below it - every row reading "NEXT - no record yet
    - nothing is blocking it", byte-identical to a genuinely never-started
    directory, in the element that occupies most of the page."""
    results.raw_state("{ this is not json")
    v = view_of(console, results, contract, log=False)
    assert v["state_status"] == "bad"
    assert v["state_known"] is False
    assert {r["state"] for r in v["rows"]} == {"none"}
    html = console.render_project(v, contract)
    assert "no dependency is blocking it" not in html
    assert "no record yet" not in html
    assert "the state file cannot be read, so whether this stage has run" in html
    # and a directory where nothing has run still reads the old way
    results.write(contract.state_name, "{}")
    fresh = view_of(console, results, contract, log=False)
    assert "no dependency is blocking it" in console.render_project(fresh,
                                                                   contract)


def test_an_output_missing_from_a_finished_stage_is_not_still_coming(console,
                                                                     results,
                                                                     contract):
    """An OK row saying "finished 2026-09-09T16:52:45" and then "not written
    yet" about its own output is a flat self-contradiction on one line, and it
    reads as "it is coming"."""
    named = next(s for s in contract.stages if s["outputs"])
    state = running_state(contract, age=5)
    state[named["name"]] = {"signature": "a", "status": "ok", "seconds": 12.0,
                            "finished": "2026-09-09T16:52:45"}
    results.state(state)
    html = console.render_project(view_of(console, results, contract,
                                          log=False), contract)
    assert "declared by this stage, but not in this directory" in html
    row = html.split(named["outputs"][0], 1)[1][:200]
    assert "not written yet" not in row


def test_a_results_directory_this_account_cannot_read_is_still_discovered(
        console, tmp_path, contract):
    """os.path.exists() is False on EACCES, so a group-unreadable results
    directory simply vanished from a --root scan with nothing said: eleven
    fixture directories, one at mode 000, and a banner reading "watching 10
    directories"."""
    for name in ("open", "shut"):
        d = tmp_path / name / "results"
        d.mkdir(parents=True)
        with open(str(d / contract.state_name), "w") as fh:
            fh.write("{}")
    shut = str(tmp_path / "shut" / "results")
    os.chmod(shut, 0o000)
    try:
        found = console.discover([str(tmp_path)], contract)
        assert sorted(found) == sorted([str(tmp_path / "open" / "results"),
                                        shut])
    finally:
        os.chmod(shut, 0o755)


def test_a_slow_rescan_parks_one_thread_not_one_every_interval(console,
                                                               contract,
                                                               tmp_path,
                                                               monkeypatch):
    """rescan() stamped _last_scan BEFORE running discover() outside the lock,
    so against a --root on an unresponsive mount one request thread was parked
    permanently every 30 s: threads 3 to 8 in 160 s, MAX_CONNS in about sixteen
    minutes, and never recovering."""
    con = console.Console(contract, [], sys.executable, METAANNOT_PY,
                          roots=[str(tmp_path)], rescan_s=0.0)
    calls, gate = [], threading.Event()

    def slow(*a, **kw):
        calls.append(1)
        gate.wait(5)
        return []

    monkeypatch.setattr(console, "discover", slow)
    threads = [threading.Thread(target=con.rescan, daemon=True)
               for _ in range(6)]
    for t in threads:
        t.start()
    time.sleep(0.4)
    assert calls == [1], "%d threads are inside discover()" % len(calls)
    gate.set()
    for t in threads:
        t.join(10)


def test_a_directory_written_on_another_host_is_named(console, results,
                                                      contract):
    """The design\'s central premise is that the console runs on the host that
    writes the results directory, and nothing said a word when it did not -
    even though describe --json supplies this console\'s own host."""
    results.state(running_state(contract, age=5))          # host lab-fedora
    html = console.render_project(view_of(console, results, contract,
                                          log=False), contract)
    assert "This directory records host lab-fedora" in html
    assert contract.host in html


# ----------------------------------------------------------------------
# the final round: six reproductions, one test each
# ----------------------------------------------------------------------

# 1. the default socket path, which creates before it refuses
# ----------------------------------------------------------------------

@pytest.mark.parametrize("how,what", [
    ("project", "watching"),
    ("root", "scanning for results directories"),
])
def test_the_default_socket_creates_nothing_under_a_watched_tree(
        console, loose_tmp, contract, monkeypatch, how, what):
    """An ordering bug, and the order is the whole guarantee.

    socket_path() called runtime_dir() BEFORE refuse_if_watched(), and
    runtime_dir -> ensure_private_dir -> os.mkdir(cand, 0o700). So with
    XDG_RUNTIME_DIR set to <results>/rt the console CREATED that directory
    inside the tree it is watching and only then refused to use it: the audit
    hook recorded the mkdir and the directory's own snapshot moved. Refusing
    after creating is not refusing.
    """
    res = os.path.join(loose_tmp, "results")
    os.mkdir(res, 0o700)
    with open(os.path.join(res, contract.state_name), "w") as fh:
        fh.write("{}")
    rt = os.path.join(res, "rt")            # not there, and must stay not there
    monkeypatch.setenv("XDG_RUNTIME_DIR", rt)
    before = (sorted(os.listdir(res)), os.stat(res).st_mtime_ns,
              sorted(os.listdir(loose_tmp)))
    argv = [res] if how == "project" else ["--root", loose_tmp]
    args = console.build_parser().parse_args(
        argv + ["--metaannot", METAANNOT_PY])
    with pytest.raises(console.Refuse) as e:
        console.socket_path(args)
    assert "is inside" in str(e.value) and what in str(e.value)
    assert not os.path.exists(rt)
    assert (sorted(os.listdir(res)), os.stat(res).st_mtime_ns,
            sorted(os.listdir(loose_tmp))) == before


# 2. one bad string, and one bad project, on a page about eight
# ----------------------------------------------------------------------

@pytest.mark.parametrize("path", ["/tmp/with\x00nul", "/tmp/lone\ud800surrogate"])
def test_a_path_the_system_will_not_accept_is_a_reason_not_a_raise(console,
                                                                   path):
    """os.stat raises ValueError("embedded null byte") for a NUL and
    UnicodeEncodeError - a ValueError subclass - for a lone surrogate, and
    neither is an OSError. stat_or_reason() caught OSError only."""
    st, why = console.stat_or_reason(path)
    assert st is None
    assert why and why != "missing"
    assert console.stat_or_none(path) is None
    assert console.regular_stat(path) is None


def test_one_projects_bad_config_path_does_not_500_the_index_for_the_others(
        console, tmp_path, contract):
    """`_run.config_path` is a string straight out of somebody's state file and
    reaches os.stat through regular_stat(). One NUL in it raised out of
    index_view, which has no per-project guard, so `/` and /api/projects were
    500 for every watched project - seven of them perfectly readable."""
    paths = []
    for i, cfg in enumerate(("/nowhere/ok.yaml", "/nowhere/b\x00d.yaml",
                             "/nowhere/also-ok.yaml")):
        d = tmp_path / ("ds%d" % i) / "results"
        d.mkdir(parents=True)
        state = running_state(contract, age=10)
        state[contract.run_key]["config_path"] = cfg
        with open(str(d / contract.state_name), "w") as fh:
            json.dump(state, fh)
        with open(str(d / contract.log_name), "w") as fh:
            fh.write("[    0.0s] INFO  x\n")
        paths.append(str(d))
    projects = [console.Project(i, p, contract) for i, p in enumerate(paths)]
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    assert [i["bucket"] for i in view["projects"]] == ["running"] * 3
    html = console.render_index(view, contract)
    for name in ("ds0", "ds1", "ds2"):
        assert name in html
    # and the project itself renders, saying only that the config was not read
    v = console.project_view(projects[1], sys.executable, METAANNOT_PY,
                             log=False)
    assert v["cfg_run_known"] is False
    assert console.render_project(v, contract)


def test_one_project_that_raises_never_takes_the_index_down(console, tmp_path,
                                                            contract):
    """The structural half, which outlives the ValueError above: whatever the
    next unanticipated shape turns out to be, the front page of a watcher must
    not be all-or-nothing. One project raises; the other two are served."""
    paths = []
    for i in range(3):
        d = tmp_path / ("ds%d" % i) / "results"
        d.mkdir(parents=True)
        with open(str(d / contract.state_name), "w") as fh:
            json.dump(running_state(contract, age=10), fh)
        with open(str(d / contract.log_name), "w") as fh:
            fh.write("[    0.0s] INFO  x\n")
        paths.append(str(d))
    projects = [console.Project(i, p, contract) for i, p in enumerate(paths)]

    def boom():
        raise RuntimeError("this project explodes")

    projects[1].read_state = boom
    view = console.index_view(projects, sys.executable, METAANNOT_PY)
    assert len(view["projects"]) == 3
    by_name = {i["name"]: i for i in view["projects"]}
    assert by_name["ds0"]["bucket"] == "running"
    assert by_name["ds2"]["bucket"] == "running"
    bad = by_name["ds1"]
    assert bad["bucket"] == "error"
    assert "RuntimeError" in bad["state_detail"]
    # nothing is invented for it: no stages, no failures, no heartbeat
    assert bad["strip"] == [] and bad["failed"] == []
    assert bad["beat_band"] == "none"
    html = console.render_index(view, contract)
    for name in ("ds0", "ds1", "ds2"):
        assert name in html
    assert "failed while reading this directory" in html
    assert "The other rows on this page are unaffected" in html
    assert "Traceback" not in html


# 3. the pinned block, and the scroll trap the last round put in it
# ----------------------------------------------------------------------

def css_rules(css):
    """(selector, declarations) per rule, with comments stripped - the fix's
    own comment quotes the declarations it removed."""
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.S)
    out = []
    for chunk in css.split("}"):
        if "{" not in chunk:
            continue
        sel, decls = chunk.split("{", 1)
        out.append((sel.strip(), decls.strip()))
    return out


def rules_for(css, classes):
    """Every rule PAGE_CSS applies to an element carrying `classes`, for the
    single-class selectors this stylesheet is written in."""
    out = []
    for sel, decls in css_rules(css):
        for one in sel.split(","):
            one = one.strip()
            if re.fullmatch(r"\.[A-Za-z0-9_-]+", one) and one[1:] in classes:
                out.append((one, decls))
    return out


def element_at(html, start):
    """The whole <div> element that begins at `start`, matched by depth."""
    depth = 0
    for m in re.finditer(r"<(/?)div\b[^>]*>", html[start:]):
        depth += -1 if m.group(1) else 1
        if depth == 0:
            return html[start:start + m.end()]
    raise AssertionError("unclosed div")


def pinned_block(html):
    """(classes, html) of the sticky vitals element."""
    m = re.search(r'<div class="([^"]*\bsticky\b[^"]*)"', html)
    assert m, "nothing is pinned"
    return m.group(1).split(), element_at(html, m.start())


def visible_text(fragment):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", fragment)).strip()


# What the pinned block may spend, in characters of visible text. The cap it
# replaces was 45vh - 324 px at 1280x720 - and at the page's 1180 px content
# width that is about fifteen lines, so about 1600 characters. Eight failures
# with 500-character errors come to 1090. The point is not the exact number: it
# is that the pinned block is short because of what goes IN it, so it never
# needs a scrollbar of its own and never covers the table.
PIN_TEXT_MAX = 1600


def test_the_pinned_block_is_never_a_nested_scroll_container(console, results,
                                                             contract):
    """Measured at 1280x720: `.vitals` was `position:sticky; top:0` AND
    `max-height:45vh; overflow:auto`, which puts everything past the cap into a
    scroller the PAGE cannot scroll. `.dead` was 262 px for three failures
    inside a `.vitals` of clientHeight 322, and the heartbeat verdict was
    unreachable at every scroll position of the document.

    So: nothing that applies to the pinned element may clip it or scroll it,
    and the three urgent facts stay in it while everything else moves below it,
    in ordinary flow, where the page can scroll to it.
    """
    results.state(many_failures(contract, 8, chars=480))
    v = view_of(console, results, contract, log=False)
    assert len(v["failures"]) == 8
    html = console.render_project(v, contract)
    classes, pin = pinned_block(html)

    for sel, decls in rules_for(console.PAGE_CSS, set(classes)):
        assert "overflow" not in decls, (sel, decls)
        assert "max-height" not in decls, (sel, decls)

    # the failure summary and the heartbeat verdict are both IN the pin
    assert "8 stages have failed in this run" in pin
    assert "Heartbeat stamped" in pin
    # and it is short because of what is in it, not because of a cap
    assert len(visible_text(pin)) < PIN_TEXT_MAX, len(visible_text(pin))
    # everything that can grow is below it, in flow, and unclipped
    assert "run 20260101T000000" not in pin
    assert "No lock file" not in pin
    assert "The log last grew" not in pin
    for line in ("run 20260101T000000", "No lock file"):
        assert line in html
    # the stage table follows the pinned block and is in no scroller either
    assert html.index("<thead>") > html.index(pin[:60])
    table_panel = element_at(html, html.rindex('<div class="panel">',
                                               0, html.index("<thead>")))
    for sel, decls in rules_for(console.PAGE_CSS, {"panel"}):
        assert "overflow" not in decls and "max-height" not in decls
    assert "<thead>" in table_panel


# 4. the heartbeat cadence the record does not give
# ----------------------------------------------------------------------

MISSING = object()


@pytest.mark.parametrize("hb", [MISSING, None, float("inf"), "thirty", 0])
def test_a_heartbeat_interval_that_is_not_recorded_is_never_invented(
        console, results, contract, hb):
    """`finite(run.get("heartbeat_s"), ...) or 30.0`, and the page then said
    "This run stamps one every 30 s, so that is late" - the same sentence for
    all of these, and the record said no such thing in any of them. 0 is the
    engine's own value for a heartbeat deliberately turned OFF, which makes the
    invented cadence a claim about a run that stamps none.
    """
    state = running_state(contract, age=4 * 3600, heartbeat=hb)
    if hb is MISSING:
        del state[contract.run_key]["heartbeat_s"]
    results.state(state)
    v = view_of(console, results, contract, log=False)
    hbv = v["heartbeat"]
    assert hbv["band"] == "unrated"
    assert hbv["heartbeat_s"] is None
    assert "4h00m" in hbv["verdict"]              # the raw age, still reported
    html = console.render_project(v, contract)
    assert "every 30" not in html
    assert "so that is late" not in html
    assert "a long way past due" not in html
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase
    if hb == 0 and hb is not MISSING:
        assert "turns the heartbeat off" in hbv["verdict"]
    else:
        assert "does not say how often" in hbv["verdict"]
    # and the index does not invent one either
    iv = console.index_view([console.Project(0, results.path, contract)],
                            sys.executable, METAANNOT_PY)
    assert iv["projects"][0]["beat_band"] == "unrated"
    index = console.render_index(iv, contract)
    assert "every 30" not in index
    assert "does not call it late" in index


# 5. a failure this run cannot be shown to own
# ----------------------------------------------------------------------

@pytest.mark.parametrize("started", [None, "yesterday", "2026-09-02 01:09:28"])
def test_no_failure_is_this_runs_when_the_run_start_cannot_be_read(
        console, results, contract, started):
    """`carried` was set only when stamp_epoch(run["started"]) returned a
    number, and stamp_epoch accepts exactly %Y-%m-%dT%H:%M:%S. So a `_run`
    carrying no readable `started` made every historical failure this run's:
    "1 stage has failed in this run - esmfold failed - 2026-09-02T01:09:28
    (8d 00h ago)", three lines above the run's own record.
    """
    now = time.time()
    state = running_state(contract, age=20, started_ago=3600)
    if started is None:
        del state[contract.run_key]["started"]
    else:
        state[contract.run_key]["started"] = started
    state["esmfold"] = {"signature": None, "status": "failed",
                        "error": "CUDA out of memory",
                        "finished": stamp(now - 8 * 86400)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert v["failures"] == []
    assert [f["name"] for f in v["undated_failures"]] == ["esmfold"]
    row = next(r for r in v["rows"] if r["name"] == "esmfold")
    assert row["era"] == "unknown" and row["carried"] is None
    html = console.render_project(v, contract)
    assert "failed in this run" not in html
    assert "from an earlier run" not in html
    # still first, still red, still the engine's own error - just not dated
    assert "cannot tell whether in this run" in html
    assert "CUDA out of memory" in html
    assert "does not say when the run started" in html
    assert "which run this record belongs to is not established" in html
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase
    # and the index makes neither claim either
    iv = console.index_view([console.Project(0, results.path, contract)],
                            sys.executable, METAANNOT_PY)
    it = iv["projects"][0]
    assert it["failed"] == [] and it["carried_failed"] == []
    assert it["undated_failed"] == ["esmfold"]
    assert "does not establish which run" in console.render_index(iv, contract)


def test_a_readable_run_start_still_dates_its_failures(console, results,
                                                       contract):
    """The other half of the same comparison: nothing above weakens the case
    the console CAN decide."""
    now = time.time()
    state = running_state(contract, age=20, started_ago=86400)
    state["esmfold"] = {"signature": None, "status": "failed",
                        "error": "CUDA out of memory",
                        "finished": stamp(now - 7200)}
    state["pfam"] = {"signature": None, "status": "failed", "error": "old",
                     "finished": stamp(now - 8 * 86400)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert [f["name"] for f in v["failures"]] == ["esmfold"]
    assert v["undated_failures"] == []
    html = console.render_project(v, contract)
    assert "1 stage has failed in this run" in html
    assert "from an earlier run" in html


# 6. thirteen rows under a sentence saying nothing is known
# ----------------------------------------------------------------------

def test_a_corrupt_state_file_does_not_deny_what_the_config_told_it(
        console, results, contract):
    """The page led with "Nothing below is known about this directory" and then
    printed thirteen rows reading "OFF - not enabled: run.topology is false",
    tallied "13 off, 8 none". The OFF verdict comes from the config file, which
    is a different source and WAS read, so the rows are right and the sentence
    was the half that was wrong: it now claims exactly what it can, which is
    that no stage's PROGRESS is known here.
    """
    effective(results, contract,
              {key: "false" for key in
               sorted({st["enabled"] for st in contract.stages
                       if st["enabled"]})})
    results.raw_state("{ this is not json")
    v = view_of(console, results, contract, log=False)
    off = [r["name"] for r in v["rows"] if r["state"] == "off"]
    assert len(off) >= 5, off
    assert v["state_known"] is False
    html = console.render_project(v, contract)
    assert "OFF" in html and "not enabled: run." in html
    assert "Nothing below is known about this directory" not in html
    assert "progress below is known" in html
    assert "come from the config file, which was read" in html
    assert "corrupt state file, not an empty directory" in html


# ======================================================================
# The v0.5.0 audit. Every test below fails against the console as it was
# tagged, and each one is named for the sentence the page was getting wrong.
# ======================================================================

# ----------------------------------------------------------------------
# 1. the console could not see `cost`, so its NEXT rows misled
# ----------------------------------------------------------------------

def queued_run(results, contract, running=("emapper", "pfam", "signalp",
                                           "tmbed")):
    """The auditor's reproduction, with the real run's config: four stages on
    the box and six ready behind them, at stage_workers 4.

    Returns the `run:` block for the caller to monkeypatch in. A config in the
    directory would put a YAML parser between the test and what it is about,
    which is the order of the six.
    """
    on = {name: True for name in
          {st["enabled"] for st in contract.stages if st["enabled"]}}
    for flag in ("structure", "smorf", "context", "jackhmmer", "hhblits",
                 "unipept", "taxonomy", "join"):
        on[flag] = False
    now = time.time()
    state = running_state(contract, age=10, started_ago=3600)
    for name in running:
        state[name] = {"signature": None, "status": "running",
                       "started": stamp(now - 1800)}
    results.state(state)
    return on


def test_a_ready_stage_says_which_other_ready_stages_the_engine_takes_first(
        console, results, contract, monkeypatch):
    """v0.4.0 made the scheduler dispatch each round's ready set longest-first
    and `describe --json` emits `cost` "so a front end can order or annotate
    the table the same way". Contract dropped the field, so the page rendered
    the queue in table order: dbcan NEXT, diamond NEXT, cluster NEXT, then
    ncbifam, kofam and interpro. The operator reads the first NEXT as the stage
    about to start; the engine starts ncbifam, and on the run this comes from
    dbcan did not get a worker for over 27 hours.
    """
    on = queued_run(results, contract)
    monkeypatch.setattr(console.Project, "config_run",
                        lambda self, py, sc, st: (on, "test"))
    v = view_of(console, results, contract, log=False)
    rows = {r["name"]: r for r in v["rows"]}
    ready = [r["name"] for r in v["rows"] if r["state"] == "next"]
    assert ready == ["dbcan", "diamond", "cluster", "ncbifam", "kofam",
                     "interpro"], ready
    # the row a reader would have taken for the next stage names the three the
    # engine ranks ahead of it
    assert "ncbifam" in rows["dbcan"]["detail"], rows["dbcan"]["detail"]
    assert rows["dbcan"]["ahead"] == ["ncbifam", "kofam", "interpro"]
    # cost 3 first, and among equals the engine's own table order, because its
    # sort is stable over exactly this sequence
    assert rows["ncbifam"]["ahead"] == []
    assert rows["kofam"]["ahead"] == ["ncbifam"]
    assert rows["cluster"]["ahead"] == ["ncbifam", "kofam", "interpro",
                                        "dbcan", "diamond"]
    assert "ranks this one first" in rows["ncbifam"]["detail"]
    # the table is still in the engine's stage order - the annotation carries
    # the queue, not the row order
    names = [r["name"] for r in v["rows"]]
    assert names.index("dbcan") < names.index("ncbifam")
    html = console.render_project(v, contract)
    assert html.index(">dbcan<") < html.index(">ncbifam<")
    assert "longest-first" in html


def test_the_ranking_claims_an_order_and_never_a_schedule(console, results,
                                                          contract,
                                                          monkeypatch):
    """The console does not know how many workers are free, and gpu_lease can
    defer a GPU stage however it is ranked. So the page says which stage is
    ranked ahead of which, and the two things it cannot see are said once,
    under the table, rather than as a caveat on ten rows."""
    on = queued_run(results, contract)
    monkeypatch.setattr(console.Project, "config_run",
                        lambda self, py, sc, st: (on, "test"))
    html = console.render_project(view_of(console, results, contract,
                                          log=False), contract)
    assert "an order and not a schedule" in html
    assert "stage_workers" in html and "gpu_workers" in html
    for claim in ("starts next", "will start next", "starts in ",
                  "about to start"):
        assert claim not in html, claim
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase


def test_an_engine_with_no_cost_field_ranks_nothing_rather_than_guessing(
        console, results, contract, monkeypatch):
    """`cost` arrived in v0.4.0. Pointed at an older engine the console gets a
    contract with none at all, and a rank invented here would be exactly the
    drift taking the stage list from `describe --json` exists to prevent."""
    older = copy.deepcopy(contract)
    for st in older.stages:
        st["cost"] = None                    # what describe --json omitted
    on = queued_run(results, older)
    monkeypatch.setattr(console.Project, "config_run",
                        lambda self, py, sc, st: (on, "test"))
    project = console.Project(0, results.path, older)
    v = console.project_view(project, sys.executable, METAANNOT_PY, log=False)
    rows = {r["name"]: r for r in v["rows"]}
    assert rows["dbcan"]["state"] == "next"
    assert rows["dbcan"]["ahead"] == []
    assert "ranks" not in rows["dbcan"]["detail"], rows["dbcan"]["detail"]
    assert "no dependency is blocking it" in rows["dbcan"]["detail"]
    html = console.render_project(v, older)
    assert "longest-first" not in html
    assert "an order and not a schedule" not in html


def test_a_ready_row_does_not_say_nothing_is_blocking_it(console, results,
                                                         contract,
                                                         monkeypatch):
    """"nothing is blocking it" is a claim about the dependency graph that a
    reader takes as a claim about the queue, and for a cost-1 stage sitting
    behind three cost-3 stages the reader's version is false. The graph half is
    now said as the graph half."""
    on = queued_run(results, contract)
    monkeypatch.setattr(console.Project, "config_run",
                        lambda self, py, sc, st: (on, "test"))
    html = console.render_project(view_of(console, results, contract,
                                          log=False), contract)
    assert "nothing is blocking it" not in html
    assert "no dependency is blocking it" in html


# ----------------------------------------------------------------------
# 2. a RUN record from a killed run rendered as this run's live stage
# ----------------------------------------------------------------------

def test_a_running_record_older_than_this_run_is_not_this_runs_live_stage(
        console, results, contract):
    """mark_running() writes {signature, status, started} and never a
    `finished`, and which run a record belongs to was decided from `finished`
    alone. So a run that was killed left `status: running` behind, the next
    run's page rendered those as ITS live stages, and the duration counted up
    from the dead run's clock - 30 hours and climbing, over a run 20 minutes
    old.
    """
    now = time.time()
    state = running_state(contract, age=10, started_ago=1200)
    state["interpro"] = {"signature": None, "status": "running",
                         "started": stamp(now - 30 * 3600)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "interpro")
    assert row["state"] == "stale" and row["label"] == "STALE"
    assert row["era"] == "earlier" and row["carried"] is not None
    # the duration that was counting from a clock that stopped is gone
    assert row["took"] == "—"
    assert "before this run started" in row["detail"]
    assert v["counts"].get("running") is None
    html = console.render_project(v, contract)
    assert "from an earlier run" in html
    # and the page still refuses the verdict it always refused
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase


def test_a_running_record_from_this_run_is_still_this_runs_live_stage(
        console, results, contract):
    """The other half of the same comparison. A record stamped after the run
    started is this run's, reads RUN, and counts up."""
    now = time.time()
    state = running_state(contract, age=10, started_ago=7200)
    state["interpro"] = {"signature": None, "status": "running",
                         "started": stamp(now - 3600)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "interpro")
    assert row["state"] == "running" and row["label"] == "RUN"
    assert row["era"] == "this" and row["carried"] is None
    assert row["took"] == "1h00m"


def test_a_stale_row_still_lists_the_outputs_its_stage_declared(
        console, results, contract, monkeypatch):
    """Naming this case STALE took the row out of the set that lists declared
    outputs, and the listing is the one piece of evidence a reader of a STALE
    row actually wants: whether the earlier run got anywhere before it stopped.
    The row was `running` before v0.5.0 and it listed them then.

    The part-file SCAN is the other half and it stays gone, because it answers
    a different question - "is a tool writing bytes into this directory right
    now" - and the row's own sentence has already said that this record is not
    evidence of that. One stat per declared output; no walk of a directory
    holding 455k structures.
    """
    now = time.time()
    os.makedirs(os.path.join(results.path, "interpro"))
    results.write(os.path.join("interpro", "interproscan.tsv"), "x" * 4096)
    scans = []
    monkeypatch.setattr(console.PartCache, "newest",
                        lambda self, out: scans.append(out))
    state = running_state(contract, age=10, started_ago=1200)
    state["interpro"] = {"signature": None, "status": "running",
                         "started": stamp(now - 30 * 3600)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "interpro")
    assert row["state"] == "stale"
    assert [o["name"] for o in row["outputs"]] == ["interpro/interproscan.tsv"]
    assert row["outputs"][0]["there"] is True
    assert row["outputs"][0]["size"] == 4096
    assert "interproscan.tsv" in console.render_project(v, contract)
    # listed, and not walked
    assert row["part"] is None
    assert scans == [], scans


def test_the_index_counts_a_stale_record_as_no_running_stage(console, results,
                                                             contract):
    """The index reads the same rows, and a killed run's leftovers were a
    running stage there too - a strip cell in the live colour and a tally
    reading "1 running" over a directory nothing is writing.

    It was called "does not put the project in the live group" and asserted
    `bucket == "done"` for it, which was true of the old console as well: this
    run record carries `final_status: "ok"` and bucket_of answers from the
    record, before any row is looked at. The group is not what a leftover can
    move. What it CAN move is everything the index builds out of the rows -
    the strip and the counts - and those are what this pins, along with the
    row-scan fallback that is the one place a row does decide a group.
    """
    now = time.time()
    state = {contract.run_key: {
        "run_id": "20260101T000000-1", "version": "0.4.0",
        "host": "lab-fedora", "pid": 4242, "started": stamp(now - 600),
        "last_seen": stamp(now - 30), "last_seen_epoch": now - 30,
        "heartbeat_s": 30, "finished": stamp(now - 60), "final_status": "ok"}}
    state["interpro"] = {"signature": None, "status": "running",
                         "started": stamp(now - 30 * 3600)}
    results.state(state)
    iv = console.index_view([console.Project(0, results.path, contract)],
                            sys.executable, METAANNOT_PY)
    it = iv["projects"][0]
    assert [c["state"] for c in it["strip"] if c["name"] == "interpro"] == [
        "stale"]
    assert it["counts"].get("running") is None
    assert it["counts"]["stale"] == 1
    # decided by the record, and the record says the run finished ok
    assert it["bucket"] == "done"
    # and the fallback under it, which is the one place a ROW picks the group:
    # a state file with no run record at all to read. One leftover is not a
    # live run there either - it is nothing this console can call.
    rows_only = [{"state": c["state"]} for c in it["strip"]]
    assert console.bucket_of({"interpro": {}}, contract,
                             rows_only) == "unknown"
    index = console.render_index(iv, contract)
    assert "running record from an earlier run" in index      # the legend
    assert ".c-stale" in console.PAGE_CSS


def test_a_run_record_that_never_says_it_ended_is_live_on_both_pages(
        console, results, contract):
    """The index badge and the project page have to answer "is this run over?"
    the same way, and for the one record shape that never says, they did not.

    run_is_over() is where that rule lives: a run record ends a run by SAYING
    so, and a missing `final_status` is not that statement - which is the
    engine's own doctrine, unprovable means alive, and why the heartbeat block
    prints "read here as still going" over exactly this shape. bucket_of
    enumerated the values the key can hold and had no branch for its absence,
    so the index answered from the stages instead: DONE once they were all ok,
    and UNCLEAR once v0.5.0 stopped calling a leftover `running`. Neither is a
    thing the page beside it was saying.
    """
    now = time.time()
    run = {"run_id": "20260101T000000-1", "version": "0.4.0",
           "host": "lab-fedora", "pid": 4242, "started": stamp(now - 600),
           "last_seen": stamp(now - 30), "last_seen_epoch": now - 30,
           "heartbeat_s": 30}                 # and no final_status, ever
    # every stage green under a record that never recorded a finish: the old
    # index called this DONE while the page under it said the run was alive
    state = {contract.run_key: dict(run)}
    for name in contract.stage_names:
        state[name] = {"signature": "a", "status": "ok", "seconds": 5.0,
                       "finished": stamp(now - 120)}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    assert v["heartbeat"]["running"] is True
    assert "read here as still going" in " ".join(v["heartbeat"]["lines"])
    assert v["bucket"] == "running"

    # and the shape this came from: the same record with a leftover `running`
    # from a run that ended 30 hours ago. The leftover is not the evidence -
    # the record is - and both pages read the record the same way.
    state = {contract.run_key: dict(run),
             "interpro": {"signature": None, "status": "running",
                          "started": stamp(now - 30 * 3600)}}
    results.state(state)
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "interpro")
    assert row["state"] == "stale"
    assert v["heartbeat"]["running"] is True
    assert v["bucket"] == "running"
    iv = console.index_view([console.Project(0, results.path, contract)],
                            sys.executable, METAANNOT_PY)
    assert iv["projects"][0]["bucket"] == "running"


# ----------------------------------------------------------------------
# 3. a disabled stage with an old record blocked nothing, and was reported
#    as blocking
# ----------------------------------------------------------------------

@pytest.mark.parametrize("status", ["failed", "running"])
def test_a_disabled_dependency_with_a_record_still_does_not_block(
        console, results, contract, monkeypatch, status):
    """A stage with a record is never called OFF - it ran, and the record is
    the evidence - and the dependency logic read that same set to decide what
    the ENGINE waits for. It is not the same question: decide() returns
    `disabled (run.X)` for any stage whose flag is falsy, whatever the state
    file holds, and finish() then adds it to `done`. So a dependency this
    config has turned off, carrying an earlier run's `failed` or a killed run's
    `running`, was reported as "waiting on dbcan (failed)" over a run the
    engine had already walked straight past.
    """
    off = {name: False for name in
           {st["enabled"] for st in contract.stages if st["enabled"]}}
    off["context"] = True
    monkeypatch.setattr(console.Project, "config_run",
                        lambda self, py, sc, st: (off, "test"))
    now = time.time()
    old = stamp(now - 5 * 86400)
    # each written the way the engine writes it: finish() stamps `finished`,
    # mark_running() stamps `started` and nothing else.
    state = running_state(contract, age=10, started_ago=600)
    state["dbcan"] = ({"signature": None, "status": "failed",
                       "error": "hmmsearch exited 1", "finished": old}
                      if status == "failed" else
                      {"signature": None, "status": "running",
                       "started": old})
    results.state(state)
    v = view_of(console, results, contract, log=False)
    rows = {r["name"]: r for r in v["rows"]}
    # context reads emapper, pfam, signalp and dbcan; all four are off
    assert rows["context"]["state"] == "next", rows["context"]["detail"]
    assert "dbcan" not in rows["context"]["detail"]
    # and the deliberate half is untouched: a stage with a record is not OFF
    assert rows["dbcan"]["state"] != "off"
    assert rows["dbcan"]["enabled_key"] == "dbcan"


def test_a_dependency_that_is_merely_unfinished_still_blocks(console, results,
                                                             contract,
                                                             monkeypatch):
    """The other half. An ENABLED dependency with no record is the case the
    WAIT row exists for, and nothing above may weaken it."""
    on = {name: True for name in
          {st["enabled"] for st in contract.stages if st["enabled"]}}
    monkeypatch.setattr(console.Project, "config_run",
                        lambda self, py, sc, st: (on, "test"))
    results.state(running_state(contract, age=10, started_ago=600))
    v = view_of(console, results, contract, log=False)
    row = next(r for r in v["rows"] if r["name"] == "context")
    assert row["state"] == "wait"
    assert "dbcan" in row["detail"]


# ----------------------------------------------------------------------
# 4. the part cache never hit for the one stage it was written for
# ----------------------------------------------------------------------

def test_the_part_cache_bounds_the_scan_of_a_directory_that_keeps_changing(
        console, results, monkeypatch):
    """esmfold's declared output is results/structures/.done and the stage
    writes one .pdb per dark protein into that same directory, so the
    directory's mtime moves between every poll and the mtime key never hits.
    On the 455,571-protein run that was a 455k-entry walk per poll per open
    tab, for the whole of a multi-day stage.
    """
    out = os.path.join(results.path, "structures", ".done")
    os.makedirs(os.path.dirname(out))
    scans = []
    real = console.newest_part_file
    monkeypatch.setattr(console, "newest_part_file",
                        lambda p: (scans.append(p), real(p))[1])
    cache = console.PartCache()
    stamp_at = time.time()
    for _ in range(20):
        # what one more .pdb landing in structures/ does, and nothing else
        stamp_at += 3
        os.utime(os.path.dirname(out), (stamp_at, stamp_at))
        assert cache.newest(out) is None
    assert len(scans) == 1, "the directory was re-scanned %d times" % len(scans)


def test_the_part_cache_still_finds_a_part_file_that_appears(console, results,
                                                             monkeypatch):
    """The floor is a rate limit, not a blindfold: a stray part file is real
    evidence and the console still finds it, on the next scan rather than on
    the next poll."""
    out = os.path.join(results.path, "structures", ".done")
    os.makedirs(os.path.dirname(out))
    cache = console.PartCache()
    assert cache.newest(out) is None
    part = os.path.join(os.path.dirname(out), "..done.4242.140.part")
    with open(part, "w") as fh:
        fh.write("z" * 900)
    monkeypatch.setattr(console, "PART_RESCAN_S", 0.0)
    found = cache.newest(out)
    assert found is not None and found["size"] == 900
    assert found["name"] == "..done.4242.140.part"


# ----------------------------------------------------------------------
# 5. the log pane grew for as long as the tab was open
# ----------------------------------------------------------------------

def test_the_log_pane_is_capped_and_says_so_when_it_drops_lines(console,
                                                                results,
                                                                contract):
    """Every line /api/log delivered became a <div> that stayed, on a page the
    run book tells you to leave open for the length of a three-day run. The cap
    is the server's, so there is one place it is decided, and it is announced:
    a pane holding the last two thousand lines looks exactly like a pane
    holding the whole log.
    """
    results.log("[    1.0s] INFO  hello\n")
    v = view_of(console, results, contract)
    html = console.render_project(v, contract)
    assert 'data-max="%d"' % console.LOG_PANE_LINES in html
    # a cap under what the first render already puts in the pane would trim on
    # the first poll, before a single line had arrived
    assert console.LOG_PANE_LINES >= console.LOG_LINES
    js = console.PAGE_JS
    assert "data-max" in js and "pane.removeChild" in js
    assert "dropped from " in js and "which keeps the last " in js
    # and the notice is a node in the pane, not a console.log nobody sees
    assert "pane.insertBefore(capline" in js


def test_the_log_pane_keeps_the_number_of_lines_its_notice_claims(console):
    """The notice the pane writes into itself says it "keeps the last 2000",
    and the trim loop has to make that sentence true. It bounded
    `pane.childElementCount`, which is not a count of lines: one marker per
    rotation and the "no log lines" placeholder from the first render live in
    that same scroller, and each of them took a slot out of the cap. The loop
    subtracted its OWN notice by hand and nothing else, which is the tell - a
    pane that had rotated twice held 1,998 lines under a sentence promising
    2,000.

    The button carried the same miscount the other way. `unseen` counted every
    line appended and nothing ever took the trimmed ones back off it, so after
    a burst larger than the cap "N new lines below" offered to scroll to lines
    the pane no longer held.

    Both are cosmetic and both are numbers this page states out loud, which is
    the whole reason the cap is announced rather than applied quietly. There is
    no JavaScript engine in this suite, so what is pinned is the source: the
    bound is a line count, and the offer is clamped to it. A DOM stub written
    here would pin the stub.
    """
    js = console.PAGE_JS
    # comments stripped, because the comment above the loop names the old
    # bound in order to explain it and would otherwise answer for it
    code = "\n".join(re.sub(r"//.*", "", line) for line in js.splitlines())
    # lines, not children - the retired bound, by name, so it cannot come back
    assert "childElementCount" not in code
    assert "while (shown > cap)" in js
    assert "var shown = pane.querySelectorAll(" in js   # what the server sent
    # `shown` moves for a line appended and for a line dropped, and for
    # nothing else: a marker is removed without freeing a slot
    assert "unseen++; shown++;" in js
    assert "{ dropped++; shown--; }" in js
    # and the offer may not name a line that is gone
    assert "if (unseen > shown) unseen = shown;" in js
    # a rotation empties the pane, so both counters go with the lines
    assert "unseen = 0; dropped = 0; capline = null; shown = 0;" in js


# ----------------------------------------------------------------------
# 6. the failure block described a scheduler that does not exist
# ----------------------------------------------------------------------

def test_the_failure_block_says_what_the_scheduler_does_at_a_failure(
        console, results, contract):
    """The block said "The engine carries on with every stage that does not
    depend on a casualty". It does not: the dispatch loop is
    `while (remaining or futures) and not failure`, so the first failure ends
    dispatch and what follows is a drain of whatever is already in flight -
    which is why a run can still look busy for hours afterwards.
    """
    results.state(failed_mid_run(contract))
    html = console.render_project(view_of(console, results, contract,
                                          log=False), contract)
    assert "Nothing new starts after the first failure" in html
    assert "cannot interrupt" in html
    # this run has NOT ended, so the past-tense half must not ship - the other
    # side of the boundary its own test pins from the ended end
    assert "This run has ended" not in html
    # and the retired sentence is guarded here, where the block that replaced
    # it is under test: unconditionally, against ever coming back
    assert "The engine carries on" not in html
    for phrase in FORBIDDEN:
        assert phrase not in html.lower(), phrase
