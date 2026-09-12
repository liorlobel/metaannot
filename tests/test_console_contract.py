"""The contract between the engine and `console/console.py`.

`docs/gui-design.md` names this file specifically, and it exists for one reason:
the console reads names and shapes out of metaannot, and if that coupling drifts
nothing at run time notices — the console simply starts describing a run that is
not happening. So every name and every shape it reads is pinned here, against
`describe --json` and against a real run's own files, not against a docstring.

Three things it asserts that no other test can:

1. **The console never writes into a results directory.** By AST scan (the only
   write-mode open in the file is the socket lock, in the console's own runtime
   directory) and empirically: a real results directory is snapshotted, every
   route is served, and the snapshot must be identical afterwards.
2. **The console never imports metaannot.** It shells out to `describe --json`.
3. **The file names come from `describe --json` by relpath**, not from string
   literals. `.metaannot_state.json` appears nowhere in the console's source.
"""
import ast
import http.client
import io
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time

import pytest

from conftest import CONSOLE_PY, METAANNOT_PY, ROOT, build_project


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def describe(*args):
    """`describe --json`, parsed from stdout alone."""
    proc = subprocess.run([sys.executable, METAANNOT_PY, "describe", "--json",
                           *args], capture_output=True, text=True, cwd=ROOT,
                          timeout=300)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def snapshot(root):
    """Every path under `root`, with everything but its access time.

    atime is excluded deliberately and is the one thing reading does change —
    `cat` changes it too. Size, mtime, ctime, inode and mode are what a write
    would move, and a new or removed file shows up as a changed key set.
    """
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        for name in list(dirnames) + list(filenames):
            path = os.path.join(dirpath, name)
            st = os.lstat(path)
            out[os.path.relpath(path, root)] = (
                st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_ino,
                st.st_mode, st.st_nlink)
    return out


class UnixConn(http.client.HTTPConnection):
    """http.client over an AF_UNIX socket — what `ssh -L` gives the browser."""

    def __init__(self, path):
        http.client.HTTPConnection.__init__(self, "localhost")
        self._path = path

    def connect(self):
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(30)
        self.sock.connect(self._path)


@pytest.fixture
def short_sock():
    """A socket path short enough for sun_path.

    pytest's tmp_path is 100+ characters deep on macOS, and AF_UNIX allows 104
    there and 108 on Linux — which is exactly why the console checks the length
    and says so instead of failing at bind(). The tests need a short directory
    of their own, so /tmp it is.
    """
    d = tempfile.mkdtemp(prefix="mac", dir="/tmp")
    try:
        yield os.path.join(d, "s.sock")
    finally:
        shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def ran(tmp_path, cli):
    """A results directory a real run actually produced."""
    project = build_project(tmp_path / "proj")
    cli("run", "--config", os.path.join(str(tmp_path / "proj"), "config.yaml"),
        cwd=str(tmp_path / "proj"))
    return os.path.join(str(tmp_path / "proj"), "results")


def source():
    with open(CONSOLE_PY, encoding="utf-8") as fh:
        return fh.read()


# ----------------------------------------------------------------------
# the document itself
# ----------------------------------------------------------------------

def test_describe_json_is_on_stdout_and_creates_nothing(tmp_path):
    """The console parses stdout only, and asking what metaannot is must not
    make a results directory."""
    project = build_project(tmp_path / "p")
    root = str(tmp_path / "p")
    before = snapshot(root)
    proc = subprocess.run([sys.executable, METAANNOT_PY, "describe", "--json",
                           "--config", os.path.join(root, "config.yaml")],
                          capture_output=True, text=True, cwd=root, timeout=300)
    assert proc.returncode == 0, proc.stderr
    json.loads(proc.stdout)                       # stdout alone must parse
    assert snapshot(root) == before
    assert not os.path.exists(os.path.join(root, "results"))


def test_the_keys_the_console_reads_are_all_present():
    doc = describe()
    for key in ("describe_version", "metaannot_version", "signature_version",
                "generated", "host", "config", "default_config", "paths",
                "run_key", "stage_names", "stages"):
        assert key in doc, key
    for key in ("results_dir", "state", "log", "lock", "effective_config"):
        assert key in doc["paths"], key
    assert isinstance(doc["config"]["run"], dict)
    assert isinstance(doc["default_config"]["heartbeat_s"], (int, float))


def test_every_stage_carries_the_fields_the_table_renders():
    doc = describe()
    names = set(doc["stage_names"])
    assert len(doc["stages"]) == len(doc["stage_names"]) == 21
    for st in doc["stages"]:
        assert set(st) >= {"name", "enabled", "deps", "keys", "gpu",
                           "empty_ok", "outputs"}
        assert st["name"] in names
        assert set(st["deps"]) <= names, st["name"]
        # `enabled` is a key into config["run"], or None for the two stages
        # that have no flag and always run. The console renders "not enabled:
        # run.<key>" off exactly this.
        assert st["enabled"] is None or st["enabled"] in doc["config"]["run"]
        assert isinstance(st["gpu"], bool) and isinstance(st["empty_ok"], bool)


def test_the_stage_graph_is_acyclic_and_in_topological_order():
    """The console renders the table in `stages` order and derives WAIT from
    `deps`; both are nonsense if a stage precedes one it depends on."""
    seen = set()
    for st in describe()["stages"]:
        assert set(st["deps"]) <= seen, st["name"]
        seen.add(st["name"])


def test_every_stage_output_lives_under_the_results_directory():
    """The console remaps outputs onto each project by relpath against
    `paths.results_dir`. An output outside it could not be remapped, and the
    console drops such a row rather than guessing at a path."""
    doc = describe()
    root = doc["paths"]["results_dir"]
    for st in doc["stages"]:
        for out in st["outputs"]:
            rel = os.path.relpath(out, root)
            assert not rel.startswith(os.pardir), (st["name"], out)


# ----------------------------------------------------------------------
# the names, derived rather than typed
# ----------------------------------------------------------------------

def test_the_four_file_names_come_from_describe_not_from_literals(console):
    doc = describe()
    contract = console.Contract(doc)
    root = doc["paths"]["results_dir"]
    assert contract.state_name == os.path.relpath(doc["paths"]["state"], root)
    assert contract.log_name == os.path.relpath(doc["paths"]["log"], root)
    assert contract.lock_name == os.path.relpath(doc["paths"]["lock"], root)
    assert contract.effective_name == os.path.relpath(
        doc["paths"]["effective_config"], root)
    # What they happen to be today. If the engine moves one, this line is the
    # only thing in the repository that has to change - the console follows.
    assert contract.state_name == ".metaannot_state.json"
    assert contract.log_name == "metaannot.log"
    assert contract.lock_name == ".metaannot.lock"
    assert contract.effective_name == "config.effective.yaml"


def test_the_console_source_hard_codes_none_of_those_names():
    """The assertion that makes the coupling real: if any of these strings were
    typed into the console, the test above would keep passing while the console
    quietly stopped following the engine."""
    # Comments and prose may name them - the module docstring does. Code may
    # not, and ast.parse leaves only the strings the code actually evaluates.
    names = {".metaannot_state.json", "metaannot.log", ".metaannot.lock",
             "config.effective.yaml"}
    for node in ast.walk(ast.parse(source())):
        if isinstance(node, ast.Constant) and node.value in names:
            raise AssertionError("console/console.py contains the literal %r; "
                                 "it must come from describe --json"
                                 % node.value)


def test_the_run_key_and_run_record_are_what_the_console_reads(console, ran):
    """A real run's own state file, against the console's Contract."""
    contract = console.Contract(describe())
    with open(os.path.join(ran, contract.state_name)) as fh:
        state = json.load(fh)
    assert contract.run_key in state
    run = state[contract.run_key]
    for field in ("run_id", "version", "config_path", "argv", "host", "pid",
                  "started", "last_seen", "last_seen_epoch", "heartbeat_s",
                  "finished", "final_status"):
        assert field in run, field
    assert isinstance(run["last_seen_epoch"], float)
    assert isinstance(run["heartbeat_s"], (int, float))
    assert run["final_status"] in ("running", "ok", "failed", "interrupted")
    # The run key is not a stage, and the console must never render it as one.
    assert contract.run_key not in contract.stage_names


def test_the_stage_record_shapes_are_what_the_table_renders(console, ran):
    """`ok` records carry seconds and finished and NO started — the console
    shows a duration from `seconds` and refuses to invent a start time.

    The exact key set moved once, deliberately, when the state file stopped
    being written from one process's snapshot: `run_id` says WHICH run wrote
    the record, which is what makes a document two runs have both written into
    auditable rather than inferred from timestamps. It is additive, it is read
    with `.get()` on the console side, and it is not a status — the vocabulary
    test below is unmoved. The assertion stays EXACT so that the next field
    cannot arrive without this being read again.
    """
    contract = console.Contract(describe())
    with open(os.path.join(ran, contract.state_name)) as fh:
        state = json.load(fh)
    done = [v for k, v in state.items() if k != contract.run_key]
    assert done, "the fixture run recorded no stages"
    for rec in done:
        assert rec["status"] in ("running", "ok", "adopted", "failed")
        if rec["status"] == "ok":
            assert set(rec) == {"signature", "status", "seconds", "finished",
                                "run_id"}
            assert "started" not in rec
            assert rec["run_id"] == state[contract.run_key]["run_id"], \
                "a record written by this run has to name this run"


def test_the_statuses_the_console_knows_are_the_statuses_the_engine_writes():
    """Grepped out of the engine, so a fifth status cannot appear without this
    test failing. The console renders an unknown status as `?`, but it should
    not have to."""
    with open(METAANNOT_PY, encoding="utf-8") as fh:
        engine = fh.read()
    written = set(re.findall(r'"status":\s*"(\w+)"', engine))
    assert written == {"running", "ok", "adopted", "failed"}, written


def test_the_lock_file_holds_what_the_heartbeat_block_prints(console, ma,
                                                             tmp_path):
    lock = ma.ResultsLock(str(tmp_path / "lk"))
    with lock:
        with open(str(tmp_path / "lk")) as fh:
            token = json.load(fh)
    assert set(token) == {"pid", "host", "started"}
    view = console.lock_view(
        type("P", (), {"lock_path": str(tmp_path / "lk")})())
    assert view["held"] is False          # released on exit; nothing holds it


def test_the_log_prefix_the_console_parses_is_the_prefix_the_engine_writes(
        console, ma, tmp_path, monkeypatch):
    """LOG_LINE_RE has to survive the real `[%7.1fs] LEVEL      tag | ` shape,
    the tagless form, and a continuation line padded with spaces."""
    path = str(tmp_path / "l.log")
    monkeypatch.setattr(ma, "_LOGFH", open(path, "w", encoding="utf-8"))
    ma.log("plain and untagged")
    ma.set_log_context("interpro")
    ma.log("two\nlines", "WARN")
    ma.set_log_context(None)
    ma._LOGFH.close()
    with open(path) as fh:
        lines = console.parse_log(fh.read())
    assert [l["level"] for l in lines] == ["INFO", "WARN", "WARN"]
    assert lines[1]["stage"] == "interpro"
    assert lines[0]["stage"] == ""
    assert lines[2].get("cont") is True    # the continuation keeps the level


def test_the_cost_the_console_ranks_by_is_the_one_the_scheduler_sorts_by(
        console, ma):
    """The console annotates its NEXT rows with the engine's dispatch order,
    and it takes that order from `cost` in `describe --json` - the field the
    engine emits "so a front end can order or annotate the table the same way".

    Pinned against stage_priority() itself rather than against the numbers,
    because the coupling is the point: the day a stage's cost changes, or the
    day the scheduler stops sorting by it, the console must move with it and
    not carry on annotating a queue that is no longer there.
    """
    contract = console.Contract(describe())
    cost = {st["name"]: st["cost"] for st in contract.stages}
    for name, value in cost.items():
        assert value == ma.STAGE_COSTS[name], name
    # every stage, so the ranking is never a partial sort over a missing key
    assert all(value is not None for value in cost.values())
    # and the console's order over one ready set is the engine's own
    ready = ["dbcan", "diamond", "cluster", "ncbifam", "kofam", "interpro"]
    rows = [{"name": n, "state": "next", "detail": "", "ahead": [],
             "cost": cost[n]} for n in ready]
    console.rank_ready(rows)
    # `ahead` IS the rank: the stage with none ahead of it is dispatched first.
    console_order = [r["name"]
                     for r in sorted(rows, key=lambda r: len(r["ahead"]))]
    engine_order = list(ready)
    engine_order.sort(key=ma.stage_priority, reverse=True)
    assert console_order == engine_order


def test_a_running_record_carries_a_start_and_never_a_finish():
    """Which run a record belongs to is one timestamp comparison, and for a
    record left behind by mark_running() the only timestamp there is is
    `started`. Comparing `finished` alone made every such record this run's,
    however old it was - so the shape mark_running writes is pinned here.
    """
    with open(METAANNOT_PY, encoding="utf-8") as fh:
        engine = fh.read()
    body = engine.split("def mark_running(", 1)[1].split("\n    def ", 1)[0]
    assert '"status": "running"' in body
    assert '"started": time.strftime' in body
    assert '"finished"' not in body


def test_heartbeat_seconds_are_read_from_the_record_not_assumed(console):
    """The engine writes `heartbeat_s` into `_run` precisely so a reader need
    not know its defaults; the console's bands are computed from it."""
    now = time.time()
    slow = {"final_status": "running", "heartbeat_s": 600,
            "last_seen_epoch": now - 900, "started": "2026-01-01T00:00:00"}
    fast = dict(slow, heartbeat_s=1)
    assert console.heartbeat_view(slow, now, None)["band"] == "fresh"
    assert console.heartbeat_view(fast, now, None)["band"] == "long"


def test_the_engine_does_not_claim_nobody_reads_heartbeat_s(console):
    """The scope of a claim is part of the claim.

    metaannot.py's comment on `heartbeat_s` said "nothing in this file reads it
    back", which was true; correcting a different error in the same comment
    widened it to "nothing ANYWHERE reads `_run.heartbeat_s` and decides
    something on it", which is false - the test above is the console doing
    exactly that, and `0`, the engine's documented way of turning the heartbeat
    off, gets a band and a sentence of its own.

    So this pins both halves against each other: the console really does decide
    on the field, and the engine's comment really is scoped to the engine. A
    reader of that comment who concluded the number is free to change is the
    person this protects.
    """
    now = time.time()
    rec = {"final_status": "running", "heartbeat_s": 0,
           "last_seen_epoch": now - 900, "started": "2026-01-01T00:00:00"}
    off = console.heartbeat_view(rec, now, None)
    rated = console.heartbeat_view(dict(rec, heartbeat_s=30), now, None)
    assert off["band"] != rated["band"] and off["verdict"] != rated["verdict"], \
        "`heartbeat_s: 0` in the record must not read the same as a cadence"

    src = io.open(METAANNOT_PY, encoding="utf-8").read()
    assert "Nothing IN THIS FILE reads" in src, \
        "the heartbeat_s comment no longer scopes its claim to this file"
    assert "Nothing anywhere reads" not in src, \
        "the heartbeat_s comment claims more than this codebase can check"


# ----------------------------------------------------------------------
# read-only, and no import of the engine
# ----------------------------------------------------------------------

def test_the_console_never_imports_metaannot():
    tree = ast.parse(source())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                assert "metaannot" not in a.name
        if isinstance(node, ast.ImportFrom):
            assert "metaannot" not in (node.module or "")
    assert "import metaannot" not in source()
    assert "importlib" not in source(), "no back door either"


def test_the_console_imports_only_the_standard_library():
    """One file, deployable by scp, with no pip install behind it."""
    tree = ast.parse(source())
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods |= {a.name.split(".")[0] for a in node.names}
        if isinstance(node, ast.ImportFrom):
            assert node.level == 0, "no relative imports: there is no package"
            mods.add((node.module or "").split(".")[0])
    stdlib = os.path.dirname(os.__file__)
    for name in sorted(mods):
        if name in sys.builtin_module_names:
            continue
        mod = __import__(name)
        path = getattr(mod, "__file__", "") or ""
        assert path.startswith(stdlib), "%s is not stdlib: %s" % (name, path)
        assert "site-packages" not in path


def test_the_only_write_mode_open_is_the_socket_lock():
    """Layer one of the read-only guarantee, where a reviewer can see it.

    Every reader in the console uses os.open(..., os.O_RDONLY). The single
    exception is take_socket_name(), which creates a flock beside the socket in
    the console's own runtime directory — never in a results directory.
    """
    tree = ast.parse(source())
    offenders = []
    for func in ast.walk(tree):
        if not isinstance(func, ast.FunctionDef):
            continue
        for node in ast.walk(func):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "id", None) or getattr(node.func, "attr",
                                                             None)
            if name not in ("open", "makedirs", "mkdir", "remove", "unlink",
                            "rename", "replace", "chmod", "chown", "utime",
                            "truncate", "mkfifo", "symlink", "link", "rmdir"):
                continue
            if name == "open":
                flags = ast.dump(node.args[1]) if len(node.args) > 1 else ""
                if "O_RDONLY" in flags and "O_CREAT" not in flags:
                    continue
            offenders.append((func.name, name))
    assert sorted(set(offenders)) == [
        # os.mkdir, and one component: os.makedirs made the parent of an
        # explicit --socket, so the console could CREATE a results directory,
        # and it applies its mode to the leaf alone, so the intermediates came
        # out 0777 & ~umask. Nothing here creates a directory for a path the
        # operator named — see socket_path().
        ("ensure_private_dir", "mkdir"),    # the console's own runtime dir
        ("serve", "unlink"),                # its own socket, on shutdown
        ("take_socket_name", "open"),       # the flock beside the socket
        ("take_socket_name", "unlink"),     # a stale socket of its own
    ], sorted(set(offenders))


def test_the_console_never_chmods_anything():
    """The socket is 0700 because bind() runs under umask 0o077, not because
    something widened it and narrowed it again. A chmod call anywhere in the
    file would mean such a window existed."""
    for node in ast.walk(ast.parse(source())):
        if isinstance(node, ast.Call):
            name = (getattr(node.func, "id", None)
                    or getattr(node.func, "attr", None))
            assert name not in ("chmod", "fchmod", "lchmod")


def test_serving_every_route_changes_nothing_in_a_results_directory(
        console, ran, short_sock):
    """Layer three, and the one that would catch a write no reviewer spotted.

    A results directory a real run produced, snapshotted, then every route the
    console serves — index, project page, both fragments, all four JSON
    endpoints, a log tail at an offset — and snapshotted again.
    """
    args = console.build_parser().parse_args(
        ["--project", ran, "--python", sys.executable,
         "--metaannot", METAANNOT_PY])
    con = console.build(args)
    assert len(con.projects) == 1

    sock = short_sock
    before = snapshot(ran)
    stop, box = threading.Event(), {}
    thread = threading.Thread(
        target=console.serve, args=(con, sock),
        kwargs={"stop": stop, "ready": lambda s: box.setdefault("srv", s)},
        daemon=True)
    thread.start()
    try:
        for _ in range(200):
            if "srv" in box:
                break
            time.sleep(0.02)
        assert "srv" in box, "the server never came up"

        conn = UnixConn(sock)
        for route in ("/", "/?f=1", "/app.css", "/app.js", "/api/contract",
                      "/api/projects", "/p/0", "/p/0?f=1", "/api/project?p=0",
                      "/api/log?p=0", "/api/log?p=0&off=10&id=1:1",
                      "/api/log?p=0&off=999999999&id=nope", "/p/9", "/nope"):
            conn.request("GET", route)
            resp = conn.getresponse()
            body = resp.read()
            assert resp.status in (200, 404), (route, resp.status)
            assert len(body) == int(resp.getheader("Content-Length"))
        conn.close()
    finally:
        stop.set()
        thread.join(10)

    after = snapshot(ran)
    assert after == before, {
        k: (before.get(k), after.get(k))
        for k in set(before) | set(after) if before.get(k) != after.get(k)}
    assert not os.path.exists(sock)         # and it cleaned up after itself


# ----------------------------------------------------------------------
# what the engine will actually run, pinned against the engine
# ----------------------------------------------------------------------

FALSY_CONFIG = """\
proteins_faa: /nowhere/proteins.faa
quant_table: /nowhere/quant.tsv
run:
  pfam: 0
  dbcan:
  interpro: ""
  kofam: []
  ncbifam: false
  eggnog: true
  cluster: yes
"""


def test_off_means_what_the_engine_means_by_off(console, tmp_path):
    """The console had a second source of truth for config semantics: it
    decided OFF with `is False` while metaannot.py decides with plain
    truthiness (`not cfg["run"].get(st["enabled"], False)`). `run.pfam: 0`,
    `: null`, `: ""` and `: []` each made the engine skip a stage the table
    called NEXT. docs/gui-design.md is arranged to prevent exactly this, so
    the agreement is asserted against `describe --json` rather than described.
    """
    cfg = tmp_path / "config.yaml"
    cfg.write_text(FALSY_CONFIG, encoding="utf-8")
    doc = describe("--config", str(cfg))
    run = doc["config"]["run"]
    # The engine's own expression, copied from cmd_describe's summary line.
    engine_on = {st["name"] for st in doc["stages"]
                 if st["enabled"] is None or run.get(st["enabled"])}

    contract = console.Contract(doc)
    rows = console.stage_rows(contract, {}, run, time.time(), False,
                              str(tmp_path))
    console_off = {r["name"] for r in rows if r["state"] == "off"}
    engine_off = {st["name"] for st in doc["stages"]
                  if st["name"] not in engine_on}
    assert console_off == engine_off, {
        "the console calls off but the engine runs":
            sorted(console_off - engine_off),
        "the engine skips but the console does not call off":
            sorted(engine_off - console_off)}
    # and the four falsy shapes are actually in this fixture
    assert {"pfam", "dbcan", "interpro", "kofam"} <= console_off
    assert "emapper" in engine_on and "cluster" in engine_on


def test_the_off_sentence_names_the_value_the_engine_saw(console, tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(FALSY_CONFIG, encoding="utf-8")
    doc = describe("--config", str(cfg))
    contract = console.Contract(doc)
    rows = console.stage_rows(contract, {}, doc["config"]["run"], time.time(),
                              False, str(tmp_path))
    by = {r["name"]: r["detail"] for r in rows}
    assert "run.pfam is 0" in by["pfam"]
    assert "run.dbcan is empty" in by["dbcan"]
    assert 'run.interpro is ""' in by["interpro"]
    assert "run.kofam is []" in by["kofam"]
    assert "run.ncbifam is false" in by["ncbifam"]


# ----------------------------------------------------------------------
# how the readers open a file, not just what they open
# ----------------------------------------------------------------------

def test_every_reader_opens_non_blocking_and_checks_it_is_a_regular_file():
    """A FIFO named .metaannot_state.json in one watched directory blocked
    every route the console serves, /app.css included, until it was restarted.
    Both halves of the fix are pinned: O_NONBLOCK so the open returns, and
    S_ISREG on the FD - never on the path, which is a race."""
    tree = ast.parse(source())
    opens = []
    for func in ast.walk(tree):
        if not isinstance(func, ast.FunctionDef):
            continue
        for node in ast.walk(func):
            if (isinstance(node, ast.Call)
                    and getattr(node.func, "attr", None) == "open"
                    and getattr(getattr(node.func, "value", None), "id", None)
                    == "os"):
                opens.append((func.name, ast.dump(node.args[1])
                              if len(node.args) > 1 else ""))
    assert opens, "no os.open in the console at all"
    for name, flags in opens:
        assert "O_NONBLOCK" in flags, "%s opens without O_NONBLOCK" % name
    # exactly one function opens files out of a watched directory
    readers = [n for n, f in opens if "O_CREAT" not in f]
    assert readers == ["open_regular"], readers
    src = source()
    body = src.split("def open_regular", 1)[1].split("\ndef ", 1)[0]
    assert "S_ISREG(st.st_mode)" in body
    assert "os.fstat(fd)" in body


def test_the_socket_lock_is_opened_with_o_nofollow():
    """One flag. Without it the call followed a planted symlink and created a
    file at the attacker's chosen path, owned by the console's uid."""
    src = source()
    body = src.split("def take_socket_name", 1)[1].split("\ndef ", 1)[0]
    assert "O_NOFOLLOW" in body
    assert "S_ISREG" in body


def test_the_runtime_directory_is_checked_with_lstat_as_the_report_claims():
    """The build report states the runtime directory "must lstat as a 0700
    directory this uid owns"; os.stat resolved a planted symlink first and
    validated the target instead."""
    src = source()
    body = src.split("def why_not_private", 1)[1].split("\ndef ", 1)[0]
    assert "lstat_or_none" in body
    assert not re.search(r"(?<!l)stat_or_none", body), "os.stat follows links"
    # and the ancestor walk, which is where the second half of that hole was
    chain = src.split("def exposed_ancestor", 1)[1].split("\ndef ", 1)[0]
    assert "lstat_or_none" in chain
    assert "os.path.realpath" in chain, "a symlink in the chain is not checked"


def test_no_command_string_is_built_from_an_unvalidated_field():
    """`"ps -p %s ..." % pid` straight from the lock file rendered
    `ps -p 21877; rm -rf / -o pid,etime,stat,args` on the served page, under a
    sentence telling the operator to paste it on the pipeline host."""
    src = source()
    body = src.split("def lock_view", 1)[1].split("\ndef ", 1)[0]
    assert "valid_pid(pid)" in body
    assert "%s -o pid" not in body
    assert "%d -o pid" in body


# ----------------------------------------------------------------------
# layer four: the socket itself is never allowed near a results directory
# ----------------------------------------------------------------------

def test_no_socket_can_be_bound_inside_a_watched_directory(console, ran,
                                                           capsys):
    """The AST scan and the snapshot are both blind to this one, which is why
    it survived all three read-only layers: the write is a bind() and an
    O_CREAT of `<socket>.lock`, both at a path the OPERATOR chose, made by a
    function whose docstring already claimed the lock lives "never in a results
    directory".

    Driven through main(), the way the reviewer drove it: the banner must not
    appear, the exit code is the refusal\'s, and the directory is untouched.
    """
    before = snapshot(ran)
    code = console.main([ran, "--socket", os.path.join(ran, "console.sock"),
                         "--metaannot", METAANNOT_PY, "--python",
                         sys.executable])
    out = capsys.readouterr()
    assert code == 2
    assert "read-only. Nothing is written" not in out.out, out.out
    assert "is inside" in out.err
    assert snapshot(ran) == before
    assert not os.path.exists(os.path.join(ran, "console.sock"))
    assert not os.path.exists(os.path.join(ran, "console.sock.lock"))


def test_no_socket_directory_is_created_by_naming_one(console, tmp_path):
    """`--socket <root>/DatasetZ/results/c.sock` created a results directory at
    0700 with a stray lock file in it, which a run started into that path
    afterwards then finds."""
    root = tmp_path / "runs"
    root.mkdir()
    code = console.main(["--root", str(root), "--socket",
                         str(root / "DatasetZ" / "results" / "c.sock"),
                         "--metaannot", METAANNOT_PY, "--python",
                         sys.executable])
    assert code == 2
    assert list(root.iterdir()) == []


UGLY = ("fifo", "dir-for-state", "unreadable", "corrupt", "failed", "empty")


def build_ugly(root, contract):
    """One results directory of every shape that has broken this console."""
    made = {}
    for kind in UGLY:
        d = os.path.join(root, kind, "results")
        os.makedirs(d)
        made[kind] = d
        with open(os.path.join(d, contract.log_name), "w") as fh:
            fh.write("[    0.0s] INFO  hello\n[    1.0s] FATAL  died\n")
    os.mkfifo(os.path.join(made["fifo"], contract.state_name))
    os.mkdir(os.path.join(made["dir-for-state"], contract.state_name))
    with open(os.path.join(made["corrupt"], contract.state_name), "w") as fh:
        fh.write("{ this is not json")
    with open(os.path.join(made["empty"], contract.state_name), "w") as fh:
        fh.write("")
    now = time.time()
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(now - 7200))
    with open(os.path.join(made["failed"], contract.state_name), "w") as fh:
        json.dump({contract.run_key: {
            "run_id": "r", "version": "0.3.0", "host": "elsewhere", "pid": 3,
            "started": time.strftime("%Y-%m-%dT%H:%M:%S",
                                     time.localtime(now - 86400)),
            "last_seen_epoch": now - 4 * 3600, "heartbeat_s": 30,
            "finished": None, "final_status": "running"},
            "esmfold": {"signature": None, "status": "failed",
                        "error": "CUDA out of memory", "finished": stamp}}, fh)
    with open(os.path.join(made["failed"], contract.lock_name), "w") as fh:
        json.dump({"pid": "0021877", "host": "elsewhere",
                   "started": stamp}, fh)
    os.chmod(made["unreadable"], 0o000)
    return made


def test_serving_every_route_changes_nothing_in_the_ugly_directories(
        console, tmp_path, short_sock):
    """The layer the reviewer called "the one that would catch a write no
    reviewer spotted", widened past the one --project and the fourteen routes
    it was written with.

    Every shape that has broken this console at once - a FIFO in place of a
    state file, a directory in place of one, a mode-000 results directory, a
    corrupt state file, an empty one, a failed stage under a stale heartbeat -
    discovered through --root so the rescan path runs too, and /api/log at a
    real prior offset rather than a made-up one.
    """
    contract = console.Contract(describe())
    root = str(tmp_path / "runs")
    os.makedirs(root)
    made = build_ugly(root, contract)
    args = console.build_parser().parse_args(
        ["--root", root, "--python", sys.executable, "--metaannot",
         METAANNOT_PY])
    con = console.build(args)
    assert len(con.projects) == len(UGLY), [p.path for p in con.projects]

    before = snapshot(root)
    stop, box, after = threading.Event(), {}, None
    thread = threading.Thread(
        target=console.serve, args=(con, short_sock),
        kwargs={"stop": stop, "ready": lambda s: box.setdefault("srv", s)},
        daemon=True)
    thread.start()
    try:
        for _ in range(300):
            if "srv" in box:
                break
            time.sleep(0.02)
        assert "srv" in box, "the server never came up"
        conn = UnixConn(short_sock)
        routes = ["/", "/?f=1", "/app.css", "/app.js", "/api/contract",
                  "/api/projects", "/nope", "/p/99"]
        for i in range(len(UGLY)):
            routes += ["/p/%d" % i, "/p/%d?f=1" % i, "/api/project?p=%d" % i,
                       "/api/log?p=%d" % i]
        offsets = {}
        for route in routes:
            conn.request("GET", route)
            resp = conn.getresponse()
            body = resp.read()
            assert resp.status in (200, 404), (route, resp.status)
            assert len(body) == int(resp.getheader("Content-Length"))
            if route.startswith("/api/log"):
                doc = json.loads(body.decode("utf-8"))
                offsets[route] = (doc.get("off"), doc.get("ident"))
        # and again from where the first read actually stopped
        for route, (off, ident) in offsets.items():
            conn.request("GET", "%s&off=%s&id=%s" % (route, off, ident or ""))
            resp = conn.getresponse()
            resp.read()
            assert resp.status == 200
        conn.request("GET", "/")            # a second index: the rescan path
        conn.getresponse().read()
        conn.close()
        # No chmod between the two snapshots: putting the mode back and
        # taking it away again moves the directory's own ctime, which is
        # exactly what this test is looking for.
        stop.set()
        thread.join(10)
        after = snapshot(root)
    finally:
        stop.set()
        thread.join(10)
        # Unconditional, including on a failure above: pytest cannot clean a
        # tmp tree with a mode-000 directory in it, and it then carries that
        # tree forward into every later session as `garbage-<uuid>`.
        os.chmod(made["unreadable"], 0o755)
    assert after == before, {
        k: (before.get(k), after.get(k))
        for k in set(before) | set(after) if before.get(k) != after.get(k)}
    assert not os.path.exists(short_sock)


def test_describe_emits_every_field_a_stage_dict_carries(ma, tmp_path):
    """A field added to STAGES must reach the console, or it is invisible.

    v0.4.0 added `cost` (1 seconds / 2 minutes / 3 hours, what the scheduler
    now sorts each round's ready set by) and `describe --json` went on emitting
    its hardcoded seven fields, so every stage came back `cost: None` and a
    front end could not have known the field existed. The projection is a list
    of names, so nothing failed - which is exactly why this compares the two
    sets rather than checking that the fields it knows about are present.
    """
    cfg = dict(ma.DEFAULT_CONFIG, results_dir=str(tmp_path / "r"))
    doc = ma.describe(cfg, ma.Paths(cfg), None)
    emitted = set(doc["stages"][0])
    # `out` and `inp` are lambdas; describe renders the first as `outputs` and
    # deliberately does not call the second, which needs a config it may not
    # have. `fn` is the stage function itself. Everything else is data and
    # belongs in the document.
    carried = set().union(*(set(st) for st in ma.STAGES)) - {"out", "inp", "fn"}
    missing = sorted(carried - emitted - {"outputs"})
    assert not missing, (
        f"STAGES carries {missing} and describe --json does not emit it, so a "
        f"console reading the contract cannot see it")


def test_a_parked_output_is_not_offered_to_the_console_as_work_in_progress(
        console, ma, tmp_path):
    """The naming rule behind `_park_superseded`, pinned against the reader it
    is a rule for.

    A superseded run's finished output is parked beside its target rather than
    renamed onto it. The console finds work in progress by matching
    `.<stem>.` in the declared output's own directory, and it consults that
    for any stage whose record says `running` - which is exactly the stage the
    REPLACEMENT is now running. A parked file named after the stem would
    therefore be offered to an operator as the live run's tool writing bytes
    right now, and, being older than the live temp, would eventually be
    reported as stalled against a run that is perfectly healthy. So the marker
    goes in FRONT of the stem, and this is what says the two halves still
    agree.
    """
    out = str(tmp_path / "hmm" / "pfam.tblout")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w", encoding="utf-8").close()
    state = str(tmp_path / ".metaannot_state.json")
    rec = ma.RunRecord(state, {}, ["metaannot", "run"], None, 999)
    ma.update_state(state, {}, (ma.RUN_KEY,), claim=rec.claim)
    ma._DECLARED_OUTPUTS.add(out)
    ma.save_state(state, {ma.RUN_KEY: {"run_id": "20260101T000000-9",
                                       "host": "elsewhere", "pid": 9,
                                       "started": "2026-01-01T00:00:00"}})
    ma._STATE_WATCH["seen"] = ma._STATE_WATCH["tried"] = 0.0
    with ma.atomic_out(out) as tmp:
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write("the superseded run's result\n")

    parked = [f for f in os.listdir(os.path.dirname(out))
              if f.startswith(ma.SUPERSEDED_SUFFIX)]
    assert parked, "nothing was parked, so there is nothing to check here"
    assert console.newest_part_file(out) is None, \
        "the console reads a parked output as the live run's work in progress"
