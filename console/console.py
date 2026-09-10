#!/usr/bin/env python3
"""metaannot console, M1 — a read-only watcher for results directories.

One stdlib-only file, deployable by `scp`. It runs on the host that writes the
results directories it watches, binds a mode-0700 UNIX socket, and serves one
page over an SSH forward:

    console$  python3 console/console.py --root /data/runs
    laptop$   ssh -N -L 8080:/run/user/1000/metaannot.sock lab-fedora
    laptop$   open http://localhost:8080/

**It never writes a single byte into a results directory.** Not a lockfile, not
a cache, not a temp file, not a log line. That is the whole reason it is safe to
point at a three-day job that is already running. Every reader below opens
O_RDONLY, the console never chdir()s into a results directory (so not even a
core dump lands there), there is no do_POST, and the one file it does create -
the flock that stops two consoles fighting over one socket - lives in the
console's own runtime directory. `tests/test_console_contract.py` asserts all of
that, by AST scan and by snapshotting a results directory around every route.

It also never imports metaannot. Stage order, dependency edges, the run-record
key and the names of the files it polls all come from `metaannot.py describe
--json`, shelled out once at startup and cached. A console that hard-codes what
a stage is named or where a state file lives drifts from the engine and starts
lying, and `docs/gui-design.md` is arranged to prevent exactly that.

What it will not do, deliberately: it will not tell you a run is dead. The
heartbeat in `_run.last_seen` is advisory - one failed write ends the heartbeat
thread while the run continues - and the engine's own rule is that unprovable
means alive. So the console reports how long it has been since the last sign of
work, names the evidence, hands over the `ps` line, and leaves the verdict to
the person reading it. There is no --force-unlock here and there is no button
that acts.

Python 3.9+, standard library only.
"""

import argparse
import fcntl
import hashlib
import html
import json
import os
import re
import signal
import socket
import socketserver
import stat
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

CONSOLE_VERSION = "0.1.0"

# The log tail is read by seeking to (size - window), never by reading the file:
# these logs run for days and reach hundreds of megabytes. 256 KB is far more
# than the 200 lines that are actually displayed, and bounds every read in the
# file to one pread-sized chunk.
LOG_WINDOW = 256 * 1024
LOG_LINES = 200

# How many lines the browser keeps in the log pane. The pane is appended to by
# byte offset for as long as the tab is open, and the whole design is that a tab
# stays open for the length of a three-day run: nothing trimmed it, so an
# InterProScan run emitting a few lines a second put a quarter of a million
# <div>s in one scroller and the tab's memory climbed until the browser killed
# it. Ten times what the first render shows, so the reader who scrolls up has
# somewhere to scroll, and the pane says out loud when it has dropped lines -
# a tail that looks like the whole log is the reason to bound it visibly rather
# than quietly.
LOG_PANE_LINES = 2000

# What a poll that is not showing the log still reads, to see the level of the
# newest line. Small on purpose: this runs every few seconds per open tab.
PEEK_WINDOW = 4096

# A state file is small - one record per stage plus `_run`. Anything past this
# is not a state file, and reading it whole would be the one unbounded read.
STATE_MAX = 8 * 1024 * 1024
LOCK_MAX = 64 * 1024

# os.replace() is atomic on POSIX, so a reader sees the old file or the new one
# and never a torn one. The retry is for everything else: an NFS revalidation, a
# Windows share, a file that genuinely does not exist yet on the first write.
STATE_RETRIES = 3
STATE_RETRY_S = 0.05

MAX_CONNS = 32
HANDLER_TIMEOUT = 30
PROBE_TIMEOUT = 0.25
DESCRIBE_TIMEOUT = 120

# The one argument this console passes to another process that a WATCHED
# directory chooses: `describe --config <config.effective.yaml>`. A 2.4 MB YAML
# made that fork take 61.6 s on a request thread holding a MAX_CONNS permit, so
# 33 concurrent requests 503'd /app.css and every healthy project - the FIFO
# wedge's exact symptom through an ordinary regular file. A 676-byte alias bomb
# does it with a tiny file, so the size cap is a courtesy and CONFIG_WAIT is
# the fix.
CONFIG_MAX = 1024 * 1024           # no config this engine writes is near this
CONFIG_TIMEOUT = 30                # the fork's own ceiling, off-thread
CONFIG_WAIT = 5.0                  # what ONE request thread will wait for it
# And a ceiling on how many of those forks may exist at once, so twenty
# projects on one index do not become twenty interpreters.
CONFIG_FORKS = threading.BoundedSemaphore(2)

# How deep --root is scanned. Deep enough for <root>/<dataset>/results, which is
# the layout examples/server-run-plan/ describes, and shallow enough that it
# cannot wander into a 200 GB foldseek tmp directory on an NFS mount.
SCAN_DEPTH = 3

# How often --root is looked at again. Not a watch: a scandir of a few
# directories, rate-limited, on the thread that serves the index.
RESCAN_S = 30.0

BUSY_RESPONSE = (
    b"HTTP/1.1 503 Service Unavailable\r\n"
    b"Content-Type: text/plain; charset=utf-8\r\n"
    b"Content-Length: 62\r\n"
    b"Connection: close\r\n\r\n"
    b"console: too many open connections; close a tab and try again\n"
)


class Refuse(Exception):
    """A startup condition the operator has to fix. Printed, not raised."""


# ======================================================================
# 1. the contract: what the engine says about itself
# ======================================================================

def describe_argv(python, script, config=None):
    """The command whose stdout is the contract. JSON on stdout, engine log
    lines on stderr - the two are never merged."""
    argv = [python, script, "describe", "--json"]
    if config:
        argv += ["--config", config]
    return argv


def run_describe(python, script, config=None, timeout=DESCRIBE_TIMEOUT):
    """`describe --json` as a dict. Raises Refuse with the command and its
    stderr, because a console that guesses the contract is the drift the
    contract exists to prevent.

    cwd is the engine's own directory: with no --config, `results_dir` is
    resolved against the process's cwd, and a console that ran this from inside
    a results directory would report that directory's own name back at itself.
    describe() creates nothing - cmd_describe says so explicitly and takes no
    lock - so this is a read even against a live run.
    """
    argv = describe_argv(python, script, config)
    try:
        proc = subprocess.run(argv, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              cwd=os.path.dirname(os.path.abspath(script)),
                              timeout=timeout)
    except OSError as e:
        raise Refuse("could not run the engine: %s\n  %s"
                     % (e, " ".join(argv)))
    except subprocess.TimeoutExpired:
        raise Refuse("the engine did not answer within %ds:\n  %s"
                     % (timeout, " ".join(argv)))
    err = proc.stderr.decode("utf-8", "replace").strip()
    if proc.returncode != 0:
        raise Refuse("the engine exited %d:\n  %s\n%s"
                     % (proc.returncode, " ".join(argv), err))
    try:
        return json.loads(proc.stdout.decode("utf-8", "replace"))
    except ValueError as e:
        raise Refuse("the engine's describe --json did not parse (%s):\n  %s"
                     % (e, " ".join(argv)))


class Contract:
    """Stage order, the run-record key, and the names of the files to poll.

    Every one of those is read out of `describe --json`; not one is a literal in
    this file. `.metaannot_state.json` in particular is derived by relpath from
    `paths.state` against `paths.results_dir`, so the day the engine moves it
    the console moves with it - and `tests/test_console_contract.py` asserts
    that derivation rather than the string, so the coupling cannot drift
    invisibly.
    """

    def __init__(self, doc):
        paths = doc["paths"]
        root = paths["results_dir"]
        self.state_name = self._rel(paths["state"], root, "state")
        self.log_name = self._rel(paths["log"], root, "log")
        self.lock_name = self._rel(paths["lock"], root, "lock")
        self.effective_name = self._rel(paths["effective_config"], root,
                                        "effective_config")
        self.run_key = doc["run_key"]
        self.stages = []
        for st in doc["stages"]:
            outs = []
            for out in st.get("outputs", ()):
                rel = os.path.relpath(out, root)
                # A stage output that is not under the results directory cannot
                # be remapped onto another project, and guessing would point the
                # page at a file belonging to a different run. Dropped, quietly:
                # no stage does this today and the contract test says so.
                if not rel.startswith(os.pardir + os.sep) and rel != os.pardir:
                    outs.append(rel)
            self.stages.append({
                "name": st["name"],
                "enabled": st["enabled"],
                "deps": list(st["deps"]),
                "gpu": bool(st.get("gpu")),
                "empty_ok": bool(st.get("empty_ok")),
                "keys": list(st.get("keys", ())),
                "outputs": outs,
                # What the scheduler sorts each round's ready set by: 1
                # seconds, 2 minutes, 3 hours, longest dispatched first.
                # describe --json says it is emitted "so a front end can order
                # or annotate the table the same way", and this file dropped
                # it, so the table listed ten NEXT rows in table order over a
                # queue the engine was taking in a different one. None for an
                # engine older than the field, and rank_ready() below declines
                # to rank at all rather than invent a number for it.
                #
                # finite(), because a cost arrives out of a subprocess's JSON
                # like every other number on this page: it refuses bools by
                # name - `true` would otherwise sort as 1 and silently rank a
                # stage - and refuses NaN and the infinities, which sort but do
                # not order. Nothing here does arithmetic on a cost, so the
                # float it hands back is only ever a sort key.
                "cost": finite(st.get("cost")),
            })
        self.stage_names = [s["name"] for s in self.stages]
        self.version = doc.get("metaannot_version", "?")
        self.describe_version = doc.get("describe_version")
        self.signature_version = doc.get("signature_version")
        self.host = doc.get("host", "")

    @staticmethod
    def _rel(path, root, what):
        rel = os.path.relpath(path, root)
        if os.sep in rel or rel.startswith(os.pardir):
            raise Refuse(
                "describe --json puts %s at %r, which is not directly inside "
                "the results directory %r. This console watches one directory "
                "per project and cannot follow that." % (what, path, root))
        return rel

    def facts(self):
        """The subset the page shows, and /api/contract returns."""
        return {
            "metaannot_version": self.version,
            "describe_version": self.describe_version,
            "signature_version": self.signature_version,
            "engine_host": self.host,
            "run_key": self.run_key,
            "state_name": self.state_name,
            "log_name": self.log_name,
            "lock_name": self.lock_name,
            "effective_name": self.effective_name,
            "stage_names": list(self.stage_names),
            "console_version": CONSOLE_VERSION,
        }


# ======================================================================
# 2. readers - pure functions, O_RDONLY, and every one of them returns
#    a value for missing / empty / truncated / invalid / raced input
# ======================================================================

def stat_or_reason(path):
    """(stat, None) or (None, why). EACCES is not absence.

    stat_or_none() folded every OSError into None, and the vitals block then
    told the operator "there is no log file in this directory" about a results
    directory that was merely group-unreadable - which on a shared lab box is
    the difference between "nothing ran" and "you cannot see what ran".

    ValueError as well as OSError, and it is not hypothetical: os.stat raises
    ValueError("embedded null byte") for a NUL in the path and
    UnicodeEncodeError - a ValueError subclass - for a lone surrogate, neither
    of which is an OSError. `_run.config_path` is a string out of somebody's
    state file and reaches here through regular_stat(), so one such string
    raised straight through this reader and 500'd the page. A path the
    operating system cannot even be asked about is not missing; it is a path
    this console cannot use, which is what it says.
    """
    try:
        return os.stat(path), None
    except FileNotFoundError:
        return None, "missing"
    except OSError as e:
        return None, "unreadable: %s" % e
    except ValueError as e:
        return None, "not a usable path: %s" % e


def stat_or_none(path):
    return stat_or_reason(path)[0]


def lstat_or_none(path):
    """os.lstat, for the checks where a symlink must not be followed.

    Every access-control check in this file is on a directory this console
    creates or is handed; os.stat resolves the link first and would validate
    the TARGET's uid and mode while binding the socket at the LINK's path.
    """
    try:
        return os.lstat(path)
    except OSError:
        return None


def open_regular(path):
    """(fd, stat, None) or (None, None, why). O_RDONLY, O_NONBLOCK, S_ISREG.

    Every file this console opens sits in a directory it does not control, and
    `os.open(path, O_RDONLY)` on a FIFO blocks until a writer appears - which
    for a FIFO named `.metaannot_state.json` is never. That thread then holds
    one of MAX_CONNS permits forever; because the index reads every project's
    state file, MAX_CONNS index loads exhaust the semaphore and every route
    after that, /app.css included, is 503 until the process is restarted.
    HANDLER_TIMEOUT does not help: it bounds socket reads, not handler work.

    O_NONBLOCK makes the open return, and S_ISREG on the fd - not on the path,
    which would be a race - refuses a pipe, a directory or a device before a
    single byte is read. It is not a cure for a hard-mounted NFS path whose
    server is unreachable; nothing in userspace is. It closes the case that a
    local file can cause.
    """
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NONBLOCK)
    except FileNotFoundError:
        return None, None, "missing"
    except OSError as e:
        return None, None, "unreadable: %s" % e
    try:
        st = os.fstat(fd)
    except OSError as e:
        os.close(fd)
        return None, None, "unreadable: %s" % e
    if not stat.S_ISREG(st.st_mode):
        os.close(fd)
        return None, None, "not a regular file"
    return fd, st, None


def read_head(path, limit):
    """At most `limit` bytes of `path`, O_RDONLY. (None, reason) on failure.

    os.open, not open(): the mode is visible to a reader of this file and to the
    AST scan in the contract test, and there is no code path here that can
    accidentally acquire O_CREAT.
    """
    fd, _, why = open_regular(path)
    if fd is None:
        return None, why
    try:
        chunks, got = [], 0
        while got <= limit:
            block = os.read(fd, 65536)
            if not block:
                break
            chunks.append(block)
            got += len(block)
        if got > limit:
            return None, "larger than %d bytes" % limit
        return b"".join(chunks), None
    except OSError as e:
        return None, "unreadable: %s" % e
    finally:
        os.close(fd)


def read_json_once(path, limit):
    """One attempt. -> (status, obj, detail), status in ok/missing/empty/bad."""
    raw, why = read_head(path, limit)
    if raw is None:
        return ("missing", None, None) if why == "missing" else ("bad", None, why)
    if not raw.strip():
        return "empty", None, None
    try:
        obj = json.loads(raw.decode("utf-8", "replace"))
    except RecursionError:
        # json.loads recurses per nesting level and raises this, which is not
        # a ValueError. A 40 KB file of `{"a":` + `[`*20000 + `]`*20000 + `}`
        # 500'd the index and every project page, for every dataset.
        return ("bad", None, "nested far too deeply to parse; this is not a "
                             "state file")
    except ValueError as e:
        return "bad", None, "not valid JSON: %s" % e
    if not isinstance(obj, dict):
        return "bad", None, "top level is %s, not an object" % type(obj).__name__
    return "ok", obj, None


def read_json_stable(path, limit, retries=STATE_RETRIES, delay=STATE_RETRY_S,
                     retry_missing=False):
    """read_json_once, retried past a rewrite.

    The engine writes state as temp-file-then-os.replace, which is atomic on
    POSIX: a reader gets the whole old file or the whole new one. The retry is
    for the cases where that guarantee is thinner - an NFS or SMB revalidation,
    a Windows share, the microseconds between the lock's O_EXCL create and its
    write - and it costs 100 ms in the worst case, once per request.

    A file that is simply NOT THERE is not retried unless the caller says it
    was there a moment ago. A results directory where nothing has run yet has
    no state file and never will until it does, and sleeping 100 ms per such
    directory would put most of a second into every poll of an eight-project
    index for nothing.
    """
    for attempt in range(retries):
        status, obj, detail = read_json_once(path, limit)
        if status == "ok" or (status == "missing" and not retry_missing):
            return status, obj, detail
        if attempt + 1 < retries:
            time.sleep(delay)
    return status, obj, detail


def join_notes(*notes):
    """One `note` string out of several, dropping the empty ones."""
    kept = [n for n in notes if n]
    return " · ".join(kept) if kept else None


def tail_log(path, want_off=None, want_ident=None, window=LOG_WINDOW):
    """The end of a log, by byte offset. Never reads the whole file.

    `want_ident` is "<dev>:<ino>" as handed out by a previous call, and
    `want_off` the byte the caller has already seen. When either fails to match
    what is on disk now - the file was rotated, truncated, or has grown by more
    than one window since - the read restarts from (size - window) and says so
    in `note` rather than guessing.
    """
    out = {"ok": False, "ident": None, "size": 0, "off": 0, "start": 0,
           "text": "", "reset": True, "note": None, "mtime": None}
    fd, st, why = open_regular(path)
    if fd is None:
        out["note"] = "no log file" if why == "missing" else "log %s" % why
        return out
    try:
        ident = "%d:%d" % (st.st_dev, st.st_ino)
        size = st.st_size
        out.update(ok=True, ident=ident, size=size, mtime=st.st_mtime)
        start, reset, note = None, True, None
        if want_ident is not None and want_off is not None:
            if want_ident != ident:
                note = ("the log was replaced - a different file now has this "
                        "name; re-reading the tail")
            elif want_off > size:
                note = ("the log is shorter than it was (rotated or "
                        "truncated); re-reading the tail")
            elif size - want_off > window:
                note = ("%s of log arrived since the last read, more than the "
                        "%s window; the middle is not shown"
                        % (human_bytes(size - want_off), human_bytes(window)))
            else:
                start, reset = want_off, False
        if start is None:
            start = max(0, size - window)
        want = size - start
        if want <= 0:
            out.update(off=size, start=start, reset=reset, note=note, text="")
            return out
        os.lseek(fd, start, os.SEEK_SET)
        chunks, got = [], 0
        while got < want:
            block = os.read(fd, min(65536, want - got))
            if not block:
                break                     # truncated under us; take what we got
            chunks.append(block)
            got += len(block)
        raw = b"".join(chunks)
        # Where the caller must resume, computed from where we actually seeked
        # and how much came back - NOT from where the text below begins. The
        # trim on the next line moves the start of the TEXT, and folding that
        # into the offset would tell the next poll to skip bytes it never saw.
        end = start + got
        # A window that starts mid-line begins with a fragment of a line nobody
        # can read. Dropped - but only when we actually seeked into the file.
        if reset and start > 0:
            cut = raw.find(b"\n")
            if cut < 0:
                # The whole window is one unterminated line: a tool writing a
                # progress bar with \r, or a single enormous record. Saying
                # nothing here renders as "no log lines", which is what an
                # EMPTY log looks like, and the two are not the same fact.
                note = join_notes(
                    note, "the last %s of this log contains no line break, so "
                          "there is no complete line to show"
                          % human_bytes(len(raw)))
            start = end if cut < 0 else start + cut + 1
            raw = b"" if cut < 0 else raw[cut + 1:]
        # errors="replace", for the same reason the engine's own opener() uses
        # it: this log carries tool stderr, and a latin-1 byte in an InterPro
        # description has already killed one stage on a real run. A console that
        # raised UnicodeDecodeError on that very log would be a bitter joke.
        out.update(off=end, start=start, reset=reset, note=note,
                   text=raw.decode("utf-8", "replace"))
        return out
    except OSError as e:
        out["note"] = "log unreadable: %s" % e
        return out
    finally:
        os.close(fd)


LOG_LINE_RE = re.compile(r"^\[\s*[\d.]+s\] (\w+) +(?:(\S+) \| )?(.*)$")


def text_lines(text):
    """`text` split on newlines, and on nothing else.

    str.splitlines() also splits on \r, \x0b, \x0c, \x1c-\x1e, \x85,
    U+2028 and U+2029, so one tool-emitted line containing U+2028 came back as
    TWO, the second parsed as a well-formed engine line the engine never wrote
    (`all stages finished ok`, forged from inside one INFO line) and the total
    that feeds `dropped` inflated with it.
    """
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()                    # the split after a trailing newline
    return [ln[:-1] if ln.endswith("\r") else ln for ln in lines]


def parse_log(text, keep=LOG_LINES):
    """Log text as [{level, stage, text}], newest last.

    The engine's prefix is fixed - `[{elapsed:7.1f}s] {level:5s} ` plus
    `{tag:>10s} | ` when a stage is set - and continuation lines are padded to
    the prefix width with spaces. A line that does not match is a continuation:
    it inherits the level above it, which is what keeps a multi-line FATAL
    traceback tinted all the way down.
    """
    lines = text_lines(text)
    if len(lines) > keep:
        lines = lines[-keep:]
    out, level = [], "INFO"
    for line in lines:
        m = LOG_LINE_RE.match(line)
        if m:
            level = m.group(1)
            out.append({"level": level, "stage": m.group(2) or "",
                        "text": line})
        else:
            out.append({"level": level, "stage": "", "text": line,
                        "cont": True})
    return out


def log_lines(text, keep=LOG_LINES):
    """(lines, dropped) — parse_log, plus how many older lines it left out.

    parse_log truncates its input to the newest `keep`, and /api/log was
    handing it a DELTA: 500 lines appended between two polls returned 200 and
    advanced the byte offset past all 500, so 300 lines the operator never saw
    were gone with nothing said. tail_log goes to some trouble to name the
    larger-than-window case; this is the same honesty one function along.
    """
    total = len(text_lines(text))
    return parse_log(text, keep), max(0, total - keep)


def dropped_note(dropped, reset, keep=LOG_LINES):
    """What to say when log_lines left some out. None when it did not."""
    if not dropped:
        return None
    if reset:
        return ("the window holds more than the %d lines this pane shows; %d "
                "older ones are above it and not in view" % (keep, dropped))
    return ("%d lines arrived since the last poll, more than the %d this pane "
            "adds at once; the %d oldest of them are not shown below"
            % (dropped + keep, keep, dropped))


def newest_part_file(out_path):
    """The in-progress file beside a declared output, if one is there.

    The engine writes through atomic_out(), which puts the work in progress in
    the output's own directory as `.<stem>.<pid>.<tid>.part<ext>`. Matched on
    the leading dot and the stem only: the suffix is the engine's business and
    is not in `describe --json`, so hard-coding it here would be exactly the
    drift this file avoids elsewhere.

    That file is the strongest evidence a console can offer that work is
    happening - it is the bytes the tool is writing right now, not a heartbeat
    thread's opinion about them.
    """
    d, base = os.path.split(out_path)
    stem = os.path.splitext(base)[0]
    if not stem:
        return None
    best = None
    try:
        with os.scandir(d) as it:
            for entry in it:
                if not entry.name.startswith("." + stem + "."):
                    continue
                try:
                    st = entry.stat(follow_symlinks=False)
                except OSError:
                    continue
                if not stat.S_ISREG(st.st_mode):
                    continue
                if best is None or st.st_mtime > best["mtime"]:
                    best = {"name": entry.name, "size": st.st_size,
                            "mtime": st.st_mtime}
    except OSError:
        return None
    return best


# A floor on how often ONE output directory is scanned for a part file,
# whatever its mtime says. The mtime key below is exactly right for a directory
# holding a handful of declared outputs, and exactly wrong for the one that
# costs the most: esmfold's declared output is results/structures/.done and the
# stage writes one .pdb per dark protein into that same directory, so the
# directory's mtime moves between every poll and the key never hits. On the
# 455,571-protein run that was a 455k-entry walk per poll per open tab, for the
# whole of a multi-day stage. A part file is evidence and not a control, so
# noticing one up to half a minute late costs a sentence on one poll; walking a
# quarter of a million directory entries every three seconds costs the box the
# run is on.
PART_RESCAN_S = 30.0


class PartCache:
    """newest_part_file, keyed on the directory's own mtime and rate-limited.

    Finding a part file is a scandir, and for esmfold that directory is
    structures/ with one file per dark protein - 455k entries on a real run,
    walked once per running stage per poll per open tab. A directory's mtime
    moves exactly when an entry is created or removed in it, which is exactly
    when the ANSWER can change; a part file merely growing changes neither. So
    the scan happens when the directory changes and every poll in between costs
    one stat of a name already known.

    And when the directory changes CONSTANTLY - which is what a stage writing
    one file per protein into its own output directory does - the mtime key
    invalidates on every poll and buys nothing at all. PART_RESCAN_S is the
    floor under that case: between scans the last answer stands, re-stat()ed so
    the size it reports is still live. The scan itself is unchanged, so a stray
    part file in any directory is still found; it is found on the next scan
    rather than on the next poll.

    monotonic() rather than time(): this is an interval between two events in
    one process, and an NTP step backwards over a lab server's first hour would
    otherwise park the floor for as long as the step.
    """

    def __init__(self):
        self._lock = threading.Lock()
        # out_path -> (dir mtime, name or None, monotonic time of that scan)
        self._by_out = {}

    def newest(self, out_path):
        dst = stat_or_none(os.path.dirname(out_path))
        if dst is None:
            return None
        with self._lock:
            hit = self._by_out.get(out_path)
        if hit is not None and (hit[0] == dst.st_mtime
                                or time.monotonic() - hit[2] < PART_RESCAN_S):
            if hit[1] is None:
                return None
            st = stat_or_none(os.path.join(os.path.dirname(out_path), hit[1]))
            if st is not None and stat.S_ISREG(st.st_mode):
                return {"name": hit[1], "size": st.st_size,
                        "mtime": st.st_mtime}
            # Gone without the directory mtime moving should not happen; fall
            # through to a real scan rather than report a file that is not
            # there.
        best = newest_part_file(out_path)
        with self._lock:
            self._by_out[out_path] = (dst.st_mtime,
                                      best["name"] if best else None,
                                      time.monotonic())
        return best


# ======================================================================
# 3. model - what a project is, and what the page says about it
# ======================================================================

# Every number and timestamp below comes out of a file the console does not
# control, and three separate permanent 500s came from trusting one. ONE gate,
# not four try blocks, because the next such value will be read by a function
# nobody has written yet. Python's json ACCEPTS NaN and Infinity, so the "not
# valid JSON" gate never sees them - `last_seen_epoch: Infinity` reached
# time.localtime, a stage's `seconds: Infinity` reached int(). A thousand-digit
# int parses and then overflows float() - `1e18` is OSError(22). time.mktime
# and time.localtime raise outside the platform's range, which for one
# directory's `_run.started: "1899-01-01T00:00:00"` was a permanent 500 on /
# and /api/projects for EVERY watched project. Out comes a real, finite,
# in-range float or None, and None is a shape every renderer here already
# handles. A page that 500s tells the operator nothing at all.
INF = float("inf")
EPOCH_MIN = -2208988800.0          # 1900-01-01Z
EPOCH_MAX = 7258118400.0           # 2200-01-01Z
SECONDS_MAX = 3153600000.0         # a hundred years, as a duration


def finite(val, lo=None, hi=None):
    """A real, finite number out of an untrusted file, or None. Never raises."""
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        return None
    try:
        num = float(val)
    except (OverflowError, ValueError):     # an int too large to be a float
        return None
    if num != num or num in (INF, -INF):    # NaN, Infinity, -Infinity
        return None
    if (lo is not None and num < lo) or (hi is not None and num > hi):
        return None
    return num


def epoch_of(val):
    """finite(), bounded to what time.localtime and time.strftime can render."""
    return finite(val, EPOCH_MIN, EPOCH_MAX)


def clip(text, limit):
    """A string out of a watched file, bounded, saying what it left out.

    The engine caps a stage's `error` at 500 characters; nothing caps what
    somebody else's file puts there, and it was printed twice - the pinned
    block and the row - into a fragment refetched every three seconds. A 4 MB
    state file made an 8 MB fragment.
    """
    text = "" if text is None else str(text)
    if len(text) <= limit:
        return text
    return "%s … (%d more characters)" % (text[:limit], len(text) - limit)


def stamp_epoch(s):
    """Seconds since the epoch for a '%Y-%m-%dT%H:%M:%S' engine stamp.

    OverflowError as well as ValueError: time.mktime raises it for a stamp
    outside the platform's range. The result goes back through epoch_of() so a
    stamp that converts here but cannot be rendered comes back None now, rather
    than as an exception inside time.localtime three functions later.
    """
    try:
        return epoch_of(time.mktime(time.strptime(str(s),
                                                  "%Y-%m-%dT%H:%M:%S")))
    except (ValueError, TypeError, OverflowError, OSError):
        return None


def dur(seconds):
    """A duration in the engine's own shapes, extended past a day.

    Below an hour `12m14s`, below a day `6h40m` - character for character what
    `_elapsed_str` puts in the log, so a duration on this page and the same
    duration in the log line beside it are not two dialects. Past a day it says
    `2d 05h`, because `53h20m` is not a number anyone reads.
    """
    seconds = finite(seconds, hi=SECONDS_MAX)
    if seconds is None:
        return "?"
    s = int(max(0, seconds))
    if s >= 86400:
        return "%dd %02dh" % (s // 86400, s % 86400 // 3600)
    if s >= 3600:
        return "%dh%02dm" % (s // 3600, s % 3600 // 60)
    return "%dm%02ds" % (s // 60, s % 60)


def human_bytes(n):
    n = finite(n, 0)
    if n is None:
        return "?"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return "%d %s" % (n, unit) if unit == "B" else "%.1f %s" % (n, unit)
        n /= 1024.0
    return "%d B" % n


def clock(epoch):
    epoch = epoch_of(epoch)
    if not epoch:
        return "?"
    try:
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(epoch))
    except (ValueError, OSError, OverflowError):
        return "?"          # a platform narrower than EPOCH_MIN..EPOCH_MAX


def collapse_home(path):
    home = os.path.expanduser("~")
    if home and path.startswith(home + os.sep):
        return "~" + path[len(home):]
    return path


def elide_path(path, keep=2):
    """A path shortened from the LEFT, because the left is the shared part.

    With the layout the run book describes - /data/runs/<dataset>/results,
    eight of them - CSS `text-overflow: ellipsis` cut off precisely the
    segment that tells one row from another and left eight identical prefixes.
    The full path is still on the row, in the title attribute.
    """
    parts = [p for p in path.split(os.sep) if p]
    if len(parts) <= keep + 1:
        return path
    return "…" + os.sep + os.sep.join(parts[-keep:])


def regular_stat(path):
    """os.stat, but only for a path that really is a regular file."""
    st = stat_or_none(path)
    return st if st is not None and stat.S_ISREG(st.st_mode) else None


def one_line(err, limit=400):
    """A multi-line Refuse squeezed onto one line, keeping the diagnosis.

    `str(e).splitlines()[0]` kept "the engine exited 1:" and threw away every
    line that says why - which for a config the engine cannot parse is the only
    sentence the operator can act on, and it governs whether nine of the
    twenty-one rows read OFF or NEXT.

    The Refuse raised by run_describe is `<what>:\n  <the command>\n<the
    engine's own stderr>`. Only the command echo is dropped, because it is
    already in the banner and it is the longest and least useful part.
    """
    raw = str(err).splitlines()
    if len(raw) > 2:
        kept = [raw[0]] + raw[2:]
    else:
        kept = raw
    kept = [ln.strip() for ln in kept if ln.strip()]
    if not kept:
        return "no message"
    return " ".join(kept)[:limit]


def valid_pid(pid):
    """The lock file's `pid` as an int, or None. Nothing else may be one.

    This value comes out of a JSON file inside a watched directory and ends up
    in a command the page tells the operator to paste on the pipeline host.
    `{"pid": "21877; rm -rf /"}` rendered, correctly escaped and perfectly
    legible, as `ps -p 21877; rm -rf / -o pid,etime,stat,args`. Escaping made
    it not XSS; it did not make it safe, because the target was the human. A
    value that is not a pid produces no command at all.
    """
    if isinstance(pid, bool):
        return None                    # bool is an int, and True is not a pid
    if isinstance(pid, float) and pid.is_integer():
        pid = int(pid)
    if isinstance(pid, str) and re.fullmatch(r"\s*[0-9]{1,10}\s*", pid):
        pid = int(pid)
    if isinstance(pid, int) and 0 < pid < 2 ** 31:
        return pid
    return None


def run_is_over(run):
    """True only when the record SAYS, in so many words, that the run ended.

    A missing or empty `final_status` is not that statement. Reading one as
    "Finished ? at ?" was the strongest possible claim in the one paragraph
    whose whole doctrine is that unprovable means alive - and it also switched
    off the staleness apparatus, so a partially written record produced a page
    asserting a finish beside a stage row counting up in RUN.
    """
    status = (run or {}).get("final_status")
    return bool(status) and status != "running"


# Buckets, in the order the index sorts them. Running first and never-run last,
# because with eight sequential datasets alphabetical order buries both the live
# one and the dead one.
# `error` is the console's own failure on one directory, which is a different
# fact from anything in that directory and sorts near the top because it is the
# one row on this page nobody can act on without being told about it.
BUCKETS = ["running", "failed", "error", "unreadable", "blank", "interrupted",
           "done", "unknown", "none"]
BUCKET_LABEL = {
    "running": "RUNNING", "failed": "FAILED", "unreadable": "UNREADABLE",
    "blank": "EMPTY STATE", "interrupted": "INTERRUPTED",
    "done": "DONE", "unknown": "UNCLEAR", "none": "NO RECORD",
    "error": "CONSOLE ERROR",
}
# What each bucket is called in the index's one-sentence summary. `unreadable`
# and `blank` exist because folding them into `none` made the front page say
# "never started" about a directory whose state file the console could not
# parse - including a job that is running right now and whose state file went
# bad before this console ever read it.
BUCKET_WORD = {
    "running": "running", "failed": "failed",
    "unreadable": "with an unreadable state file",
    "blank": "with an empty state file",
    "interrupted": "stopped without finishing", "done": "finished",
    "unknown": "unclear", "none": "never started",
    "error": "that this console could not read",
}


class Project:
    """One results directory, and the little memory the console keeps for it.

    That memory - the last state file that parsed, and the last `describe
    --config` for its effective config - lives here, in this process, and
    nowhere else. Nothing about a project is ever written down: closing the tab,
    restarting the console and opening cold on day three are one code path, and
    that property is free only as long as nothing is stored.
    """

    def __init__(self, index, path, contract):
        self.index = index
        self.path = os.path.abspath(path)
        self.contract = contract
        self.state_path = os.path.join(self.path, contract.state_name)
        self.log_path = os.path.join(self.path, contract.log_name)
        self.lock_path = os.path.join(self.path, contract.lock_name)
        self.effective_path = os.path.join(self.path, contract.effective_name)
        self._lock = threading.Lock()
        self._good = None            # last state file that parsed
        self._good_at = None
        self._bad_key = None         # identity of a state file already known bad
        self._bad_status = None      # and what it was found to be
        self._cfg_key = None         # (path, mtime, size) of the config probed
        self._cfg_run = None         # its config["run"], or None for "cannot tell"
        self._cfg_note = "not looked up yet"
        self._cfg_done = None        # set while a probe thread is running
        self.parts = PartCache()

    @property
    def name(self):
        """What the operator calls this dataset.

        The basename of eight results directories is `results` eight times over,
        which is the failure this avoids: the parent directory is the name a
        person actually uses.
        """
        base = os.path.basename(self.path)
        if base in ("results", "result", ""):
            parent = os.path.basename(os.path.dirname(self.path))
            if parent:
                return parent
        return base

    def read_state(self):
        """One snapshot of the state file, and the honest story about it.

        Read once, into bytes, and everything on the page derives from that one
        parse: reading it twice inside one request is how a header and the table
        under it end up describing two different moments.

        `cached` is True when this read failed and the snapshot is the last one
        that parsed. It was called `stale` until v0.5.0, and that word now
        belongs to a stage row: a STALE row is a `running` record left behind
        by a run that ended, which is the sense this file already used for a
        stale socket and a stale lock - something a dead process left. A
        snapshot the console is re-serving because the file would not parse
        THIS second is not that. It is a cache, and it says so.
        """
        # A state file that was there a moment ago and is not there now is
        # worth another look; one that has never been there is not.
        st = stat_or_none(self.state_path)
        ident = None if st is None else (st.st_dev, st.st_ino, st.st_mtime,
                                         st.st_size)
        with self._lock:
            known = self._bad_status if ident == self._bad_key else None
        if known is not None:
            # Not read at all, not merely read once: the short-circuit
            # shortened the retry loop and left the READ, so an 8 MB state
            # file already known corrupt was read whole on every poll of every
            # affected project. The identity has not moved, so neither have the
            # bytes: os.replace() gives a new inode, and a rewrite in place
            # moves the mtime or the size.
            status, obj, detail = known
        else:
            status, obj, detail = read_json_stable(
                self.state_path, STATE_MAX,
                retry_missing=self._good is not None)
        now = time.time()
        with self._lock:
            bad = status in ("bad", "empty")
            self._bad_key = ident if bad else None
            self._bad_status = (status, None, detail) if bad else None
            if status == "ok":
                self._good, self._good_at = obj, now
                return {"state": obj, "status": "ok", "detail": None,
                        "cached": False, "good_at": now}
            if self._good is not None:
                # Never blank the table over a transient read. The page says
                # which moment it is showing, in the footer, always.
                return {"state": self._good, "status": status,
                        "detail": detail, "cached": True,
                        "good_at": self._good_at}
        return {"state": {}, "status": status, "detail": detail,
                "cached": False, "good_at": None}

    def config_run(self, python, script, state):
        """`config["run"]` for this project, or None for "cannot tell".

        This is the difference between a stage that is OFF and one that simply
        has not been reached, and there are nine of the former in a default
        config - a table that called them all "not reached" would be wrong about
        nearly half its rows.

        The engine is the parser, always: the console shells out `describe
        --json --config <config.effective.yaml>`, which is the merged config
        this run actually used, defaults included, sitting inside the results
        directory. It is cached on (path, mtime, size), so a run that is not
        editing its config costs one fork for the life of the console. If the
        file is absent or the engine cannot read it, this returns None and every
        stage falls back to WAIT/NEXT - the console does not guess that a stage
        is disabled.

        THE FORK DOES NOT RUN ON THIS THREAD. Both candidates are files
        inside a watched directory handed to a YAML parser, and S_ISREG does
        not make a regular file quick: a 2.4 MB config.effective.yaml took the
        engine 61.6 s, and 33 concurrent requests against a page that waited
        for it 503'd every other project and /app.css. The regular-file check
        closes `--config /dev/zero` and `--config <a FIFO>`, which park a fork
        forever, and closes nothing else - the docstring used to claim
        otherwise. So the probe runs on a thread of its own under a global
        two-fork ceiling, one request thread waits CONFIG_WAIT for it, and
        until it answers this returns None, which every caller already renders
        as "cannot tell". The next poll has the answer.
        """
        st = regular_stat(self.effective_path)
        path, note = self.effective_path, None
        if st is None:
            # An older run, or one killed before the effective config was
            # written. `_run.config_path` names the file it was launched with -
            # today's contents, not March's, and the page says so.
            rec = state.get(self.contract.run_key) or {}
            cand = rec.get("config_path") if isinstance(rec, dict) else None
            if isinstance(cand, str) and cand:
                st, path = regular_stat(cand), cand
                note = ("read from the config file the run names, %s - which is "
                        "what it says today, not necessarily what it said when "
                        "the run started" % collapse_home(cand))
        if st is None:
            return None, ("no %s in this directory and no readable config "
                          "beside it, so enabled/off cannot be told apart from "
                          "not-yet-reached" % self.contract.effective_name)
        if st.st_size > CONFIG_MAX:
            return None, ("%s is %s, which is not a config this engine wrote; "
                          "it is not handed to the parser, so enabled/off "
                          "cannot be told apart from not-yet-reached"
                          % (collapse_home(path), human_bytes(st.st_size)))
        key = (path, st.st_mtime, st.st_size)
        with self._lock:
            if key == self._cfg_key:
                return self._cfg_run, self._cfg_note
            done = self._cfg_done
            first = done is None
            if first:
                done = self._cfg_done = threading.Event()
        if first:
            # One fork per project, not one per concurrent first request:
            # without this, N tabs opening at once fork N identical
            # `describe --config`.
            threading.Thread(target=self._probe_config, daemon=True,
                             name="cfg-%d" % self.index,
                             args=(python, script, path, key, note,
                                   done)).start()
            if done.wait(CONFIG_WAIT):
                with self._lock:
                    if key == self._cfg_key:
                        return self._cfg_run, self._cfg_note
        return None, ("the effective config (%s) is still being read; until "
                      "the engine answers, no stage below is called OFF"
                      % collapse_home(path))

    def _probe_config(self, python, script, path, key, note, done):
        """The fork, on a thread of its own. Never raises into the caller."""
        run = None
        try:
            with CONFIG_FORKS:
                doc = run_describe(python, script, config=path,
                                   timeout=CONFIG_TIMEOUT)
            run = doc.get("config", {}).get("run")
            if not isinstance(run, dict):
                run, note = None, "the engine reported no run: block"
            elif note is None:
                note = "read from %s" % collapse_home(path)
        except Refuse as e:
            run, note = None, "the engine could not read %s (%s)" % (
                collapse_home(path), one_line(e))
        except Exception as e:                       # never a silent thread
            run, note = None, "reading %s raised %r" % (collapse_home(path), e)
        finally:
            with self._lock:
                self._cfg_key, self._cfg_run, self._cfg_note = key, run, note
                self._cfg_done = None
            done.set()


def run_record(state, contract):
    rec = state.get(contract.run_key)
    return rec if isinstance(rec, dict) else None


def off_reason(key, val):
    """Why the engine will skip this stage, in the value's own terms.

    metaannot.py skips on `not cfg["run"].get(st["enabled"], False)` - plain
    truthiness - while this file used `is False`, so `run.pfam: 0`, `: null`,
    `: ""` and `: []` each made the engine skip a stage the console announced
    as NEXT. `run:\\n  pfam:` - a value forgotten in YAML - is the realistic
    form of that, and the console must not contradict the engine about what
    will run.
    """
    if val is False:
        shown = "false"
    elif val is None:
        shown = "empty (`run.%s:` with no value after it)" % key
    else:
        try:
            shown = json.dumps(val)
        except (TypeError, ValueError):
            shown = repr(val)
    return ("not enabled: run.%s is %s, and the engine runs a stage only when "
            "its flag is true" % (key, shown))


# What one stage's error string is allowed to cost the page. The engine caps
# its own at 500 characters; a state file need not have been written by this
# engine, and this string is rendered twice - the row and the pinned block.
ROW_ERROR_MAX = 500

# How many of the stages ranked ahead of a NEXT row are named in its sentence.
# The rest are counted, by name_list(): on a fresh run every dependency-free
# stage is ready at once, and a row that names nine of them is a directory
# listing where a sentence was wanted.
AHEAD_NAMED = 4


def ready_rows(rows):
    """The rows that are ready to start, when they can be RANKED at all.

    One predicate, because rank_ready() and the note under the table have to
    agree about it: a table that ranks its NEXT rows and a note that does not
    explain the ranking, or the other way round, is worse than neither.

    Empty when any ready row has no cost - which is what an engine older than
    the field gives for every stage, and what a stage added to a newer engine
    without one would give for itself. A partial sort over a missing key would
    rank a stage by a number this file made up, and the whole point of taking
    the order from `describe --json` is that it is the engine's order and not
    this file's guess.
    """
    ready = [r for r in rows if r["state"] == "next"]
    return [] if any(r["cost"] is None for r in ready) else ready


def next_detail(ahead):
    """What a NEXT row says, given what the scheduler ranks ahead of it.

    `ahead` is None when the contract cannot be ranked, [] when this is the
    stage the engine ranks first among those that are ready, and otherwise the
    names it ranks ahead of this one.

    "nothing is blocking it" was the whole sentence, and it is a claim about
    the dependency graph that a reader takes as a claim about the queue. On the
    455,571-protein run the page listed dbcan NEXT above ncbifam NEXT, the
    engine dispatched ncbifam, and dbcan did not get a worker for 27 hours. So
    the graph half is now said as the graph half, and the queue half is said
    separately - as a RANKING and nothing more. This console does not know how
    many workers are free or whether gpu_lease will defer a GPU stage, so it
    says which stage is ranked ahead of which and stops there; the note under
    the table carries the two unknowns once, rather than on every row.
    """
    if ahead is None:
        rank = ""
    elif not ahead:
        rank = (", and of the stages that are ready the engine ranks this one "
                "first")
    elif len(ahead) == 1:
        rank = ", but the engine ranks one other ready stage ahead of it: %s" \
               % ahead[0]
    else:
        rank = (", but the engine ranks %d of the stages that are ready ahead "
                "of it: %s" % (len(ahead), name_list(ahead, AHEAD_NAMED)))
    return ("no record yet — no dependency is blocking it%s. It has not been "
            "reached, or this run did not select it." % rank)


def rank_ready(rows):
    """Fill in which ready stages the engine ranks ahead of which.

    v0.4.0 made the scheduler dispatch each round's ready set longest-first -
    `run_now.sort(key=stage_priority, reverse=True)`, over cost 3 hours, 2
    minutes, 1 seconds - and describe --json emits `cost` for exactly this,
    "so a front end can order or annotate the table the same way". This file
    dropped the field, so the table went on listing its NEXT rows in stage
    order over a queue the engine was taking in a different one.

    The rows are annotated rather than reordered. The table is in the engine's
    own stage order - the strip on the index reads the same way, and a
    dependency graph read out of order is harder, not easier - so the order
    goes into the sentence instead.

    Stable, like the engine's own sort, and over the same sequence: stages of
    equal cost keep table order there, so they keep it here.
    """
    ahead = []
    for r in sorted(ready_rows(rows), key=lambda row: -row["cost"]):
        r["ahead"] = list(ahead)
        r["detail"] = next_detail(r["ahead"])
        ahead.append(r["name"])


def stage_rows(contract, state, cfg_run, now, with_outputs, project_path,
               parts=None, known=True):
    """All 21 stages, always.

    The state file records only what was reached, so a table built from it alone
    is mostly blank and says nothing about why. `describe --json` supplies the
    other rows, their dependency edges and their outputs, which turns "no
    record" into "waiting on interpro" or "not enabled: run.topology".

    `known=False` means the state file could not be read and there is no
    earlier good one. The vitals block said so at the top and then every row
    below it read "NEXT - no record yet - nothing is blocking it": a claim
    about progress derived from an empty dict, filling most of the screen,
    byte-identical to a directory nothing has ever run in. OFF survives it -
    that one comes from the config file, which WAS read.
    """
    rec_of, odd = {}, {}
    for st in contract.stages:
        r = state.get(st["name"])
        rec_of[st["name"]] = r if isinstance(r, dict) else None
        if r is not None and not isinstance(r, dict):
            # `"dbcan": "a string"` used to be indistinguishable from a stage
            # that was never reached, which reads as "nothing is wrong here".
            # The `bad` row type exists for exactly this.
            odd[st["name"]] = r
    run = run_record(state, contract)
    run_started = stamp_epoch(run.get("started")) if run else None
    # Which run a record belongs to is decided by ONE comparison, against the
    # run's start - so without a start there is no comparison and no answer.
    # `carried` was set only when there WAS one, and everything else therefore
    # defaulted to this run: a `_run` with no readable `started` announced an
    # eight-day-old failure as "1 stage has failed in this run" three lines
    # above the run record that says nothing of the sort. Undecidable is its
    # own answer here, and it is not "this run".
    era_known = run_started is not None
    over = run_is_over(run)
    over_when = (run or {}).get("finished") or "?"
    # Which stages this run turned off, decided once - and in two forms,
    # because two different questions are asked of it below.
    #
    # `off_set` answers "may this stage's OWN row read OFF?", and a stage with
    # a record is never called off whatever the config says now: it ran, and
    # the record is the evidence.
    #
    # `skip_set` answers the dependency question - "will the engine WAIT for
    # this one?" - and there the record is beside the point. decide() returns
    # `disabled (run.X)` for any stage whose flag is falsy, whatever the state
    # file holds, finish() then counts it as skipped and adds it to `done`, and
    # unmet_deps() lets it pass with "a disabled dependency is fine". Reading
    # `d in off_set` for that answered the wrong question by a record: a
    # dependency disabled in this config but carrying last week's `failed` or a
    # killed run's `running` fell through to `named` and its dependents were
    # reported WAIT, "waiting on foldseek (failed)", over a run the engine had
    # already walked straight past.
    off_set, off_why, skip_set = set(), {}, set()
    if isinstance(cfg_run, dict):
        for st in contract.stages:
            if st["enabled"] is None:
                continue
            val = cfg_run.get(st["enabled"], False)
            if val:
                continue
            skip_set.add(st["name"])
            if rec_of.get(st["name"]) is None:
                off_set.add(st["name"])
                off_why[st["name"]] = off_reason(st["enabled"], val)
    rows = []
    for st in contract.stages:
        name = st["name"]
        rec = rec_of[name]
        row = {"name": name, "gpu": st["gpu"], "empty_ok": st["empty_ok"],
               "deps": st["deps"], "enabled_key": st["enabled"],
               "state": "none", "label": "—", "took": "—", "detail": "",
               "carried": None, "era": None, "outputs": [], "part": None,
               "started": None, "error": None, "finished": None,
               # The scheduler's own rank, and which ready stages it puts
               # ahead of this one. rank_ready() fills `ahead` in after the
               # loop, when it is known which rows ended up ready at all.
               "cost": st["cost"], "ahead": []}
        status = rec.get("status") if rec else None
        if name in odd:
            row.update(state="bad", label="?",
                       detail="the state file's record for this stage is a %s, "
                              "not an object: %s"
                              % (type(odd[name]).__name__,
                                 repr(odd[name])[:120]))
        elif rec is not None and status not in ("running", "ok", "adopted",
                                                "failed"):
            row.update(state="bad", label="?",
                       detail="state record unreadable: %r" % (status,))
        elif status == "running":
            started = stamp_epoch(rec.get("started"))
            row.update(state="running", label="RUN", started=started,
                       took=dur(now - started) if started else "?",
                       detail="started %s" % clock(started))
        elif status == "ok":
            secs = finite(rec.get("seconds"), 0, SECONDS_MAX)
            row.update(state="ok", label="OK",
                       took=dur(secs) if secs is not None else "—",
                       detail="finished %s" % clip(rec.get("finished", "?"),
                                                   60))
        elif status == "adopted":
            row.update(state="adopted", label="ADOPT", took="—",
                       detail="signature matched, the existing output was "
                              "reused. A reuse records no duration, which is "
                              "why the time is a dash. Verified %s"
                              % rec.get("finished", "?"))
        elif status == "failed":
            err = clip(rec.get("error") or "no error recorded", ROW_ERROR_MAX)
            row.update(state="failed", label="FAIL", took="—", error=err,
                       finished=rec.get("finished") or None,
                       detail="%s — finished %s"
                              % (err, clip(rec.get("finished", "?"), 60)))
        elif name in off_set:
            row.update(state="off", label="OFF", detail=off_why[name])
        else:
            # A DISABLED dependency does not block anything, and saying it does
            # would be wrong about half the table on a default config. The
            # engine's own unmet_deps() puts it plainly - "a disabled
            # dependency is fine, its evidence is legitimately absent" - and
            # finish() adds a skipped stage to `done`, so the scheduler treats
            # it as satisfied. This follows the engine rather than the graph.
            named, silent, doneish, offish = [], [], 0, []
            for d in st["deps"]:
                dstat = (rec_of.get(d) or {}).get("status")
                if dstat in ("ok", "adopted"):
                    doneish += 1
                elif d in skip_set:
                    # skip_set, not off_set: whether the ENGINE waits for this
                    # dependency is decided by the config alone. A record only
                    # decides whether the dependency's own row may read OFF.
                    offish.append(d)
                elif dstat:
                    named.append("%s (%s)" % (d, dstat))
                else:
                    silent.append(d)
            if named or silent:
                # Past tense once the run is over, here too: "(no record yet)"
                # is the same present-tense claim one clause deeper.
                who = named + ([", ".join(silent)
                                + (" (no record)" if over else
                                   " (no record yet)")]
                               if silent else [])
                extra = []
                if doneish:
                    extra.append("%d of its %d dependencies %s done"
                                 % (doneish, len(st["deps"]),
                                    "is" if doneish == 1 else "are"))
                if offish:
                    extra.append("%s %s disabled and do%s not block it"
                                 % (", ".join(offish),
                                    "is" if len(offish) == 1 else "are",
                                    "es" if len(offish) == 1 else ""))
                waiting = "%s%s" % (", ".join(who),
                                    (". " + "; ".join(extra) + ".")
                                    if extra else "")
                # Once a run is over nothing is waiting on anything. The
                # heartbeat block already suppresses its whole staleness
                # apparatus for a finished run; a table underneath it still
                # saying "waiting on pfam, signalp" is the same bug one
                # function along, and it is the sentence a reader believes.
                if over:
                    row.update(state="none", label="—",
                               detail="no record. This run ended (%s) without "
                                      "reaching it; it was still short of %s"
                                      % (over_when, waiting))
                else:
                    row.update(state="wait", label="WAIT",
                               detail="no record yet — waiting on %s" % waiting)
            elif over:
                row.update(state="none", label="—",
                           detail="no record. This run ended (%s) without "
                                  "reaching it, and nothing was blocking it — "
                                  "it was not selected, or the run stopped "
                                  "first." % over_when)
            else:
                # Not "nothing is blocking it": nothing in the DEPENDENCY graph
                # is, which is a narrower claim than the one a reader takes
                # from a row labelled NEXT. rank_ready() below replaces this
                # sentence with the ranked one wherever the contract carries
                # enough to rank it.
                row.update(state="next", label="NEXT",
                           detail=next_detail(None))
        # A state file is cumulative across runs, so a stage disabled this time
        # keeps last week's green record. Presenting that as this run's success
        # is a lie by omission, and one timestamp comparison avoids it.
        if rec is not None:
            if not era_known:
                row["era"] = "unknown"
            else:
                row["era"] = "this"
                # `finished` for a record that reached an end; `started` for the
                # one kind that never does. mark_running() writes {signature,
                # status, started} and no `finished` at all, so a RUN record a
                # killed run left behind compared None against the run's start,
                # fell through to "this run", and rendered as THIS run's live
                # stage - with a duration counting up, for ever, from a clock
                # that stopped days ago. Every record carries one stamp or the
                # other, and the one it carries is the one to date it by.
                when = stamp_epoch(rec.get("finished"))
                if when is None and status == "running":
                    when = stamp_epoch(rec.get("started"))
                if when is not None and when < run_started:
                    row["carried"], row["era"] = clock(when), "earlier"
        # A running record that predates this run is not this run's live stage,
        # and it is not this console's business to say what it is instead. The
        # engine's own reading is in decide(): "finish() never ran, so the
        # writer was killed", and it recomputes the stage. This page may not go
        # that far - it never asserts that anything is dead - so it says the
        # one thing the two timestamps establish, drops the duration that was
        # counting from the older run's clock, and leaves the verdict where
        # every other verdict on this page is left.
        if row["state"] == "running" and row["era"] == "earlier":
            row.update(state="stale", label="STALE", took="—",
                       detail="recorded as running since %s, which is before "
                              "this run started (%s). The record was written "
                              "by an earlier run, so it is not evidence that "
                              "this stage is running now. When a run reaches a "
                              "stage in this state the engine reads it as an "
                              "interrupted write and recomputes it."
                              % (clock(row["started"]), clock(run_started)))
        if not known and row["state"] != "off":
            row.update(state="none", label="—", carried=None, era=None,
                       detail="the state file cannot be read, so whether this "
                              "stage has run is not known from here")
        # STALE is in the listing set and not in the scanning one, and the two
        # halves are decided by different questions. What a stale row leaves a
        # reader with is "did the killed run get anywhere before it died?", and
        # a stat of each declared output is the only evidence on this page that
        # answers it - the row was `running` before v0.5.0 named this case, it
        # listed its outputs then, and dropping the listing would take away the
        # one thing worth looking at. The part-file SCAN is the other question,
        # "is a tool writing bytes right now", which for a record from a dead
        # run is settled: nothing is. So a stale row is stat()ed and never
        # walked, which is also what keeps a killed esmfold from costing a
        # 455k-entry scandir on every poll for ever.
        if with_outputs and row["state"] in ("running", "stale", "ok",
                                             "adopted"):
            for rel in st["outputs"]:
                full = os.path.join(project_path, rel)
                fst = stat_or_none(full)
                row["outputs"].append({
                    "name": rel, "size": fst.st_size if fst else None,
                    "mtime": fst.st_mtime if fst else None,
                    "there": fst is not None})
                if row["state"] == "running" and row["part"] is None:
                    row["part"] = (parts.newest(full) if parts is not None
                                   else newest_part_file(full))
        rows.append(row)
    # After the loop, because which rows ended up ready is not known until
    # every row has been decided.
    rank_ready(rows)
    return rows


def stage_failures(rows, now, era="this"):
    """The stages of one era that have already failed, newest first.

    `carried` is the row's own record that the failure finished before this run
    started - stage_rows computes it for this reason and the table renders it
    as "from an earlier run". Without the filter, a failure from eight days ago
    was announced in red at the top of a healthy page as "1 stage has failed in
    this run", with the run's own start time on the next line. The state file
    is cumulative and adopted reuse is a headline feature, so a reused
    directory is metaannot's normal mode: a permanent red alarm over a healthy
    multi-day run is how a reader is trained past the real one.

    The index computed this and printed it in red; the project page - the one
    an operator leaves open - dropped it, so a run whose esmfold died two hours
    ago with `CUDA out of memory` presented a pinned block about run ids, a
    fresh heartbeat and a held lock, and said nothing at all about the failure.
    The FAIL row was three hundred pixels below the fold behind fifteen green
    ones. This is the fact the page exists to deliver, so it is lifted out
    here and rendered first.

    `era` is which of the answers the caller wants - "this" for the failures
    this run owns, "unknown" for the ones no run can be attached to because the
    run's start could not be read. They are rendered differently for the reason
    the filter exists at all: one is a claim about this run, and the other is a
    refusal to make one.
    """
    out = []
    for r in rows:
        if r["state"] != "failed" or r.get("era") != era:
            continue
        when = stamp_epoch(r["finished"])
        out.append({"name": r["name"], "error": r["error"] or "",
                    "finished": r["finished"], "epoch": when,
                    "age": None if when is None else max(0.0, now - when)})
    out.sort(key=lambda f: -(f["epoch"] or 0))
    return out


def era_note(run):
    """Why a stage record cannot be placed in or before this run, or None.

    Two shapes reach here, and both are real: a state file whose engine predates
    the run record, and a `_run` whose `started` is missing or is not the stamp
    this console parses.
    """
    if run is None:
        return ("There is no run record in this state file, so there is no run "
                "start to measure a stage's finish against.")
    if stamp_epoch(run.get("started")) is not None:
        return None
    raw = run.get("started")
    return ("This run record does not say when the run started (%s), so a "
            "stage's finish cannot be placed before or after it."
            % ("`started` is not in it" if raw is None else
               "`started` is %s, which is not a stamp this console can read"
               % clip(repr(raw), 40)))


def beat_of(run, now):
    """The advisory heartbeat as (age, heartbeat_s, band).

    One function, because the project page's staleness apparatus and the
    index's badge must not disagree about what "late" means - and the index did
    not read `last_seen` at all, so a project whose heartbeat died five hours
    ago was badged RUNNING with an empty annotation, byte-identical to one at
    13 s. `band` is fresh/late/long/skew/unrated/unknown/none, and not one of
    them is a verdict: `long` is "a long way past due", never "dead".

    `hb` is None when the record does not say, and then the band is `unrated`:
    an age, and no judgement on it. `finite(...) or 30.0` invented one instead,
    and the page stated the invention as this run's own fact - "This run stamps
    one every 30 s, so that is late" - identically for heartbeat_s of None, of
    Infinity, of "thirty", and of 0, which is the engine's own value for a
    heartbeat deliberately turned OFF. Three of those records say nothing about
    a cadence and the fourth says there is none; none of them says 30. The
    late/long bands are arithmetic on `hb`, so without one there is no late.
    """
    if not run or run_is_over(run):
        return None, None, "none"
    hb = finite(run.get("heartbeat_s"), 0.001, SECONDS_MAX)
    seen = epoch_of(run.get("last_seen_epoch"))
    if seen is None:
        return None, hb, "unknown"
    age = now - seen
    if age < -5:
        return age, hb, "skew"
    if hb is None:
        return age, None, "unrated"
    if age <= 3 * hb:
        return age, hb, "fresh"
    return age, hb, "late" if age <= 20 * hb else "long"


def heartbeat_view(run, now, newest_level, failures=(), state_ok=True,
                   state_status="ok"):
    """What the page says about whether anything is still happening.

    Three independent clocks exist and showing one of them would be dishonest:
    `last_seen` proves the heartbeat thread ran and the filesystem was
    writable, the log's mtime proves the run emitted output, and a part file's
    mtime proves the tool is writing bytes. This one reads `last_seen`;
    render_vitals puts the other two beside it, and part_sentence weighs them
    against each other rather than reciting them.

    When the run is over the whole staleness apparatus is suppressed. A finished
    run whose log last moved a day ago is not late; rendering an alarm there is
    the easy bug in this design and it is foreclosed rather than remembered.

    `failures` is passed in rather than derived because none of what is below
    is reassurance while a stage is dead, and a fresh heartbeat over a
    casualty is the exact shape that reads as "all well" to someone scanning.
    """
    v = {"running": False, "band": "none", "lines": [], "verdict": None,
         "beat_age": None, "heartbeat_s": None, "failures": list(failures)}
    if not run:
        if state_status == "missing" and state_ok:
            # state_trouble() returns None for a file that is not there, so
            # this branch led with "There is no run record in the state file"
            # - about a file that does not exist, whose second clause ("or it
            # was written by an older metaannot") is impossible, since an
            # older metaannot still writes a state file - while the accurate
            # sentence sat last, in muted small text.
            v["verdict"] = ("There is no state file in this directory, so "
                            "nothing has been recorded here. That is what a "
                            "directory nothing has run in looks like, and "
                            "also what one looks like in the seconds before "
                            "the first stage writes.")
            return v
        if not state_ok:
            # The state file could not be read at all. "Either nothing has run
            # here, or it was written by an older metaannot" are both wrong,
            # and they were the first thing on the page while the true reason
            # sat four lines down in muted small text.
            v["band"] = "trouble"
            v["verdict"] = ("The state file cannot be read, so there is no run "
                            "record to report — the reason is above. This is "
                            "not evidence that nothing has run here.")
            return v
        v["verdict"] = ("There is no run record in the state file. Either "
                        "nothing has run here, or it was written by a "
                        "metaannot older than the one that added the record.")
        return v
    final = run.get("finished")
    status = run.get("final_status")
    started = stamp_epoch(run.get("started"))
    if run_is_over(run):
        v["verdict"] = "Finished %s at %s%s." % (
            status, final or "?",
            " · ran %s" % dur(stamp_epoch(final) - started)
            if started and stamp_epoch(final) else "")
        return v
    v["running"] = True
    if not status:
        # Neither "running" nor a finish: an interrupted write, or an engine
        # that predates the key. Read as still going, and said out loud.
        v["lines"].append(
            "This run record carries no final_status, so whether the run ended "
            "was never written down. It is read here as still going, which is "
            "the engine's own rule — unprovable means alive — and it is why "
            "the ages below are still being reported.")
    age, hb, band = beat_of(run, now)
    seen = run.get("last_seen_epoch")
    v["heartbeat_s"] = hb
    if band == "unknown":
        v["band"] = "unknown"
        v["verdict"] = (
            "This run records no machine-readable heartbeat, so there is no "
            "age to report. The log's own mtime is the evidence below."
            if seen is None else
            "The heartbeat stamp in this run record is not a number this "
            "console can use (%s), so there is no age to report. The log's "
            "own mtime is the evidence below." % clip(repr(seen), 40))
        return v
    if band == "unrated":
        # An age with no cadence to judge it by. The two readings differ enough
        # to be worth telling apart: heartbeat_s 0 is the engine's documented
        # way of turning the heartbeat off, and last_seen is then the stamp
        # from when the run started and nothing more.
        v["beat_age"] = age
        v["band"] = "unrated"
        raw = run.get("heartbeat_s")
        if finite(raw, 0, SECONDS_MAX) == 0:
            v["verdict"] = (
                "Last heartbeat stamp %s ago. This run records heartbeat_s 0, "
                "which turns the heartbeat off, so that stamp is from when the "
                "run started and its age says nothing about whether the run is "
                "alive." % dur(age))
        else:
            v["verdict"] = (
                "Last heartbeat stamp %s ago. This run record does not say how "
                "often it stamps one (heartbeat_s is %s), so this console has "
                "no interval to call that late or on time — it will not invent "
                "one." % (dur(age), clip(repr(raw), 40)))
        v["lines"].append(
            "That leaves the log's own mtime below, and any part file, as the "
            "only clocks on this page with a meaning of their own — which is "
            "what to read here instead of a heartbeat age.")
        return v
    if band == "skew":
        # Two hosts sharing one filesystem across two clocks. Say so; do not
        # render a negative age and do not pretend it is fresh.
        v["band"] = "skew"
        v["verdict"] = ("The last heartbeat is stamped %s, which is in the "
                        "future for this machine. The run's clock and this "
                        "one disagree, so the age below cannot be trusted."
                        % clock(epoch_of(seen)))
        return v
    v["beat_age"] = age
    v["band"] = band
    if v["band"] == "fresh":
        v["verdict"] = ("Heartbeat stamped %s ago; this run stamps one every "
                        "%g s." % (dur(age), hb))
        if failures:
            # A fresh heartbeat directly above a dead stage is the shape that
            # reads as "all well" to someone scanning. It is not.
            v["lines"].append(
                "A fresh heartbeat means the process is alive and the "
                "filesystem is writable. It says nothing about the %s above."
                % ("failure" if len(failures) == 1 else "failures"))
        return v
    v["verdict"] = ("No heartbeat for %s. This run stamps one every %g s, so "
                    "that is %s. Last stamp %s."
                    % (dur(age), hb,
                       "late" if v["band"] == "late" else "a long way past due",
                       clock(seen)))
    v["lines"].append(
        "That is not evidence the run is dead, and this console will not say "
        "it is. One failed write ends the heartbeat thread while the run "
        "continues, and metaannot itself never treats a missing heartbeat as "
        "proof — its rule is that unprovable means alive.")
    if newest_level == "FATAL":
        v["lines"].append(
            "The newest log line in view is FATAL while the run record still "
            "says running. That is the killed-mid-write shape, and it is the "
            "one case where the silence is informative.")
    return v


def lock_view(project):
    status, obj, detail = read_json_stable(project.lock_path, LOCK_MAX)
    if status == "missing":
        return {"held": False,
                "text": "No lock file. Nothing holds this directory."}
    if status != "ok":
        return {"held": None, "text": "Lock file present but unreadable (%s)."
                                      % (detail or status)}
    pid, host = obj.get("pid"), obj.get("host") or "?"
    started = obj.get("started") or "?"
    # Nothing out of a file becomes part of a command string until it has been
    # validated as the type it claims to be. See valid_pid().
    num = valid_pid(pid)
    out = {
        "held": True, "pid": pid, "pid_ok": num is not None, "host": host,
        "started": started, "why_no_ps": None,
        "ps": "ps -p %d -o pid,etime,stat,args" % num if num else None,
        # The VALIDATED pid where there is one: `{"pid": "0021877"}` printed
        # "Held by pid 0021877" one clause before "ps -p 21877", which is two
        # different numbers on one line, in the block whose entire purpose is
        # handing over an exact command. Every field here is bounded, because
        # all three come out of a file in a watched directory.
        "text": "Held by pid %s on %s since %s."
                % (num if num is not None else clip(pid, 120),
                   clip(host, 120), clip(started, 60)),
    }
    if num is None:
        out["why_no_ps"] = (
            "The pid recorded in the lock file is not a number, so this "
            "console will not build a command out of it — the page would be "
            "handing you something to paste on the pipeline host, and it "
            "cannot vouch for what is in that field.")
    return out


def project_view(project, python, script, with_outputs=True, log=True):
    """Everything one project page shows, as plain data.

    The HTML renderer and /api/project both consume this dict, so there is one
    place where a fact is decided and two places that display it.
    """
    now = time.time()
    snap = project.read_state()
    state = snap["state"]
    run = run_record(state, project.contract)
    cfg_run, cfg_note = project.config_run(python, script, state)
    # Nothing parsed and there is no earlier read to fall back on: the table
    # below must not describe a directory out of an empty dict.
    known = not (snap["status"] in ("bad", "empty") and not snap["cached"])
    rows = stage_rows(project.contract, state, cfg_run, now, with_outputs,
                      project.path, project.parts, known=known)
    failures = stage_failures(rows, now)
    undated = stage_failures(rows, now, era="unknown")
    log_stat, log_why = stat_or_reason(project.log_path)
    state_stat = stat_or_none(project.state_path)
    trouble = state_trouble(snap, rows)
    unknown = sorted(k for k in state
                     if k != project.contract.run_key
                     and k not in project.contract.stage_names)
    tail, lines, dropped = None, [], 0
    if log:
        tail = tail_log(project.log_path)
        lines, dropped = log_lines(tail["text"])
        peek = lines
    else:
        # Still the last few kilobytes, and only those: the heartbeat block
        # says so when a FATAL line sits under a run record that still reads
        # "running", and a fragment that quietly dropped that sentence three
        # seconds after the page loaded would be worse than never showing it.
        peek = parse_log(tail_log(project.log_path,
                                  window=PEEK_WINDOW)["text"])
    newest_level = peek[-1]["level"] if peek else None
    view = {
        "index": project.index,
        "name": project.name,
        "path": project.path,
        "path_short": collapse_home(project.path),
        "taken": now,
        "state_status": snap["status"],
        "state_detail": snap["detail"],
        "state_cached": snap["cached"],
        "state_good_at": snap["good_at"],
        "state_mtime": state_stat.st_mtime if state_stat else None,
        "state_trouble": trouble,
        "unknown_stages": unknown,
        "log_mtime": log_stat.st_mtime if log_stat else None,
        "log_size": log_stat.st_size if log_stat else None,
        "log_why": log_why,
        "run": run,
        "cfg_run_known": isinstance(cfg_run, dict),
        "cfg_note": cfg_note,
        "rows": rows,
        "failures": failures,
        "undated_failures": undated,
        "era_note": era_note(run),
        "bucket": bucket_of(state, project.contract, rows, snap["status"],
                            snap["cached"]),
        "heartbeat": heartbeat_view(run, now, newest_level, failures,
                                    state_ok=trouble is None,
                                    state_status=snap["status"]),
        "state_known": known,
        "lock": lock_view(project),
        "counts": count_states(rows),
    }
    if log:
        view["log"] = {"ident": tail["ident"], "off": tail["off"],
                       "size": tail["size"],
                       "note": join_notes(tail["note"],
                                          dropped_note(dropped,
                                                       tail["reset"])),
                       "reset": tail["reset"], "lines": lines,
                       "dropped": dropped,
                       "warn": sum(1 for l in lines if l["level"] == "WARN"),
                       "fatal": sum(1 for l in lines
                                    if l["level"] in ("FATAL", "ERROR"))}
    return view


def state_trouble(snap, rows=()):
    """The one sentence that has to lead when the state file cannot be read.

    None when there is nothing wrong. When there IS, this is the strongest
    fact on the page and belongs at the top of the pinned block, not fourth in
    muted small text under two explanations that are both wrong.

    `rows` is passed because the sentence has to describe the table that is
    actually below it. An earlier version named the OFF rows unconditionally
    and so claimed a config had been read on a directory that has none, which
    is the same species of confident wrongness this function exists to remove.
    """
    if snap["cached"]:
        return ("The state file could not be read just now (%s), so this page "
                "is showing the last good read, from %s — not this moment. "
                "The table below is not blanked over a transient read."
                % (snap["detail"] or snap["status"], clock(snap["good_at"])))
    if snap["status"] == "bad":
        # "Nothing below is known about this directory" was not true of the
        # whole table, and the rows proved it: thirteen of them read OFF, from
        # the config file, which was read perfectly well - a different source
        # from the unreadable state file. Of the two, the ROWS are right, so
        # the sentence is the half that changes: it now claims exactly what it
        # can, which is that no stage's PROGRESS is known. Suppressing the OFF
        # rows to match the old sentence would have thrown away the one thing
        # about this directory the console can still establish.
        off = ("the rows marked OFF come from the config file, which was "
               "read, and nothing else in the table is evidence about this "
               "directory"
               if any(r.get("state") == "off" for r in rows) else
               "nothing in the table below is evidence about this directory")
        return ("The state file is unreadable — %s — and there is no earlier "
                "good read to fall back on. No stage's progress below is "
                "known: %s. This is a corrupt state file, not an empty "
                "directory." % (snap["detail"] or "no detail", off))
    if snap["status"] == "empty":
        return ("The state file is there but empty, and there is no earlier "
                "good read to fall back on. The engine writes it whole, by "
                "rename, so an empty one is a write that did not complete — "
                "not a directory where nothing has run.")
    return None


def count_states(rows):
    counts = {}
    for r in rows:
        counts[r["state"]] = counts.get(r["state"], 0) + 1
    return counts


def bucket_of(state, contract, rows, state_status="ok", cached=False):
    """Which group this project sorts into.

    Liveness first: a run that is still going and has already lost a stage is
    still running. The engine dispatches nothing new after the first failure,
    but it cannot interrupt the stages already in flight and waits for them,
    and an InterProScan that started an hour before the failure has hours left
    in it. The failure is not hidden by that - it rides along as a separate
    mark on the row - but burying a run that is still writing under FAILED
    would put the one directory that still needs watching in the wrong group.

    `state_status` is here because without it a state file that would not parse
    on the first read came out as NO RECORD and was counted under "never
    started" - the same badge, strip and sentence as a directory nothing has
    ever run in. The project page already worded it correctly; the index
    asserted the opposite of what the console knew, on the front page, for a
    job that might be running right now.
    """
    if not cached and state_status in ("bad", "empty"):
        return "unreadable" if state_status == "bad" else "blank"
    if not state:
        return "none"
    run = run_record(state, contract)
    if run:
        # run_is_over(), and not a fourth reading of `final_status`, because
        # the two pages have to answer this the same way. That function is
        # where the rule lives - a run record only ends a run by SAYING so -
        # and the project page's heartbeat block reads a record with no
        # final_status as still going, out loud, in the engine's own terms:
        # unprovable means alive. This branch enumerated the values the key can
        # HOLD and had no branch for its absence, so such a record fell through
        # to the row scan below and the index answered from whatever the stages
        # happened to look like: RUNNING while a leftover row still read
        # running, DONE once every stage was ok, UNCLEAR once v0.5.0 named the
        # leftover STALE. Three different groups for one shape, none of them
        # the sentence the project page was printing about it.
        if not run_is_over(run):
            return "running"
        final = run.get("final_status")
        if final == "failed":
            return "failed"
        if final == "ok":
            return "done"
        return "interrupted"
    # No run record at all - an engine older than the key, or a directory whose
    # first stage has not written one yet. Only here are the rows the evidence,
    # because only here is there nothing better.
    if any(r["state"] == "running" for r in rows):
        return "running"
    if any(r["state"] == "failed" for r in rows):
        return "failed"
    if any(r["state"] in ("ok", "adopted") for r in rows):
        return "done"
    return "unknown"


def index_item(p, python, script, now):
    """One project's row on the index. Two stats and one small JSON parse.

    Split out of index_view for the guard around it: everything that can go
    wrong for ONE project goes wrong inside this function, where it can be
    caught per project rather than taking the page down for all of them.
    """
    snap = p.read_state()
    state = snap["state"]
    cfg_run, _ = p.config_run(python, script, state)
    known = not (snap["status"] in ("bad", "empty") and not snap["cached"])
    rows = stage_rows(p.contract, state, cfg_run, now, False, p.path,
                      known=known)
    log_stat = stat_or_none(p.log_path)
    state_stat = stat_or_none(p.state_path)
    moved, source, moved_size = None, "nothing here", 0
    for st, what in ((log_stat, p.contract.log_name),
                     (state_stat, p.contract.state_name)):
        if st is not None and (moved is None or st.st_mtime > moved):
            moved, source, moved_size = st.st_mtime, what, st.st_size
    run = run_record(state, p.contract)
    bucket = bucket_of(state, p.contract, rows, snap["status"],
                       snap["cached"])
    beat_age, heartbeat_s, beat_band = beat_of(run, now)
    return {
        "index": p.index, "name": p.name, "path": p.path,
        "path_short": collapse_home(p.path),
        "path_elided": elide_path(collapse_home(p.path)),
        "bucket": bucket,
        "moved": moved, "moved_age": None if moved is None else max(0.0, now - moved),
        "moved_source": source,
        # "Newest activity" naming a directory somebody merely touched is
        # not activity in any useful sense.
        "evidence": bucket != "none" or moved_size > 0,
        "strip": [{"name": r["name"], "state": r["state"],
                   "label": r["label"]} for r in rows],
        "counts": count_states(rows),
        # This run's casualties and an earlier run's are two different
        # sentences, on this page as on the other one.
        "failed": [r["name"] for r in rows
                   if r["state"] == "failed" and r.get("era") == "this"],
        "carried_failed": [r["name"] for r in rows
                           if r["state"] == "failed"
                           and r.get("era") == "earlier"],
        # Neither this run's nor an earlier one's, because the record does not
        # establish which. The index says so rather than picking.
        "undated_failed": [r["name"] for r in rows
                           if r["state"] == "failed"
                           and r.get("era") == "unknown"],
        "beat_age": beat_age, "beat_band": beat_band,
        "heartbeat_s": heartbeat_s,
        "state_status": snap["status"], "state_cached": snap["cached"],
        "state_detail": snap["detail"],
        "run_status": (run or {}).get("final_status"),
        "run_started": (run or {}).get("started"),
    }


def broken_item(p, exc):
    """The row for a project whose facts could not be assembled at all.

    Every field the sort, the summary and the renderer read, and not one of
    them a claim: the bucket is its own, the strip is empty, and the detail is
    the exception. Deriving anything else here would be inventing it.

    Nothing in here may raise, because it runs INSIDE the handler that exists so
    that nothing raises: `name` is a property, and a property that throws in the
    recovery path is the same 500 one function later.
    """
    try:
        path, name = p.path, p.name
    except Exception:
        path = name = "?"
    return {"index": getattr(p, "index", -1), "name": name,
            "path": path, "path_short": collapse_home(path),
            "path_elided": elide_path(collapse_home(path)),
            "bucket": "error", "moved": None, "moved_age": None,
            "moved_source": "not read", "evidence": False,
            "strip": [], "counts": {}, "failed": [], "carried_failed": [],
            "undated_failed": [],
            "beat_age": None, "beat_band": "none", "heartbeat_s": None,
            "state_status": "error", "state_cached": False,
            "state_detail": one_line(repr(exc), 200),
            "run_status": None, "run_started": None}


def index_view(projects, python, script, scan=None):
    """The project list, and the guarantee that one project cannot take it away.

    The loop had no per-project guard, so anything raised for ONE directory -
    a path string os.stat will not accept, a shape nothing here anticipated -
    was a 500 on `/` and on /api/projects for EVERY watched project, including
    the seven that were perfectly readable. The front page of a watcher is the
    last thing that should be all-or-nothing, so a project that cannot be read
    becomes a row that says so and the other rows are served.
    """
    now = time.time()
    items = []
    for p in projects:
        try:
            items.append(index_item(p, python, script, now))
        except Exception as e:                       # never a silent swallow
            sys.stderr.write("console: %r while reading %s for the index\n"
                             % (e, getattr(p, "path", "?")))
            items.append(broken_item(p, e))
    # Bucket, then whatever is WRONG inside the bucket, then recency. The key
    # was (bucket, -moved), which sorts a running project by how recently its
    # files moved - so the one that stopped moving five hours ago, the only one
    # in the group that needs anything, sorted to the BOTTOM of it.
    items.sort(key=index_order)
    view = {"taken": now, "projects": items, "summary": summarise(items, now)}
    view.update(scan or {"roots": [], "scanned_at": None, "rescan_s": 0})
    return view


STALE_BANDS = ("late", "long")


def index_order(item):
    hurt = 1 if item["beat_band"] in STALE_BANDS + ("skew",) else 0
    return (BUCKETS.index(item["bucket"]), -hurt,
            -(item["beat_age"] or 0) if hurt else -(item["moved"] or 0))


def name_list(names, cap=3):
    """Names for a one-line summary, capped. Uncapped, summarise() degenerated
    into "Fail8, Fail5, OneBigErr, Messy, Cohort_B, Reused_dir have a failed
    stage but are still going" - a directory listing with a verb in it."""
    if len(names) <= cap:
        return ", ".join(names)
    return "%s and %d more" % (", ".join(names[:cap]), len(names) - cap)


def summarise(items, now):
    """The tmux replacement, in one sentence.

    It earns its place because `02_run.sh` pipes each dataset through `tail
    -80`: a failure in dataset three at 02:00 scrolls off the screen and is
    invisible until someone thinks to scroll back. Here it is permanent.
    """
    if not items:
        return "No results directories are being watched."
    counts = {}
    for i in items:
        counts[i["bucket"]] = counts.get(i["bucket"], 0) + 1
    parts = ["%d %s" % (counts[b], BUCKET_WORD[b])
             for b in BUCKETS if b in counts]
    hurt = [i["name"] for i in items if i["failed"] and i["bucket"] != "failed"]
    if hurt:
        parts.append("%s %s a failed stage but %s still going"
                     % (name_list(hurt), "has" if len(hurt) == 1 else "have",
                        "is" if len(hurt) == 1 else "are"))
    # The sentence the docstring calls "the tmux replacement" said "10 running"
    # over a project whose heartbeat had been silent for 5h18m. Staleness, in
    # the summary, in the same words the project page uses - and still not a
    # verdict: "has not stamped a heartbeat", never "has died".
    quiet = sorted((i for i in items if i["beat_band"] in STALE_BANDS),
                   key=lambda i: -(i["beat_age"] or 0))
    if quiet:
        parts.append("%s %s not stamped a heartbeat for %s%s"
                     % (name_list([i["name"] for i in quiet]),
                        "has" if len(quiet) == 1 else "have",
                        "" if len(quiet) == 1 else "up to ",
                        dur(quiet[0]["beat_age"])))
    newest = max((i for i in items if i["moved"] and i["evidence"]),
                 default=None, key=lambda i: i["moved"])
    tail = ""
    if newest:
        tail = " — newest activity %s ago (%s)" % (
            dur(now - newest["moved"]), newest["name"])
    return " · ".join(parts) + tail


# ======================================================================
# 4. the page - strings in this file, because there is no second asset
#    to scp and no CDN a lab server can reach
# ======================================================================

PAGE_CSS = """
:root {
  --bg:#f6f6f4; --panel:#ffffff; --ink:#16181d; --dim:#5d6470; --line:#dfe1e6;
  --run:#a05a00; --run-bg:#fdf3e2; --ok:#1d6b3a; --ok-bg:#e8f4ec;
  --fail:#a52020; --fail-bg:#fbeaea; --wait:#5d6470; --off:#9aa1ad;
  --warn:#8a6100; --accent:#2b4d8a;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg:#14161a; --panel:#1b1e24; --ink:#e6e8ec; --dim:#98a0ad; --line:#2c313a;
    --run:#e8a33d; --run-bg:#2e2410; --ok:#5fca8a; --ok-bg:#12281c;
    --fail:#ff8080; --fail-bg:#2e1414; --wait:#98a0ad; --off:#5b6270;
    --warn:#e0b34a; --accent:#89b4f8;
  }
}
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }
a { color:var(--accent); }
.wrap { max-width:1180px; margin:0 auto; padding:18px 20px 60px; }
header.top { display:flex; align-items:baseline; gap:14px; flex-wrap:wrap;
  border-bottom:1px solid var(--line); padding-bottom:10px; margin-bottom:16px; }
header.top h1 { font-size:16px; margin:0; font-weight:600; letter-spacing:.2px; }
.ro { font-size:11px; text-transform:uppercase; letter-spacing:.8px;
  color:var(--ok); border:1px solid var(--ok); border-radius:3px; padding:1px 6px; }
.muted { color:var(--dim); }
.mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
.small { font-size:12px; }
.summary { font-size:15px; margin:0 0 16px; }
.panel { background:var(--panel); border:1px solid var(--line); border-radius:6px;
  padding:14px 16px; margin-bottom:16px; }
table { border-collapse:collapse; width:100%; table-layout:fixed; }
th { text-align:left; font-size:11px; text-transform:uppercase;
  letter-spacing:.6px; color:var(--dim); font-weight:600;
  border-bottom:1px solid var(--line); padding:6px 8px; }
td { padding:7px 8px; border-bottom:1px solid var(--line);
  vertical-align:top; overflow-wrap:anywhere; }
tr:last-child td { border-bottom:none; }
.badge { display:inline-block; font-size:11px; font-weight:700;
  letter-spacing:.5px; padding:1px 6px; border-radius:3px; border:1px solid; }
.b-running { color:var(--run); border-color:var(--run); background:var(--run-bg); }
.b-done { color:var(--ok); border-color:var(--ok); background:var(--ok-bg); }
.b-failed { color:var(--fail); border-color:var(--fail); background:var(--fail-bg); }
.b-interrupted, .b-unknown, .b-blank { color:var(--warn); border-color:var(--warn); }
.b-unreadable { color:var(--fail); border-color:var(--fail); background:var(--fail-bg); }
.b-none { color:var(--off); border-color:var(--off); }
.b-error { color:var(--fail); border-color:var(--fail); background:var(--fail-bg); }
.strip { display:flex; gap:2px; }
.cell { width:6px; height:15px; border-radius:1px; background:var(--off);
  opacity:.35; }
.c-ok, .c-adopted { background:var(--ok); opacity:1; }
.c-running { background:var(--run); opacity:1; }
/* A RUN record left behind by an earlier run. Deliberately the running colour
   at the not-reached opacity: it is what a running record looks like, and it
   is not this run's. Colour is never the only channel here - the cell's title
   and the State column both read STALE. */
.c-stale { background:var(--run); opacity:.3; }
.c-failed { background:var(--fail); opacity:1; }
.c-wait, .c-next, .c-none { background:var(--wait); opacity:.3; }
.c-bad { background:var(--warn); opacity:1; }
.c-off { background:none; opacity:1; display:flex; align-items:center; }
.c-off::after { content:""; display:block; width:6px; height:1px;
  background:var(--off); }
.st { font-weight:700; font-size:11px; letter-spacing:.4px; }
.s-running { color:var(--run); } .s-ok, .s-adopted { color:var(--ok); }
.s-failed { color:var(--fail); } .s-off { color:var(--off); }
.s-wait, .s-next, .s-none { color:var(--dim); } .s-bad { color:var(--warn); }
.s-stale { color:var(--warn); }
tr.r-failed { background:var(--fail-bg); }
tr.r-running { background:var(--run-bg); }
/* Not r-running's background: the row is not this run's live stage, and the
   one thing it must not do is look like one. */
tr.r-stale td { color:var(--dim); }
tr.r-off td { color:var(--off); }
.chip { font-size:10px; border:1px solid var(--line); border-radius:3px;
  padding:0 4px; color:var(--dim); margin-left:5px; white-space:nowrap; }
.vitals h2 { font-size:15px; margin:0 0 6px; }
.vitals p { margin:6px 0; }
/* No max-height and no overflow on the pinned block, and that is the fix, not
   an omission. `position:sticky; top:0` plus `max-height:45vh; overflow:auto`
   put everything past the cap into a nested scroller the PAGE cannot scroll:
   measured at 1280x720 with three failures, .dead was 262 px inside a .vitals
   of clientHeight 322, and the heartbeat verdict was off the bottom of that
   box at every scroll position of the document. A cap that hides the answer is
   worse than the tall block it was capping.
   What keeps this short is what goes IN it: only the three urgent facts - the
   failures, anything wrong with the state file, and the verdict - each already
   bounded in Python (BLOCK_FAILS, BLOCK_ERROR_MAX, clip()). Everything else
   the block used to carry - the run identifiers, the explanations, the log,
   part and lock lines - is rendered below it in ordinary flow, where it may be
   as tall as it likes and the page can scroll to it. See render_vitals. */
.vitals-more p { margin:6px 0; }
.sticky { position:sticky; top:0; z-index:5; }
.band-late, .band-trouble { border-left:4px solid var(--warn); }
.band-long, .band-skew, .band-fail { border-left:4px solid var(--fail); }
.band-fresh { border-left:4px solid var(--ok); }
/* The failure block, first in the pinned vitals. A run whose esmfold died two
   hours ago must not present as a fresh heartbeat and a held lock. */
.dead { background:var(--fail-bg); border:1px solid var(--fail);
  border-radius:4px; padding:8px 10px; margin:0 0 10px; }
.dead h3 { margin:0 0 4px; font-size:14px; color:var(--fail); }
.dead p { margin:3px 0; }
.dead .err { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:12px; overflow-wrap:anywhere; }
.trouble { background:var(--run-bg); border:1px solid var(--warn);
  border-radius:4px; padding:8px 10px; margin:0 0 10px; }
.legend { display:flex; gap:12px; align-items:center; flex-wrap:wrap;
  font-size:12px; color:var(--dim); margin-top:10px; }
.legend span.cell { margin-right:4px; }
.legend i { font-style:normal; display:inline-flex; align-items:center; }
.logpane { background:var(--panel); border:1px solid var(--line);
  border-radius:6px; }
.loghead { display:flex; gap:12px; align-items:center; flex-wrap:wrap;
  padding:8px 12px; border-bottom:1px solid var(--line); font-size:12px; }
#logbody { max-height:44vh; overflow:auto; padding:6px 0;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:12px; line-height:1.45; }
.ln { padding:0 12px 0 9px; border-left:3px solid transparent;
  white-space:pre-wrap; overflow-wrap:anywhere; }
.ln.l-WARN { border-left-color:var(--warn); background:var(--run-bg); }
.ln.l-FATAL, .ln.l-ERROR { border-left-color:var(--fail); background:var(--fail-bg); }
.marker { padding:4px 12px; color:var(--dim); font-style:italic; }
[hidden] { display:none !important; }
.newbtn { position:sticky; bottom:6px; margin:0 12px 6px; display:block;
  background:var(--accent); color:#fff; border:none; border-radius:4px;
  padding:3px 10px; font-size:12px; cursor:pointer; }
footer { color:var(--dim); font-size:12px; margin-top:22px;
  border-top:1px solid var(--line); padding-top:10px; }
.row-name { font-weight:600; }
/* Column widths live here, not in a style attribute: the page is served under
   Content-Security-Policy default-src 'self', which blocks inline styles - and
   blocks them silently, so the table simply came out with five equal columns.
   Nothing in this console emits a style attribute; test_console pins that. */
.c-proj { width:32%; } .c-badge { width:16%; } .c-strip { width:22%; }
.c-stage { width:12%; } .c-state { width:7%; } .c-took { width:9%; }
.c-out { width:24%; }
.note { margin-top:8px; }
.failmark { color:var(--fail); }
.latemark { color:var(--warn); }
/* No text-overflow here: the paths share their whole left-hand side, so the
   ellipsis cut off the only segment that told one row from another. They are
   shortened from the left, in Python, by elide_path(). */
.path { display:block; overflow:hidden; white-space:nowrap; }
code { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:12px; }
"""

PAGE_JS = """
// Refresh, without ever reloading the page.
//
// Two independent loops, and they do not overlap. The server renders every
// panel, so no fact is decided twice: the first loop swaps the already-rendered
// HTML of the vitals-and-table region, and the second appends to the log pane by
// byte offset. The log pane is not inside the swapped region, because two
// writers over one pane would show every new line twice and throw away the
// operator's scroll position each time - and keeping that scroll is the one
// interaction detail in M1 that decides whether this is usable.
(function () {
  var box = document.getElementById('live');
  if (!box) return;
  var frag = box.getAttribute('data-frag');
  var floor = parseFloat(box.getAttribute('data-interval') || '3');
  var quiet = 0, timer = null;

  // Both of these live on #livestate, which is INSIDE the fragment, because
  // the fragment is #live's innerHTML and so can never replace an attribute of
  // #live itself. They used to be on #live, and the page therefore re-read its
  // load-time values forever: a tab opened before a run started polled every
  // 60 s through the whole run, and one opened during it kept polling every
  // 3 s hours after it ended.
  function marker(attr) {
    var m = box.querySelector('#livestate');
    return m ? m.getAttribute(attr) : null;
  }
  var running = marker('data-running') !== '0';
  var alarm = marker('data-alert') === '1';
  // The back-off compared the fetched text against the DOM's re-serialised
  // innerHTML, which differs on every poll because the fragment carries
  // relative ages ("last grew 7m49s"). So `quiet` never incremented and the
  // 15 s branch below was unreachable. data-sig is the server's digest of the
  // FACTS - mtimes, sizes, stage states, bands - with no clock in it.
  var sig = marker('data-sig');

  // Slow down on a run that is over, and on one that has gone quiet - but
  // never while anything is wrong. What the engine actually does is rewrite
  // the state file for every 30 s heartbeat, so data-sig flips every 30 s and
  // `quiet` cannot climb past about ten on a healthy live run: the 15 s branch
  // is really the branch for a run whose heartbeat has stopped, which is
  // exactly when you are watching for it to come back. data-alert pins the
  // page at the floor whenever the server says something is wrong.
  function period() {
    if (!running) return Math.max(floor, 60);
    if (alarm) return floor;
    return quiet > (120 / floor) ? Math.max(floor, 15) : floor;
  }

  // Do not fetch at all while the reader has something selected in here: the
  // region is replaced by innerHTML every three seconds, so selecting a CUDA
  // error string to paste into a search was wiped on the next poll. Ages
  // frozen while a selection is held is the cheaper of the two costs.
  function selecting() {
    var s = window.getSelection && window.getSelection();
    return !!(s && !s.isCollapsed && s.rangeCount &&
              box.contains(s.getRangeAt(0).commonAncestorContainer));
  }

  function schedule(fn) {
    clearTimeout(timer);
    if (document.visibilityState === 'hidden') return;   // a lidded laptop
    timer = setTimeout(fn, period() * 1000);             // fires no requests
  }

  function pull() {
    if (selecting()) { schedule(pull); return; }
    fetch(frag, {cache: 'no-store'}).then(function (r) { return r.text(); })
      .then(function (t) {
        box.innerHTML = t;              // the ages are live; always swap
        var s = marker('data-sig');
        if (s !== null && s === sig) { quiet++; } else { quiet = 0; }
        sig = s;
        running = marker('data-running') !== '0';
        alarm = marker('data-alert') === '1';
      }).catch(function () {}).then(function () { schedule(pull); });
  }

  document.addEventListener('visibilitychange', function () {
    if (document.visibilityState === 'visible') pull(); else clearTimeout(timer);
  });
  schedule(pull);

  // ---- the log, appended by byte offset -----------------------------
  var pane = document.getElementById('logbody');
  if (!pane) return;
  var btn = document.getElementById('lognew'), unseen = 0, logTimer = null;

  function atBottom() {
    return pane.scrollHeight - pane.scrollTop - pane.clientHeight < 40;
  }
  function el(cls, text) {
    var d = document.createElement('div');
    d.className = cls; d.textContent = text; return d;
  }
  // `hidden`, not style.display: the page is served under a strict CSP, and
  // keeping every style decision in app.css is what makes that survivable.
  if (btn) btn.onclick = function () {
    pane.scrollTop = pane.scrollHeight; unseen = 0; btn.hidden = true;
  };

  // ---- and bounded, because the tab is meant to stay open for days --
  // Nothing trimmed this pane. Every line /api/log delivered became a <div>
  // that stayed, on a page the run book tells you to leave open for the length
  // of a three-day run, so a stage emitting a few lines a second put a quarter
  // of a million nodes in one scroller. The cap is the server's (data-max), and
  // it is announced rather than applied quietly: a pane holding the last two
  // thousand lines looks exactly like a pane holding the whole log, and a
  // reader who believes the second one will conclude a stage never logged
  // something it logged four hours ago.
  //
  // The || is for a page cached from a console that served no data-max; a
  // bound this script picked is still a bound, and an unbounded pane is the
  // one outcome that is not acceptable.
  var cap = parseInt(pane.getAttribute('data-max'), 10) || 2000;
  var dropped = 0, capline = null;
  // LINES, counted, and not the pane's children. Markers share the scroller
  // with the lines - one per rotation, and the "no log lines" placeholder the
  // first render leaves behind - and bounding childElementCount gave each of
  // them a slot out of the cap. Subtracting the cap's own notice by hand was
  // the tell: it fixed the one marker whose arithmetic was obvious and left
  // the rest, so a pane that had rotated twice held 1,998 lines under a
  // sentence saying it keeps the last 2,000. The number the page states has
  // to be the number the page keeps. It starts at what the server sent.
  var shown = pane.querySelectorAll('div[class^="ln"]').length;

  function trim() {
    while (shown > cap) {
      var first = pane.firstElementChild;
      if (first === capline) first = capline.nextElementSibling;
      if (!first) break;
      // A marker is not a log line: dropping one is counted as no line and
      // frees no slot. It goes because the lines it was about have gone.
      if (first.className.indexOf('ln') === 0) { dropped++; shown--; }
      pane.removeChild(first);
    }
    // "N new lines below" is an offer to scroll down to them, so it may not
    // name a line the pane no longer holds. A burst bigger than the cap
    // between two polls trims the very lines it just counted, and the button
    // was offering 2,300 over a pane holding 2,000.
    if (unseen > shown) unseen = shown;
    if (!dropped) return;
    if (!capline) {
      capline = el('marker', '');
      pane.insertBefore(capline, pane.firstChild);
    }
    capline.textContent = '— ' + dropped + ' earlier line(s) dropped from ' +
      'this pane, which keeps the last ' + cap + '. The log file named above ' +
      'is complete; this pane is not —';
  }

  function poll() {
    var u = '/api/log?p=' + encodeURIComponent(pane.getAttribute('data-project')) +
            '&off=' + pane.getAttribute('data-off') +
            '&id=' + encodeURIComponent(pane.getAttribute('data-ident') || '');
    fetch(u, {cache: 'no-store'}).then(function (r) { return r.json(); })
      .then(function (d) {
        if (!d.ok) return;
        // Stick to the bottom only if you were already there. Scroll up to read
        // something and the pane freezes until you ask for the new lines.
        var stick = atBottom();
        // The pane was emptied because the log rotated or was truncated, and
        // d.note says so. Nothing was dropped BY THE CAP, so the cap's own
        // notice goes with the lines it was counting.
        if (d.reset) {
          pane.innerHTML = '';
          unseen = 0; dropped = 0; capline = null; shown = 0;
        }
        if (d.note) pane.appendChild(el('marker', '— ' + d.note + ' —'));
        d.lines.forEach(function (l) {
          pane.appendChild(el('ln l-' + l.level, l.text)); unseen++; shown++;
        });
        if (d.lines.length) {
          var tall = pane.scrollHeight;
          trim();
          if (stick) {
            pane.scrollTop = pane.scrollHeight; unseen = 0;
            if (btn) btn.hidden = true;
          } else {
            // Trimming takes lines off the TOP, above the reader, so the
            // document shortens under them and the line they were reading
            // would slide up the pane. Give back exactly what was taken.
            pane.scrollTop -= (tall - pane.scrollHeight);
            if (btn) {
              btn.hidden = false;
              btn.textContent = unseen + ' new lines below';
            }
          }
        }
        pane.setAttribute('data-off', d.off);
        pane.setAttribute('data-ident', d.ident || '');
      }).catch(function () {})
      .then(function () {
        clearTimeout(logTimer);
        logTimer = setTimeout(poll, (document.visibilityState === 'hidden'
                                     ? 60 : period()) * 1000);
      });
  }
  pane.scrollTop = pane.scrollHeight;
  logTimer = setTimeout(poll, floor * 1000);
})();
"""


# C0 and C1 controls, minus tab and newline, plus the two line separators that
# are not line separators here.
CONTROL_RE = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f\x7f-\x9f"
                        "\u2028\u2029]")


def esc(text):
    """HTML-escape, and neutralise what is underneath the markup.

    Escaping is correct for markup and does nothing about the bytes that are
    not: NUL and the ANSI escapes a tool writes into its own stderr went into
    the served body raw, invisible in the page and live for whatever
    `curl | less` hands them to. U+2028 and U+2029 for text_lines()' reason.
    """
    return CONTROL_RE.sub("\ufffd",
                          html.escape("" if text is None else str(text),
                                      quote=True))


def page(title, body):
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>%s</title><link rel="stylesheet" href="/app.css"></head>'
        '<body><div class="wrap">%s</div>'
        '<script src="/app.js"></script></body></html>' % (esc(title), body))


def top_bar(contract, extra=""):
    return (
        '<header class="top"><h1>metaannot console</h1>'
        '<span class="ro">read-only</span>'
        '<span class="muted small">console %s · engine %s · describe schema %s '
        '· signatures v%s</span>%s</header>'
        % (esc(CONSOLE_VERSION), esc(contract.version),
           esc(contract.describe_version), esc(contract.signature_version),
           extra))


def render_strip(strip):
    """One 6x15 cell per stage, in the engine's own order.

    Colour is never the only channel: the State column carries a word, every
    cell carries a title, and the strip carries the whole reading as an
    aria-label. `off` is drawn as a hairline rather than a grey box, so
    "excluded by config" is visually a different kind of thing from "not yet
    reached".
    """
    cells = "".join(
        '<span class="cell c-%s" title="%s: %s"></span>'
        % (esc(c["state"]), esc(c["name"]), esc(c["label"])) for c in strip)
    words = ", ".join("%s %s" % (c["name"], c["label"]) for c in strip)
    return '<span class="strip" role="img" aria-label="%s">%s</span>' % (
        esc(words), cells)


LEGEND = (("ok", "finished"), ("running", "running"), ("failed", "failed"),
          ("adopted", "reused"), ("bad", "unreadable record"),
          ("stale", "running record from an earlier run"),
          ("wait", "not reached"), ("off", "off in the config"))


def render_legend():
    """What the 21 cells mean, in words, without hovering for it."""
    return ('<div class="legend">%s</div>' % "".join(
        '<i><span class="cell c-%s"></span>%s</i>' % (esc(k), esc(word))
        for k, word in LEGEND))


def live_marker(running, sig, alert=False):
    """The three facts the poller needs, INSIDE the swapped fragment.

    They were attributes of #live, which innerHTML replacement never touches,
    so neither the cadence nor the back-off in app.js did anything at all.

    `alert` exists because the back-off could not fire on a live run and could
    fire only on a dead one: the engine rewrites the state file for every 30 s
    heartbeat, so data-sig flips every 30 s and `quiet` never reaches the 40
    unchanged polls the 15 s branch needs - while the moment the heartbeat
    STOPS, quiet climbs freely and the page slows down, which is the inverse of
    what you want. Anything wrong pins the page at the floor.
    """
    return ('<span id="livestate" hidden data-running="%d" data-alert="%d" '
            'data-sig="%s"></span>'
            % (1 if running else 0, 1 if alert else 0, esc(sig)))


def index_sig(view):
    return digest("%s:%s:%s:%s" % (i["index"], i["bucket"], i["moved"],
                                   i["beat_band"])
                  for i in view["projects"])


def digest(parts):
    """A short, stable hash of the facts a fragment carries.

    A hash rather than the facts themselves, because some of those facts are
    strings out of a file in a watched directory - the lock's own `Held by pid
    ... on ...` line among them - and there is no reason to echo attacker-
    chosen bytes into an attribute nobody reads. Not a security control: the
    escaping is. This just keeps a page free of content it has no use for, and
    it costs 16 bytes instead of 400 on every poll.
    """
    joined = "|".join(parts).encode("utf-8", "replace")
    return hashlib.sha256(joined).hexdigest()[:16]


def render_index_body(view):
    rows = []
    for it in view["projects"]:
        moved = ("%s ago" % dur(it["moved_age"]) if it["moved"]
                 else "nothing here yet")
        marks = ""
        if it["failed"] and it["bucket"] != "failed":
            # A live run that has already lost a stage sorts as running, so the
            # casualty needs saying here or it is invisible until you open it.
            marks += ('<div class="small failmark">%s failed</div>'
                      % esc(", ".join(it["failed"])))
        if it["beat_band"] in STALE_BANDS:
            # Never a verdict, in the index any more than on the project page:
            # the age and the cadence, and the reader decides.
            marks += ('<div class="small %s">no heartbeat for %s; this run '
                      'stamps one every %s s</div>'
                      % ("failmark" if it["beat_band"] == "long"
                         else "latemark", esc(dur(it["beat_age"])),
                         esc("%g" % it["heartbeat_s"])))
        elif it["beat_band"] == "unrated":
            # No cadence in the record, so no "late" here either - the age, and
            # the reason there is no judgement on it.
            marks += ('<div class="small muted">last heartbeat stamp %s ago; '
                      'this run record gives no heartbeat interval to judge '
                      'that by, so this console does not call it late</div>'
                      % esc(dur(it["beat_age"])))
        elif it["beat_band"] == "skew":
            marks += ('<div class="small latemark">the last heartbeat is '
                      'stamped in the future for this machine; the ages in '
                      'this row cannot be trusted</div>')
        if it["carried_failed"]:
            marks += ('<div class="muted small">%s failed in an earlier run, '
                      'before this one started</div>'
                      % esc(", ".join(it["carried_failed"])))
        if it.get("undated_failed"):
            # Not "failed" and not "in an earlier run": with no readable run
            # start there is nothing to compare a finish against, and the index
            # picks neither answer rather than the flattering one.
            marks += ('<div class="small failmark">%s failed; this record does '
                      'not establish which run, so this console does not say '
                      'it was this one</div>'
                      % esc(", ".join(it["undated_failed"])))
        if it["state_cached"]:
            marks += ('<div class="muted small">state file unreadable just '
                      'now; showing the last good read</div>')
        elif it["bucket"] == "unreadable":
            marks += ('<div class="small failmark">the state file will not '
                      'parse (%s) — this is not an empty directory</div>'
                      % esc(it["state_detail"] or "no detail"))
        elif it["bucket"] == "blank":
            marks += ('<div class="small failmark">the state file is there and '
                      'empty — a write that did not complete</div>')
        elif it["bucket"] == "error":
            marks += ('<div class="small failmark">this console failed while '
                      'reading this directory (%s), so nothing below is known '
                      'about it. The other rows on this page are unaffected; '
                      'the error is on the console\'s stderr</div>'
                      % esc(it["state_detail"] or "no detail"))
        rows.append(
            '<tr><td><a class="row-name" href="/p/%d">%s</a>'
            '<div class="muted small mono path" title="%s">%s</div></td>'
            '<td><span class="badge b-%s">%s</span>%s</td>'
            '<td>%s</td><td>%s<div class="muted small">%s</div></td></tr>'
            % (it["index"], esc(it["name"]), esc(it["path_short"]),
               esc(it["path_elided"]),
               esc(it["bucket"]), esc(BUCKET_LABEL[it["bucket"]]), marks,
               render_strip(it["strip"]), esc(moved), esc(it["moved_source"])))
    empty = ('<tr><td colspan="4" class="muted">No results directories are '
             'being watched.</td></tr>')
    running = any(i["bucket"] == "running" for i in view["projects"])
    alert = any(i["failed"] or i.get("undated_failed")
                or i["beat_band"] in STALE_BANDS + ("skew",)
                or i["bucket"] in ("failed", "error", "unreadable", "blank")
                for i in view["projects"])
    return (
        live_marker(running, index_sig(view), alert)
        + '<p class="summary">%s</p>'
          '<div class="panel"><table><thead><tr>'
          '<th class="c-proj">Project</th><th class="c-badge">State</th>'
          '<th class="c-strip">Stages</th><th>Last moved</th>'
          '</tr></thead><tbody>%s</tbody></table>%s</div>'
          '<footer>Sorted %s — and inside each group, whatever has stopped '
          'stamping a heartbeat first, then most recently moved. '
          '&ldquo;Last moved&rdquo; is the newer of the state file and the log, '
          'and the line under it names which of the two it came from. Snapshot '
          '%s. %s</footer>'
        % (esc(view["summary"]), "".join(rows) or empty, render_legend(),
           esc(", ".join(BUCKET_WORD[b] for b in BUCKETS)),
           esc(clock(view["taken"])), scan_sentence(view)))


def scan_sentence(view):
    """Where the watched list comes from, and whether it can still change.

    The footer said "Snapshot <fresh time>" beside a list that was fixed at
    startup, so a dataset launched an hour into the session was silently
    missing from the page that exists to replace looking. Either the list
    changes or the page says it does not; this says which.
    """
    roots = view.get("roots") or []
    if not roots:
        return ('These directories were named on the command line; the list '
                'does not change while the console runs.')
    return ('Found by scanning %s, re-scanned about every %s; last scan %s. A '
            'directory found later is added at the end and keeps its number, '
            'so a bookmarked link does not move.'
            % (esc(", ".join(collapse_home(r) for r in roots)),
               esc(dur(view.get("rescan_s") or 0)),
               esc(clock(view.get("scanned_at")))))


def render_index(view, contract, interval=15.0):
    return page("metaannot console", top_bar(contract) +
                '<div id="live" data-frag="/?f=1" data-interval="%g">%s</div>'
                % (interval, render_index_body(view)))


def render_row(r):
    bits = []
    if r["carried"]:
        bits.append('<span class="chip">from an earlier run (%s)</span>'
                    % esc(r["carried"]))
    elif r.get("era") == "unknown":
        # The row's half of the same refusal: with no readable run start, this
        # record is not evidence about THIS run either way.
        bits.append('<span class="chip">which run this record belongs to is '
                    'not established</span>')
    if r["gpu"]:
        bits.append('<span class="chip">GPU</span>')
    if r["empty_ok"] and r["state"] in ("ok", "adopted"):
        bits.append('<span class="chip">an empty output is a valid result'
                    '</span>')
    out = []
    for o in r["outputs"]:
        if o["there"]:
            out.append('<div class="mono small">%s<br><span class="muted">'
                       '%s · %s</span></div>'
                       % (esc(o["name"]), esc(human_bytes(o["size"])),
                          esc(clock(o["mtime"]))))
        else:
            # "not written yet" under a row that says "finished 16:52:45" is a
            # flat self-contradiction on one line, and it reads as "it is
            # coming". The state file and the directory disagreeing is a real
            # fact worth saying: somebody cleaned up intermediates, or the run
            # wrote somewhere else.
            #
            # It is the right sentence for a STALE row too, and for the same
            # reason: nothing is coming, because whatever was writing this
            # stopped when the earlier run did. The line claims only what a
            # stat established - declared here, absent here - and leaves why
            # to the reader, which is where this page leaves every verdict.
            out.append('<div class="muted small mono">%s<br>%s</div>'
                       % (esc(o["name"]),
                          "not written yet" if r["state"] == "running"
                          else "declared by this stage, but not in this "
                               "directory"))
    if r["part"]:
        out.append('<div class="mono small">%s<br><span class="muted">%s · '
                   'grew %s ago</span></div>'
                   % (esc(r["part"]["name"]),
                      esc(human_bytes(r["part"]["size"])),
                      esc(dur(time.time() - r["part"]["mtime"]))))
    return ('<tr class="r-%s"><td><span class="row-name">%s</span>%s</td>'
            '<td><span class="st s-%s">%s</span></td><td>%s</td>'
            '<td>%s</td><td>%s</td></tr>'
            % (esc(r["state"]), esc(r["name"]), "".join(bits), esc(r["state"]),
               esc(r["label"]), esc(r["took"]), esc(r["detail"]),
               "".join(out) or '<span class="muted">—</span>'))


# How much of the pinned block one failure may spend. The full error is on the
# FAIL row below, which is not pinned and can be as tall as it likes.
BLOCK_ERROR_MAX = 200
BLOCK_FAILS = 3


def render_failures(fails, now, over=False):
    """The first thing in the pinned block, when there is one.

    Stage, time, and the engine's own error string. Nothing else on this page
    is more urgent, and until this existed the answer was three hundred pixels
    below the fold behind fifteen green OK rows.

    Bounded, because the fix for "the failure is below the fold" must not push
    everything else below the fold: this sits at the top of a `position:
    sticky` block, and uncapped, five failures made that block 799 px tall at
    1440x860 and eight made it 1123 px, with the heartbeat line and the whole
    stage table off the screen. One wrong `db.*` root fails five stages in the
    same run, so five is not a stress case.

    `over` because a run that has ENDED is not carrying on with anything - the
    present-tense bug stage_rows had, printed over a run that finished 26 h ago
    with final_status failed.
    """
    if not fails:
        return ""
    n = len(fails)
    return ('<div class="dead"><h3>%d stage%s %sfailed in this run</h3>%s'
            '<p class="small">%s</p></div>'
            % (n, "" if n == 1 else "s",
               "" if over else ("has " if n == 1 else "have "),
               failure_lines(fails),
               "This run has ended. What is above is what it lost on the way; "
               "its own outcome is on the next line." if over else
               "Nothing new starts after the first failure: the engine stops "
               "dispatching there and waits for the stages already running, "
               "which it cannot interrupt. Those can take hours, so a run can "
               "still look busy long after this."))


def failure_lines(fails):
    """Stage, time and the engine's own error string, capped. Shared by the two
    blocks below, because the difference between them is what may be CLAIMED
    about the failures, never how a failure is shown."""
    lines = []
    for f in fails[:BLOCK_FAILS]:
        when = ("%s (%s ago)" % (clip(f["finished"], 60), dur(f["age"]))
                if f["age"] is not None
                else (clip(f["finished"], 60) or "at an unrecorded time"))
        lines.append('<p><strong>%s</strong> failed — %s<br>'
                     '<span class="err">%s</span></p>'
                     % (esc(f["name"]), esc(when),
                        esc(clip(f["error"], BLOCK_ERROR_MAX))))
    if len(fails) > BLOCK_FAILS:
        lines.append('<p class="small">…and %d more, in the table below.</p>'
                     % (len(fails) - BLOCK_FAILS))
    return "".join(lines)


def render_undated(fails, note):
    """Failures the console cannot attach to any run, and will not.

    The same block, minus the one clause that was not true: `carried` was set
    only when the run's start parsed, so a `_run` with no readable `started`
    put every historical failure under "N stages have failed in this run" -
    with the run's own record, saying no such thing, three lines below. The
    failure is still first and still red; what it is NOT is dated.
    """
    if not fails:
        return ""
    n = len(fails)
    return ('<div class="dead"><h3>%d stage%s failed, and this console cannot '
            'tell whether in this run</h3>%s<p class="small">%s The state file '
            'is cumulative across runs, so a record here may be this run\'s or '
            'an earlier one\'s; the times above are what it says, and they are '
            'all this page will claim.</p></div>'
            % (n, "" if n == 1 else "s", failure_lines(fails), esc(note or "")))


def part_sentence(r, now, beat_age, log_age):
    """What a .part file is evidence OF, which depends on how old it is.

    The claim "a more direct sign of work than any heartbeat" is true only
    while the file is actually growing. It was emitted unconditionally, so a
    part file that had not moved for 26 minutes was offered as stronger
    evidence of work than a heartbeat that had been silent for 2 - and on a
    real InterProScan run a part file stalled for four hours printed the same
    reassuring clause. When it is the stalest clock on the page, that is
    itself the interesting fact.
    """
    age = now - r["part"]["mtime"]
    others = [a for a in (beat_age, log_age) if a is not None]
    head = ('<strong>%s</strong> is writing <span class="mono">%s</span> — %s, '
            'last grew %s ago.'
            % (esc(r["name"]), esc(r["part"]["name"]),
               esc(human_bytes(r["part"]["size"])), esc(dur(age))))
    if others and age > min(others) + 60:
        return ('<p class="small failmark">%s That file is STALER than every '
                'other clock on this page (%s), so it is not evidence of work '
                '— it is the strongest sign here that this stage has stopped '
                'producing bytes.</p>'
                % (head, esc(dur(min(others)))))
    return ('<p class="small">%s That is a more direct sign of work than any '
            'heartbeat.</p>' % head)


def render_vitals(v):
    """The block that answers the only urgent question, pinned to the top.

    Scrolling down twenty-one stage rows must not cost you the answer, so it is
    position: sticky. Order is by urgency, not by chronology: a failed stage,
    then anything wrong with the state file, then whether work is happening.

    TWO panels, and the split is the whole point. What is pinned is only what
    is urgent and short: the failures, the state-file trouble, the verdict. The
    rest - the run identifiers, the explanatory lines, the log, part and lock
    sentences - goes in a second panel directly below, in ordinary flow.
    A sticky block that grows without limit takes the viewport with it, and the
    cap that used to answer that made its own tail unreachable (see app.css).
    Keeping the pinned half small by CONSTRUCTION is the only version of this
    that needs neither a cap nor a scrollbar.
    """
    hb, run = v["heartbeat"], v["run"]
    now = time.time()
    head = ""
    if run:
        head = ('<p class="small muted mono">run %s · metaannot %s · pid %s on '
                '%s · started %s</p>'
                % (esc(run.get("run_id")), esc(run.get("version")),
                   esc(run.get("pid")), esc(run.get("host")),
                   esc(run.get("started"))))
        # The console asks whatever --metaannot points at for the stage graph,
        # and that need not be the engine that made this directory. Where the
        # two disagree, the table below may be describing a different pipeline
        # from the one that ran, and only the reader can judge how much.
        if run.get("version") and run["version"] != v.get("engine_version"):
            head += ('<p class="small">This directory was written by metaannot '
                     '%s, and the stage list below comes from metaannot %s. '
                     'Where the two differ, this table is the newer one\'s '
                     'idea of the pipeline.</p>'
                     % (esc(run["version"]), esc(v.get("engine_version"))))
        # The design's central premise is that the console runs on the host
        # that writes the directory, and nothing was said when it did not.
        # `describe --json` reports this machine's host, so this is free.
        if (run.get("host") and v.get("engine_host")
                and run["host"] != v["engine_host"]):
            head += ('<p class="small">This directory records host %s and this '
                     'console is running on %s. Every age here is measured '
                     'against THIS machine\'s clock, the ps line below is for '
                     'that other one, and a file mtime crossing a shared mount '
                     'is only as good as the two clocks.</p>'
                     % (esc(clip(run["host"], 120)), esc(v["engine_host"])))
    # Ahead of the run id, the version and the pid: those are identifiers, and
    # this is the answer. Nothing in `head` is more urgent than a dead stage or
    # a state file that will not parse.
    lead = render_failures(hb["failures"], now, over=run_is_over(run))
    lead += render_undated(v.get("undated_failures") or (), v.get("era_note"))
    if v.get("state_trouble"):
        lead += ('<div class="trouble"><p><strong>%s</strong></p></div>'
                 % esc(v["state_trouble"]))
    # The verdict goes ABOVE the run id, not below it. For a run that has
    # finished failed, "Finished failed at 2026-09-08T21:54:42 · ran 4h00m" is
    # the headline and `run 20260908T175442-1 · metaannot 0.3.0 · pid 4242` is
    # a set of identifiers; they were the other way round.
    lead += '<p><strong>%s</strong></p>' % esc(hb["verdict"])
    body = ['<p>%s</p>' % esc(line) for line in hb["lines"]]
    log_age = None
    if v["log_mtime"]:
        log_age = now - v["log_mtime"]
        body.append('<p class="small">The log last grew %s ago (%s, now %s).'
                    '</p>' % (esc(dur(log_age)), esc(clock(v["log_mtime"])),
                              esc(human_bytes(v["log_size"]))))
    elif v.get("log_why") and v["log_why"] != "missing":
        # stat_or_none folded EACCES into None, so a results directory this
        # account cannot read was reported as one with no log in it.
        body.append('<p class="small failmark">The log file cannot be read '
                    '(%s). That is not the same as there being no log.</p>'
                    % esc(v["log_why"]))
    else:
        body.append('<p class="small">There is no log file in this directory.'
                    '</p>')
    for r in (r for r in v["rows"] if r["part"]):
        body.append(part_sentence(r, now, hb.get("beat_age"), log_age))
    lock = v["lock"]
    settle = ""
    if lock.get("ps"):
        settle = (' To settle it yourself, on that host: <code>%s</code>. This '
                  'console offers no --force-unlock — two runs writing one '
                  'results directory corrupt it silently, and that decision is '
                  'yours.' % esc(lock["ps"]))
    elif lock.get("why_no_ps"):
        settle = " " + esc(lock["why_no_ps"])
    body.append('<p class="small">%s%s</p>' % (esc(lock["text"]), settle))
    band = "fail" if hb["failures"] or v.get("undated_failures") else hb["band"]
    return ('<div class="panel vitals sticky band-%s"><h2>%s</h2>%s</div>'
            '<div class="panel vitals-more">%s%s</div>'
            % (esc(band), esc(v["name"]), lead, head, "".join(body)))


def vitals_alert(v):
    """Whether anything on this page is wrong enough to pin the poll at the
    floor. See live_marker()."""
    return bool(v["heartbeat"]["failures"] or v.get("undated_failures")
                or v.get("state_trouble")
                or v["heartbeat"]["band"] in ("late", "long", "skew"))


def project_sig(v):
    """A digest of the facts on this page, with no clock in it.

    The poller compares this across fetches to decide whether anything has
    actually changed. Relative ages must stay out of it: they move every second
    on a page where nothing is happening, which is exactly the case the
    back-off exists for.
    """
    parts = [str(v["state_status"]), str(v["state_mtime"]), str(v["log_mtime"]),
             str(v["log_size"]), str(v["bucket"]), str(v["heartbeat"]["band"]),
             str(v["heartbeat"]["running"]), str(v["lock"]["text"])]
    for r in v["rows"]:
        parts.append("%s=%s" % (r["name"], r["state"]))
        if r["part"]:
            parts.append("%s@%s" % (r["part"]["size"], r["part"]["mtime"]))
    return digest(parts)


def queue_note(rows):
    """What the ranking on the NEXT rows means, said once under the table.

    The rows themselves claim only that the engine ranks one ahead of another,
    which is all this console can establish. The two things it cannot see go
    here, once, rather than as a caveat on ten rows: how many stages start
    together is stage_workers against what is still running, and gpu_lease can
    hold a GPU stage back however it is ranked - a deferred stage stays in
    `remaining` and is reconsidered next round.

    Empty when nothing is ranked, so an engine older than `cost` gets neither
    the ranking nor an explanation of a ranking that is not there.
    """
    if not ready_rows(rows):
        return ""
    return ('<div class="muted small note">The NEXT rows carry the order the '
            'engine considers them in: it sorts each round\'s ready stages '
            'longest-first — hours, then minutes, then seconds — and '
            'dispatches from the top of that sort. That is an order and not a '
            'schedule. How many start together depends on stage_workers and '
            'on what is still running, and a stage that needs the GPU waits '
            'for gpu_workers however it is ranked.</div>')


def render_project_body(v):
    """The live region: vitals, the state-file note, the stage table.

    The log pane is deliberately NOT in here. This fragment is re-fetched and
    swapped whole every few seconds; the log pane is appended to by byte offset
    on its own clock. Two writers over one pane would show every new line twice
    and lose the operator's scroll each time, so they are kept apart.
    """
    counts = v["counts"]
    order = ["running", "stale", "failed", "ok", "adopted", "bad", "wait",
             "next", "off", "none"]
    tally = ", ".join("%d %s" % (counts[k], k) for k in order if k in counts)
    if v["cfg_run_known"]:
        cfg = ('<div class="muted small note">enabled/off %s.</div>'
               % esc(v["cfg_note"]))
    else:
        cfg = ('<div class="muted small note">enabled/off: cannot tell — %s. '
               'No stage below is reported OFF on a guess; they fall back to '
               'waiting or not-yet-reached.</div>' % esc(v["cfg_note"]))
    cfg += queue_note(v["rows"])
    # Anything actually WRONG with the state file is now said at the top of the
    # vitals block, by state_trouble(); what is left here is the quiet cases.
    if v.get("state_trouble"):
        note = ""
    elif not v["run"] and v["state_status"] != "missing":
        note = ('<p class="muted small">The state file carries no run record — '
                'written by a metaannot older than the one that added it, or '
                'by a run that never got that far.</p>')
    else:
        note = ""
    if v.get("unknown_stages"):
        # A directory written by a newer engine. These records were simply
        # dropped, with no mention, whenever _run.version happened to match.
        note += ('<p class="small failmark">The state file also has records '
                 'for %s, which this contract has no stage for. They are not '
                 'in the table below — the engine that wrote this directory '
                 'knows stages the one at --metaannot does not.</p>'
                 % esc(", ".join(v["unknown_stages"])))
    return (
        live_marker(v["heartbeat"]["running"], project_sig(v), vitals_alert(v))
        + render_vitals(v) + note
        + '<div class="panel"><table><thead><tr>'
          '<th class="c-stage">Stage</th><th class="c-state">State</th>'
          '<th class="c-took">Took / for</th>'
          '<th>What it is doing, or waiting on</th>'
          '<th class="c-out">Newest output</th></tr></thead><tbody>%s'
          '</tbody></table>%s'
          '<div class="muted small note">%d stages as '
          'metaannot %s defines them (describe --json, schema %s): %s. '
          'Snapshot %s; %s.</div></div>'
          % ("".join(render_row(r) for r in v["rows"]), cfg, len(v["rows"]),
             esc(v.get("engine_version")), esc(v.get("describe_version")),
             esc(tally), esc(clock(v["taken"])),
             esc("the state file itself last changed %s"
                 % clock(v["state_mtime"]) if v["state_mtime"]
                 else "there is no state file here")))


def render_log_pane(v):
    """Owned by the log poller alone, and never re-rendered by the fragment.

    `data-max` is the pane's cap, decided here rather than in the script for
    the same reason LOG_LINES is decided here: how much of a three-day log this
    console will hold is a property of the console, and one place to change it
    is one place to read it. app.js trims to it and says so when it has.
    """
    log = v.get("log") or {}
    lines = "".join('<div class="ln l-%s">%s</div>'
                    % (esc(l["level"]), esc(l["text"]))
                    for l in log.get("lines", ()))
    if log.get("note"):
        lines = ('<div class="marker">— %s —</div>' % esc(log["note"])) + lines
    return (
        '<div class="logpane"><div class="loghead"><strong>%s</strong>'
        '<span class="muted">last %d lines</span>'
        '<span class="muted">%d WARN in view</span>'
        '<span class="muted">%d FATAL in view</span></div>'
        '<div id="logbody" data-off="%s" data-ident="%s" data-project="%d" '
        'data-max="%d">'
        '%s</div><button class="newbtn" id="lognew" hidden>new lines below</button>'
        '</div>'
        % (esc(v["log_name"]), len(log.get("lines", ())), log.get("warn", 0),
           log.get("fatal", 0), esc(log.get("off", 0)),
           esc(log.get("ident") or ""), v["index"], LOG_PANE_LINES,
           lines or '<div class="marker">— no log lines —</div>')
        + '<footer>The counts above are for the lines shown, not for the whole '
          '%s log — a three-day log is never scanned end to end. Failure is not '
          'detected here: a stage that fails is recorded in the state file and '
          'appears in the table above whether or not its FATAL line is still '
          'inside this window.</footer>' % esc(human_bytes(v["log_size"])))


def decorate(v, contract):
    """The three contract facts the renderers quote, added to a view."""
    return dict(v, engine_version=contract.version,
                describe_version=contract.describe_version,
                engine_host=contract.host, log_name=contract.log_name)


def render_project(v, contract, interval=3.0):
    v = decorate(v, contract)
    body = (top_bar(contract, '<span class="small"><a href="/">&larr; all '
                              'projects</a></span>')
            + '<p class="muted small mono">%s</p>' % esc(v["path_short"])
            + '<div id="live" data-frag="/p/%d?f=1" data-interval="%g">%s</div>'
              % (v["index"], interval, render_project_body(v))
            + render_log_pane(v))
    return page("%s — metaannot console" % v["name"], body)


# ======================================================================
# 5. http - GET and HEAD, and nothing else
# ======================================================================

class Console:
    """The watched set, and the one place a route turns an index into a path.

    No endpoint accepts a filesystem path. Projects are numbered 0..N-1 from
    what the operator named on the command line, and the API takes `?p=3`; a
    path parameter would turn a console anyone on the box could reach into a
    whole-filesystem read oracle.
    """

    def __init__(self, contract, projects, python, script, interval=3.0,
                 roots=(), rescan_s=RESCAN_S):
        self.contract = contract
        self.projects = list(projects)
        self.python = python
        self.script = script
        self.interval = interval
        self.roots = [os.path.abspath(r) for r in roots]
        self.rescan_s = rescan_s
        self.scanned_at = time.time()
        self._scan_lock = threading.Lock()
        self._seen = {p.path for p in self.projects}
        self._last_scan = self.scanned_at
        self._scanning = False

    def project(self, raw):
        try:
            idx = int(raw)
        except (TypeError, ValueError):
            return None
        projects = self.projects
        if 0 <= idx < len(projects):
            return projects[idx]
        return None

    def rescan(self, force=False):
        """Pick up a dataset launched after this console started.

        --root was scanned once, inside build(), before the server existed, and
        nothing on the page said so while the footer's "Snapshot <time>"
        implied the list was current: on the run book's own deployment you kick
        off dataset five, refresh, and it is simply not there.

        Appending only, and never renumbering: a project keeps the index it was
        given, so a bookmarked /p/3 is still that directory after the ninth
        dataset appears. The page orders by bucket rather than by index, so a
        late arrival is not stuck at the bottom. Rate-limited because this is a
        scandir on what may be an NFS mount, and it runs on a request thread.
        """
        if not self.roots:
            return
        now = time.time()
        with self._scan_lock:
            # _last_scan was stamped HERE, before the scan, and nothing
            # stopped a second thread starting another: against a --root on an
            # unresponsive mount, where os.path.exists blocks in the kernel,
            # that parked one request thread every 30 s - threads 3 to 8 in
            # 160 s and MAX_CONNS in sixteen minutes. One at a time, stamped
            # when it finishes, so a slow scan costs one thread and not one an
            # interval.
            if self._scanning:
                return
            if not force and now - self._last_scan < self.rescan_s:
                return
            self._scanning = True
        try:
            found = discover(self.roots, self.contract)
        finally:
            with self._scan_lock:
                self._scanning = False
                self._last_scan = time.time()
        with self._scan_lock:
            fresh = [p for p in found if p not in self._seen]
            for path in fresh:
                self._seen.add(path)
            n = len(self.projects)
            # Rebound, not appended in place: `self.projects` is read without
            # the lock by every other request thread, and a list that is only
            # ever replaced is one such a reader cannot catch mid-write.
            self.projects = self.projects + [
                Project(n + i, p, self.contract) for i, p in enumerate(fresh)]
            self.scanned_at = time.time()

    def scan_facts(self):
        return {"roots": list(self.roots), "scanned_at": self.scanned_at,
                "rescan_s": self.rescan_s}

    def index(self):
        self.rescan()
        return index_view(self.projects, self.python, self.script,
                          self.scan_facts())

    def view(self, project, **kw):
        return project_view(project, self.python, self.script, **kw)


class Handler(BaseHTTPRequestHandler):
    """GET and HEAD. There is no do_POST, so the stdlib answers 501 to
    everything that could change anything, and CSRF has nothing to attack.

    Two overrides that are not cosmetic. On an AF_UNIX socket `client_address`
    is the empty string, so the inherited `address_string()` evaluates ''[0] and
    raises IndexError from inside send_response - after the socket is closed,
    which delivers the browser zero bytes and a blank tab. That is precisely the
    failure the design doc warns about, so both are replaced.
    """

    protocol_version = "HTTP/1.1"      # so every reply carries Content-Length
    timeout = HANDLER_TIMEOUT          # reaps parked keep-alive threads
    server_version = "metaannot-console"
    sys_version = ""

    def address_string(self):
        return "unix"

    def log_message(self, fmt, *args):
        pass                           # one line per poll, forever, otherwise

    # -- plumbing ------------------------------------------------------
    def send(self, code, ctype, body, with_body=True):
        if isinstance(body, str):
            # errors="replace", because every string on the page came out of a
            # file this console does not own. A lone surrogate anywhere in
            # _run.host, a stage name or a log line raises UnicodeEncodeError
            # here - AFTER the whole page has rendered - and the project is
            # then a 500 for a byte it could have shown as U+FFFD. The engine
            # made the same call in opener() for the same reason: a bad byte
            # should cost one character, not the page.
            body = body.encode("utf-8", "replace")
        # Recorded before the status line goes out, so the handler below can
        # tell "the render blew up" from "the render was already on the wire".
        # A second response written into the same stream is a corrupt page,
        # which is worse than the error it was trying to report.
        self.sent = True
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        # No CORS headers at all: once the forward is up, a page open in the
        # operator's own browser can issue requests to localhost:8080, and the
        # absence of CORS is what stops it reading the answers.
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy",
                         "default-src 'self'; base-uri 'none'; form-action 'none'")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        if with_body:
            self.wfile.write(body)

    def send_json(self, obj, with_body=True, code=200):
        self.send(code, "application/json; charset=utf-8",
                  json.dumps(obj, default=str), with_body)

    def do_GET(self):
        self.route(True)

    def do_HEAD(self):
        self.route(False)

    def route(self, with_body):
        self.sent = False
        try:
            self.dispatch(with_body)
        except (BrokenPipeError, ConnectionResetError):
            raise
        except Exception as e:                       # never a blank tab
            sys.stderr.write("console: %s while serving %s\n"
                             % (repr(e), self.path))
            if self.sent:
                # Half a page already went out. Close the connection rather
                # than append a second status line to it.
                self.close_connection = True
                return
            self.send(500, "text/html; charset=utf-8",
                      page("console error",
                           '<div class="panel"><h2>The console failed to '
                           'render this page.</h2><p class="mono">%s</p>'
                           '<p>Nothing was written to any results directory. '
                           'The error is on the console\'s stderr.</p>'
                           '<p><a href="/">All projects</a></p></div>'
                           % esc(repr(e))), with_body)

    def dispatch(self, with_body):
        con = self.server.console
        url = urlparse(self.path)
        q = parse_qs(url.query)
        path = url.path

        if path == "/app.css":
            return self.send(200, "text/css; charset=utf-8", PAGE_CSS, with_body)
        if path == "/app.js":
            return self.send(200, "application/javascript; charset=utf-8",
                             PAGE_JS, with_body)
        if path == "/":
            view = con.index()
            if q.get("f"):
                return self.send(200, "text/html; charset=utf-8",
                                 render_index_body(view), with_body)
            return self.send(200, "text/html; charset=utf-8",
                             render_index(view, con.contract,
                                          max(con.interval, 15.0)), with_body)
        if path.startswith("/p/"):
            project = con.project(path[3:])
            if project is None:
                return self.not_found(with_body)
            # The fragment carries the live region alone, so it reads the last
            # PEEK_WINDOW bytes of the log rather than the 256 KB tail - the
            # vitals block's FATAL sentence needs the newest line's level, and
            # a fragment that quietly dropped that sentence three seconds after
            # the page loaded would be worse than never showing it.
            wants_fragment = bool(q.get("f"))
            view = con.view(project, log=not wants_fragment)
            if wants_fragment:
                return self.send(200, "text/html; charset=utf-8",
                                 render_project_body(decorate(view,
                                                              con.contract)),
                                 with_body)
            return self.send(200, "text/html; charset=utf-8",
                             render_project(view, con.contract, con.interval),
                             with_body)
        if path == "/api/contract":
            return self.send_json(con.contract.facts(), with_body)
        if path == "/api/projects":
            return self.send_json(con.index(), with_body)
        if path == "/api/project":
            project = con.project((q.get("p") or [None])[0])
            if project is None:
                return self.not_found(with_body)
            return self.send_json(con.view(project, log=False), with_body)
        if path == "/api/log":
            project = con.project((q.get("p") or [None])[0])
            if project is None:
                return self.not_found(with_body)
            off = (q.get("off") or [None])[0]
            ident = (q.get("id") or [None])[0] or None
            try:
                off = int(off) if off is not None else None
            except ValueError:
                off = None
            tail = tail_log(project.log_path, off, ident)
            # The offset advances past every byte read, so any line log_lines
            # leaves out is gone from the live pane for good. It gets said.
            lines, dropped = log_lines(tail["text"])
            return self.send_json(
                {"ok": tail["ok"], "ident": tail["ident"], "off": tail["off"],
                 "size": tail["size"], "reset": tail["reset"],
                 "dropped": dropped, "lines": lines,
                 "note": join_notes(tail["note"],
                                    dropped_note(dropped, tail["reset"]))},
                with_body)
        return self.not_found(with_body)

    def not_found(self, with_body):
        self.send(404, "text/html; charset=utf-8",
                  page("not found",
                       '<div class="panel"><h2>No such page.</h2>'
                       '<p><a href="/">All projects</a></p></div>'), with_body)


# ======================================================================
# 6. server - a 0700 UNIX socket, and nothing listening on a TCP port
# ======================================================================

class ConsoleServer(socketserver.ThreadingMixIn, socketserver.UnixStreamServer):
    """UnixStreamServer, deliberately not HTTPServer.

    `HTTPServer.server_bind` does `host, port = self.server_address[:2]`, which
    on an AF_UNIX address slices two CHARACTERS off the path (and raises
    ValueError outright on a one-character one), and its
    `allow_reuse_address = 1` leaves the socket at mode 0755 - world-connectable
    on the shared lab box this exists to stay off.
    """

    daemon_threads = True
    block_on_close = False             # 3.12 semantics, stated for 3.9
    allow_reuse_address = False        # SO_REUSEADDR means nothing on AF_UNIX
    request_queue_size = 16
    # BaseHTTPRequestHandler reads these; on AF_UNIX there is no host or port to
    # take them from.
    server_name, server_port = "metaannot-console", 0

    def __init__(self, path, handler, console):
        self.console = console
        self._slots = threading.BoundedSemaphore(MAX_CONNS)
        socketserver.UnixStreamServer.__init__(self, path, handler)

    def server_bind(self):
        # bind() creates the inode with 0777 & ~umask, so the socket is 0700 at
        # creation. bind()-then-chmod() leaves a window in which any account on
        # a shared box can connect, and this closes it. umask is process-global
        # and not thread-safe, which is fine here and only here: this runs once,
        # at startup, before a single serving thread exists.
        old = os.umask(0o077)
        try:
            socketserver.UnixStreamServer.server_bind(self)
        finally:
            os.umask(old)

    def process_request(self, request, client_address):
        # A hard cap, refused without parsing a byte. This machine is running
        # twenty-one stages under a CPU and RAM budget; the console must be
        # incapable of thread-bombing it, whatever a browser does.
        if not self._slots.acquire(blocking=False):
            try:
                request.sendall(BUSY_RESPONSE)
            except OSError:
                pass
            # The BASE close, not self.shutdown_request: this connection never
            # took a permit, and returning one it never held would hand out a
            # slot per refusal and raise the cap the refusal exists to enforce.
            socketserver.UnixStreamServer.shutdown_request(self, request)
            return
        socketserver.ThreadingMixIn.process_request(self, request,
                                                    client_address)

    def shutdown_request(self, request):
        # Called exactly once per ACCEPTED connection, by
        # ThreadingMixIn.process_request_thread, so it is where the permit goes
        # back.
        socketserver.UnixStreamServer.shutdown_request(self, request)
        try:
            self._slots.release()
        except ValueError:
            pass

    def handle_error(self, request, client_address):
        exc = sys.exc_info()[1]
        if isinstance(exc, (BrokenPipeError, ConnectionResetError,
                            socket.timeout)):
            return                     # a tab was closed; not news
        sys.stderr.write("console: %r\n" % (exc,))


def path_chain(path):
    """`path` and every directory above it, innermost first."""
    out, cur = [], os.path.abspath(path)
    while True:
        out.append(cur)
        parent = os.path.dirname(cur)
        if parent == cur:
            return out
        cur = parent


def why_not_private(path, uid):
    """Why `path` is not a private directory this uid owns, or None.

    lstat, not stat. The last runtime candidate is /tmp/metaannot-console-<uid>
    and the bind-failure message names exactly that; a symlink planted there by
    another account resolved through os.stat, whose uid and mode checks then
    described the TARGET while the socket was bound at the LINK.
    """
    st = lstat_or_none(path)
    if st is None:
        return "does not exist"
    if stat.S_ISLNK(st.st_mode):
        return "is a symlink"
    if not stat.S_ISDIR(st.st_mode):
        return "is not a directory"
    if st.st_uid != uid:
        return "is owned by uid %d, not by this account" % st.st_uid
    if st.st_mode & 0o077:
        return ("is reachable by other accounts (mode %04o)"
                % (st.st_mode & 0o7777))
    return None


def exposed_ancestor(path, uid, mask=0o002):
    """(directory, why) for the first directory ABOVE `path` that another
    account could interfere with, or (None, None).

    The half that was never checked. os.makedirs applies its `mode` to the LEAF
    only, so every intermediate it invents comes out 0777 & ~umask: measured,
    `--socket /tmp/x/a/b/m.sock` under umask 0 left /tmp/x/a at 0777 with only
    b at 0700, and /tmp/x/a is not sticky, so any account can rename b away and
    put its own directory there - after which the operator's `ssh -L` reaches a
    socket somebody else is serving. The same hole is open on an ancestor that
    merely already exists that way, so this walks the whole chain rather than
    what the console created. Group- or world-writable is tolerated WITH the
    sticky bit (/tmp has it), because sticky is exactly what stops another
    account renaming an entry it does not own.

    `mask` is 0o002 - ANY account - because that is the measured attack and
    the only one worth a refusal. GROUP-writable is asked for separately, by
    the callers, and only warned about: pass 0o020. A home directory at 0775
    is the default on any Linux with umask 002 and per-user groups, where the
    group has exactly one member, and refusing there would send a diligent
    operator round the loop the old "Try: mkdir -m 700 -p /tmp" message sent
    them round. Whether that group is one person or forty is the site's fact,
    not this console's, so it is reported rather than judged.

    The chain walked is the REAL one, so no component is a symlink and no check
    can be pointed at the wrong inode. That leaves a symlink component whose
    own directory is outside the real chain - reachable only from a directory
    an attacker already controls.
    """
    for anc in path_chain(os.path.realpath(path))[1:]:
        st = lstat_or_none(anc)
        if st is None:
            return anc, "does not exist"
        if st.st_uid not in (0, uid):
            return anc, ("is owned by uid %d, which is neither root nor this "
                         "account" % st.st_uid)
        if st.st_mode & mask and not st.st_mode & stat.S_ISVTX:
            return anc, ("is writable by %s (mode %04o) and is not sticky, so "
                         "they can replace what is inside it"
                         % ("any account on this machine" if mask & 0o002
                            else "gid %d" % st.st_gid, st.st_mode & 0o7777))
    return None, None


def ensure_private_dir(path, uid):
    """Create `path` at 0700 if it is not there. Returns why_not_private().

    os.mkdir, not os.makedirs: makedirs invents intermediates at 0777 & ~umask
    (see exposed_ancestor), and a console that creates a chain of directories
    is a console that can create a RESULTS directory - which is what an
    explicit --socket made it do. A candidate whose parent is not already there
    is simply skipped; the next one down the list will do.
    """
    try:
        os.mkdir(path, 0o700)
    except FileExistsError:
        pass
    except OSError as e:
        return "cannot be created (%s)" % e
    return why_not_private(path, uid)


def runtime_dir(uid=None, candidates=None, why=None, veto=None):
    """A directory this uid owns at 0700. Its search bit is the access control.

    The socket's own mode is belt. This is braces, and on the platforms that
    ignore the mode bits of an AF_UNIX inode entirely it is the only thing
    holding. $XDG_RUNTIME_DIR already is such a directory by construction, set
    by pam_systemd on any interactive SSH login, which is why it is first.

    This creates only the last component, and only for a candidate of its OWN
    choosing. An explicit --socket goes through socket_path(), which creates
    nothing at all.

    `veto` is consulted on each candidate BEFORE that candidate is brought into
    being, and raises to reject it. Order is the whole point: ensure_private_dir
    creates, so a veto applied to the RESULT of this function has already lost -
    see socket_path().
    """
    uid = os.getuid() if uid is None else uid
    if candidates is None:
        candidates = (os.environ.get("XDG_RUNTIME_DIR"),
                      "/run/user/%d" % uid,
                      os.path.expanduser("~/.cache/metaannot-console"),
                      "/tmp/metaannot-console-%d" % uid)
    for cand in candidates:
        if not cand:
            continue
        if veto is not None:
            # Nothing has been created yet, and if this raises nothing will be.
            veto(cand)
        bad = ensure_private_dir(cand, uid)
        if bad is None:
            anc, anc_why = exposed_ancestor(cand, uid)
            if anc is None:
                return cand
            bad = "sits under %s, which %s" % (anc, anc_why)
        if lstat_or_none(cand) is not None:
            # A candidate that is simply not there, on a machine that has no
            # /run/user, is not news; one that EXISTS and is reachable by
            # somebody else is.
            sys.stderr.write("console: %s is not a private directory: it %s; "
                             "trying the next one\n" % (cand, bad))
    raise Refuse(why or "no private runtime directory found; pass --socket PATH")


def take_socket_name(path):
    """Own `path` exclusively, then leave it free for bind(). Returns the lock fd.

    flock first, so the stale-socket probe below cannot race a second console
    into unlinking a live one. The lock is on the inode and the kernel drops it
    when this process dies, so there is no stale lock and no --force-unlock to
    write. It lives beside the socket, in the console's own runtime directory,
    and never in a results directory.

    A stale socket file is not self-healing: server_close() closes the fd and
    leaves the inode, and so does a SIGKILL. Recovery is a connect probe -
    ECONNREFUSED means nothing is listening, so unlink; a successful connect
    means a live console, so refuse. And a path that is not a socket is never
    unlinked, because `--socket ~/thesis.docx` must not delete the thesis.

    O_NOFOLLOW on the lock, and S_ISREG on the fd. Without them this call
    followed a symlink planted at `<sock>.lock`: with an existing target it
    opened and flocked someone else's file, and with a missing one it CREATED a
    file at the attacker's chosen path, owned by this uid. Unreachable at the
    default socket path, whose runtime directory is validated 0700 - and
    perfectly reachable the moment --socket points into a world-writable
    directory, which is what the bind-failure message below suggests.
    """
    if len(os.fsencode(path)) > 100:
        # sun_path is 108 bytes on Linux and 104 on macOS, and the failure is an
        # unhelpful OSError at bind time.
        raise Refuse("socket path is too long for AF_UNIX (%d bytes): %s\n"
                     "Pass a shorter --socket PATH." % (len(path), path))
    try:
        fd = os.open(path + ".lock",
                     os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK,
                     0o600)
    except OSError as e:
        raise Refuse("cannot create the console's own lock beside %s: %s\n"
                     "Pass --socket PATH somewhere this account can write."
                     % (path, e))
    try:
        if not stat.S_ISREG(os.fstat(fd).st_mode):
            os.close(fd)
            raise Refuse("%s.lock exists and is not a regular file; refusing "
                         "to use it as this console's lock" % path)
    except OSError as e:
        os.close(fd)
        raise Refuse("cannot inspect %s.lock: %s" % (path, e))
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        os.close(fd)
        raise Refuse("another console already holds %s" % path)
    st = None
    try:
        st = os.lstat(path)
    except FileNotFoundError:
        return fd                                    # clean slate
    except OSError as e:
        os.close(fd)
        raise Refuse("cannot inspect %s: %s" % (path, e))
    if not stat.S_ISSOCK(st.st_mode):
        os.close(fd)
        raise Refuse("%s exists and is not a socket; refusing to remove it"
                     % path)
    if st.st_uid != os.getuid():
        # S_ISSOCK alone was the whole test before the unlink below. Reachable
        # only through a --socket in a directory another account can write,
        # which is now refused - but "the way in is refused" is how the last
        # three holes were argued.
        os.close(fd)
        raise Refuse("%s is a socket owned by uid %d, not by this account; "
                     "refusing to remove it" % (path, st.st_uid))
    probe = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    probe.settimeout(PROBE_TIMEOUT)
    try:
        probe.connect(path)
        os.close(fd)
        raise Refuse("a console is already serving %s" % path)
    except (ConnectionRefusedError, FileNotFoundError):
        try:
            os.unlink(path)                          # stale: a SIGKILL, or a reboot
        except OSError as e:
            os.close(fd)
            raise Refuse("cannot remove the stale socket %s: %s" % (path, e))
    finally:
        probe.close()
    return fd


def banner(console, sock_path, out=sys.stdout):
    uid = os.getuid()
    lines = [
        "metaannot console %s — read-only. Nothing is written to any results "
        "directory." % CONSOLE_VERSION,
        "  socket    %s   (mode 0700, uid %d)" % (sock_path, uid),
        "  engine    %s %s" % (console.python, console.script),
        "            metaannot %s, describe schema %s, signatures v%s"
        % (console.contract.version, console.contract.describe_version,
           console.contract.signature_version),
    ]
    if console.projects:
        lines.append("  watching  %d director%s:"
                     % (len(console.projects),
                        "y" if len(console.projects) == 1 else "ies"))
        for p in console.projects:
            lines.append("            [%d] %s" % (p.index, p.path))
    else:
        lines.append("  watching  nothing yet — pass --root DIR or a results "
                     "directory")
    if console.roots:
        lines.append("  scanning  %s every %s; a directory found later is "
                     "added and keeps its number"
                     % (", ".join(console.roots), dur(console.rescan_s)))
    else:
        lines.append("  scanning  nothing — this list is fixed for the life "
                     "of the process (pass --root DIR to have it re-scanned)")
    lines += [
        "",
        "From your workstation:",
        "    ssh -N -L 8080:%s %s" % (sock_path, socket.gethostname()),
        "then open http://localhost:8080/    (any free local port; -L 8081:… "
        "if 8080 is taken)",
        "",
        "File permissions are the whole access control, by decision: there is "
        "no token,",
        "because a token leaks into ps, shell history and the URL bar. A "
        "UNIX-socket forward",
        "target needs OpenSSH 6.7 or newer on your side.",
        "Ctrl-C to stop. To outlive your shell: tmux, setsid, or nohup — a "
        "SIGHUP already",
        "ignored on the way in is left ignored, so nohup keeps working.",
    ]
    out.write("\n".join(lines) + "\n")
    out.flush()


def serve(console, sock_path, ready=None, stop=None, announce=None):
    """Bind, serve until a signal, unlink. Returns the exit code.

    serve_forever() runs in a daemon thread and the MAIN thread parks: calling
    shutdown() from inside the serving thread deadlocks, and the signal handler
    therefore only sets an event. SIGHUP is handled beside SIGINT and SIGTERM
    because the run book uses tmux over ssh, where a dropped connection is a
    SIGHUP and a console that ignored it would leave the socket behind.

    `stop` and `ready` exist so the tests can drive the real lifecycle - the
    flock, the 0700 bind, the shutdown and the unlink - rather than a
    look-alike assembled beside it.
    """
    lock_fd = take_socket_name(sock_path)
    try:
        srv = ConsoleServer(sock_path, Handler, console)
    except OSError as e:
        os.close(lock_fd)
        # The one that actually happens: an NFS or CIFS home that will not hold
        # an AF_UNIX socket at all. Nothing about that is the operator's fault
        # and nothing about it is fixable except by binding elsewhere.
        raise Refuse("could not bind %s: %s\nIf that filesystem cannot hold a "
                     "socket (NFS homes often cannot), pass --socket PATH "
                     "somewhere local, such as /tmp/metaannot-console-%d/"
                     "metaannot.sock." % (sock_path, e, os.getuid()))
    stop = threading.Event() if stop is None else stop
    for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        try:
            if (sig == signal.SIGHUP
                    and signal.getsignal(sig) == signal.SIG_IGN):
                # nohup and setsid install SIG_IGN before exec. Replacing it
                # meant `nohup python3 console/console.py … &` exited cleanly
                # the moment the operator logged out - the one command they
                # would reach for on a box without tmux, silently defeated.
                # The socket then outlives the process, and take_socket_name's
                # stale-socket probe is what recovers from that.
                continue
            signal.signal(sig, lambda n, f: stop.set())
        except (ValueError, OSError):
            pass                       # not the main thread: a test drove us
    thread = threading.Thread(target=srv.serve_forever,
                              kwargs={"poll_interval": 0.2}, daemon=True)
    thread.start()
    # After the flock, the stale-socket probe and the bind, so that a refusal
    # like "thesis.docx exists and is not a socket" cannot appear underneath a
    # complete success banner ending "Ctrl-C to stop."
    if announce is not None:
        announce()
    if ready is not None:
        ready(srv)
    try:
        while not stop.wait(1.0):
            pass
    finally:
        srv.shutdown()
        srv.server_close()             # closes the fd; it does not unlink
        try:
            os.unlink(sock_path)
        except OSError:
            pass
        os.close(lock_fd)
    return 0


# ======================================================================
# 7. cli
# ======================================================================

def looks_like_results(path, contract):
    """True when this directory holds a state file or a log — or when this
    account cannot tell.

    os.path.exists() is False on EACCES, so a group-unreadable results
    directory vanished from a --root scan with nothing said: eleven fixtures,
    one at mode 000, and a banner reading "watching 10 directories". Same
    "you cannot see what ran" against "nothing ran" distinction stat_or_reason
    keeps one layer up, so an unreadable candidate is watched and the page says
    what it could not read.
    """
    blocked = False
    for name in (contract.state_name, contract.log_name):
        st, why = stat_or_reason(os.path.join(path, name))
        if st is not None:
            return True
        if why != "missing":
            blocked = True
    return blocked


def discover(roots, contract, depth=SCAN_DEPTH):
    """Results directories under each --root, breadth first.

    A directory qualifies when it holds a state file or a log. Its subtree is
    never descended into: below a results directory sit eggnog/, hmm/,
    structures/ and a foldseek tmp that can reach hundreds of gigabytes, and
    walking that on an NFS mount is exactly the thing a watcher must not do.

    Symlinks ARE followed. `<root>/DatasetA -> /mnt/nvme/DatasetA` is an
    ordinary way to keep one dataset on another mount, and refusing to see it
    made the index quietly short with nothing said. Loop protection does not
    need the exclusion: `seen` already dedupes on realpath, before descending.
    """
    found, seen = [], set()
    for root in roots:
        root = os.path.abspath(root)
        queue = [(root, 0)]
        while queue:
            path, level = queue.pop(0)
            real = os.path.realpath(path)
            if real in seen:
                continue
            seen.add(real)
            if looks_like_results(path, contract):
                found.append(path)
                continue               # do not descend into a results directory
            if level >= depth:
                continue
            try:
                with os.scandir(path) as it:
                    kids = sorted(e.path for e in it
                                  if not e.name.startswith(".")
                                  and e.is_dir())
            except OSError:
                continue
            queue += [(k, level + 1) for k in kids]
    return found


def default_engine():
    """metaannot.py beside this file's parent, which is where it lives."""
    here = os.path.dirname(os.path.abspath(__file__))
    cand = os.path.join(os.path.dirname(here), "metaannot.py")
    return cand if os.path.exists(cand) else None


MIN_INTERVAL = 0.5


def poll_interval(text):
    """--interval, checked. argparse's float accepts 0, -1, nan and inf.

    period() returns the floor unchanged, so `--interval 0` became
    `setTimeout(pull, 0)` - a busy fetch loop against the NFS mount the
    comments in this file worry about, one request per turn of the event loop,
    per open tab.
    """
    try:
        val = float(text)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError("not a number: %r" % text)
    if val != val or val in (float("inf"), float("-inf")):
        raise argparse.ArgumentTypeError("not a usable interval: %r" % text)
    if val < MIN_INTERVAL:
        raise argparse.ArgumentTypeError(
            "must be at least %g seconds; %s would have the page poll faster "
            "than any file it reads can change" % (MIN_INTERVAL, text))
    return val


def build_parser():
    ap = argparse.ArgumentParser(
        prog="console.py",
        description="A read-only watcher for metaannot results directories. "
                    "It never writes into one.")
    ap.add_argument("dirs", nargs="*", metavar="DIR",
                    help="a results directory to watch")
    ap.add_argument("--project", action="append", default=[], metavar="DIR",
                    help="a results directory to watch (repeatable)")
    ap.add_argument("--root", action="append", default=[], metavar="DIR",
                    help="scan DIR for results directories (repeatable)")
    ap.add_argument("--socket", metavar="PATH",
                    default=os.environ.get("METAANNOT_CONSOLE_SOCK"),
                    help="where to bind; default $XDG_RUNTIME_DIR/"
                         "metaannot.sock or the first private directory found")
    ap.add_argument("--metaannot", metavar="PATH", default=default_engine(),
                    help="the engine to ask for the contract")
    ap.add_argument("--python", metavar="PATH", default=sys.executable,
                    help="the interpreter that runs it")
    ap.add_argument("--interval", type=poll_interval, default=3.0, metavar="S",
                    help="the fastest the page will poll, at least 0.5 "
                         "(default 3)")
    ap.add_argument("--version", action="version",
                    version="metaannot console " + CONSOLE_VERSION)
    return ap


def build(args):
    """Contract, then projects. Returns a Console, having served nothing."""
    if not args.metaannot:
        raise Refuse("cannot find metaannot.py; pass --metaannot PATH")
    if not os.path.exists(args.metaannot):
        raise Refuse("no such engine: %s" % args.metaannot)
    contract = Contract(run_describe(args.python, args.metaannot))
    paths = []
    for d in list(args.dirs) + list(args.project):
        d = os.path.abspath(d)
        if not os.path.isdir(d):
            raise Refuse("not a directory: %s" % d)
        if not looks_like_results(d, contract):
            sys.stderr.write(
                "console: %s holds neither %s nor %s; watching it anyway, "
                "because a run that has not started yet looks exactly like "
                "this\n" % (d, contract.state_name, contract.log_name))
        paths.append(d)
    paths += [d for d in discover(args.root, contract) if d not in paths]
    projects = [Project(i, p, contract) for i, p in enumerate(paths)]
    return Console(contract, projects, args.python, args.metaannot,
                   args.interval, roots=args.root)


def inside(path, other):
    """True when `path` is `other`, or lies inside it. Real paths both, so a
    symlink is not a way round it."""
    path, other = os.path.realpath(path), os.path.realpath(other)
    return path == other or path.startswith(other.rstrip(os.sep) + os.sep)


# What the operator is told instead of `mkdir -m 700 -p <the parent they
# passed>`, which for `--socket /tmp/pwn.sock` read "mkdir -m 700 -p /tmp":
# advice they must not take, a silent no-op on an existing directory, and for a
# --socket inside a results directory, the console explaining how to make the
# violation happen. This names a directory of the console's own.
SOCKET_ADVICE = (
    "Leave --socket off entirely and the console picks a private directory of "
    "its own ($XDG_RUNTIME_DIR, /run/user/<uid>, ~/.cache/metaannot-console).\n"
    "To choose one yourself, make it first — this console will not create a "
    "directory for you:\n"
    "    mkdir -p ~/.cache/metaannot-console\n"
    "    chmod 700 ~/.cache/metaannot-console\n"
    "    --socket ~/.cache/metaannot-console/metaannot.sock")


def refuse_if_watched(sock, args):
    """`sock` back, unless it is in or under something this console watches.

    The watched set comes off `args` rather than off the built Console, so this
    cannot be defeated by a caller that forgets to pass it, and it covers a
    --root subtree whose results directories do not exist yet - precisely the
    case that let the console CREATE one. The default path goes through here
    too: $METAANNOT_CONSOLE_SOCK and $XDG_RUNTIME_DIR are both somebody's
    environment, and the guarantee is not "unless you asked for it".
    """
    where = os.path.dirname(sock) or "."
    for watched, what in ([(d, "watching") for d in
                           list(args.dirs) + list(args.project)]
                          + [(r, "scanning for results directories")
                             for r in args.root]):
        if inside(sock, watched) or inside(where, watched):
            raise Refuse(
                "--socket %s is inside %s, which this console is %s.\n"
                "This console does not write into a directory it watches, and "
                "binding there writes two things: the socket itself and "
                "%s.lock, which outlives the process.\n%s"
                % (sock, os.path.abspath(watched), what,
                   os.path.basename(sock), SOCKET_ADVICE))
    return sock


def socket_path(args):
    """Where to bind. Never inside a watched directory, never one this console
    would have to create, never one another account can reach.

    The banner's first line is "read-only. Nothing is written to any results
    directory", and `--socket <results>/console.sock` made it false in the most
    direct way available: the console bound a socket and created
    `console.sock.lock` inside the directory it was watching, moving that
    directory's mtime, ctime, nlink and size, and leaving the lock behind after
    a clean SIGTERM. take_socket_name's docstring already claimed the lock
    lives "never in a results directory"; nothing enforced it.

    Worse than accepting one: MAKING one. runtime_dir's os.makedirs created the
    parent of an explicit --socket, so `--socket <root>/DatasetZ/results/c.sock`
    produced a results directory at 0700 with a stray lock file in it, which a
    run started into that path afterwards then finds. So nothing here creates a
    directory at all: a --socket whose directory is not already there is
    refused, and the operator makes it where they choose.

    The DEFAULT path had the same defect one level down, in the order of two
    calls. `refuse_if_watched(os.path.join(runtime_dir(), ...))` runs
    runtime_dir FIRST, and runtime_dir -> ensure_private_dir -> os.mkdir: with
    XDG_RUNTIME_DIR set to <results>/rt the console CREATED that directory
    inside the tree it watches and only then refused to use it. Refusing after
    creating is not refusing. So the veto goes INSIDE the candidate loop, ahead
    of the mkdir: a candidate under a watched tree is refused while it is still
    only a string.
    """
    if not args.socket:
        return refuse_if_watched(
            os.path.join(runtime_dir(veto=lambda cand: refuse_if_watched(
                os.path.join(cand, "metaannot.sock"), args)),
                "metaannot.sock"), args)
    sock = refuse_if_watched(
        os.path.abspath(os.path.expanduser(args.socket)), args)
    parent = os.path.dirname(sock) or "."
    # realpath, then check: what the kernel binds in is the target, so that is
    # what has to be private, and `/tmp` reported as "is a symlink" is a true
    # sentence about the wrong thing. Nothing here creates it, so resolving a
    # link the operator named cannot create one either.
    real = os.path.realpath(parent)
    named = parent if real == parent else "%s (%s)" % (parent, real)
    uid = os.getuid()
    why = why_not_private(real, uid)
    if why is None:
        anc, anc_why = exposed_ancestor(real, uid)
        if anc is not None:
            why = "sits under %s, which %s" % (anc, anc_why)
        else:
            soft, soft_why = exposed_ancestor(real, uid, 0o020)
            if soft is not None:
                sys.stderr.write(
                    "console: %s is 0700, but %s %s. If that group is more "
                    "than you, they can replace the directory the socket is "
                    "in; this console cannot tell how many members it has.\n"
                    % (real, soft, soft_why))
    if why is not None:
        raise Refuse(
            "--socket %s will not be used: its directory %s %s.\n"
            "File permissions are this console's whole access control, so the "
            "socket has to sit in a private directory this account owns.\n%s"
            % (sock, named, why, SOCKET_ADVICE))
    return sock


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        console = build(args)
        sock = socket_path(args)
        return serve(console, sock,
                     announce=lambda: banner(console, sock))
    except Refuse as e:
        sys.stderr.write("console: %s\n" % e)
        return 2
    except KeyboardInterrupt:
        return 0


if __name__ == "__main__":
    sys.exit(main())
