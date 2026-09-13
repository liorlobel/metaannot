"""Shared fixtures: the module under test, temporary projects, tool stubs.

The whole suite runs offline. Where a stage shells out to hmmsearch, DIAMOND,
MMseqs2 or Foldseek, the binary is replaced by a script on PATH that writes a
canned file, so the stage's own plumbing (atomic writes, adoption, signatures,
argument assembly) is exercised without the tool.
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import stat
import subprocess
import sys
import textwrap

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METAANNOT_PY = os.path.join(ROOT, "metaannot.py")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: excluded from the default run")
    config.addinivalue_line("markers", "R: needs Rscript and R packages")


def _load():
    spec = importlib.util.spec_from_file_location("metaannot", METAANNOT_PY)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["metaannot"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def ma():
    """The tool, imported as a module."""
    return _load()


@pytest.fixture(autouse=True)
def _fresh_ownership_watch(ma):
    """Start every test with the module's ownership watch disarmed.

    `_STATE_WATCH`, `_DECLARED_OUTPUTS` and `_STOPPING` are process-wide
    because a process runs exactly one metaannot run; a test session runs
    hundreds in the same interpreter. Without this, one test that drove a run
    to the point of being superseded would leave `lost` set for the rest of
    the session, and every
    later atomic_out() in the process would decline its rename - a green suite
    turning red three files away, for a reason nothing in the failing test
    mentions. Reset on the way in as well as out, so a test that fails half
    way through does not take the next one with it.
    """
    def clear():
        ma._watch_state(None, None, ())
        ma._DECLARED_OUTPUTS.clear()
        # The stopping latch, for the same reason and with a sharper edge: one
        # test that drives an in-process run to an interrupt would leave it
        # set, and every _tool_process() in the session afterwards would
        # refuse to launch anything - a green suite turning red three files
        # away, for a reason nothing in the failing test mentions.
        ma._STOPPING = False

    clear()
    yield
    clear()


CONSOLE_PY = os.path.join(ROOT, "console", "console.py")


@pytest.fixture(scope="session")
def console():
    """The console, imported as a module.

    It is a plain module with pure readers and a Console object, so the tests
    drive it directly rather than through the CLI. Registered under a name of
    its own: nothing in the console imports metaannot, and nothing here should
    make it look as though it did.
    """
    spec = importlib.util.spec_from_file_location("metaannot_console",
                                                  CONSOLE_PY)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["metaannot_console"] = mod
    spec.loader.exec_module(mod)
    return mod


# ----------------------------------------------------------------------
# R availability
# ----------------------------------------------------------------------
R_PACKAGES_CACHE = {}


def r_has(*packages):
    """True when Rscript exists and every named package is installed."""
    if not shutil.which("Rscript"):
        return False
    want = [p for p in packages if p not in R_PACKAGES_CACHE]
    if want:
        expr = ('cat(paste(vapply(c(%s), function(p) paste0(p, "=", '
                'requireNamespace(p, quietly=TRUE)), character(1)), '
                'collapse=" "))' % ",".join(f'"{p}"' for p in want))
        try:
            r = subprocess.run(["Rscript", "-e", expr], capture_output=True,
                               text=True, timeout=300)
        except (OSError, subprocess.SubprocessError):
            for p in want:
                R_PACKAGES_CACHE[p] = False
            return False
        for tok in (r.stdout or "").split():
            name, _, val = tok.partition("=")
            R_PACKAGES_CACHE[name] = (val == "TRUE")
        for p in want:
            R_PACKAGES_CACHE.setdefault(p, False)
    return all(R_PACKAGES_CACHE.get(p, False) for p in packages)


def needs_r(*packages):
    """Mark a test as needing R, and skip it cleanly when R is not there.

    The `R` mark is applied ALWAYS, so `pytest -m R` selects the same set of
    tests on every machine and a missing package shows up as a skip with its
    name rather than as a test that quietly disappeared.
    """
    reason = None
    if not shutil.which("Rscript"):
        reason = "Rscript not on PATH"
    else:
        missing = [p for p in packages if not r_has(p)]
        if missing:
            reason = "missing R package(s): " + ", ".join(missing)

    def decorate(fn):
        fn = pytest.mark.R(fn)
        if reason:
            fn = pytest.mark.skip(reason=reason)(fn)
        return fn
    return decorate


# ----------------------------------------------------------------------
# tool stubs
# ----------------------------------------------------------------------
STUBS = {
    # hmmsearch writes whatever --tblout / --domtblout names.
    "hmmsearch": r"""
import subprocess, sys, os, time
a = sys.argv[1:]
def val(flag):
    return a[a.index(flag) + 1] if flag in a else None
# STUB_PIDFILE records every process this stub is responsible for, one
# "<role> <pid>" line each, so a test can ask the PROCESS TABLE whether a real
# child of a real run is still alive after that run has been killed. Nothing
# else in the suite could answer that question.
pidfile = os.environ.get("STUB_PIDFILE")
def note(role, pid):
    if pidfile:
        with open(pidfile, "a", encoding="utf-8") as fh:
            fh.write(role + " " + str(pid) + "\n")
            fh.flush()
note("tool", os.getpid())
# STUB_FORK_CHILD gives the stub a GRANDCHILD of the run, which is the shape
# proc.kill() never reached: interproscan.sh -> java, emapper -> its children,
# torch -> its dataloader workers. It outlives this stub on purpose.
if os.environ.get("STUB_FORK_CHILD"):
    subprocess.Popen([sys.executable, "-c",
                      "import os,sys,time\n"
                      "p = sys.argv[1]\n"
                      "if p:\n"
                      "    fh = open(p, 'a')\n"
                      "    fh.write('child ' + str(os.getpid()) + chr(10))\n"
                      "    fh.close()\n"
                      "time.sleep(float(sys.argv[2]))\n",
                      pidfile or "",
                      os.environ.get("STUB_CHILD_SLEEP", "30")])
# STUB_GROW appends to the output the run handed us, over and over, so a test
# can produce the thing the leftover census is about: a real orphaned tool
# writing into a real results directory, its `.part` file growing, after the
# run that launched it is gone.
if os.environ.get("STUB_GROW") and val("--tblout"):
    grow = open(val("--tblout"), "a", encoding="utf-8")
# STUB_WAIT_FOR blocks the stage until the test creates that path, which is
# what turns "a run that unwinds slowly" from a sleep race into something a
# test can time exactly: the run stays inside the stage, and so inside the
# executor's shutdown wait, until the test says otherwise. STUB_GATE_MAX_S is
# a backstop and not a feature: a regression that stops the gate from ever
# being reached must FAIL a test rather than wedge the suite for ever.
gate = os.environ.get("STUB_WAIT_FOR")
deadline = time.time() + float(os.environ.get("STUB_GATE_MAX_S", "120"))
while gate and not os.path.exists(gate) and time.time() < deadline:
    if os.environ.get("STUB_GROW") and val("--tblout"):
        grow.write("# still searching\n")
        grow.flush()
    time.sleep(0.02)
time.sleep(float(os.environ.get("STUB_SLEEP", "0")))
out = val("--tblout")
dom = val("--domtblout")
# The last two positionals are <hmm library> <query fasta>.
faa = a[-1]
ids = [l[1:].split()[0] for l in open(faa, encoding="utf-8") if l.startswith(">")]
canned = os.environ.get("STUB_HMM_HITS", "")
want = set(canned.split(",")) if canned else set(ids[:2])
if out:
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("# target name        accession  query name\n")
        for i in ids:
            if i in want:
                fh.write(" ".join([i, "-", "Peptidase_S8", "PF00082.1",
                                   "1e-30", "100.0", "0.0", "1e-30", "100.0",
                                   "0.0", "1.0", "1", "0", "0", "1", "1", "1",
                                   "1", "-"]) + "\n")
if dom:
    with open(dom, "w", encoding="utf-8") as fh:
        fh.write("# target name accession tlen query name\n")
        for i in ids:
            if i in want:
                fh.write(" ".join([i, "-", "300", "GH13.hmm", "-", "100",
                                   "1e-20", "80.0", "0.0", "1", "1", "1e-20",
                                   "1e-20", "80.0", "0.0", "1", "80", "10",
                                   "200", "10", "200", "0.95", "-"]) + "\n")
""",
    "diamond": r"""
import sys, os, re
a = sys.argv[1:]
def val(flag):
    return a[a.index(flag) + 1] if flag in a else None
if a and a[0] == "dbinfo":
    # Real dbinfo reports Sequences and Letters; F.write_dmnd puts them in the
    # first line of the stand-in file so the shape of a database is a test
    # parameter rather than something diamond has to be installed to know.
    txt = open(val("-d"), encoding="utf-8", errors="replace").readline()
    m = re.search(r"sequences=(\d+) letters=(\d+)", txt)
    n, L = (m.group(1), m.group(2)) if m else ("5000", "1750000")
    print("          Database type  Diamond database")
    print("              Sequences  " + n)
    print("                Letters  " + L)
    sys.exit(0)
out = val("-o")
faa = val("-q")
ids = [l[1:].split()[0] for l in open(faa, encoding="utf-8") if l.startswith(">")]
with open(out, "w", encoding="utf-8") as fh:
    for i in ids[:2]:
        fh.write("\t".join([i, "VFG0001", "88.0", "150", "1e-40", "300.0",
                            "90", "80", "hemolysin BL binding component"])
                 + "\n")
""",
    "mmseqs": r"""
import sys, os
a = sys.argv[1:]
faa, pref = a[1], a[2]
ids = [l[1:].split()[0] for l in open(faa, encoding="utf-8") if l.startswith(">")]
with open(pref + "_cluster.tsv", "w", encoding="utf-8") as fh:
    rep = ids[0] if ids else "none"
    for i in ids:
        fh.write(f"{rep if i in ids[:2] else i}\t{i}\n")
open(pref + "_rep_seq.fasta", "w", encoding="utf-8").close()
""",
    # A tool that exits 0 and writes nothing; used to prove atomic_out refuses.
    "silent_tool": r"""
import sys
sys.exit(0)
""",
}


@pytest.fixture
def stub_bin(tmp_path, monkeypatch):
    """Put fake search tools on PATH. Returns the directory."""
    d = tmp_path / "stubbin"
    d.mkdir()
    for name, body in STUBS.items():
        path = d / name
        path.write_text("#!/usr/bin/env python3\n" + textwrap.dedent(body),
                        encoding="utf-8")
        path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP
                   | stat.S_IXOTH)
    # The stubs are python scripts; make sure the interpreter running the
    # tests is the one they get.
    py = d / "python3"
    if not py.exists():
        try:
            os.symlink(sys.executable, py)
        except OSError:
            pass
    if os.name == "nt":
        # Windows will not run an extension-less file no matter where it sits
        # on PATH: only the suffixes in PATHEXT are executable. Without a .cmd
        # beside each stub the whole fixture is inert here, and every test that
        # depends on it silently reaches for the real tool instead. That makes
        # the suite pass or fail according to what happens to be installed on
        # the machine, which is worse than failing outright, so give each stub
        # a shim rather than skipping these tests on Windows.
        for name in STUBS:
            (d / (name + ".cmd")).write_text(
                "@echo off\r\n"
                f'"{sys.executable}" "%~dp0{name}" %*\r\n',
                encoding="utf-8")
    monkeypatch.setenv("PATH", str(d) + os.pathsep + os.environ["PATH"])
    return d


# ----------------------------------------------------------------------
# projects
# ----------------------------------------------------------------------
class Project:
    """A directory holding config.yaml, inputs and results."""

    def __init__(self, root, cfg):
        self.root = str(root)
        self.cfg = cfg
        self.config_path = os.path.join(self.root, "config.yaml")

    @property
    def results(self):
        rd = self.cfg.get("results_dir", "results")
        return rd if os.path.isabs(rd) else os.path.join(self.root, rd)

    def path(self, *parts):
        return os.path.join(self.root, *parts)

    def rpath(self, *parts):
        return os.path.join(self.results, *parts)

    def write_config(self, **overrides):
        import yaml
        cfg = dict(self.cfg)
        cfg.update(overrides)
        self.cfg = cfg
        with open(self.config_path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(cfg, fh, sort_keys=False)
        return self.config_path

    def run(self, *args, expect=0, env=None, cwd=None, timeout=300):
        return run_metaannot("run", "--config", self.config_path, *args,
                             expect=expect, env=env, cwd=cwd or self.root,
                             timeout=timeout)

    def state(self):
        with open(os.path.join(self.results, ".metaannot_state.json"),
                  encoding="utf-8") as fh:
            return json.load(fh)


def run_metaannot(*args, expect=0, env=None, cwd=None, timeout=300):
    """Invoke the CLI in a subprocess.

    A subprocess, not an in-process call, because cmd_run holds the results
    lock until the process exits: two sequential in-process runs against one
    results directory would refuse each other.
    """
    e = dict(os.environ)
    e.setdefault("PYTHONHASHSEED", "0")
    if env:
        e.update(env)
    proc = subprocess.run([sys.executable, METAANNOT_PY, *[str(a) for a in args]],
                          capture_output=True, text=True, env=e, cwd=cwd,
                          timeout=timeout)
    if expect is not None and proc.returncode != expect:
        raise AssertionError(
            f"metaannot {' '.join(str(a) for a in args)} exited "
            f"{proc.returncode}, expected {expect}\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}")
    return proc


@pytest.fixture
def cli():
    return run_metaannot


def build_project(root, proteins=None, samples=("A_1", "A_2", "B_1", "B_2"),
                  fractions=1, n_extra=0, **cfg_overrides):
    """A complete, runnable project: fasta, eggNOG table, peptides, manifest."""
    import fixtures as F

    proteins = proteins if proteins is not None else F.protein_set(n_extra=n_extra)
    root = str(root)
    os.makedirs(os.path.join(root, "input"), exist_ok=True)
    faa = F.write_fasta(os.path.join(root, "input", "proteins.faa"), proteins)
    emp = F.write_emapper(os.path.join(root, "input", "cat.emapper.annotations"),
                          proteins)
    quant = os.path.join(root, "input", "combined_peptide.tsv")
    F.write_peptide_table(quant, proteins, list(samples))
    groups = {s: s.split("_")[0] for s in samples}
    man = F.simple_manifest(os.path.join(root, "input", "experiment.fp-manifest"),
                            list(samples), groups, fractions=fractions)
    cfg = {
        "proteins_faa": faa,
        "quant_table": quant,
        "quant_format": "fragpipe_peptide",
        "manifest": man,
        "emapper_precomputed": [emp],
        "results_dir": "results",
        "threads": 2,
        "stage_workers": 2,
        "ram_gb": 4,
        "run": {"eggnog": True, "pfam": False, "dbcan": False,
                "diamond": False, "cluster": False, "join": True,
                "topology": False, "structure": False, "context": False,
                "unipept": False, "taxonomy": False, "ncbifam": False,
                "kofam": False, "interpro": False, "hhblits": False,
                "jackhmmer": False, "smorf": False, "effectors": False},
    }
    cfg.update(cfg_overrides)
    p = Project(root, cfg)
    p.write_config()
    p.proteins = proteins
    p.samples = list(samples)
    return p


@pytest.fixture
def project(tmp_path):
    """A minimal project whose default run needs no external tool."""
    return build_project(tmp_path / "proj")


@pytest.fixture
def paths_for(ma, tmp_path):
    """A Paths object over a fresh results directory."""
    def make(name="results"):
        cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
        cfg["results_dir"] = str(tmp_path / name)
        p = ma.Paths(cfg)
        p.mkdirs()
        return cfg, p
    return make
