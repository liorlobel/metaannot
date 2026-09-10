"""Config loading, path resolution, cache keys and CLI argument handling."""
from __future__ import annotations

import json
import io
import os

import pytest
import yaml

from conftest import build_project, run_metaannot


# --- finding 1 --------------------------------------------------------
def test_misspelled_run_key_is_reported_with_a_suggestion(ma, tmp_path, capsys):
    # symptom: `run: {unipep: true}` left the stage disabled with nothing said.
    cfg = tmp_path / "c.yaml"
    cfg.write_text("run:\n  unipep: true\n", encoding="utf-8")
    ma.load_config(str(cfg))
    err = capsys.readouterr().err
    assert "unrecognised key 'run.unipep'" in err
    assert "did you mean 'run.unipept'" in err
    assert "NOT in effect" in err


def test_misspelled_analysis_key_is_reported(ma, tmp_path, capsys):
    # symptom: analysis: keys are the Rmd's params, so `fdrr: 0.01` used to be
    # written into the header as a spurious param while `fdr` stayed default.
    cfg = tmp_path / "c.yaml"
    cfg.write_text("analysis:\n  fdrr: 0.01\n", encoding="utf-8")
    ma.load_config(str(cfg))
    err = capsys.readouterr().err
    assert "unrecognised key 'analysis.fdrr'" in err


@pytest.mark.parametrize("block,key", [
    ("tool_args", "hmmsearch"),
    ("diamond_weights", "my_own_db"),
])
def test_freeform_blocks_do_not_warn_about_their_own_keys(ma, tmp_path,
                                                          block, key):
    # symptom: the key check must stop at blocks whose keys the user chooses,
    # or every legitimate tool/database/predictor name is a "typo".
    body = {block: {key: (["--flag"] if block == "tool_args"
                          else 3 if block == "diamond_weights"
                          else {"file": "x.tsv", "score_col": "s"})}}
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump(body), encoding="utf-8")
    assert ma.unknown_keys(body, ma.DEFAULT_CONFIG) == []


@pytest.mark.parametrize("block", ["db.diamond", "sources.diamond"])
def test_nested_freeform_database_blocks_do_not_warn(ma, block):
    top, sub = block.split(".")
    body = {top: {sub: {"my_db": "/data/db/mine.dmnd"}}}
    assert ma.unknown_keys(body, ma.DEFAULT_CONFIG) == []


def test_unknown_key_inside_a_checked_block_is_still_found(ma):
    body = {"thresholds": {"diamond_evalu": 1e-5}}
    assert ma.unknown_keys(body, ma.DEFAULT_CONFIG) == ["thresholds.diamond_evalu"]


# --- finding 2 --------------------------------------------------------
def test_relative_config_paths_resolve_against_the_config_file(tmp_path):
    # symptom: a config that only works from the directory it was written in.
    proj = build_project(tmp_path / "proj")
    # rewrite every path as relative to the config file
    rel = {
        "proteins_faa": "input/proteins.faa",
        "quant_table": "input/combined_peptide.tsv",
        "manifest": "input/experiment.fp-manifest",
        "emapper_precomputed": ["input/cat.emapper.annotations"],
        "results_dir": "results",
    }
    proj.write_config(**rel)

    hashes = []
    for cwd in (str(tmp_path), proj.root, os.path.expanduser("~")):
        # a fresh results dir each time, so this measures resolution, not cache
        import shutil
        shutil.rmtree(proj.results, ignore_errors=True)
        run_metaannot("run", "--config", proj.config_path, cwd=cwd)
        with open(proj.rpath("annotation_final.tsv"), encoding="utf-8") as fh:
            hashes.append(fh.read())
    assert hashes[0] == hashes[1] == hashes[2]
    assert os.path.isdir(proj.results), "results landed beside the config"


def test_home_and_env_vars_in_a_config_path_are_expanded(ma, tmp_path,
                                                         monkeypatch):
    # symptom: `~/db/Pfam-A.hmm` became <configdir>/~/db/Pfam-A.hmm and doctor
    # reported a database the user does have as missing.
    # Compared through abspath/normcase because resolve_paths makes every
    # db path absolute: on Windows "/somewhere/db" acquires the current drive
    # and "~/x" comes back from expanduser with a mixed separator, and
    # neither of those is the bug this test is about.
    monkeypatch.setenv("MA_TEST_DB", str(tmp_path / "db"))
    cfg = tmp_path / "c.yaml"
    cfg.write_text("db:\n  pfam_hmm: '$MA_TEST_DB/Pfam-A.hmm'\n"
                   "  dbcan_hmm: '~/dbcan.txt'\n", encoding="utf-8")
    c = ma.load_config(str(cfg))

    def same(a, b):
        return os.path.normcase(os.path.abspath(a)) == \
            os.path.normcase(os.path.abspath(b))
    assert same(c["db"]["pfam_hmm"], str(tmp_path / "db" / "Pfam-A.hmm"))
    assert same(c["db"]["dbcan_hmm"], os.path.expanduser("~/dbcan.txt"))
    assert "~" not in c["db"]["dbcan_hmm"]
    assert "$MA_TEST_DB" not in c["db"]["pfam_hmm"]


# --- finding 3 --------------------------------------------------------
def test_tool_version_bump_does_not_invalidate_the_stage_cache(ma, project,
                                                               monkeypatch):
    # symptom: the cache was keyed to __version__, so a patch release that
    # fixed a log message discarded days of InterProScan and Foldseek compute.
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    st = [s for s in ma.STAGES if s["name"] == "emapper"][0]
    before = ma.signature(st, cfg, p)
    monkeypatch.setattr(ma, "__version__", "99.9.9")
    assert ma.signature(st, cfg, p) == before


def test_signature_version_bump_does_invalidate_the_stage_cache(ma, project,
                                                                monkeypatch):
    # the other half: SIGNATURE_VERSION exists precisely so a real change of
    # meaning CAN discard cached results.
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    st = [s for s in ma.STAGES if s["name"] == "emapper"][0]
    before = ma.signature(st, cfg, p)
    monkeypatch.setattr(ma, "SIGNATURE_VERSION", ma.SIGNATURE_VERSION + 1)
    assert ma.signature(st, cfg, p) != before


# --- finding 4 --------------------------------------------------------
@pytest.mark.parametrize("given,want", [
    ("512M", 1), ("64G", 64), ("64GB", 64), ("64", 64), (64, 64),
    ("65536M", 64), ("1T", 1024), ("1500M", 1),
])
def test_a_small_explicit_ram_budget_never_collapses_to_auto_detect(ma, given,
                                                                    want):
    # symptom: 512M truncated to 0, and 0 means "auto-detect", so an explicit
    # budget was silently replaced by 80% of the machine.
    assert ma.parse_ram(given) == want


@pytest.mark.parametrize("given", ["", None, 0, "0"])
def test_an_absent_ram_budget_still_means_auto_detect(ma, given):
    assert ma.parse_ram(given) == 0


def test_an_unreadable_ram_budget_is_a_message_not_a_traceback(ma):
    with pytest.raises(ma.StageError) as e:
        ma.parse_ram("sixty-four")
    assert "could not read a memory size" in str(e.value)


# --- finding 5 --------------------------------------------------------
@pytest.mark.parametrize("total,workers", [(2, 4), (1, 8), (3, 4)])
def test_per_stage_ram_floors_at_one_gb(ma, project, total, workers):
    # symptom: integer division to 0 GB disabled every memory flag (diamond
    # -b, mmseqs --split-memory-limit, hhblits -maxmem) on a small machine.
    cfg = ma.load_config(project.config_path)
    per = max(1, total // workers)
    c = ma.stage_cfg(cfg, 1, per)
    assert c["ram_gb"] >= 1
    # and a genuinely absent budget stays absent rather than becoming 1
    assert ma.stage_cfg(cfg, 1, 0)["ram_gb"] == 0


# --- finding 6 --------------------------------------------------------
def test_config_pointing_at_a_directory_is_one_line_not_a_traceback(tmp_path):
    proc = run_metaannot("run", "--config", str(tmp_path), expect=1)
    assert "Traceback" not in proc.stderr
    assert "is a directory" in proc.stderr + proc.stdout


def test_results_dir_pointing_at_a_file_is_one_line_not_a_traceback(tmp_path,
                                                                    project):
    f = tmp_path / "not_a_dir"
    f.write_text("x", encoding="utf-8")
    proc = run_metaannot("run", "--config", project.config_path,
                         "--results-dir", str(f), expect=1)
    assert "Traceback" not in proc.stderr
    assert "exists and is not a directory" in proc.stderr


def test_unknown_stage_name_is_one_line_not_a_traceback(project):
    proc = run_metaannot("run", "--config", project.config_path,
                         "--only", "pfamm", expect=1)
    assert "Traceback" not in proc.stderr
    assert "unknown stage 'pfamm'" in proc.stderr


# --- finding 7 --------------------------------------------------------
def test_duplicate_top_level_run_block_is_refused(tmp_path):
    # symptom: yaml.safe_load keeps only the LAST mapping with a given key, so
    # the first run: block was discarded, every setting in it reverted to the
    # default, and the key check could not see it — the dict had collapsed.
    cfg = tmp_path / "c.yaml"
    cfg.write_text("run:\n  structure: false\n"
                   "threads: 4\n"
                   "run:\n  unipept: true\n", encoding="utf-8")
    proc = run_metaannot("doctor", "--config", str(cfg), expect=1)
    out = proc.stdout + proc.stderr
    assert "duplicate key 'run'" in out
    assert "Traceback" not in proc.stderr


def test_duplicate_nested_key_is_also_refused(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text("db:\n  pfam_hmm: /a\n  pfam_hmm: /b\n", encoding="utf-8")
    proc = run_metaannot("doctor", "--config", str(cfg), expect=1)
    assert "duplicate key 'pfam_hmm'" in proc.stdout + proc.stderr


def test_an_empty_block_keeps_the_defaults_and_says_so(ma, tmp_path, capsys):
    # `db:` with every child commented out parses as None; letting that
    # replace the block turned an ordinary edit into a NoneType traceback.
    cfg = tmp_path / "c.yaml"
    cfg.write_text("db:\n", encoding="utf-8")
    c = ma.load_config(str(cfg))
    assert c["db"]["pfam_hmm"]
    assert "keeping the built-in defaults" in capsys.readouterr().err


# --- finding 24 -------------------------------------------------------
def test_no_stage_config_key_resolves_to_none(ma):
    # symptom: a key whose value is None hashes identically for ever, so that
    # stage can never be invalidated by a change to it.
    bad = []
    for st in ma.STAGES:
        for k in st["keys"]:
            if ma._dig(ma.DEFAULT_CONFIG, k) is None:
                bad.append((st["name"], k))
    assert bad == [], f"stage keys missing from DEFAULT_CONFIG: {bad}"


def test_every_stage_enabled_flag_exists_in_the_default_run_block(ma):
    missing = [st["name"] for st in ma.STAGES
               if st["enabled"] and st["enabled"] not in ma.DEFAULT_CONFIG["run"]]
    assert missing == []


def test_deep_merge_keeps_defaults_the_user_did_not_override(ma):
    # symptom: a shallow update dropped every default inside a nested block a
    # user partially overrode.
    merged = ma.deep_merge(ma.DEFAULT_CONFIG, {"thresholds": {"smorf_max_len": 50}})
    assert merged["thresholds"]["smorf_max_len"] == 50
    assert merged["thresholds"]["diamond_evalue"] == \
        ma.DEFAULT_CONFIG["thresholds"]["diamond_evalue"]


def test_diamond_database_block_replaces_rather_than_merges(ma, capsys):
    # symptom: listing only vfdb still left four /data/db/*.dmnd phantoms in
    # doctor, in the signature and in the stage's "missing" warnings.
    merged = ma.deep_merge(ma.DEFAULT_CONFIG,
                           {"db": {"diamond": {"vfdb": "/x/vfdb.dmnd"}}})
    assert set(merged["db"]["diamond"]) == {"vfdb"}
    assert "are NOT used" in capsys.readouterr().err


def test_init_writes_a_config_that_loads_back_without_warnings(tmp_path, ma,
                                                               capsys):
    out = tmp_path / "config.yaml"
    run_metaannot("init", "--out", str(out))
    assert out.exists()
    ma.load_config(str(out))
    err = capsys.readouterr().err
    assert "unrecognised key" not in err


# ======================================================================
# doctor and install
# ======================================================================
import shutil                                                    # noqa: E402
import stat                                                      # noqa: E402
import subprocess                                                # noqa: E402
import sys                                                       # noqa: E402

from conftest import METAANNOT_PY                                # noqa: E402


def _doctor_project(tmp_path, **cfg):
    """A config whose only requirement is one HMM database, so doctor's
    output is short and deterministic."""
    proj = build_project(tmp_path / "doc")
    base = {"run": {"eggnog": True, "pfam": True, "dbcan": False,
                    "diamond": False, "cluster": False, "join": True,
                    "topology": False, "structure": False, "context": False,
                    "unipept": False, "taxonomy": False, "ncbifam": False,
                    "kofam": False, "interpro": False, "hhblits": False,
                    "jackhmmer": False, "smorf": False, "effectors": False}}
    base.update(cfg)
    proj.write_config(**base)
    return proj


def _no_r_env():
    """PATH without Rscript, so doctor's verdict does not depend on which R
    packages happen to be installed on the machine running the tests."""
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep)
             if p and not os.path.exists(os.path.join(p, "Rscript"))]
    return {"PATH": os.pathsep.join(parts)}


# --- finding 8 --------------------------------------------------------
def test_a_database_path_containing_a_space_survives_the_install_plan(
        tmp_path):
    # symptom: everything spliced into these commands comes from the config —
    # paths and URLs the user wrote. Unquoted, a path with a space silently
    # breaks every command.
    db = tmp_path / "my databases"
    proj = _doctor_project(tmp_path, db={"pfam_hmm": str(db / "Pfam-A.hmm")})
    plan = tmp_path / "install.sh"
    run_metaannot("doctor", "--config", proj.config_path,
                  "--install-plan", str(plan), expect=None, env=_no_r_env())
    text = plan.read_text(encoding="utf-8")
    assert "'" + str(db / "Pfam-A.hmm") + "'" in text or \
        str(db / "Pfam-A.hmm").replace(" ", "\\ ") in text
    # and the plan is still a valid script
    assert subprocess.run(["bash", "-n", str(plan)]).returncode == 0
    mkdirs = [l for l in text.splitlines() if l.startswith("mkdir -p")]
    assert mkdirs
    assert subprocess.run(["bash", "-c", mkdirs[0]]).returncode == 0
    assert db.is_dir(), "a path with a space did not survive quoting"


def test_a_database_path_containing_shell_metacharacters_is_a_literal_name(
        tmp_path):
    # symptom: unquoted, a path containing `;` or `$(...)` EXECUTES.
    evil = tmp_path / "we ; touch EVIL"
    proj = _doctor_project(tmp_path,
                           db={"pfam_hmm": str(evil / "Pfam-A.hmm")})
    plan = tmp_path / "install.sh"
    run_metaannot("doctor", "--config", proj.config_path,
                  "--install-plan", str(plan), expect=None, env=_no_r_env())
    text = plan.read_text(encoding="utf-8")
    mkdirs = [l for l in text.splitlines() if l.startswith("mkdir -p")]
    assert mkdirs
    subprocess.run(["bash", "-c", mkdirs[0]], cwd=str(tmp_path))
    assert evil.is_dir(), "the metacharacters were not treated as a name"
    assert not (tmp_path / "EVIL").exists(), "the config executed a command"


def test_a_source_url_is_quoted_too(ma, tmp_path):
    proj = _doctor_project(
        tmp_path,
        db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")},
        sources={"pfam": "http://x/y;touch URLEVIL"})
    plan = tmp_path / "install.sh"
    run_metaannot("doctor", "--config", proj.config_path,
                  "--install-plan", str(plan), expect=None, env=_no_r_env())
    text = plan.read_text(encoding="utf-8")
    assert "'http://x/y;touch URLEVIL'" in text
    assert subprocess.run(["bash", "-n", str(plan)]).returncode == 0


def test_q_quotes_everything_it_is_given(ma):
    assert ma.q("/a b/c") == "'/a b/c'"
    assert ma.q("x;rm -rf /") == "'x;rm -rf /'"
    assert ma.q("plain") == "plain"


# --- finding 9 --------------------------------------------------------
def test_an_unset_database_path_is_reported_not_guessed(tmp_path):
    # symptom: an unset path used to produce `curl -o .gz <url>` and
    # `hmmpress -f ` with no argument, writing a file called `.gz` into the
    # working directory.
    proj = _doctor_project(tmp_path, db={"pfam_hmm": ""})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=None,
                         env=_no_r_env())
    assert "MANUAL" in proc.stdout
    assert "no path configured" in proc.stdout
    assert "set db.pfam_hmm" in proc.stdout
    plan = tmp_path / "i.sh"
    run_metaannot("doctor", "--config", proj.config_path, "--install-plan",
                  str(plan), expect=None, env=_no_r_env())
    text = plan.read_text(encoding="utf-8")
    assert "-o .gz" not in text
    assert "curl" not in text, "an unconfigured database produced a download"


def test_a_diamond_database_with_no_url_is_reported_as_manual(tmp_path):
    # the tag is configured and the file is absent, but sources.diamond has no
    # URL for it: doctor says where to get it rather than inventing a download.
    proj = _doctor_project(tmp_path,
                           run=dict(_doctor_project(tmp_path).cfg["run"],
                                    diamond=True),
                           db={"pfam_hmm": str(tmp_path / "p.hmm"),
                               "diamond": {"merops": str(tmp_path / "m.dmnd")}})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "DIAMOND merops" in proc.stdout
    assert "no URL in sources.diamond" in proc.stdout


def test_an_empty_diamond_entry_disables_that_database(tmp_path):
    # the documented way to drop one of the five stock databases.
    proj = _doctor_project(tmp_path,
                           run=dict(_doctor_project(tmp_path).cfg["run"],
                                    diamond=True),
                           db={"pfam_hmm": str(tmp_path / "p.hmm"),
                               "diamond": {"vfdb": "", "merops": ""}})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=None,
                         env=_no_r_env())
    assert "DIAMOND vfdb" not in proc.stdout
    assert "DIAMOND merops" not in proc.stdout


@pytest.mark.xfail(reason="LIVE DOC DRIFT: README says a DIAMOND database "
                          "disabled with an empty string is 'still listed by "
                          "doctor as MANUAL: no path configured'. "
                          "resolve_paths() drops empty db.diamond entries "
                          "before doctor ever sees them, so it is listed "
                          "nowhere and there is no record that the user "
                          "turned it off.",
                   strict=True)
def test_a_disabled_diamond_database_is_still_listed_as_manual(tmp_path):
    proj = _doctor_project(tmp_path,
                           run=dict(_doctor_project(tmp_path).cfg["run"],
                                    diamond=True),
                           db={"pfam_hmm": str(tmp_path / "p.hmm"),
                               "diamond": {"vfdb": ""}})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=None,
                         env=_no_r_env())
    assert "no path configured" in proc.stdout


# --- finding 10 -------------------------------------------------------
def test_the_install_plan_is_a_script_bash_accepts(tmp_path):
    proj = _doctor_project(
        tmp_path, db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")},
        run={"eggnog": True, "pfam": True, "dbcan": True, "diamond": True,
             "cluster": True, "join": True, "topology": True,
             "structure": True, "context": False, "unipept": False,
             "taxonomy": True, "ncbifam": True, "kofam": True,
             "interpro": True, "hhblits": True, "jackhmmer": True,
             "smorf": True, "effectors": False})
    plan = tmp_path / "install.sh"
    run_metaannot("doctor", "--config", proj.config_path, "--install-plan",
                  str(plan), expect=None, env=_no_r_env())
    r = subprocess.run(["bash", "-n", str(plan)], capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stderr
    assert plan.read_text(encoding="utf-8").startswith("#!/usr/bin/env bash")
    assert "set -euo pipefail" in plan.read_text(encoding="utf-8")
    assert os.access(plan, os.X_OK)


def test_the_install_plan_is_written_even_when_nothing_is_missing(tmp_path):
    # symptom: `doctor --install-plan i.sh && bash i.sh` failed with "No such
    # file" on a machine that had everything.
    proj = _doctor_project(tmp_path, run={"eggnog": True, "pfam": False,
                                          "dbcan": False, "diamond": False,
                                          "cluster": False, "join": True,
                                          "topology": False,
                                          "structure": False,
                                          "context": False, "unipept": False,
                                          "taxonomy": False, "ncbifam": False,
                                          "kofam": False, "interpro": False,
                                          "hhblits": False, "jackhmmer": False,
                                          "smorf": False, "effectors": False})
    plan = tmp_path / "install.sh"
    run_metaannot("doctor", "--config", proj.config_path, "--install-plan",
                  str(plan), expect=None, env=_no_r_env())
    assert plan.exists()
    assert "nothing to install" in plan.read_text(encoding="utf-8")
    assert subprocess.run(["bash", "-n", str(plan)]).returncode == 0


# --- finding 11 -------------------------------------------------------
def _fake_installers(tmp_path):
    """curl/gunzip/hmmpress that all exit 0 and produce nothing — the shape
    of a download that wrote a 404 page, or of a mirror that has moved."""
    d = tmp_path / "fakebin"
    d.mkdir(exist_ok=True)
    for name in ("curl", "gunzip", "hmmpress", "tar", "gzip", "diamond",
                 "foldseek", "download_eggnog_data.py"):
        f = d / name
        f.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        f.chmod(f.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return d


@pytest.mark.skipif(
    os.name == "nt",
    reason="doctor --fix is refused on Windows: its install commands are\n           POSIX shell and cmd.exe mis-executes them. See\n           test_fix_is_refused_on_windows_rather_than_half_working.")
def test_fix_verifies_afterwards_instead_of_trusting_exit_codes(tmp_path):
    # symptom: a download that wrote a 404 page exits 0. doctor re-runs its
    # checks after installing and reports UNVERIFIED rather than passing.
    proj = _doctor_project(tmp_path,
                           db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")})
    fake = _fake_installers(tmp_path)
    env = dict(_no_r_env())
    env["PATH"] = str(fake) + os.pathsep + env["PATH"]
    proc = run_metaannot("doctor", "--config", proj.config_path, "--fix",
                         "--yes", expect=1, env=env)
    assert "UNVERIFIED" in proc.stdout
    assert "commands succeeded but the file is still not what the config "\
           "expects" in proc.stdout


@pytest.mark.skipif(
    os.name == "nt",
    reason="doctor --fix is refused on Windows: its install commands are\n           POSIX shell and cmd.exe mis-executes them. See\n           test_fix_is_refused_on_windows_rather_than_half_working.")
def test_fix_refuses_without_confirmation(tmp_path):
    proj = _doctor_project(tmp_path,
                           db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")})
    proc = subprocess.run(
        [sys.executable, METAANNOT_PY, "doctor", "--config",
         proj.config_path, "--fix"],
        input="", capture_output=True, text=True,
        env=dict(os.environ, **_no_r_env()))
    assert "aborted; nothing was downloaded" in proc.stdout
    assert proc.returncode == 1


@pytest.mark.skipif(
    os.name == "nt",
    reason="doctor --fix is refused on Windows: its install commands are\n           POSIX shell and cmd.exe mis-executes them. See\n           test_fix_is_refused_on_windows_rather_than_half_working.")
def test_fix_prints_both_totals_before_doing_anything(tmp_path):
    # the number the user approves is the one that decides whether the volume
    # survives the night, so download size and on-disk size are separate.
    proj = _doctor_project(
        tmp_path, db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")},
        run=dict(_doctor_project(tmp_path).cfg["run"], structure=True))
    proc = subprocess.run(
        [sys.executable, METAANNOT_PY, "doctor", "--config",
         proj.config_path, "--fix"],
        input="n\n", capture_output=True, text=True,
        env=dict(os.environ, **_no_r_env()))
    assert "GB to download" in proc.stdout
    assert "GB on disk" in proc.stdout


def test_fix_is_refused_on_windows_rather_than_half_working(tmp_path):
    """`mkdir -p C:\\db` under cmd.exe makes a directory called -p.

    Every command doctor generates is POSIX shell and every one of them can
    exit 0 having done nothing there, which is exactly the failure the
    post-install verification exists to catch -- so the answer is to refuse
    and say where to run the plan, not to run it and report UNVERIFIED for
    everything.
    """
    if os.name != "nt":
        pytest.skip("the refusal is Windows-only")
    proj = _doctor_project(tmp_path,
                           db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")})
    proc = run_metaannot("doctor", "--config", proj.config_path, "--fix",
                         "--yes", expect=1, env=dict(_no_r_env()))
    assert "doctor --fix cannot run on Windows" in proc.stderr
    assert "--install-plan" in proc.stderr
    assert "wsl bash install.sh" in proc.stderr


def test_the_install_plan_is_still_written_on_windows(tmp_path):
    # only --fix is refused; writing the script the user then runs elsewhere
    # is the whole point of having an alternative to offer.
    proj = _doctor_project(tmp_path,
                           db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")})
    out = str(tmp_path / "install.sh")
    run_metaannot("doctor", "--config", proj.config_path,
                  "--install-plan", out, expect=None, env=dict(_no_r_env()))
    body = io.open(out, encoding="utf-8").read()
    assert body.startswith("#!/usr/bin/env bash")
    assert "Pfam" in body


def test_a_manual_item_is_never_attempted(tmp_path):
    # SignalP 6.0 needs an academic licence; InterProScan is version-specific.
    proj = _doctor_project(
        tmp_path, db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")},
        run=dict(_doctor_project(tmp_path).cfg["run"], topology=True,
                 interpro=True))
    plan = tmp_path / "i.sh"
    proc = run_metaannot("doctor", "--config", proj.config_path,
                         "--install-plan", str(plan), expect=None,
                         env=_no_r_env())
    assert "need a manual download" in proc.stdout
    assert "SignalP 6.0" in proc.stdout
    assert "SignalP" not in plan.read_text(encoding="utf-8")


def test_the_download_guard_rejects_an_html_landing_page(ma, tmp_path):
    # symptom: `curl -fL` follows redirects, so a moved URL answers 200 with a
    # landing page, curl exits 0, and the web page lands on disk under the
    # database's name — passing the existence check for ever.
    page = tmp_path / "Pfam-A.hmm.gz"
    page.write_text("<!DOCTYPE html>\n<html><body>moved</body></html>\n",
                    encoding="utf-8")
    guard = ma._reject_html(str(page))
    r = subprocess.run(["bash", "-c", guard], capture_output=True, text=True)
    assert r.returncode != 0
    assert "is an HTML page, not the database" in r.stderr
    assert not page.exists(), "the web page was left on disk"


def test_the_download_guard_passes_a_real_binary_file(ma, tmp_path):
    real = tmp_path / "Pfam-A.hmm"
    real.write_text("HMMER3/f [3.4]\nNAME  x\n", encoding="utf-8")
    r = subprocess.run(["bash", "-c", ma._reject_html(str(real))],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert real.exists()


def test_an_hmm_library_must_be_pressed_to_count_as_installed(ma, tmp_path):
    # requiring the four hmmpress siblings is also what rejects an HTML page
    # and a truncated download, which existence alone cannot tell apart.
    hmm = tmp_path / "Pfam-A.hmm"
    hmm.write_text("HMMER3/f\n", encoding="utf-8")
    assert not ma._pressed(str(hmm))
    for suf in (".h3f", ".h3i", ".h3m", ".h3p"):
        (tmp_path / ("Pfam-A.hmm" + suf)).write_text("x", encoding="utf-8")
    assert ma._pressed(str(hmm))


def test_a_partial_download_does_not_count_as_installed(ma, tmp_path):
    # symptom: the check used to be `exists(path) or glob(path + "*")`, so a
    # curl killed halfway left the partial .gz and counted for ever.
    target = tmp_path / "uniref50.fasta"
    (tmp_path / "uniref50.fasta.gz").write_text("partial", encoding="utf-8")
    assert not ma._exists(str(target))
    assert not ma._exists(str(tmp_path / "empty"))
    (tmp_path / "empty").write_text("", encoding="utf-8")
    assert not ma._exists(str(tmp_path / "empty"))


# --- finding 12 -------------------------------------------------------
def test_doctor_exits_non_zero_when_something_is_missing(tmp_path):
    proj = _doctor_project(tmp_path,
                           db={"pfam_hmm": str(tmp_path / "absent.hmm")})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "problems found" in proc.stdout


def test_doctor_exits_zero_when_every_requirement_is_satisfied(tmp_path):
    proj = _doctor_project(tmp_path, run={
        k: (k in ("eggnog", "join")) for k in
        ["eggnog", "pfam", "dbcan", "diamond", "cluster", "join", "topology",
         "structure", "context", "unipept", "taxonomy", "ncbifam", "kofam",
         "interpro", "hhblits", "jackhmmer", "smorf", "effectors"]})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=0,
                         env=_no_r_env())
    assert "all checks passed" in proc.stdout


def test_a_manual_item_does_not_count_as_satisfied(tmp_path):
    # symptom: doctor said "all checks passed" with no InterProScan and no
    # SignalP 6, and the run died at that stage hours later.
    proj = _doctor_project(tmp_path, run=dict(
        _doctor_project(tmp_path).cfg["run"], pfam=False, topology=True))
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "MANUAL" in proc.stdout
    assert "all checks passed" not in proc.stdout


def test_a_misspelled_config_key_makes_doctor_fail(tmp_path):
    proj = _doctor_project(tmp_path)
    with open(proj.config_path, "a", encoding="utf-8") as fh:
        fh.write("\nemapper_min_coverag: 0.4\n")
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "== config ==" in proc.stdout
    assert "emapper_min_coverag" in proc.stdout


def test_doctor_survives_a_manifest_it_cannot_read(tmp_path):
    # doctor is the one command that has to survive every other failure.
    proj = _doctor_project(tmp_path)
    with open(proj.cfg["manifest"], "w", encoding="utf-8") as fh:
        fh.write("just one field\n")
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "== resources ==" in proc.stdout, \
        "a bad manifest took the rest of doctor down with it"


def test_doctor_reports_the_manifest_to_column_mapping(tmp_path):
    proj = _doctor_project(tmp_path)
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=None,
                         env=_no_r_env())
    assert "every run maps to a quant column" in proc.stdout


def test_doctor_names_unmatched_manifest_runs(tmp_path):
    proj = _doctor_project(tmp_path)
    import fixtures as _F
    _F.simple_manifest(proj.cfg["manifest"], ["Z_9", "Z_8"])
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "match no column" in proc.stdout


# --- doctor's `== tmt ==` block ---------------------------------------
# Every check below is really one question: does doctor's verdict match what
# read_fragpipe_tmt will do with the same config? doctor is what a TMT user
# reads before committing to a long run.
def _tmt_doctor_project(tmp_path, plexes, no_annotation=(), **tmt):
    """A fragpipe_tmt project whose only enabled stages are eggnog and join,
    so doctor's verdict is about the TMT layout and nothing else.

    `plexes` is {plex: [(channel, sample), ...]} in annotation order.
    """
    import fixtures as F
    proj = _doctor_project(tmp_path, run={
        k: (k in ("eggnog", "join")) for k in
        ["eggnog", "pfam", "dbcan", "diamond", "cluster", "join", "topology",
         "structure", "context", "unipept", "taxonomy", "ncbifam", "kofam",
         "interpro", "hhblits", "jackhmmer", "smorf", "effectors"]})
    root = str(tmp_path / "tmtrun")
    for plex, channels in plexes.items():
        F.write_tmt_plex(root, plex, list(channels),
                         [{"peptide": "PEPTIDEK", "razor": "P_ko_path"}],
                         write_annotation=plex not in no_annotation)
    proj.write_config(quant_table=root, quant_format="fragpipe_tmt",
                      manifest="", tmt=dict(tmt))
    return proj


def test_a_sample_name_in_two_plexes_is_a_doctor_failure_not_a_warning(
        tmp_path):
    # symptom: doctor WARNed that the plexes "will be treated as one sample
    # measured in each" and exited 0, while read_fragpipe_tmt refuses the same
    # config outright -- so the long run died on what doctor had blessed.
    proj = _tmt_doctor_project(
        tmp_path, {"TMT1": [("126", "A1"), ("127N", "A2")],
                   "TMT2": [("126", "A1"), ("127N", "B2")]})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "MISS   1 sample name(s) appear in more than one plex" in proc.stdout
    assert "'A1'" in proc.stdout
    assert "treated as one sample" not in proc.stdout


def test_a_bridge_named_in_every_plex_is_not_a_duplicate_sample_name(
        tmp_path):
    # the one case doctor must NOT flag: under both reference treatments the
    # reader drops the reference before its own collision check, so a pool
    # carrying one name in every plex is the design, not a name clash.
    proj = _tmt_doctor_project(
        tmp_path, {"TMT1": [("126", "A1"), ("131C", "Pool01")],
                   "TMT2": [("126", "B1"), ("131N", "Pool01")]},
        reference_name="Pool*")
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=0,
                         env=_no_r_env())
    assert "more than one plex" not in proc.stdout
    assert "all checks passed" in proc.stdout


def test_doctor_reads_tmt_reference_channel_as_a_glob(tmp_path):
    # symptom: the reader and the README both glob the channel, doctor
    # compared it as a literal, so '131*' -- a valid config -- was reported as
    # matching nothing in any plex and doctor exited 1 on a run that works.
    proj = _tmt_doctor_project(
        tmp_path, {"TMT1": [("126", "A1"), ("131C", "Pool01")],
                   "TMT2": [("126", "B1"), ("131N", "Pool02")]},
        reference_channel="131*")
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=0,
                         env=_no_r_env())
    assert "tmt.reference_channel '131*' resolves in every plex" in proc.stdout


def test_a_reference_pattern_matching_two_channels_of_a_plex_is_reported(
        tmp_path):
    # symptom: doctor said the pattern "resolves in every plex" and exited 0
    # whenever it matched at least once, but the reader requires EXACTLY one
    # reference per plex and dies on "matches 2 of its channels".
    proj = _tmt_doctor_project(
        tmp_path, {"TMT1": [("126", "A1"), ("131N", "Pool01"),
                            ("131C", "Pool02")],
                   "TMT2": [("126", "B1"), ("131N", "Pool03"),
                            ("131C", "Pool04")]},
        reference_name="Pool*")
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "matches more than one channel in 2 plex(es)" in proc.stdout
    assert "resolves in every plex" not in proc.stdout


def test_doctor_does_not_say_a_reference_resolves_when_it_found_no_plexes(
        tmp_path):
    # symptom: with no plex directory at all the "matches nothing" set was
    # empty, so doctor printed OK for a reference it had never looked for --
    # an OK about a run it never opened.
    proj = _tmt_doctor_project(
        tmp_path, {"TMT1": [("126", "A1"), ("131C", "Pool01")]},
        plex_glob="PLEX*", reference_name="Pool*")
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "no plex directory matches" in proc.stdout
    assert "resolves in every plex" not in proc.stdout


def test_the_channel_count_is_not_claimed_for_plexes_doctor_could_not_read(
        tmp_path):
    # symptom: "N channels in every plex" was printed from the plexes whose
    # annotation could be READ, so a run whose second plex had no annotation
    # file was still described as uniform.
    proj = _tmt_doctor_project(
        tmp_path, {"TMT1": [("126", "A1"), ("127N", "A2")],
                   "TMT2": [("126", "B1"), ("127N", "B2")]},
        no_annotation=("TMT2",))
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "channels in every plex" not in proc.stdout
    assert "2 channels in each of the 1 plex(es) read" in proc.stdout


# --- finding 13 -------------------------------------------------------
def test_doctor_output_is_byte_stable_across_two_invocations(tmp_path):
    # symptom: an unstable report cannot be diffed between machines or across
    # a change, which is the only way to see what a config edit did.
    proj = _doctor_project(
        tmp_path, db={"pfam_hmm": str(tmp_path / "db" / "Pfam-A.hmm")},
        run=dict(_doctor_project(tmp_path).cfg["run"], dbcan=True,
                 diamond=True, ncbifam=True, kofam=True, taxonomy=True))
    a = run_metaannot("doctor", "--config", proj.config_path, expect=None,
                      env=_no_r_env())
    b = run_metaannot("doctor", "--config", proj.config_path, expect=None,
                      env=_no_r_env())
    assert a.stdout == b.stdout
    assert a.returncode == b.returncode


def test_doctor_output_is_stable_with_the_r_block_too(tmp_path):
    proj = _doctor_project(tmp_path)
    a = run_metaannot("doctor", "--config", proj.config_path, expect=None)
    b = run_metaannot("doctor", "--config", proj.config_path, expect=None)
    assert a.stdout == b.stdout


def test_doctor_prints_the_per_stage_resource_split(tmp_path):
    proj = _doctor_project(tmp_path)
    proj.write_config(threads=32, ram_gb=128, stage_workers=4)
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=None,
                         env=_no_r_env())
    assert "32 cpu, 128 GB budget, up to 4 stage(s) at once" in proc.stdout
    assert "per stage: 8 cpu, 32 GB" in proc.stdout


def test_doctor_reads_a_suffixed_ram_budget_without_a_traceback(tmp_path):
    # symptom: `ram_gb: 64G` is the form --ram advertises, and int() turned it
    # into an unhandled ValueError halfway through doctor's output.
    proj = _doctor_project(tmp_path)
    proj.write_config(ram_gb="not-a-size")
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1,
                         env=_no_r_env())
    assert "Traceback" not in proc.stderr
    assert "could not read a memory size" in proc.stdout
    assert "== resources ==" in proc.stdout


# ----------------------------------------------------------------------
# RAM detection has to work on every platform the tool runs any part of
# ----------------------------------------------------------------------
def test_ram_detection_answers_on_this_platform(ma):
    """It used to return 0 anywhere without /proc: macOS defines
    _SC_PHYS_PAGES but sysconf returns EINVAL for it, and Darwin has no /proc,
    so both probes fell through and the memory budget silently became 0 — every
    stage then ran with no allocation at all. Windows has neither and runs
    doctor, report and object directly."""
    gb = ma.detect_ram_gb()
    assert gb > 0, "no probe answered on this platform"
    assert gb < 100_000, f"implausible: {gb} GB"


def test_ram_detection_falls_back_rather_than_raising(ma, monkeypatch):
    """Every probe failing must give 0, not a traceback: the caller treats 0 as
    'could not detect' and asks the user for --ram."""
    import builtins
    # raising=False because os.sysconf does not exist on Windows at all, which
    # is half the reason this function needed a third probe.
    monkeypatch.setattr(ma.os, "sysconf",
                        lambda *_: (_ for _ in ()).throw(OSError("nope")),
                        raising=False)
    monkeypatch.setattr(ma.sys, "platform", "sunos5")
    monkeypatch.setattr(ma.os, "name", "posix")
    real_open = builtins.open

    def no_proc(path, *a, **k):
        if str(path).startswith("/proc"):
            raise OSError("no /proc here")
        return real_open(path, *a, **k)

    monkeypatch.setattr(builtins, "open", no_proc)
    assert ma.detect_ram_gb() == 0
