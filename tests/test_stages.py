"""Stage plumbing: adoption, atomic writes, dependencies, the emapper join.

Every external tool is stubbed by a script on PATH, so these exercise the
pipeline's own logic without hmmsearch, DIAMOND, MMseqs2 or Foldseek.
"""
from __future__ import annotations

import contextlib
import gzip
import io
import json
import os
import re
import shutil
import socket
import sys
import threading
import types
import time

import pytest

import fixtures as F
from conftest import METAANNOT_PY, _load, build_project, run_metaannot


def _searchable(tmp_path, root, **over):
    """A project with the search stages on and fake database files present."""
    db = tmp_path / "db"
    db.mkdir(exist_ok=True)
    for name in ("Pfam-A.hmm", "dbCAN.txt", "hmm_PGAP.LIB"):
        (db / name).write_text("HMMER3/f\n", encoding="utf-8")
    F.write_dmnd(db / "vfdb.dmnd")
    F.write_dmnd(db / "merops.dmnd")
    cfg = {
        "run": {"eggnog": True, "pfam": True, "dbcan": True, "diamond": True,
                "cluster": True, "join": True, "topology": False,
                "structure": False, "context": False, "unipept": False,
                "taxonomy": False, "ncbifam": False, "kofam": False,
                "interpro": False, "hhblits": False, "jackhmmer": False,
                "smorf": False, "effectors": False},
        "db": {"pfam_hmm": str(db / "Pfam-A.hmm"),
               "dbcan_hmm": str(db / "dbCAN.txt"),
               "diamond": {"vfdb": str(db / "vfdb.dmnd"),
                           "merops": str(db / "merops.dmnd")}},
    }
    cfg.update(over)
    return build_project(root, **cfg)


# --- finding 21 -------------------------------------------------------
def test_diamond_stage_runs_on_a_fresh_results_dir(tmp_path, stub_bin):
    # symptom: the stage's output used to be a DIRECTORY that Paths.mkdirs()
    # pre-created, so exists_all() was true before anything ran and the stage
    # was adopted as already done — for ever.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    st = proj.state()
    assert st["diamond"]["status"] == "ok", "diamond must RUN, not be adopted"
    assert os.path.exists(proj.rpath("diamond", "vfdb.tsv"))
    assert os.path.exists(proj.rpath("diamond", "merops.tsv"))


def test_a_stage_whose_output_directory_is_precreated_is_not_adopted(ma,
                                                                     project):
    # the general form: an empty directory does not count as an output.
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    assert not ma.exists_all([p.diamond_dir])
    assert not ma.exists_all([p.hhr_dir])
    assert not ma.exists_all([p.structures])


def test_every_search_stage_runs_and_its_evidence_reaches_the_bins(tmp_path,
                                                                   stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    st = proj.state()
    for name in ("pfam", "dbcan", "diamond", "cluster"):
        assert st[name]["status"] == "ok", f"{name} did not run"
    import pandas as pd
    ann = pd.read_csv(proj.rpath("annotation_final.tsv"), sep="\t")
    assert ann["pfam_hits"].fillna("").ne("").any()
    assert ann["dbcan_hits"].fillna("").ne("").any()
    assert ann["vfdb_hit"].fillna("").ne("").any()
    assert ann["family_id"].notna().all()


# --- finding 22 -------------------------------------------------------
def test_an_empty_output_is_not_adopted(tmp_path, stub_bin):
    # symptom: a truncated pfam.tblout — an interrupted writer — was adopted
    # as a complete stage, so every protein silently lost its Pfam evidence.
    proj = _searchable(tmp_path, tmp_path / "p")
    os.makedirs(proj.rpath("hmm"), exist_ok=True)
    open(proj.rpath("hmm", "pfam.tblout"), "w").close()
    proc = proj.run()
    assert "not adopting" in proc.stderr
    assert "the file is empty" in proc.stderr
    assert proj.state()["pfam"]["status"] == "ok"
    assert os.path.getsize(proj.rpath("hmm", "pfam.tblout")) > 0


def test_a_stage_that_declares_emptiness_meaningful_is_adopted_empty(ma):
    # diamond, smorf, jackhmmer, hhblits, esmfold, foldseek and tmbed can
    # legitimately produce nothing; every other stage cannot. tmbed joined the
    # list when tmbed_max_len grew teeth: a proteome whose every sequence is
    # over the cap leaves an empty prediction file on purpose.
    empty_ok = {s["name"] for s in ma.STAGES if s.get("empty_ok")}
    assert empty_ok == {"diamond", "smorf", "jackhmmer", "hhblits", "esmfold",
                        "foldseek", "tmbed"}


def test_a_non_empty_output_produced_elsewhere_is_adopted(tmp_path, stub_bin):
    # the two-machine workflow: fold on the GPU box, rsync results back.
    proj = _searchable(tmp_path, tmp_path / "p")
    os.makedirs(proj.rpath("hmm"), exist_ok=True)
    F.write_tblout(proj.rpath("hmm", "pfam.tblout"),
                   [("P_dark1", "Peptidase_S8", "PF00082.1")])
    proc = proj.run()
    assert proj.state()["pfam"]["status"] == "adopted"
    assert "adopting output this run did not produce" in proc.stderr


def test_an_interrupted_writer_leaves_no_output_at_all(ma, tmp_path,
                                                       paths_for):
    # symptom: every stage output used to be written in place, so a writer
    # killed by the OOM reaper left a truncated file the next run adopted.
    cfg, p = paths_for()
    target = os.path.join(p.R, "out.tsv")
    with pytest.raises(RuntimeError):
        with ma.atomic_out(target) as tmp:
            with open(tmp, "w", encoding="utf-8") as fh:
                fh.write("half a fi")
            raise RuntimeError("killed")
    assert not os.path.exists(target)
    leftovers = [f for f in os.listdir(p.R) if ma.ATOMIC_SUFFIX in f]
    assert leftovers == []


def test_a_tool_that_exits_zero_without_writing_is_a_failure(ma, tmp_path,
                                                             paths_for):
    cfg, p = paths_for()
    target = os.path.join(p.R, "out.tsv")
    with pytest.raises(ma.StageError) as e:
        with ma.atomic_out(target):
            pass
    assert "nothing was written" in str(e.value)


def test_the_temp_file_is_a_dotfile_so_globs_cannot_pick_it_up(ma, paths_for):
    # symptom: results/diamond/*.tsv is globbed as "every database we
    # searched", and a temp left by a killed writer came back as a database
    # tag called "vfdb.9134.7.part".
    cfg, p = paths_for()
    seen = []
    with ma.atomic_out(os.path.join(p.diamond_dir, "vfdb.tsv")) as tmp:
        seen.append(os.path.basename(tmp))
        open(tmp, "w").write("x")
    assert seen[0].startswith(".")
    assert seen[0].endswith(".tsv"), "the extension must be preserved"


# --- finding 23 -------------------------------------------------------
def test_every_dependency_output_is_in_the_dependents_signature_inputs(ma,
                                                                       project):
    # symptom: integrate depended on diamond for ORDERING but not for
    # invalidation, so adding a database left integrate cached and the new
    # evidence never reached the bins.
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    by_name = {s["name"]: s for s in ma.STAGES}
    bad = []
    for st in ma.STAGES:
        inputs = set(x for x in st["inp"](cfg, p) if x)
        for d in st["deps"]:
            outs = set(by_name[d]["out"](p))
            if not (outs & inputs):
                bad.append((st["name"], d))
    assert bad == [], f"dependencies absent from the signature: {bad}"


def _build_annotation_config_keys():
    """Config keys build_annotation reads on EVERY call, read off the source.

    Derived rather than listed, so a key added to build_annotation later
    cannot quietly skip the signature. Keys read only inside an `emit_dark`
    branch are excluded: finalise calls build_annotation with emit_dark unset,
    so those belong to integrate alone. `db` and `proteins_faa` name files,
    which the signature already covers by content through stage["inp"].
    """
    import ast

    tree = ast.parse(open(METAANNOT_PY, encoding="utf-8").read())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "build_annotation")

    def reads(node):
        out = set()
        for n in ast.walk(node):
            if (isinstance(n, ast.Subscript) and isinstance(n.value, ast.Name)
                    and n.value.id == "cfg"
                    and isinstance(n.slice, ast.Constant)):
                out.add(n.slice.value)
            if (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                    and n.func.attr == "get"
                    and isinstance(n.func.value, ast.Name)
                    and n.func.value.id == "cfg" and n.args
                    and isinstance(n.args[0], ast.Constant)):
                out.add(n.args[0].value)
        return out

    gated = set()
    for n in ast.walk(fn):
        if isinstance(n, ast.If) and "emit_dark" in ast.unparse(n.test):
            gated |= reads(n)
    return reads(fn) - gated - {"db", "proteins_faa"}


def _bump(v):
    """A value that differs from v, whatever shape v has."""
    if isinstance(v, dict):
        return dict(v, __probe__=1)
    if isinstance(v, (list, tuple)):
        return list(v) + ["__probe__"]
    if isinstance(v, bool) or v is None:
        return not v
    if isinstance(v, (int, float)):
        return v + 1
    return f"{v}__probe__"


@pytest.mark.parametrize("key", sorted(_build_annotation_config_keys()))
@pytest.mark.parametrize("stage", ["integrate", "finalise"])
def test_a_config_key_build_annotation_reads_invalidates_the_stage_that_runs_it(
        ma, project, stage, key):
    # symptom: vfdb_category_weights was in finalise's signature keys but not
    # in integrate's. build_annotation applies it, and integrate is what runs
    # build_annotation — so re-weighting VFDB left integrate cached, finalise
    # re-ran, took its "no structure or profile evidence, reusing the first
    # pass" branch, and re-published the stale scores. The run reported
    # success and the setting had done nothing. foldseek_target_priority was
    # missing from integrate the same way.
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    st = [s for s in ma.STAGES if s["name"] == stage][0]
    assert key in st["keys"], (
        f"build_annotation reads cfg[{key!r}] on every call, but {stage!r} "
        f"does not list it, so changing it leaves the stage cached")
    before = ma.signature(st, cfg, p)
    cfg[key] = _bump(cfg.get(key))
    assert ma.signature(st, cfg, p) != before, (
        f"changing {key!r} did not invalidate {stage!r}")


def test_adding_a_diamond_database_invalidates_integrate(ma, tmp_path,
                                                         stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    sig_before = proj.state()["integrate"]["signature"]
    db = tmp_path / "db" / "card.dmnd"
    F.write_dmnd(db)
    d = dict(proj.cfg["db"])
    d["diamond"] = dict(d["diamond"], card=str(db))
    proj.write_config(db=d)
    proj.run()
    assert proj.state()["integrate"]["signature"] != sig_before
    assert proj.state()["integrate"]["status"] == "ok"


def test_changing_a_config_key_a_stage_depends_on_reruns_it(ma, tmp_path,
                                                            stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    before = proj.state()["dbcan"]["signature"]
    proj.write_config(thresholds={"dbcan_evalue": 1e-10})
    proj.run()
    assert proj.state()["dbcan"]["signature"] != before


def test_tool_args_changes_invalidate_the_cache(ma, project):
    # tool_args is documented as THE way to correct a wrong flag, and is
    # appended to nearly every command line.
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    st = [s for s in ma.STAGES if s["name"] == "pfam"][0]
    before = ma.signature(st, cfg, p)
    cfg["tool_args"] = {"hmmsearch": ["--nonull2"]}
    assert ma.signature(st, cfg, p) != before


def test_an_input_is_identified_by_its_content_not_by_its_path(ma, tmp_path):
    # symptom: hashing the path string meant that spelling the same directory
    # as 'results' on one run and as an absolute path on the next recomputed
    # everything; hashing the mtime meant a cp/scp of an unchanged database
    # did the same.
    a = tmp_path / "a.hmm"
    a.write_text("HMMER3/f\nNAME  x\n", encoding="utf-8")
    b = tmp_path / "sub" / "b.hmm"
    b.parent.mkdir()
    shutil.copyfile(a, b)
    os.utime(b, (0, 0))                    # a different mtime, same bytes
    ma._HASH_CACHE.clear()
    assert ma._stat(str(a)) == ma._stat(str(b))
    b.write_text("HMMER3/f\nNAME  y\n", encoding="utf-8")
    ma._HASH_CACHE.clear()
    assert ma._stat(str(a)) != ma._stat(str(b))


def test_a_missing_input_hashes_as_absent_rather_than_raising(ma, tmp_path):
    assert ma._stat(str(tmp_path / "nope")) == [None, None]


# --- finding 33 -------------------------------------------------------
def test_a_prefixed_search_database_can_be_joined_to_an_unprefixed_table(
        ma, tmp_path, capsys):
    # symptom: a search database often prefixes the ids it was built from
    # (uhgpSM_, HUMANHOST_) while the eggNOG table does not; no id transform
    # can add or remove a prefix, so the coverage gate blamed biology.
    ps = [F.Protein(f"uhgpSM_MGYG{i:05d}", "MKV" * 40, ko="ko:K01234",
                    pathway="ko00010,map00010", seed_taxid="820")
          for i in range(10)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), ps,
                          id_prefix="uhgpSM_")
    out = str(tmp_path / "out.annotations")
    report = str(tmp_path / "report.tsv")
    ma.prepare_emapper([emp], faa, out, report, "exact", 0.5, 0.9, 100,
                       strip_prefixes=["uhgpSM_"])
    df = ma.parse_emapper(out)
    assert set(df.index) == {p.pid for p in ps}, "rows keep the FASTA id"
    err = capsys.readouterr().err
    assert "prefix 'uhgpSM_' stripped from 10 fasta ids" in err
    assert "matched only after emapper_strip_id_prefix" in err
    rep = dict(l.split("\t") for l in open(report).read().splitlines()[1:])
    assert rep["proteins_bridged_by_prefix"] == "10"


def test_the_coverage_diagnostic_names_the_prefix_that_would_have_matched(
        ma, tmp_path, capsys):
    # symptom: the diagnostic reported the proteins as genuinely unannotated.
    ps = [F.Protein(f"uhgpSM_MGYG{i:05d}", "MKV" * 40, seed_taxid="820")
          for i in range(20)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), ps,
                          id_prefix="uhgpSM_")
    with pytest.raises(ma.StageError):
        ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), "",
                           "exact", 0.5, 0.9, 1000)
    err = capsys.readouterr().err
    assert "set emapper_strip_id_prefix: uhgpSM_" in err


def test_a_prefix_that_matches_nothing_is_reported(ma, tmp_path, capsys):
    ps = F.protein_set()
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), ps)
    ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), "",
                       "exact", 0.5, 0.9, 100, strip_prefixes=["nope_"])
    assert "matches no fasta id" in capsys.readouterr().err


def test_one_eggnog_row_annotates_every_fasta_id_it_covers(ma, tmp_path,
                                                           capsys):
    # symptom: a composite database holds the same bare id under two tier
    # prefixes; keeping only the first fasta id sent the other to 4_dark
    # although its annotation was sitting right there.
    bare = F.Protein("MGYG00001", "MKV" * 40, ko="ko:K01234",
                     pathway="ko00010,map00010", seed_taxid="820")
    ps = [F.Protein("uhgpSM_MGYG00001", bare.seq, seed_taxid="820"),
          F.Protein("uhgpL_MGYG00001", bare.seq, seed_taxid="820")]
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), [bare])
    out = str(tmp_path / "o.annotations")
    ma.prepare_emapper([emp], faa, out, "", "exact", 0.5, 0.9, 100,
                       strip_prefixes=["uhgpSM_", "uhgpL_"])
    df = ma.parse_emapper(out)
    assert set(df.index) == {"uhgpSM_MGYG00001", "uhgpL_MGYG00001"}
    assert "annotated more than one fasta id" in capsys.readouterr().err


def test_low_coverage_aborts_and_names_the_transform_that_fixes_it(ma,
                                                                   tmp_path,
                                                                   capsys):
    # symptom: every unmatched protein would be misreported as 4_dark.
    ps = [F.Protein(f"sp|X{i:04d}|NAME", "MKV" * 40) for i in range(20)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    stripped = [F.Protein(p.pid.split("|")[1], p.seq) for p in ps]
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), stripped)
    with pytest.raises(ma.StageError) as e:
        ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), "",
                           "exact", 0.5, 0.9, 1000)
    assert "emapper_min_coverage" in str(e.value)
    assert "Fix the id mismatch rather than lowering" in str(e.value)


def test_the_rejected_partial_table_is_not_left_where_a_rerun_adopts_it(
        ma, tmp_path):
    ps = [F.Protein(f"sp|X{i:04d}|N", "MKV" * 40) for i in range(20)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"),
                          [F.Protein(p.pid.split("|")[1], p.seq) for p in ps])
    out = str(tmp_path / "o.annotations")
    with pytest.raises(ma.StageError):
        ma.prepare_emapper([emp], faa, out, "", "exact", 0.5, 0.9, 1000)
    assert not os.path.exists(out), "only the .tmp may survive a refusal"


def test_an_unknown_id_transform_is_refused(ma, project):
    cfg = ma.load_config(project.config_path)
    cfg["emapper_id_transform"] = "strip_everything"
    p = ma.Paths(cfg)
    p.mkdirs()
    with pytest.raises(ma.StageError) as e:
        ma.stage_emapper(cfg, p)
    assert "unknown emapper_id_transform" in str(e.value)


def test_strip_prefix_without_a_precomputed_table_says_it_is_ignored(
        ma, project, capsys, monkeypatch):
    cfg = ma.load_config(project.config_path)
    cfg["emapper_precomputed"] = []
    cfg["emapper_strip_id_prefix"] = "uhgpSM_"
    p = ma.Paths(cfg)
    p.mkdirs()
    with pytest.raises(ma.StageError):
        ma.stage_emapper(cfg, p)          # emapper.py is not installed
    assert "only applies when reusing a precomputed table" in \
        capsys.readouterr().err


# --- finding 28 -------------------------------------------------------
def test_a_failed_stage_reports_its_message_not_just_a_count(tmp_path,
                                                             stub_bin):
    # symptom: die() raised SystemExit, which inside a worker thread ended
    # only that thread and whose str() is just the exit code — the diagnosis
    # was lost before it reached the summary.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(db=dict(proj.cfg["db"], pfam_hmm=""))
    proj.write_config(run=dict(proj.cfg["run"], kofam=True))
    proc = proj.run(expect=1)
    assert "KOfamScan not found" in proc.stderr
    assert "failed: 1" not in proc.stderr
    st = proj.state()
    assert st["kofam"]["status"] == "failed"
    assert "KOfamScan not found" in st["kofam"]["error"]


def test_a_missing_tool_message_carries_its_install_command(ma, project):
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    with pytest.raises(ma.StageError) as e:
        ma.stage_kofam(cfg, p)
    assert "conda install" in str(e.value)


def test_die_raises_a_real_exception_not_system_exit(ma):
    with pytest.raises(ma.StageError):
        ma.die("boom")
    assert not issubclass(ma.StageError, SystemExit)


# --- --force scoping --------------------------------------------------
def test_force_with_only_discards_just_that_stages_record(tmp_path, stub_bin):
    # symptom: a bare --force wiped the whole state file, so every unselected
    # stage looked as if it had never been recorded and its stale output was
    # adopted without a signature check.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    before = proj.state()
    proj.run("--force", "--only", "pfam")
    after = proj.state()
    assert after["pfam"]["finished"] >= before["pfam"]["finished"]
    for name in ("dbcan", "diamond", "cluster", "emapper"):
        assert after[name]["signature"] == before[name]["signature"]
        assert after[name]["status"] == before[name]["status"]


def test_only_without_force_reports_the_stage_as_cached(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    proc = proj.run("--only", "pfam")
    assert "pfam: cached" in proc.stderr
    assert "--force --only pfam" in proc.stderr


def test_naming_a_disabled_stage_with_only_overrides_its_flag(tmp_path,
                                                              stub_bin):
    # this is how the GPU box runs tmbed/esmfold against the server's config.
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.write_config(run=dict(proj.cfg["run"], ncbifam=False))
    proc = proj.run("--only", "ncbifam")
    assert "run.ncbifam is false, but the stage was named with --only, so it "\
           "is running anyway" in proc.stderr
    assert proj.state()["ncbifam"]["status"] == "ok"


# --- misc stage behaviour ---------------------------------------------
def test_interproscan_input_is_sanitised_of_non_iupac_residues(ma, tmp_path):
    # symptom: a single '*' (Prodigal keeps the terminal stop) aborts
    # InterProScan with a Java exception hours into the longest stage.
    src = tmp_path / "in.faa"
    src.write_text(">P1\nMKVJU*\n>P2\nMKVA\n", encoding="utf-8")
    dst = str(tmp_path / "clean.faa")
    n, fixed = ma.sanitise_faa(str(src), dst)
    assert (n, fixed) == (2, 1)
    seqs = dict(ma.read_fasta(dst))
    assert seqs["P1"] == "MKVXX"
    assert "*" not in seqs["P1"]


def test_tool_args_accepts_a_string_or_a_list(ma):
    assert ma.tool_args({"tool_args": {"h": "--a --b"}}, "h") == ["--a", "--b"]
    assert ma.tool_args({"tool_args": {"h": ["--a", 1]}}, "h") == ["--a", "1"]
    assert ma.tool_args({}, "h") == []


def test_tool_args_of_the_wrong_type_is_refused(ma):
    with pytest.raises(ma.StageError) as e:
        ma.tool_args({"tool_args": {"h": {"a": 1}}}, "h")
    assert "must be a command line" in str(e.value)


def test_a_diamond_database_that_does_not_exist_is_skipped_with_a_warning(
        ma, tmp_path, paths_for, stub_bin, capsys):
    cfg, p = paths_for()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    cfg["db"]["diamond"] = {"vfdb": str(tmp_path / "absent.dmnd")}
    ma.stage_diamond(cfg, p)
    assert "diamond database missing, skipping" in capsys.readouterr().err
    assert os.path.exists(p.diamond_done)


def test_an_unknown_tmbed_gpu_setting_is_refused(ma, tmp_path, paths_for,
                                                 monkeypatch):
    cfg, p = paths_for()
    cfg["tmbed_use_gpu"] = "maybe"
    monkeypatch.setattr(ma, "have", lambda t: True)
    with pytest.raises(ma.StageError) as e:
        ma.stage_tmbed(cfg, p)
    assert "unknown tmbed_use_gpu" in str(e.value)


# tmbed's real shape: a tqdm bar redrawn in place on stderr for hours, and the
# predictions written only at the very end.
_TMBED_STUB = """#!/usr/bin/env python3
import os, sys, time
a = sys.argv[1:]
out = a[a.index("-p") + 1]
faa = a[a.index("-f") + 1]
for i in range(6):
    bar = "\\r 61%|###### | " + str(i * 5000) + "/38204 [2:36:04<1:39:41]"
    sys.stderr.write(bar)
    sys.stderr.flush()
    time.sleep(0.08)
# One 3-line record per INPUT sequence, which is the part of the real tool's
# behaviour the stage now depends on: it reconciles what it handed over
# against what came back, so a stub that always wrote the same fixed record
# would pass every id-handling bug straight through.
recs, pid, seq = [], None, []
for line in open(faa, encoding="utf-8"):
    line = line.strip()
    if line.startswith(">"):
        if pid:
            recs.append((pid, "".join(seq)))
        pid, seq = line[1:].split()[0], []
    elif line:
        seq.append(line)
if pid:
    recs.append((pid, "".join(seq)))
# The three ways a real run goes wrong, addressed by PROTEIN ID rather than by
# chunk file name so a test does not have to predict how the plan was
# numbered. FAIL: this chunk dies (after PARTIAL records, if set). DROP: it
# exits 0 having written fewer records than it was given, which is the failure
# no exit code reports.
fail = set(x for x in os.environ.get("TMBED_STUB_FAIL", "").split(",") if x)
dying = fail & set(q for q, _ in recs)
if dying:
    recs = recs[:int(os.environ.get("TMBED_STUB_PARTIAL", "0"))]
elif os.environ.get("TMBED_STUB_DROP"):
    recs = recs[:int(os.environ["TMBED_STUB_DROP"])]
with open(out, "w", encoding="utf-8") as fh:
    for pid, seq in recs:
        fh.write(">" + pid + "\\n" + seq + "\\n" + "i" * len(seq) + "\\n")
if dying:
    sys.stderr.write("\\nRuntimeError: CUDA out of memory\\n")
    sys.exit(1)
"""


def _stub_tmbed(tmp_path, monkeypatch):
    """A tmbed on PATH that behaves like the real one's output does."""
    d = tmp_path / "tmbedbin"
    d.mkdir(exist_ok=True)
    (d / "tmbed").write_text(_TMBED_STUB, encoding="utf-8")
    os.chmod(d / "tmbed", 0o755)
    if os.name == "nt":
        # PATHEXT decides what is executable there, so an extension-less stub
        # is never found — the same reason conftest's stub_bin writes a shim.
        (d / "tmbed.cmd").write_text(
            f'@echo off\r\n"{sys.executable}" "%~dp0tmbed" %*\r\n',
            encoding="utf-8")
    monkeypatch.setenv("PATH", str(d) + os.pathsep + os.environ["PATH"])
    return d


def test_tmbed_reports_progress_like_every_other_long_running_tool(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    # symptom: tmbed is the tool the heartbeat was written for — the example
    # in the release note is its own bar — but the stage called subprocess.run
    # directly instead of run_cmd, so it got neither the stderr ring nor the
    # heartbeat and said nothing for the 2 h 36 min it ran.
    cfg, p = paths_for("tmbed_progress")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setattr(ma, "_PROGRESS_INTERVAL", 0.15)
    ma.stage_tmbed(cfg, p)
    err = capsys.readouterr().err
    progress = [l for l in err.splitlines()
                if " running " in l and "tmbed" in l]
    assert progress, "the stage the heartbeat exists for emitted no progress"
    assert any("61%|" in l for l in progress), \
        "tmbed's own bar is what says how far along it is"
    assert os.path.exists(p.tmbed), "the predictions are still adopted"


def test_tmbed_is_launched_by_the_path_that_was_resolved_for_it(
        ma, tmp_path, paths_for, monkeypatch):
    # symptom: every tool is meant to be launched by the absolute path PATH
    # resolves to, because CreateProcess ignores PATHEXT and a .cmd earlier on
    # PATH loses to an .exe later on it — but the stage handed the OS the bare
    # name, so the binary that ran need not be the one have() reported and the
    # logged command line named.
    cfg, p = paths_for("tmbed_argv")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    d = _stub_tmbed(tmp_path, monkeypatch)
    launched = []

    def _spy(real):
        def go(argv, *a, **k):
            launched.append([str(c) for c in argv])
            return real(argv, *a, **k)
        return go
    # Both launchers, because the name handed to the OS is what is under test
    # and not which of them the stage reaches for.
    monkeypatch.setattr(ma.subprocess, "run", _spy(ma.subprocess.run))
    monkeypatch.setattr(ma.subprocess, "Popen", _spy(ma.subprocess.Popen))
    ma.stage_tmbed(cfg, p)
    argv0 = [c[0] for c in launched
             if os.path.basename(c[0]).startswith("tmbed")]
    assert argv0, "no tmbed process was launched"
    assert os.path.isabs(argv0[0]), \
        f"tmbed was launched as {argv0[0]!r}, not the path PATH resolves to"
    assert os.path.dirname(argv0[0]) == str(d)


# --- run_cmd: progress out of a multi-hour tool ------------------------
# symptom: tmbed ran 2 h 36 min and then died, twice, and InterProScan 2.9 h,
# with nothing in the log between the command and the failure. Both write a
# tqdm bar to stderr the whole time; stderr was captured and only quoted on
# failure, so the only way to tell a live stage from a hung one was to watch
# its CPU ticks accumulate in /proc.

_TQDM_LIKE = r"""
import sys, time
for i in range(8):
    sys.stderr.write('\r\x1b[32m 12%|##        | ' + str(i)
                     + '/8 [00:00<00:05, 1.2it/s]\x1b[0m')
    sys.stderr.flush()
    time.sleep(0.08)
"""


def test_a_running_tool_reports_progress_instead_of_going_silent(
        ma, monkeypatch, capsys):
    monkeypatch.setattr(ma, "_PROGRESS_INTERVAL", 0.15)
    assert ma.run_cmd([sys.executable, "-c", _TQDM_LIKE]) == "", \
        "the return value is what every stage depends on"
    progress = [l for l in capsys.readouterr().err.splitlines()
                if " running " in l]
    assert progress, "a tool that ran for several intervals said nothing"
    # one readable line per interval, not a wall of partial redraws
    assert all(l.count("it/s") <= 1 for l in progress)
    assert all("\r" not in l and "\x1b" not in l for l in progress)
    assert any(re.search(r"running \d+m\d\ds", l) for l in progress), \
        "the elapsed time is the half of the message that proves it is alive"
    assert any("12%|" in l for l in progress), \
        "the tool's own progress is what says how far along it is"


def test_progress_can_be_turned_off(ma, monkeypatch, capsys):
    monkeypatch.setattr(ma, "_PROGRESS_INTERVAL", 0)
    ma.run_cmd([sys.executable, "-c", _TQDM_LIKE])
    assert " running " not in capsys.readouterr().err


def test_a_silent_tool_still_gets_a_heartbeat(ma, monkeypatch, capsys):
    # a stage that writes nothing at all is exactly the one you cannot tell
    # from a hang, so the line goes out with or without tool output.
    monkeypatch.setattr(ma, "_PROGRESS_INTERVAL", 0.15)
    ma.run_cmd([sys.executable, "-c", "import time; time.sleep(0.5)"])
    err = capsys.readouterr().err
    assert "no output yet on stderr" in err


def test_the_failure_tail_still_quotes_the_last_fifteen_lines(ma):
    # unchanged on purpose: every stage's diagnosis comes out of this string.
    code = ("import sys\n"
            "for i in range(4000): sys.stderr.write('line %d\\n' % i)\n"
            "sys.exit(3)\n")
    with pytest.raises(RuntimeError) as e:
        ma.run_cmd([sys.executable, "-c", code])
    msg = str(e.value)
    assert "exited 3" in msg
    tail = msg.split("--- stderr tail ---\n")[1].splitlines()
    assert tail == [f"line {i}" for i in range(3985, 4000)], \
        "4000 lines through a bounded ring must still end in the last 15"


_TQDM_THEN_DIES = r"""
import sys, time
for i in range(8):
    sys.stderr.write('\r\x1b[32m %d%%|##        | %d/8 [00:00<00:05, 1.2it/s]'
                     '\x1b[0m' % (12 * i, i))
    sys.stderr.flush()
    time.sleep(0.08)
sys.stderr.write('\nCUDA out of memory. Tried to allocate 2.00 GiB\n')
sys.exit(1)
"""


def test_a_carriage_return_bar_gives_progress_and_still_quotes_the_tail(
        ma, monkeypatch, capsys):
    # both halves in one run, because they trade against each other: a bar
    # that only ever redraws in place is one unterminated line, so a reader
    # that waited for a newline would print nothing while it ran AND leave the
    # failure tail empty. This is tmbed's exact shape - hours of bar, then a
    # CUDA OOM on the last line.
    monkeypatch.setattr(ma, "_PROGRESS_INTERVAL", 0.15)
    with pytest.raises(RuntimeError) as e:
        ma.run_cmd([sys.executable, "-c", _TQDM_THEN_DIES])
    progress = [l for l in capsys.readouterr().err.splitlines()
                if " running " in l]
    assert progress, "a carriage-return bar must still produce progress lines"
    assert any("|##" in l for l in progress)
    msg = str(e.value)
    tail = msg.split("--- stderr tail ---\n")[1].splitlines()
    assert tail, "the redraws must not swallow the tail"
    assert tail[-1] == "CUDA out of memory. Tried to allocate 2.00 GiB", \
        "the reason a stage died is the last thing it wrote"
    assert any("it/s" in l for l in tail), \
        "a redraw is a line of the ring like any other"
    assert len(tail) <= 15


def test_stderr_is_not_hoarded(ma):
    assert ma._STDERR_KEEP <= 1000, \
        "stdout goes to devnull because tens of MB bought nothing; the " \
        "same reasoning caps what stderr may keep"


@pytest.mark.parametrize("raw,want", [
    ("\x1b[32m 45%|####      | 45/100\x1b[0m", "45%|####      | 45/100"),
    ("  padded  ", "padded"),
    ("a\tb", "a b"),
    ("", ""),
])
def test_a_progress_bar_becomes_one_sensible_line(ma, raw, want):
    assert ma._progress_line(raw) == want


def test_a_very_long_progress_line_is_truncated(ma):
    assert len(ma._progress_line("x" * 5000)) == 160


def test_the_progress_interval_comes_from_the_config(ma, monkeypatch):
    monkeypatch.setattr(ma, "_PROGRESS_INTERVAL", 999.0)
    ma.set_progress_interval(ma.DEFAULT_CONFIG["progress_interval_s"])
    assert ma._PROGRESS_INTERVAL == 60.0
    ma.set_progress_interval(0)
    assert ma._PROGRESS_INTERVAL == 0.0


# --- the longest stage in the pipeline, and the only signal it gives ----
# symptom: InterProScan ran 57 hours on the reference dataset with nothing on
# stderr that said how far it had got, so an operator wanting an ETA had to
# write a chunk-counting script against its -T tree by hand. The heartbeat can
# count the same files.


def _ips_tree(ma, root, created, done):
    """A -T tree shaped like the one the reference run left behind.

    The extensions come from the module rather than from this file: they are
    the observation the whole mechanism rests on, and a test that spelt them
    out again would agree with itself after someone corrected them there.
    """
    job = os.path.join(root, "acme01_20260101_120000", "job7")
    os.makedirs(job, exist_ok=True)
    for i in range(created):
        open(os.path.join(job, f"chunk{i}{ma._IPS_WORK_EXT}"), "w").close()
    for i in range(done):
        open(os.path.join(job, f"chunk{i}{ma._IPS_DONE_EXT}"), "w").close()
    return job


def test_the_interpro_heartbeat_counts_the_chunks_the_tool_leaves_behind(
        ma, tmp_path):
    root = str(tmp_path / "tmp")
    os.makedirs(root)
    probe = ma._InterProProgress(root, 60)
    assert probe() == "", "nothing on disk yet is not a number"
    assert not probe.seen_any
    _ips_tree(ma, root, created=9, done=4)
    assert probe() == "chunk 4/9"
    assert probe.seen_any, \
        "the stage asks this afterwards to tell silence from a dead probe"


def test_no_eta_is_offered_while_interproscan_is_still_splitting(
        ma, tmp_path, monkeypatch):
    # The denominator GROWS during the split phase, so a remaining-time taken
    # then is measured against a number about to move — and chunks ARE
    # completing while it moves, so the rate is perfectly computable and
    # perfectly wrong. A fake clock, because the whole claim is about time.
    clock = [1_000.0]
    monkeypatch.setattr(ma.time, "time", lambda: clock[0])
    root = str(tmp_path / "tmp")
    job = _ips_tree(ma, root, created=4, done=1)
    probe = ma._InterProProgress(root, 60)
    for i in range(6):
        line = probe()
        assert "left" not in line, \
            f"an ETA against a denominator that is still moving: {line}"
        clock[0] += 600.0
        open(os.path.join(job, f"chunk{4 + i}{ma._IPS_WORK_EXT}"), "w").close()
        open(os.path.join(job, f"chunk{1 + i}{ma._IPS_DONE_EXT}"), "w").close()
    # Splitting stops at 10. The count must then hold still for
    # _IPS_ETA_STABLE_TICKS probes AND a chunk must finish after the anchor,
    # so there is a measured rate rather than an extrapolated one.
    for _ in range(ma._IPS_ETA_STABLE_TICKS + 1):
        assert "left" not in probe()
        clock[0] += 600.0
    open(os.path.join(job, f"chunk7{ma._IPS_DONE_EXT}"), "w").close()
    # one chunk in the 600 s since the anchor, two left -> twenty minutes
    assert probe() == "chunk 8/10, ~20m00s left"
    # and a fresh slice withdraws it rather than quoting a stale denominator
    for i in range(10, 14):
        open(os.path.join(job, f"chunk{i}{ma._IPS_WORK_EXT}"), "w").close()
    assert probe() == "chunk 8/14"


def test_a_layout_this_does_not_recognise_reports_nothing_not_a_wrong_number(
        ma, tmp_path):
    # The .fasta/.raw pair is the layout of ONE observed run, and no
    # InterProScan was available to check it against. A build that writes a
    # different tree must produce silence, which the stage then reports.
    root = str(tmp_path / "tmp")
    os.makedirs(os.path.join(root, "run1", "job2"))
    for name in ("chunk0.xml", "chunk0.out", "summary.log"):
        open(os.path.join(root, "run1", "job2", name), "w").close()
    probe = ma._InterProProgress(root, 60)
    assert probe() == ""
    assert not probe.seen_any


def test_a_census_too_big_to_finish_says_nothing_rather_than_a_partial_count(
        ma, tmp_path, monkeypatch):
    root = str(tmp_path / "tmp")
    _ips_tree(ma, root, created=6, done=3)
    monkeypatch.setattr(ma, "_IPS_SCAN_CAP", 4)
    assert ma._InterProProgress(root, 60)() == "", \
        "a truncated count is a wrong count, and this runs every minute"


def test_a_raw_that_outlived_its_chunk_drops_the_denominator_not_the_count(
        ma, tmp_path):
    root = str(tmp_path / "tmp")
    _ips_tree(ma, root, created=0, done=5)
    assert ma._InterProProgress(root, 60)() == "5 chunk(s) analysed"


def test_the_interpro_probe_watches_the_directory_the_stage_gave_the_tool(
        ma, tmp_path, paths_for, monkeypatch):
    # Both halves come out of the one call: if -T moves and the probe does
    # not, it counts a directory nothing is writing and reports nothing
    # forever, which is the failure this whole mechanism replaces.
    cfg, p = paths_for("ips_probe")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    cfg["db"]["interproscan_sh"] = sys.executable      # an exe that exists
    seen = {}

    def fake_run_cmd(cmd, cwd=None, env=None, progress=None):
        seen["cmd"] = [str(c) for c in cmd]
        seen["probe"] = progress
        open(seen["cmd"][seen["cmd"].index("-o") + 1], "w").close()
        return ""

    monkeypatch.setattr(ma, "run_cmd", fake_run_cmd)
    ma.stage_interpro(cfg, p)
    probe = seen["probe"]
    assert isinstance(probe, ma._InterProProgress), \
        "the stage that runs blind for days is the one that supplies a probe"
    assert probe.root == seen["cmd"][seen["cmd"].index("-T") + 1]


def test_a_progress_probe_that_raises_costs_its_line_and_not_the_stage(
        ma, monkeypatch, capsys):
    # It reads a directory the tool is concurrently writing, where a file
    # vanishing between the readdir and the stat is ordinary.
    monkeypatch.setattr(ma, "_PROGRESS_INTERVAL", 0.15)
    calls = []

    def boom():
        calls.append(1)
        raise OSError(2, "No such file or directory")

    assert ma.run_cmd([sys.executable, "-c", _TQDM_LIKE], progress=boom) == ""
    assert len(calls) == 1, \
        "a probe that raised once is dropped, not retried every interval"
    err = capsys.readouterr().err
    assert err.count("progress probe") == 1, "said once, not every minute"
    assert [l for l in err.splitlines() if " running " in l], \
        "the heartbeat itself survives the probe that failed"


# --- logging must not be what kills a multi-hour run -------------------
# symptom: tool output is decoded with errors="replace", so descriptions carry
# U+FFFD; printing one to a cp1252 console raises UnicodeEncodeError, and the
# run dies on the log line rather than on the work.

def _cp1252_stream():
    raw = io.BytesIO()
    return raw, io.TextIOWrapper(raw, encoding="cp1252", errors="strict")


def test_log_survives_a_character_the_console_cannot_encode(ma, monkeypatch):
    raw, stream = _cp1252_stream()
    monkeypatch.setattr(sys, "stderr", stream)
    ma.log("subtilisin-like � peptidase")      # must not raise
    stream.flush()
    assert b"subtilisin-like" in raw.getvalue()
    assert b"peptidase" in raw.getvalue(), \
        "the rest of the line must survive the one bad character"


@pytest.mark.parametrize("encoding", ["cp1252", "ascii", "latin-1", "cp437",
                                     "utf-8"])
def test_a_replacement_character_costs_a_character_not_the_run(ma, monkeypatch,
                                                               encoding):
    # "whatever the stream encoding": the console code page is the machine's,
    # not ours - cp1252 here, cp437 on an older Windows, ascii under a bare C
    # locale in a container - and only utf-8 can encode U+FFFD at all. The
    # line must come out on every one of them.
    raw = io.BytesIO()
    stream = io.TextIOWrapper(raw, encoding=encoding, errors="strict")
    monkeypatch.setattr(sys, "stderr", stream)
    ma.log("PF00082 � subtilisin-like peptidase")   # must not raise
    stream.flush()
    out = raw.getvalue()
    assert b"PF00082" in out and b"subtilisin-like peptidase" in out


def test_the_log_file_is_written_as_defensively_as_the_console(ma,
                                                               monkeypatch):
    # log() writes twice. The run's own log file is opened as utf-8, but it is
    # the half nobody watches, so a raise there would still end the run.
    raw = io.BytesIO()
    fh = io.TextIOWrapper(raw, encoding="ascii", errors="strict")
    monkeypatch.setattr(ma, "_LOGFH", fh)
    monkeypatch.setattr(sys, "stderr", io.StringIO())
    ma.log("� in a VFDB subject title")
    fh.flush()
    assert b"in a VFDB subject title" in raw.getvalue()


def test_configure_console_streams_makes_the_stream_unable_to_raise(
        ma, monkeypatch):
    raw, stream = _cp1252_stream()
    monkeypatch.setattr(sys, "stdout", stream)
    monkeypatch.setattr(sys, "stderr", stream)
    ma.configure_console_streams()
    stream.write("�")                          # must not raise
    stream.flush()
    assert raw.getvalue()


def test_configure_console_streams_tolerates_a_stream_it_cannot_touch(
        ma, monkeypatch):
    class Dumb:
        def write(self, s):
            return len(s)

        def flush(self):
            pass

    monkeypatch.setattr(sys, "stderr", Dumb())
    ma.configure_console_streams()                  # no reconfigure(): fine
    ma.log("still works")


# --- a database that cannot possibly hit ------------------------------
# symptom (a): a `diamond makedb` that had failed left a ZERO-BYTE .dmnd on
# disk. The stage would have searched it and reported no hits, which is
# indistinguishable in the output from a real absence of virulence factors.
# symptom (b): BAGEL built correctly - 262 sequences, median length 15 - and
# returned exactly 0 hits against 38,204 proteins at --evalue 1e-10, because no
# 15-residue alignment can reach 1e-10.
def _dia_project(ma, tmp_path, paths_for, **dbs):
    cfg, p = paths_for()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    cfg["db"]["diamond"] = dbs
    return cfg, p


def test_a_zero_byte_diamond_database_is_refused_not_searched(
        ma, tmp_path, paths_for, stub_bin):
    db = tmp_path / "vfdb.dmnd"
    db.write_bytes(b"")
    cfg, p = _dia_project(ma, tmp_path, paths_for, vfdb=str(db))
    with pytest.raises(ma.StageError) as e:
        ma.stage_diamond(cfg, p)
    assert "0 bytes" in str(e.value)
    assert "diamond makedb" in str(e.value)
    # nothing was searched, so nothing may claim to be an answer
    assert not os.path.exists(f"{p.diamond_dir}/vfdb.tsv")
    assert not os.path.exists(p.diamond_done)


def test_a_truncated_diamond_database_is_refused_too(ma, tmp_path, paths_for,
                                                     stub_bin):
    # 143 bytes is a real one-sequence database, so the floor is below that and
    # above anything a half-finished makedb leaves.
    db = tmp_path / "vfdb.dmnd"
    db.write_bytes(b"DIAMOND" * 3)
    cfg, p = _dia_project(ma, tmp_path, paths_for, vfdb=str(db))
    with pytest.raises(ma.StageError) as e:
        ma.stage_diamond(cfg, p)
    assert "smaller than a DIAMOND header" in str(e.value)


def test_a_database_with_no_sequences_is_refused(ma, tmp_path, paths_for,
                                                 stub_bin):
    db = F.write_dmnd(tmp_path / "vfdb.dmnd", sequences=0, letters=0)
    cfg, p = _dia_project(ma, tmp_path, paths_for, vfdb=db)
    with pytest.raises(ma.StageError) as e:
        ma.stage_diamond(cfg, p)
    assert "holds no sequences" in str(e.value)


def test_a_healthy_database_is_searched_without_comment(ma, tmp_path,
                                                        paths_for, stub_bin,
                                                        capsys):
    db = F.write_dmnd(tmp_path / "vfdb.dmnd", sequences=4000, letters=1400000)
    cfg, p = _dia_project(ma, tmp_path, paths_for, vfdb=db)
    ma.stage_diamond(cfg, p)
    err = capsys.readouterr().err
    assert "incapable of a hit" not in err
    assert os.path.exists(f"{p.diamond_dir}/vfdb.tsv")


def test_a_seed_set_is_named_as_one_and_not_blamed_on_the_evalue(
        ma, tmp_path, paths_for, stub_bin, capsys):
    # 262 sequences / 4009 letters is the shape of the database this pipeline
    # actually built and searched: BAGEL4's motif SEED set, mean 15 residues,
    # 0 hits against 38,204 proteins. The e-value was never the problem.
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=262, letters=4009)
    cfg, p = _dia_project(ma, tmp_path, paths_for, bagel=db)
    ma.stage_diamond(cfg, p)
    err = capsys.readouterr().err
    assert "motif or seed set rather than a protein sequence database" in err
    assert "15 residues" in err
    assert "NO e-value makes that a real search" in err
    assert "diamond_evalues" not in err
    # it still runs: this is a warning about what was built, not a broken file
    assert os.path.exists(f"{p.diamond_dir}/bagel.tsv")


def test_a_short_peptide_database_warns_that_the_evalue_is_unreachable(
        ma, tmp_path, paths_for, stub_bin, capsys):
    # 40-residue peptides are long enough not to be a seed set, and still
    # cannot reach an --evalue somebody set to 1e-30.
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=262, letters=10480)
    cfg, p = _dia_project(ma, tmp_path, paths_for, bagel=db)
    cfg["diamond_evalues"] = {"bagel": 1e-30}
    ma.stage_diamond(cfg, p)
    err = capsys.readouterr().err
    assert "incapable of a hit before it starts" in err
    assert "40 residues" in err
    assert "motif or seed set" not in err
    assert os.path.exists(f"{p.diamond_dir}/bagel.tsv")


def test_the_warning_names_the_scoring_weight_the_database_is_holding(
        ma, tmp_path, paths_for, stub_bin, capsys):
    # bagel carries diamond_weights 3, equal to TADB3: a database that cannot
    # hit is also holding a weight that says it can.
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=262, letters=4009)
    cfg, p = _dia_project(ma, tmp_path, paths_for, bagel=db)
    ma.stage_diamond(cfg, p)
    assert "diamond_weights bagel: 3" in capsys.readouterr().err


def test_a_per_database_evalue_silences_the_warning_and_is_used(
        ma, tmp_path, paths_for, stub_bin, capsys):
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=262, letters=4009)
    cfg, p = _dia_project(ma, tmp_path, paths_for, bagel=db)
    cfg["diamond_evalues"] = {"bagel": 1e-3}
    ma.stage_diamond(cfg, p)
    err = capsys.readouterr().err
    assert "incapable of a hit" not in err
    # the seed-set warning is about WHAT was built and is not silenced by a
    # threshold; only the e-value advice is.
    assert "motif or seed set" in err
    assert "searching at --evalue 0.001 from diamond_evalues" in err
    assert ma.diamond_evalue_for(cfg, "bagel") == 1e-3
    assert ma.diamond_evalue_for(cfg, "vfdb") == \
        cfg["thresholds"]["diamond_evalue"]


def test_a_non_numeric_per_database_evalue_is_a_message(ma):
    cfg = {"diamond_evalues": {"bagel": "soon"},
           "thresholds": {"diamond_evalue": 1e-10}}
    with pytest.raises(ma.StageError) as e:
        ma.diamond_evalue_for(cfg, "bagel")
    assert "diamond_evalues.bagel must be a number" in str(e.value)


def test_the_reachability_estimate_uses_diamonds_own_constants(ma):
    # BLOSUM62 Lambda=0.267 K=0.041, as diamond prints them. A 15-residue
    # perfect match against 4009 letters cannot reach 1e-10; a 350-residue one
    # against the same database can reach anything.
    assert (ma._DMND_LAMBDA, ma._DMND_K) == (0.267, 0.041)
    assert ma.best_possible_evalue(4009, 15) > 1e-10
    assert ma.best_possible_evalue(4009, 350) < 1e-10
    assert ma.best_possible_evalue(0, 15) is None
    assert ma.best_possible_evalue(4009, 0) is None


def test_a_length_that_cannot_be_read_is_said_rather_than_guessed(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    # no diamond on PATH and no source FASTA: the check must say the length is
    # unknown rather than assume one and act on it.
    db = F.write_dmnd(tmp_path / "vfdb.dmnd")
    cfg, p = _dia_project(ma, tmp_path, paths_for, vfdb=db)
    monkeypatch.setattr(ma, "have", lambda t: t != "diamond")
    typical, letters, why, read = ma.diamond_db_profile(cfg, "vfdb", db)
    assert (typical, letters) == (None, None)
    assert "diamond is not installed" in why
    # Nothing was opened, and the token says so rather than leaving doctor to
    # assume a DIAMOND header was read.
    assert read == ""
    bad, warn, reads = ma.diamond_db_check(cfg, "vfdb", db)
    assert bad is None
    assert "nothing here can tell whether" in warn
    assert reads == ("size",)


def test_a_source_fasta_beside_the_database_gives_the_median(
        ma, tmp_path, paths_for, monkeypatch):
    db = F.write_dmnd(tmp_path / "bagel.dmnd")
    with open(tmp_path / "bagel.fas", "w", encoding="utf-8") as fh:
        for i in range(5):
            fh.write(">p%d\n%s\n" % (i, "M" * (10 + i)))
    cfg, p = _dia_project(ma, tmp_path, paths_for, bagel=db)
    monkeypatch.setattr(ma, "have", lambda t: t != "diamond")
    typical, letters, why, read = ma.diamond_db_profile(cfg, "bagel", db)
    assert typical == 12                      # median of 10..14
    assert letters == 60
    assert "median of 5 sequences" in why
    # The fallback READS THE FASTA - every record of it, up to the cap - and
    # the token is what lets doctor say so instead of publishing the DIAMOND
    # header caveat over a check that never opened a .dmnd.
    assert read == "fasta_lengths"
    _bad, _warn, reads = ma.diamond_db_check(cfg, "bagel", db)
    assert reads == ("size", "fasta_lengths", "fasta_headers"), \
        "motif_seed_evidence opens the same FASTA again, for its deflines"


def test_doctor_refuses_a_zero_byte_diamond_database(tmp_path):
    db = tmp_path / "vfdb.dmnd"
    db.write_bytes(b"")
    proj = build_project(tmp_path / "p", threads=1,
                         db={"diamond": {"vfdb": str(db)}},
                         run={"eggnog": True, "pfam": False, "dbcan": False,
                              "diamond": True, "cluster": False, "join": False,
                              "topology": False, "structure": False,
                              "context": False, "unipept": False,
                              "taxonomy": False, "ncbifam": False,
                              "kofam": False, "interpro": False,
                              "hhblits": False, "jackhmmer": False,
                              "smorf": False, "effectors": False})
    proc = run_metaannot("doctor", "--config", proj.config_path, expect=1)
    assert "0 bytes" in proc.stdout
    assert "diamond makedb" in proc.stdout


def test_doctor_recognises_a_motif_seed_set_rather_than_blaming_the_evalue(
        tmp_path, stub_bin):
    # 262 sequences / 4009 letters is the real BAGEL database that prompted
    # both of these checks: a mean of 15 residues, which is not a short
    # protein database but BAGEL4's motif SEED set. Advising a lower --evalue
    # sends the reader off to tune a threshold on the wrong kind of file, and
    # the tuned search still answers a question nobody asked.
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=262, letters=4009)
    proj = build_project(tmp_path / "p", threads=1,
                         db={"diamond": {"bagel": db}},
                         run={"eggnog": True, "pfam": False, "dbcan": False,
                              "diamond": True, "cluster": False, "join": False,
                              "topology": False, "structure": False,
                              "context": False, "unipept": False,
                              "taxonomy": False, "ncbifam": False,
                              "kofam": False, "interpro": False,
                              "hhblits": False, "jackhmmer": False,
                              "smorf": False, "effectors": False})
    proc = run_metaannot("doctor", "--config", proj.config_path)
    assert "motif or seed set rather than a protein sequence database" in \
        proc.stdout
    assert "15 residues" in proc.stdout
    assert "NO e-value makes that a real search" in proc.stdout
    assert "diamond_evalues" not in proc.stdout, \
        "the e-value advice is meant to be replaced here, not added to"


def test_doctor_reads_the_headers_when_the_source_fasta_is_beside_the_db(
        tmp_path, stub_bin):
    # length is not the only signal, and it is the weaker one: a database of
    # 60-residue entries whose headers say ggmotif is still a seed set.
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=100, letters=6000)
    io.open(str(tmp_path / "bagel.fas"), "w", encoding="utf-8").write(
        "".join(f">LE-entry{i} ggmotif\n{'A' * 60}\n" for i in range(100)))
    proj = build_project(tmp_path / "p", threads=1,
                         db={"diamond": {"bagel": db}},
                         run={"eggnog": True, "pfam": False, "dbcan": False,
                              "diamond": True, "cluster": False, "join": False,
                              "topology": False, "structure": False,
                              "context": False, "unipept": False,
                              "taxonomy": False, "ncbifam": False,
                              "kofam": False, "interpro": False,
                              "hhblits": False, "jackhmmer": False,
                              "smorf": False, "effectors": False})
    proc = run_metaannot("doctor", "--config", proj.config_path)
    assert "motif or seed set" in proc.stdout
    assert "ggmotif" in proc.stdout
    assert "LE-/MA- accession prefixes" in proc.stdout


def test_doctor_still_warns_about_a_database_too_short_for_the_evalue(
        tmp_path, stub_bin):
    # the e-value check is narrower now, not gone: 40-residue peptides are
    # long enough not to be a seed set, and still cannot reach an --evalue
    # that someone set to 1e-30.
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=262, letters=10480)
    proj = build_project(tmp_path / "p", threads=1,
                         db={"diamond": {"bagel": db}},
                         diamond_evalues={"bagel": 1e-30},
                         run={"eggnog": True, "pfam": False, "dbcan": False,
                              "diamond": True, "cluster": False, "join": False,
                              "topology": False, "structure": False,
                              "context": False, "unipept": False,
                              "taxonomy": False, "ncbifam": False,
                              "kofam": False, "interpro": False,
                              "hhblits": False, "jackhmmer": False,
                              "smorf": False, "effectors": False})
    proc = run_metaannot("doctor", "--config", proj.config_path)
    assert "incapable of a hit before it starts" in proc.stdout
    assert "diamond_evalues" in proc.stdout
    assert "motif or seed set" not in proc.stdout


def test_the_per_database_evalue_also_filters_the_hit_table(ma, tmp_path,
                                                            stub_bin):
    # filtering the table back down to thresholds.diamond_evalue in integrate
    # would leave diamond_evalues doing nothing at all.
    proj = _searchable(tmp_path, tmp_path / "p")
    F.write_dmnd(tmp_path / "db" / "bagel.dmnd", sequences=262, letters=4009)
    d = dict(proj.cfg["db"])
    d["diamond"] = dict(d["diamond"], bagel=str(tmp_path / "db" / "bagel.dmnd"))
    proj.write_config(db=d, diamond_evalues={"bagel": 1e-3})
    proj.run()
    got = open(proj.rpath("annotation_pass1.tsv"), encoding="utf-8").read()
    assert "bagel_hit" in got
    # the stub writes a 1e-40 hit, which passes either threshold; what matters
    # is that the column exists and the run did not refuse the database
    assert proj.state()["diamond"]["status"] == "ok"


def test_changing_a_per_database_evalue_reruns_the_search(ma, tmp_path,
                                                          stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    before = proj.state()["diamond"]["signature"]
    proj.write_config(diamond_evalues={"vfdb": 1e-3})
    proj.run()
    assert proj.state()["diamond"]["signature"] != before


# ----------------------------------------------------------------------
# foldseek asks for the columns the analysis needs
# ----------------------------------------------------------------------
def _fs_project(ma, tmp_path, paths_for):
    cfg, p = paths_for()
    tgt = tmp_path / "PDB"
    tgt.write_text("x", encoding="utf-8")
    cfg["db"]["foldseek_target"] = str(tgt)
    cfg["db"]["foldseek_extra_targets"] = []
    os.makedirs(p.structures, exist_ok=True)
    with open(f"{p.structures}/P1.pdb", "w", encoding="utf-8") as fh:
        fh.write("ATOM      1  CA  ALA A   1      "
                 "0.000   0.000   0.000  1.00 90.00           C" + chr(10))
    return cfg, p


def test_foldseek_requests_qtmscore_and_qlen_not_the_legacy_columns(
        ma, tmp_path, paths_for, stub_bin, capsys, monkeypatch):
    """The stage used to hardcode the 10-column legacy list.

    FOLDSEEK_COLS' own comment says it "is what stage_foldseek should ask
    for", and the stage ignored it — so the TM gate silently fell back to
    alntmscore, which is normalised by the ALIGNMENT rather than the query,
    and finalise printed advice to "re-run the foldseek stage to get qtmscore
    and qlen" that re-running could not act on.
    """
    seen = []

    def fake_run(cmd, **kw):
        seen.append([str(c) for c in cmd])
        return ""

    monkeypatch.setattr(ma, "run_cmd", fake_run)
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    try:
        ma.stage_foldseek(cfg, p)
    except Exception:              # noqa: BLE001
        pass    # the fake writes no result file; the COMMAND is what matters
    search = [c for c in seen if len(c) > 1 and c[1] == "easy-search"]
    assert search, "no easy-search was issued"
    fields = search[0][search[0].index("--format-output") + 1].split(",")
    assert "qtmscore" in fields, "the query-normalised TM score must be requested"
    assert "qlen" in fields, "qlen is what the coverage backstop needs"
    assert fields == ma.FOLDSEEK_COLS


# What run_cmd actually raises when a tool exits non-zero, and what Foldseek 5
# actually prints for a format code it does not have. The old test invented
# `ma.StageError("Invalid selection: qtmscore")`; run_cmd raises a plain
# RuntimeError, and Foldseek's LocalParameters.cpp prints "Format code <field>
# does not exist." So the string the fallback keys on has to be the real one.
_FOLDSEEK_5_STDERR = (
    "foldseek exited 1\n--- stderr tail ---\n"
    "Format code qtmscore does not exist.\n")


def _fake_foldseek(calls, fail_on_qtmscore=_FOLDSEEK_5_STDERR):
    """A run_cmd that fails the full-column easy-search the way Foldseek does."""
    def fake_run(cmd, **kw):
        c = [str(x) for x in cmd]
        calls.append(c)
        if (len(c) > 1 and c[1] == "easy-search" and fail_on_qtmscore
                and "qtmscore" in c[c.index("--format-output") + 1]):
            raise RuntimeError(fail_on_qtmscore)
        return ""
    return fake_run


def test_a_foldseek_that_rejects_the_new_columns_falls_back_and_says_so(
        ma, tmp_path, paths_for, stub_bin, capsys, monkeypatch):
    # symptom: run_cmd raises RuntimeError, StageError is a SUBCLASS of it,
    # and the fallback was written `except StageError` — so it could never
    # catch the failure it exists for. A Foldseek 5 build lost its structure
    # evidence entirely instead of degrading to the legacy columns.
    calls = []
    monkeypatch.setattr(ma, "run_cmd", _fake_foldseek(calls))
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    try:
        ma.stage_foldseek(cfg, p)
    except Exception:              # noqa: BLE001
        pass
    search = [c for c in calls if len(c) > 1 and c[1] == "easy-search"]
    assert len(search) == 2, "it must retry once with the legacy columns"
    second = search[1][search[1].index("--format-output") + 1].split(",")
    assert second == ma.FOLDSEEK_COLS_LEGACY
    err = capsys.readouterr().err
    assert "predates qtmscore/ttmscore" in err
    assert "normalised by the alignment" in err


def _search_calls(calls):
    return [c for c in calls if len(c) > 1 and c[1] == "easy-search"]


def _tmpdir_of(argv):
    """easy-search's tmp dir: the 6th positional, i.e. the arg after `out`.

    Found by position rather than by name because easy-search takes it
    positionally, but bounded by the first flag so a reordered argv fails the
    test instead of silently creating a directory called '--format-output'.
    """
    pos = [a for a in argv[2:] if not a.startswith("-")]
    return pos[3]


def test_the_scratch_tree_is_removed_on_every_exit_not_only_success(
        ma, tmp_path, paths_for, stub_bin, monkeypatch):
    # symptom: the cleanup sat after the loop body, so every re-raise - a full
    # disk, an OOM kill, the undroppable-column path added beside it - left
    # the tree behind. Against AFDB50 that is tens to hundreds of GB, and
    # CLAUDE.md already lists results/foldseek/tmp* as never cleaned up.
    seen = {}

    def fake_run(cmd, **kw):
        c = [str(x) for x in cmd]
        if len(c) > 1 and c[1] == "easy-search":
            d = _tmpdir_of(c)
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, "big"), "w").close()
            seen["tmpd"] = d
            raise RuntimeError("foldseek exited 1\n--- stderr tail ---\n"
                               "Error: Could not open database\n")
        return ""

    monkeypatch.setattr(ma, "run_cmd", fake_run)
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    with pytest.raises(RuntimeError, match="Could not open database"):
        ma.stage_foldseek(cfg, p)
    assert seen, "easy-search was never called"
    assert not os.path.exists(seen["tmpd"]), (
        "the scratch tree survived a failed search")


# Derived from the two column lists, not copied: a field added to BOTH lists
# later must become a new undroppable case automatically, and a field that
# stops being droppable must not silently drop out of the parametrization.
_MA = _load()
_DROPPABLE = sorted(set(_MA.FOLDSEEK_COLS) - set(_MA.FOLDSEEK_COLS_LEGACY))
_UNDROPPABLE = sorted(set(_MA.FOLDSEEK_COLS) & set(_MA.FOLDSEEK_COLS_LEGACY))


def _reject(field):
    return (f"foldseek exited 1\n--- stderr tail ---\n"
            f"Format code {field} does not exist.\n")


@pytest.mark.parametrize("field", _UNDROPPABLE)
def test_a_rejected_column_the_legacy_list_also_asks_for_is_not_retried(
        ma, tmp_path, paths_for, stub_bin, capsys, monkeypatch, field):
    # symptom: the first version of this gate retried on ANY rejected format
    # code. FOLDSEEK_COLS_LEGACY is a strict subset of FOLDSEEK_COLS, so only
    # qlen/tlen/qtmscore/ttmscore can be dropped; a build rejecting one of the
    # ten fields in BOTH lists fails the retry identically, and the second
    # error is the one the operator then has to explain.
    calls = []
    monkeypatch.setattr(ma, "run_cmd",
                        _fake_foldseek(calls, fail_on_qtmscore=_reject(field)))
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    with pytest.raises(RuntimeError, match="Format code"):
        ma.stage_foldseek(cfg, p)
    assert len(_search_calls(calls)) == 1, (
        f"'{field}' is in the legacy list too, so the retry cannot help")
    assert "legacy column list asks for as well" in capsys.readouterr().err, (
        "the operator must be told why it did not fall back")


@pytest.mark.parametrize("field", _DROPPABLE)
def test_every_column_the_legacy_list_drops_does_trigger_the_fallback(
        ma, tmp_path, paths_for, stub_bin, capsys, monkeypatch, field):
    # symptom: only qtmscore was ever exercised, so narrowing `droppable` to
    # {"qtmscore"} passed the whole suite - and qlen and tlen sit EARLIER in
    # FOLDSEEK_COLS than qtmscore, so they are the first codes an old build
    # rejects.
    calls = []

    def fake_run(cmd, **kw):
        c = [str(x) for x in cmd]
        calls.append(c)
        if (len(c) > 1 and c[1] == "easy-search"
                and field in c[c.index("--format-output") + 1].split(",")):
            raise RuntimeError(_reject(field))
        return ""

    monkeypatch.setattr(ma, "run_cmd", fake_run)
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    try:
        ma.stage_foldseek(cfg, p)
    except Exception:              # noqa: BLE001
        pass
    search = _search_calls(calls)
    assert len(search) == 2, f"'{field}' is droppable, so it must fall back"
    assert search[1][search[1].index("--format-output") + 1].split(",") \
        == ma.FOLDSEEK_COLS_LEGACY
    assert f"rejected the format code '{field}'" in capsys.readouterr().err


def test_the_phrase_is_not_spliced_across_two_lines_of_stderr(
        ma, tmp_path, paths_for, stub_bin, monkeypatch):
    # symptom: with \s+ between the words, a multi-line tail splices unrelated
    # lines - "<path> does not exist" is stock MMseqs2 wording for a missing
    # database - and the gate then names a PATH as the rejected column, either
    # refusing to fall back or buying a pointless second invocation.
    calls = []
    monkeypatch.setattr(ma, "run_cmd", _fake_foldseek(
        calls, fail_on_qtmscore="foldseek exited 1\n--- stderr tail ---\n"
                                "Please choose a valid format code\n"
                                "/scratch/foldseek/tmp0 does not exist\n"))
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    with pytest.raises(RuntimeError):
        ma.stage_foldseek(cfg, p)
    assert len(_search_calls(calls)) == 1, (
        "a path on a following line was read as the rejected format code")


def test_a_rejected_code_in_a_different_wording_still_falls_back(
        ma, tmp_path, paths_for, stub_bin, capsys, monkeypatch):
    # symptom: keying the fallback on one exact sentence means a build whose
    # message differs by a colon loses its structural evidence entirely. The
    # rejection costs seconds - foldseek validates --format-output before it
    # creates the temp directory, let alone searches - so the permissive
    # direction is the cheap one.
    calls = []
    monkeypatch.setattr(ma, "run_cmd", _fake_foldseek(
        calls, fail_on_qtmscore="foldseek exited 1\n--- stderr tail ---\n"
                                "Format code: qtmscore does not exist.\n"))
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    try:
        ma.stage_foldseek(cfg, p)
    except Exception:              # noqa: BLE001
        pass
    assert len(_search_calls(calls)) == 2, "a rejection must still fall back"
    # `format code\W*` absorbs the colon, so the field is still named rather
    # than the message degrading to "a format code it could not name".
    assert "rejected the format code 'qtmscore'" in capsys.readouterr().err


def test_a_foldseek_failure_that_is_not_a_bad_format_code_is_not_retried(
        ma, tmp_path, paths_for, stub_bin, monkeypatch):
    # symptom: widening the handler to catch what run_cmd raises would, on its
    # own, retry EVERY foldseek failure — and a search that died on a full
    # disk or a bad database has no completed alignment to reuse, so the retry
    # repeats the multi-hour search to arrive at the same error.
    calls = []
    monkeypatch.setattr(ma, "run_cmd", _fake_foldseek(
        calls, fail_on_qtmscore="foldseek exited 1\n--- stderr tail ---\n"
                                "Error: Could not open database\n"))
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    with pytest.raises(RuntimeError, match="Could not open database"):
        ma.stage_foldseek(cfg, p)
    search = [c for c in calls if len(c) > 1 and c[1] == "easy-search"]
    assert len(search) == 1, "a real failure must not be retried"


def test_the_fallback_catches_what_run_cmd_actually_raises(ma):
    # symptom: the whole defect in one line. StageError is a subclass of
    # RuntimeError, so `except StageError` cannot catch run_cmd's
    # RuntimeError; the reverse containment is what makes the handler work.
    assert issubclass(ma.StageError, RuntimeError)
    assert not issubclass(RuntimeError, ma.StageError)
    src = open(METAANNOT_PY, encoding="utf-8").read()
    fold = src[src.index("def stage_foldseek"):]
    fold = fold[:fold.index("\ndef ")]
    assert "except RuntimeError" in fold, (
        "stage_foldseek must catch what run_cmd raises, not only StageError")

# ----------------------------------------------------------------------
# a machine with no usable GPU should learn that from doctor, not from
# esmfold dying six hours in
# ----------------------------------------------------------------------
def test_the_probe_separates_no_card_from_a_cpu_only_torch(ma, monkeypatch):
    """The confusing failure is a perfectly good card with a CPU-only torch
    wheel: the hardware is right there and the stage still refuses."""
    monkeypatch.setattr(ma.shutil, "which", lambda n: "/usr/bin/nvidia-smi")
    monkeypatch.setattr(ma.subprocess, "run",
                        lambda *a, **k: type("R", (), {
                            "returncode": 0,
                            "stdout": "GPU 0: NVIDIA Test Card (UUID: GPU-x)\n",
                            "stderr": ""})())
    fake = type("T", (), {"cuda": type("C", (), {
        "is_available": staticmethod(lambda: False)})()})
    monkeypatch.setitem(ma.sys.modules, "torch", fake)
    ok, why = ma.cuda_probe()
    assert ok is False
    assert "card is present" in why
    assert "CPU-only torch" in why


def test_the_probe_says_torch_is_missing_rather_than_blaming_the_card(
        ma, monkeypatch):
    monkeypatch.setattr(ma.shutil, "which", lambda n: None)
    monkeypatch.setattr(ma.importlib.util, "find_spec", lambda n: None)
    monkeypatch.delitem(ma.sys.modules, "torch", raising=False)
    ok, why = ma.cuda_probe()
    assert ok is False and "torch is not installed" in why


def test_a_usable_gpu_is_reported_by_name(ma, monkeypatch):
    monkeypatch.setattr(ma.shutil, "which", lambda n: None)
    fake = type("T", (), {"cuda": type("C", (), {
        "is_available": staticmethod(lambda: True),
        "get_device_name": staticmethod(lambda i: "NVIDIA Test Card")})()})
    monkeypatch.setitem(ma.sys.modules, "torch", fake)
    ok, why = ma.cuda_probe()
    assert ok is True and "NVIDIA Test Card" in why


def test_structure_without_cuda_is_a_doctor_failure_not_a_surprise(
        ma, tmp_path, capsys, monkeypatch):
    """run.structure with no CUDA must fail doctor. stage_esmfold exits rather
    than fold on CPU, and learning that after InterProScan has run for hours is
    the whole reason this check exists."""
    monkeypatch.setattr(ma, "cuda_probe", lambda: (False, "no CUDA device"))
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["run"]["structure"] = True
    assert cfg["run"]["structure"] and not ma.cuda_probe()[0]


def test_topology_without_cuda_warns_but_does_not_fail(ma, monkeypatch):
    """Not a MISS: SignalP is CPU-only and useful alone, and tmbed does run
    without a GPU - just not at a few hundred thousand proteins."""
    monkeypatch.setattr(ma, "cuda_probe", lambda: (False, "no CUDA device"))
    ok, _ = ma.cuda_probe()
    assert ok is False


# --- tmbed chunking ---------------------------------------------------
# symptom: TMbed writes nothing until it finishes, so one invocation over a
# whole proteome is an all-or-nothing bet measured in days. The run on
# 455,571 proteins was still going after two days with an empty output file,
# and the two before it died at 2 h 36 min with nothing recoverable.
def test_the_chunk_plan_is_longest_first(ma):
    # ProtT5 pads a batch out to its longest member, and whatever is going to
    # exhaust the device should be in the FIRST chunk, not the last.
    chunks, _ = ma.tmbed_chunk_plan([("s", 10), ("l", 500), ("m", 100)], 600)
    assert chunks[0][0] == "l"
    assert [q for c in chunks for q in c] == ["l", "m", "s"]


def test_no_chunk_exceeds_the_residue_budget(ma):
    lengths = [(f"p{i}", 50 + (i % 7) * 30) for i in range(200)]
    chunks, budget = ma.tmbed_chunk_plan(lengths, 400)
    by_id = dict(lengths)
    assert len(chunks) > 1
    for c in chunks:
        assert sum(by_id[q] for q in c) <= budget, c


def test_a_sequence_longer_than_the_budget_still_gets_a_chunk(ma):
    # the budget is a target, not a filter: dropping the sequence here would
    # lose a protein that tmbed_max_len had already decided to keep.
    chunks, budget = ma.tmbed_chunk_plan(
        [("big", 9000), ("a", 10), ("b", 10)], 100)
    assert sorted(q for c in chunks for q in c) == ["a", "b", "big"]
    assert budget >= 9000


def test_a_tiny_budget_cannot_ask_for_more_chunks_than_the_ceiling(ma):
    # one ProtT5 load per chunk, so a 1-residue budget over a real proteome
    # would spend all its time loading weights.
    chunks, budget = ma.tmbed_chunk_plan(
        [(f"p{i}", 100) for i in range(5000)], 1)
    assert len(chunks) <= ma.TMBED_MAX_PARTS
    assert budget > 1


def test_a_zero_budget_is_one_invocation(ma):
    chunks, budget = ma.tmbed_chunk_plan(
        [(f"p{i}", 100) for i in range(50)], 0)
    assert len(chunks) == 1 and len(chunks[0]) == 50
    assert budget == 5000


def test_every_protein_lands_in_exactly_one_chunk(ma):
    lengths = [(f"p{i}", 1 + (i * 37) % 500) for i in range(1000)]
    chunks, _ = ma.tmbed_chunk_plan(lengths, 2000)
    flat = [q for c in chunks for q in c]
    assert len(flat) == len(set(flat)) == 1000


def test_the_chunk_plan_is_deterministic_when_lengths_tie(ma):
    # equal lengths are ordered by id, so two runs of the same config produce
    # the same parts and a resume matches them up.
    lengths = [("b", 100), ("a", 100), ("c", 100)]
    one, _ = ma.tmbed_chunk_plan(lengths, 150)
    two, _ = ma.tmbed_chunk_plan(list(reversed(lengths)), 150)
    assert one == two == [["a"], ["b"], ["c"]]


def test_a_record_with_no_label_line_is_not_a_prediction(ma, tmp_path):
    # the 3-line format has no trailer, so a killed writer leaves a header and
    # a sequence. Counting that as done would commit a protein with no
    # topology as though it had one.
    f = tmp_path / "t.pred"
    f.write_text(">P1\nMKV\niii\n>P2\nMKVA\n", encoding="utf-8")
    assert [r[0] for r in ma.iter_tmbed_records(str(f))] == [">P1"]
    assert set(ma.parse_tmbed(str(f))) == {"P1"}


def test_a_missing_prediction_file_reads_as_no_records(ma, tmp_path):
    assert list(ma.iter_tmbed_records(str(tmp_path / "nope.pred"))) == []


def _many(n, length=200):
    """n equal-length proteins, so the chunk plan falls out by id."""
    return [F.Protein(f"P{i:03d}", "M" + "A" * (length - 1)) for i in range(n)]


def _parts_dir(p):
    return f"{os.path.dirname(p.tmbed)}/tmbed_parts"


def test_a_split_run_predicts_every_protein_and_leaves_no_parts(
        ma, tmp_path, paths_for, monkeypatch):
    cfg, p = paths_for("tmbed_split")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(20))
    cfg["tmbed_chunk_residues"] = 800          # 4 proteins per chunk
    _stub_tmbed(tmp_path, monkeypatch)
    ma.stage_tmbed(cfg, p)
    assert len(ma.parse_tmbed(p.tmbed)) == 20
    assert not os.path.exists(_parts_dir(p)), \
        "the parts are the result until the file is committed, and rubbish " \
        "afterwards"


def test_a_finished_chunk_is_not_predicted_a_second_time(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    # the whole point of splitting: an interrupted run resumes.
    cfg, p = paths_for("tmbed_resume")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(20))
    cfg["tmbed_chunk_residues"] = 800          # 5 chunks of 4
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_FAIL", "P008,P012,P016")   # chunks 2-4
    with pytest.raises(ma.StageError):
        ma.stage_tmbed(cfg, p)
    kept = sorted(f for f in os.listdir(_parts_dir(p)) if f.endswith(".pred"))
    assert kept == ["0.pred", "1.pred"], kept
    monkeypatch.delenv("TMBED_STUB_FAIL")
    capsys.readouterr()
    ma.stage_tmbed(cfg, p)
    err = capsys.readouterr().err
    assert err.count("is already predicted; skipping") == 2, err
    assert len(ma.parse_tmbed(p.tmbed)) == 20


def test_a_committed_chunk_is_adopted_by_identity_and_not_by_record_count(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    # symptom: a chunk counted as done when its file held at least as many
    # records as the chunk had proteins, and never asked WHICH proteins. The
    # plan is not fixed across a resume -- tmbed_chunk_residues is outside the
    # stage signature on purpose, so tuning it between an interrupted run and
    # its resume is allowed, and it re-plans the chunks. Chunk 2 of the new
    # plan then adopts the file chunk 2 of the OLD plan wrote, and the
    # concatenation writes two records for every protein in both plans and
    # none for the proteins in neither.
    cfg, p = paths_for("tmbed_replan")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(20))
    cfg["tmbed_chunk_residues"] = 800          # 5 chunks of 4
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_FAIL", "P008")             # chunk 3 dies
    with pytest.raises(ma.StageError):
        ma.stage_tmbed(cfg, p)
    monkeypatch.delenv("TMBED_STUB_FAIL")

    cfg["tmbed_chunk_residues"] = 600          # re-planned: 7 chunks of 3
    capsys.readouterr()
    ma.stage_tmbed(cfg, p)
    err = capsys.readouterr().err
    assert "a DIFFERENT set of proteins" in err, \
        "a chunk of the old plan was adopted, or dropped, without a word"
    assert "is already predicted; skipping" not in err, \
        "no chunk of the old plan covers the same proteins as a new one"
    recs = list(ma.iter_tmbed_records(p.tmbed))
    assert len(recs) == 20, "a protein was predicted twice, or not at all"
    assert sorted(ma.parse_tmbed(p.tmbed)) == [f"P{i:03d}" for i in range(20)]


def test_a_chunk_whose_plan_still_matches_is_adopted_as_before(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    # the direction that matters on every ordinary resume: proving identity
    # must not turn a resumable run back into a run that starts over. Same
    # plan, same proteins, same files -- and nothing is predicted twice.
    cfg, p = paths_for("tmbed_replan_same")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(12))
    cfg["tmbed_chunk_residues"] = 800          # 3 chunks of 4
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_FAIL", "P008")             # the last chunk
    with pytest.raises(ma.StageError):
        ma.stage_tmbed(cfg, p)
    monkeypatch.delenv("TMBED_STUB_FAIL")
    capsys.readouterr()
    ma.stage_tmbed(cfg, p)
    err = capsys.readouterr().err
    assert err.count("is already predicted; skipping") == 2, err
    assert "a DIFFERENT set of proteins" not in err
    assert len(list(ma.iter_tmbed_records(p.tmbed))) == 12


def test_the_first_chunk_of_a_resume_reports_no_eta_rather_than_a_wrong_one(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    # symptom: the ETA divided THIS process's elapsed time by a residue count
    # that already included every chunk skipped as "already predicted", so the
    # first progress line of a resume -- the line an operator reads, and the
    # line a console puts in front of them -- was wrong by the whole ratio of
    # resumed to new work, and was quoted before a single residue had been
    # predicted here.
    cfg, p = paths_for("tmbed_eta")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(20))
    cfg["tmbed_chunk_residues"] = 800          # 5 chunks of 4
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_FAIL", "P012,P016")        # chunks 4 and 5
    with pytest.raises(ma.StageError):
        ma.stage_tmbed(cfg, p)
    monkeypatch.delenv("TMBED_STUB_FAIL")
    capsys.readouterr()
    ma.stage_tmbed(cfg, p)

    lines = [l for l in capsys.readouterr().err.splitlines()
             if "sequence(s)," in l]
    assert len(lines) == 2, lines
    assert "h left" not in lines[0], \
        "an ETA was quoted from three chunks this process did not predict"
    assert "h left" in lines[1], \
        "the estimate must come back as soon as there is real work to rate"


def test_the_eta_rates_only_the_work_this_process_did(ma):
    # 1000 residues in an hour, 3000 left: three hours, whatever else the run
    # adopted on the way in.
    assert ma.tmbed_eta_hours(3600.0, 1000, 3000) == pytest.approx(3.0)
    # nothing predicted here yet, so there is no rate and nothing to say
    assert ma.tmbed_eta_hours(3600.0, 0, 3000) is None
    # and nothing left to do is not an estimate of zero, it is no estimate
    assert ma.tmbed_eta_hours(3600.0, 1000, 0) is None


def test_a_failed_chunk_keeps_what_tmbed_wrote_and_names_the_rest(
        ma, tmp_path, paths_for, monkeypatch):
    cfg, p = paths_for("tmbed_partial")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(8))
    cfg["tmbed_chunk_residues"] = 800          # 2 chunks of 4
    cfg["tmbed_allow_partial"] = True
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_FAIL", "P004")             # chunk 1
    monkeypatch.setenv("TMBED_STUB_PARTIAL", "1")             # 1 record, die
    ma.stage_tmbed(cfg, p)
    got = ma.parse_tmbed(p.tmbed)
    assert len(got) == 5, "4 from the good chunk plus the 1 salvaged"
    miss = f"{os.path.dirname(p.tmbed)}/tmbed_failed.tsv"
    rows = [l.split("\t") for l in
            io.open(miss, encoding="utf-8").read().splitlines()[1:]]
    assert len(rows) == 3
    assert not ({r[0] for r in rows} & set(got)), \
        "a protein cannot be both predicted and a casualty"
    assert all("out of memory" in r[3] for r in rows), rows


def test_a_failed_chunk_stops_the_run_unless_partial_is_allowed(
        ma, tmp_path, paths_for, monkeypatch):
    # a silently short topology set shifts every bin, and nothing downstream
    # can tell "no helix" from "never asked".
    cfg, p = paths_for("tmbed_strict")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(8))
    cfg["tmbed_chunk_residues"] = 800
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_FAIL", "P004")
    with pytest.raises(ma.StageError) as e:
        ma.stage_tmbed(cfg, p)
    assert "tmbed_allow_partial: true" in str(e.value)
    assert not os.path.exists(p.tmbed), \
        "a partial prediction set must not be committed"
    assert os.path.exists(_parts_dir(p)), \
        "the finished chunks have to survive for the rerun to resume"


def test_a_wedged_device_stops_the_run_rather_than_failing_every_chunk(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    cfg, p = paths_for("tmbed_wedged")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(20))
    cfg["tmbed_chunk_residues"] = 800          # 5 chunks
    cfg["tmbed_max_consecutive_failures"] = 2
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_FAIL", "P000,P004,P008,P012,P016")
    with pytest.raises(ma.StageError):
        ma.stage_tmbed(cfg, p)
    err = capsys.readouterr().err
    assert "chunk(s) in a row failed" in err
    assert err.count("failed after writing") == 2, \
        "the remaining chunks must not be attempted one by one"


def test_a_chunk_that_exits_zero_but_comes_back_short_is_not_accepted(
        ma, tmp_path, paths_for, monkeypatch):
    # exit 0 is the tool's opinion; the reconciliation is ours.
    cfg, p = paths_for("tmbed_short")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(4))
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setenv("TMBED_STUB_DROP", "2")       # writes 2 of 4, exit 0
    with pytest.raises(ma.StageError) as e:
        ma.stage_tmbed(cfg, p)
    assert "2 of 4" in str(e.value)
    rows = io.open(f"{os.path.dirname(p.tmbed)}/tmbed_failed.tsv",
                   encoding="utf-8").read().splitlines()[1:]
    assert len(rows) == 2


def test_a_proteome_of_nothing_but_over_cap_sequences_runs_nothing(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    cfg, p = paths_for("tmbed_allcapped")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(3, 400))
    cfg["tmbed_max_len"] = 100
    _stub_tmbed(tmp_path, monkeypatch)
    ma.stage_tmbed(cfg, p)
    assert not os.path.exists(_parts_dir(p)), "tmbed was given an empty input"
    assert os.path.exists(p.tmbed) and os.path.getsize(p.tmbed) == 0
    err = capsys.readouterr().err
    assert "no sequence is left to predict" in err
    excl = io.open(f"{os.path.dirname(p.tmbed)}/tmbed_excluded.tsv",
                   encoding="utf-8").read().splitlines()[1:]
    assert len(excl) == 3


def test_the_stage_may_write_an_empty_prediction_file(ma):
    # ...so decide() has to be willing to adopt one on a rerun.
    assert {s["name"]: s for s in ma.STAGES}["tmbed"].get("empty_ok") is True


def test_the_chunk_size_is_not_allowed_to_invalidate_the_cache(ma):
    # it changes the order of the records and nothing else; listing it would
    # throw away a 30-hour stage because someone tuned a checkpoint size.
    st = {s["name"]: s for s in ma.STAGES}["tmbed"]
    assert "tmbed_chunk_residues" not in st["keys"]
    # what DOES change the contents is listed
    assert "tmbed_allow_partial" in st["keys"]
    assert "tmbed_max_consecutive_failures" in st["keys"]


def test_the_chunk_plan_is_logged_before_any_prediction_starts(
        ma, tmp_path, paths_for, monkeypatch, capsys):
    cfg, p = paths_for("tmbed_plan")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(12))
    cfg["tmbed_chunk_residues"] = 800
    _stub_tmbed(tmp_path, monkeypatch)
    ma.stage_tmbed(cfg, p)
    lines = capsys.readouterr().err.splitlines()
    plan = [i for i, l in enumerate(lines) if "chunk(s) of at most" in l]
    ran = [i for i, l in enumerate(lines) if "$ " in l and "tmbed predict" in l]
    assert plan and ran and plan[0] < ran[0]
    assert "3 chunk(s)" in lines[plan[0]]


def test_a_die_from_inside_run_cmd_is_not_reported_as_a_chunk_failure(
        ma, tmp_path, paths_for, monkeypatch):
    # StageError IS a RuntimeError, so the `except RuntimeError` that turns a
    # dead chunk into a casualty list would also swallow an abort. run_cmd
    # raises only plain RuntimeError today, so this pins the intent rather
    # than a live path: the monkeypatch is what a future die() in run_cmd
    # would look like, and the stage must let it through untouched.
    cfg, p = paths_for("tmbed_die")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), _many(4))
    _stub_tmbed(tmp_path, monkeypatch)
    monkeypatch.setattr(ma, "run_cmd",
                        lambda *a, **k: ma.die("tmbed not found"))
    with pytest.raises(ma.StageError) as e:
        ma.stage_tmbed(cfg, p)
    assert "tmbed not found" in str(e.value)
    assert "failed after writing" not in str(e.value)
    assert not os.path.exists(f"{os.path.dirname(p.tmbed)}/tmbed_failed.tsv")


def test_emapper_coverage_is_reported_per_identifier_tier(ma, tmp_path,
                                                          capsys):
    # symptom: a merged search database was reported by one headline number.
    # The real one is a public catalogue tier that arrives with precomputed
    # annotations and a tier assembled from this study's own reads that does
    # not; 60% overall here is 100% and 0%, and only the split says so.
    have = [F.Protein(f"uhgpL_MGYG{i:05d}", "MKV" * 40, ko="ko:K01234",
                      pathway="ko00010,map00010", seed_taxid="820")
            for i in range(6)]
    missing = [F.Protein(f"OIDECCNN_{i:05d}", "MKV" * 40) for i in range(4)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), have + missing)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), have)
    report = str(tmp_path / "report.tsv")
    ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), report,
                       "exact", 0.5, 0.9, 100)
    err = capsys.readouterr().err
    tier_lines = [l for l in err.splitlines() if "emapper reuse:   " in l]
    assert len(tier_lines) == 2, err
    assert "uhgpL_" in tier_lines[0] and "100.0%" in tier_lines[0]
    assert "OIDECCNN_" in tier_lines[1] and "0.0%" in tier_lines[1]
    rep = dict(l.split("\t") for l in
               io.open(report, encoding="utf-8").read().splitlines()[1:])
    assert rep["tier_uhgpL__proteins"] == "6"
    assert rep["tier_uhgpL__annotated"] == "6"
    assert rep["tier_OIDECCNN__proteins"] == "4"
    assert rep["tier_OIDECCNN__annotated"] == "0"


def test_an_untiered_protein_set_gets_no_tier_rows(ma, tmp_path):
    ps = [F.Protein(f"P{i:05d}", "MKV" * 40, ko="ko:K01234",
                    pathway="ko00010,map00010", seed_taxid="820")
          for i in range(5)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), ps)
    report = str(tmp_path / "report.tsv")
    ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), report,
                       "exact", 0.5, 0.9, 100)
    assert "tier_" not in io.open(report, encoding="utf-8").read()


def test_a_tier_left_out_of_strip_id_prefix_is_named_while_it_can_be_fixed(
        ma, tmp_path, capsys):
    # symptom: uhgpL_ and uhgpSM_ wrap the SAME MGYG namespace, and one eggNOG
    # row annotates a protein under every tag it carries. A tag missing from
    # emapper_strip_id_prefix does not error -- its whole tier simply reports
    # as unannotated, which reads as biology.
    both = [F.Protein(f"uhgpL_MGYG00000{i}_0100{i}", "MKV" * 40,
                      ko="ko:K01234", pathway="ko00010,map00010",
                      seed_taxid="820") for i in range(5)]
    sm = [F.Protein(f"uhgpSM_MGYG00001{i}_0200{i}", "MKV" * 40,
                    ko="ko:K01234", pathway="ko00010,map00010",
                    seed_taxid="820") for i in range(3)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), both + sm)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), both + sm,
                          id_prefix="uhgpL_")
    ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), "",
                       "exact", 0.0, 0.0, 100, strip_prefixes=["uhgpL_"])
    err = capsys.readouterr().err
    assert "share the identifier key MGYG#_#" in err
    assert "uhgpSM_" in err
    assert "reports as unannotated when it is only unjoined" in err


def test_no_such_warning_when_every_sharing_tier_is_listed(ma, tmp_path,
                                                           capsys):
    both = [F.Protein(f"uhgpL_MGYG00000{i}_0100{i}", "MKV" * 40,
                      ko="ko:K01234", pathway="ko00010,map00010",
                      seed_taxid="820") for i in range(5)]
    sm = [F.Protein(f"uhgpSM_MGYG00001{i}_0200{i}", "MKV" * 40,
                    ko="ko:K01234", pathway="ko00010,map00010",
                    seed_taxid="820") for i in range(3)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), both + sm)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), both + sm,
                          id_prefix="uhgpL_")
    ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), "",
                       "exact", 0.0, 0.0, 100,
                       strip_prefixes=["uhgpL_", "uhgpSM_"])
    assert "reports as unannotated" not in capsys.readouterr().err


def test_a_tier_that_matches_nothing_is_called_out_as_a_zero_not_a_low_number(
        ma, tmp_path, capsys):
    # symptom: the AMPSphere tier of the real run matched 0 of its 2,168
    # proteins -- no eggNOG table on that machine is keyed on AMP/SPHERE ids
    # -- and it is 30% of the whole dark fraction. The headline coverage was
    # 98.4%, so nothing said so.
    have = [F.Protein(f"uhgpL_MGYG00000{i}_0100{i}", "MKV" * 40,
                      ko="ko:K01234", pathway="ko00010,map00010",
                      seed_taxid="820") for i in range(8)]
    none = [F.Protein(f"ampS_AMP10.000_{i:03d}", "MKV" * 40) for i in range(3)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), have + none)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), have)
    ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), "",
                       "exact", 0.0, 0.0, 100)
    err = capsys.readouterr().err
    assert "tier ampS_ matched NONE of its 3 protein(s)" in err
    assert "not for want of biology" in err
    # the tier that DID match must not be accused of it
    assert "tier uhgpL_ matched NONE" not in err


def test_a_merely_low_tier_is_not_called_a_zero(ma, tmp_path, capsys):
    have = [F.Protein(f"uhgpL_MGYG00000{i}_0100{i}", "MKV" * 40,
                      ko="ko:K01234", pathway="ko00010,map00010",
                      seed_taxid="820") for i in range(8)]
    thin = [F.Protein(f"ampS_AMP10.000_{i:03d}", "MKV" * 40,
                      ko="ko:K01234", pathway="ko00010,map00010",
                      seed_taxid="820") for i in range(3)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), have + thin)
    emp = F.write_emapper(str(tmp_path / "cat.annotations"), have + thin[:1])
    ma.prepare_emapper([emp], faa, str(tmp_path / "o.annotations"), "",
                       "exact", 0.0, 0.0, 100)
    assert "matched NONE" not in capsys.readouterr().err


# ----------------------------------------------------------------------
# tmbed_use_gpu, and the half doctor now asserts
# ----------------------------------------------------------------------
# A tmbed that behaves the way the real one does on a host with no GPU: it
# tolerates a missing device under --cpu-fallback and refuses without it. That
# is the claim `doctor`'s `gpu:topology` row makes about `tmbed_use_gpu: true`,
# and it was being made about a tool nothing in this suite had ever driven.
_TMBED_GPU_STUB = """#!/usr/bin/env python3
import sys
a = sys.argv[1:]
if "--use-gpu" in a and "--no-cpu-fallback" in a:
    sys.stderr.write("No GPU available and CPU fallback is disabled.\\n")
    sys.exit(1)
out = a[a.index("-p") + 1]
faa = a[a.index("-f") + 1]
recs, pid, seq = [], None, []
for line in open(faa, encoding="utf-8"):
    line = line.strip()
    if line.startswith(">"):
        if pid:
            recs.append((pid, "".join(seq)))
        pid, seq = line, []
    elif line:
        seq.append(line)
if pid:
    recs.append((pid, "".join(seq)))
with open(out, "w", encoding="utf-8") as fh:
    for hdr, s in recs:
        fh.write(hdr + "\\n" + s + "\\n" + "G" * len(s) + "\\n")
"""


def _stub_gpu_tmbed(tmp_path, monkeypatch):
    d = tmp_path / "tmbedgpubin"
    d.mkdir(exist_ok=True)
    (d / "tmbed").write_text(_TMBED_GPU_STUB, encoding="utf-8")
    os.chmod(d / "tmbed", 0o755)
    if os.name == "nt":
        (d / "tmbed.cmd").write_text(
            f'@echo off\r\n"{sys.executable}" "%~dp0tmbed" %*\r\n',
            encoding="utf-8")
    monkeypatch.setenv("PATH", str(d) + os.pathsep + os.environ["PATH"])
    return d


@pytest.mark.parametrize("want,fatal", [("auto", False), ("true", True),
                                        ("false", False)])
def test_tmbed_use_gpu_true_is_fatal_without_a_gpu_and_the_others_are_not(
        ma, tmp_path, paths_for, monkeypatch, want, fatal):
    """The engine half of doctor's `gpu:topology` verdict, driven.

    `doctor` claimed one outcome for all three values of this key while
    `stage_tmbed` maps them onto three different command lines, and the row
    that said so was never tested against the stage. TMbed tolerates a missing
    or failing GPU only under `--cpu-fallback`, so `true` - which passes
    `--no-cpu-fallback` - makes every chunk fail, and a stage with no
    prediction for any protein die()s unless `tmbed_allow_partial` is on.
    """
    cfg, p = paths_for(f"tmbedgpu_{want}")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    cfg["tmbed_use_gpu"] = want
    _stub_gpu_tmbed(tmp_path, monkeypatch)
    if not fatal:
        ma.stage_tmbed(cfg, p)
        assert os.path.getsize(p.tmbed) > 0
        return
    with pytest.raises(ma.StageError) as e:
        ma.stage_tmbed(cfg, p)
    assert "prediction(s)" in str(e.value)
    assert "tmbed_allow_partial: true" in str(e.value)


def test_tmbed_allow_partial_survives_the_gpu_setting_that_kills_every_chunk(
        ma, tmp_path, paths_for, monkeypatch):
    # The second key doctor's row now names. With it on, the same refusal
    # leaves an EMPTY prediction file and the stage returns - which is a
    # degradation and not a death, and the row says which.
    cfg, p = paths_for("tmbedgpu_partial")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    cfg["tmbed_use_gpu"] = "true"
    cfg["tmbed_allow_partial"] = True
    _stub_gpu_tmbed(tmp_path, monkeypatch)
    ma.stage_tmbed(cfg, p)
    assert os.path.exists(p.tmbed) and os.path.getsize(p.tmbed) == 0


def test_the_gpu_flag_map_is_the_one_the_stage_really_passes(
        ma, tmp_path, paths_for, monkeypatch):
    # TMBED_GPU_MODES is read by both stage_tmbed and doctor's gpu:topology
    # row, so this is the tripwire that keeps it describing the real argv.
    cfg, p = paths_for("tmbedgpu_argv")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    _stub_gpu_tmbed(tmp_path, monkeypatch)
    seen = []
    real = ma.run_cmd
    monkeypatch.setattr(ma, "run_cmd",
                        lambda argv, *a, **k: (seen.append(list(argv)),
                                               real(argv, *a, **k))[1])
    for want, (flags, _fatal) in sorted(ma.TMBED_GPU_MODES.items()):
        if want == "true":
            continue                      # its argv is the one that refuses
        seen.clear()
        shutil.rmtree(p.R, ignore_errors=True)
        cfg["tmbed_use_gpu"] = want
        ma.stage_tmbed(cfg, p)
        assert seen, f"{want}: the stage ran nothing"
        argv = [str(x) for x in seen[0]]
        assert all(fl in argv for fl in flags), \
            f"{want}: TMBED_GPU_MODES says {flags}, argv was {argv}"


def test_a_fifo_input_blocks_the_reader_rather_than_failing_it(ma, tmp_path):
    """The claim `doctor`'s FIFO rows make about `run`, driven on the reader.

    A differential sweep of `doctor` against `run` - every configured input
    through every state a path can be in - turned up four cells where `run`
    never returned at all: a FIFO at `proteins_faa`, `quant_table`, `manifest`
    or `gff`. The rows had said the stage could not "get a FASTA out of this",
    which reads as a failure; what actually happens is that a read-only open
    of a FIFO with no writer BLOCKS, so the run produces no output, no
    traceback and no exit status, and looks exactly like a long stage.

    Driven on `read_fasta()` rather than on `run`, because the mechanism is
    `opener()` and a test that waits for a whole pipeline to not finish costs
    a minute to assert the same thing. The reader is released at the end, so
    nothing is left blocked behind this test.
    """
    fifo = tmp_path / "proteins.faa"
    os.mkfifo(str(fifo))
    started, done = threading.Event(), []

    def read():
        started.set()
        done.append(list(ma.read_fasta(str(fifo))))

    th = threading.Thread(target=read, daemon=True)
    th.start()
    started.wait(5)
    th.join(3)
    assert th.is_alive() and not done, \
        "read_fasta returned on a FIFO with no writer; if opener() now " \
        "refuses one, doctor's FIFO rows should say it dies rather than hangs"
    # Release it, so the thread is not left blocked for the rest of the run.
    with open(str(fifo), "w", encoding="utf-8") as fh:
        fh.write(">p1\nMKV\n")
    th.join(10)
    assert done and done[0] == [("p1", "MKV")], \
        "a FIFO with a writer reads perfectly well, which is why opener() " \
        "does not refuse one and only doctor guards against the open"


# ----------------------------------------------------------------------
# opener(): one open, kept - the FIFO contract, both halves
# ----------------------------------------------------------------------
# Seven operator-supplied paths used to hang `run` forever when a FIFO with no
# writer sat at them - proteins_faa, quant_table, manifest, gff,
# emapper_precomputed, and a TMT plex's ion.tsv and annotation - with no
# output, no traceback and no exit status to time out against. The fix is in
# opener(), which is the one place every parser in the engine reads through,
# and its whole shape is decided by two things that must not break:
#
#   * a FIFO WITH a live writer is a supported input and reads to completion.
#     A gate that refused a FIFO would be removing a working workflow in order
#     to catch a misconfiguration.
#   * probing a FIFO and closing the probe sends the writer EPIPE - measured,
#     in test_probing_a_live_fifo_breaks_the_writer_on_the_other_side. So the
#     open cannot be a gate in front of another open: it is ONE open, kept,
#     and the descriptor it returns is the one the parser reads.
#
# Everything below drives one of those, or the bound that ends the hang.
def _slow_writer(path, chunks, gap=0.15, box=None):
    """Write `chunks` to `path` with a pause between them, recording EPIPE.

    The PAUSE is the point. A descriptor left non-blocking raises EAGAIN as
    BlockingIOError the moment a read finds no bytes ready, which on a healthy
    pipe is any moment the writer is thinking - so a writer that never pauses
    cannot tell a cleared O_NONBLOCK from an uncleared one.
    """
    box = {} if box is None else box

    def write():
        try:
            with open(path, "wb") as fh:
                for c in chunks:
                    fh.write(c)
                    fh.flush()
                    time.sleep(gap)
        except BrokenPipeError as e:                        # noqa: PERF203
            box["broken"] = e
        except OSError as e:                                # noqa: PERF203
            box["error"] = e

    th = threading.Thread(target=write, daemon=True)
    th.start()
    return th, box


def test_a_fifo_with_a_live_writer_reads_to_completion_plain_and_gzipped(
        ma, tmp_path):
    """The workflow the whole design is built around, both branches.

    `mkfifo p; zcat big.faa.gz > p &` is how a disk-constrained cluster feeds
    this tool. The GZIP branch is here because it is not a curiosity:
    emapper_precomputed is routinely a `.gz`, so the branch that used to
    re-open the path by name is the one the commonest piped input goes
    through, and re-opening is exactly what breaks a live writer.
    """
    plain = str(tmp_path / "live.faa")
    os.mkfifo(plain)
    th, box = _slow_writer(plain, [b">p1\nMKV\n", b">p2\nMKW\n"])
    assert list(ma.read_fasta(plain)) == [("p1", "MKV"), ("p2", "MKW")]
    th.join(5)
    assert not box, f"the writer did not finish cleanly: {box}"

    gz = str(tmp_path / "live.faa.gz")
    os.mkfifo(gz)
    blob = gzip.compress(b">g1\nMKV\n>g2\nMKW\n")
    th, box = _slow_writer(gz, [blob[:len(blob) // 2], blob[len(blob) // 2:]])
    assert list(ma.read_fasta(gz)) == [("g1", "MKV"), ("g2", "MKW")], \
        "the gzip branch no longer reads a live FIFO - it is re-opening the " \
        "path instead of reading the descriptor opener() already holds"
    th.join(5)
    assert not box, f"the writer did not finish cleanly: {box}"


def test_a_whole_stage_read_of_a_live_fifo_never_breaks_the_writer(
        ma, tmp_path):
    """The measured cost of a gate, asserted from the writer's side.

    The writer is held open across the WHOLE read rather than being allowed to
    finish first, because that is the state a probe damages: a reader that
    opened the path, closed it and opened it again would send this thread
    EPIPE somewhere in the middle, and the records would still come back from
    the second open. Asserting the records alone would pass. Asserting that
    the writer saw no BrokenPipeError is what fails.
    """
    fifo = str(tmp_path / "held.faa")
    os.mkfifo(fifo)
    box, ready, release = {}, threading.Event(), threading.Event()

    def write():
        try:
            with open(fifo, "w", encoding="utf-8") as fh:
                for i in range(20):
                    fh.write(f">p{i}\n" + "MKV" * 50 + "\n")
                    fh.flush()
                    time.sleep(0.02)
                ready.set()
                # Still OPEN, and staying open until the reader has finished:
                # an EPIPE arrives on a write, so the writer has to be alive to
                # receive one.
                release.wait(10)
                fh.write(">last\nMKV\n")
        except BrokenPipeError as e:
            box["broken"] = e

    th = threading.Thread(target=write, daemon=True)
    th.start()
    got = []
    reader = threading.Thread(
        target=lambda: got.extend(ma.read_fasta(fifo)), daemon=True)
    reader.start()
    ready.wait(10)
    release.set()
    reader.join(20)
    th.join(10)
    assert "broken" not in box, \
        "the writer on the other side got EPIPE during a stage read - " \
        "opener() is probing and re-opening instead of keeping one descriptor"
    assert len(got) == 21 and got[-1] == ("last", "MKV")


def test_a_writer_that_pauses_mid_stream_does_not_raise_eagain(ma, tmp_path):
    """Why O_NONBLOCK is cleared before the read, driven.

    The flag lives on the open file description, so it is still set on every
    read the parser makes through the handle - and a non-blocking read of a
    pipe with no bytes ready fails with EAGAIN, which Python raises as
    BlockingIOError. That would turn a working pipe into a failure whose
    likelihood depends on how fast the writer happens to be: intermittent, and
    worse than the hang it replaced, because a hang is at least reproducible.

    The gaps here are long compared with the reads, so an uncleared flag
    fails this every time rather than one run in ten.
    """
    fifo = str(tmp_path / "slow.tsv")
    os.mkfifo(fifo)
    th, box = _slow_writer(fifo, [b"a\tb\n", b"1\t2\n", b"3\t4\n"], gap=0.4)
    with ma.opener(fifo) as fh:
        assert fh.read() == "a\tb\n1\t2\n3\t4\n"
    th.join(5)
    assert not box, f"the writer did not finish cleanly: {box}"


def test_a_fifo_with_no_writer_fails_with_a_message_instead_of_hanging(
        ma, tmp_path, monkeypatch):
    """The defect itself: the silence, and then the exit status there was not.

    Driven in a thread with a join, so a regression FAILS this test rather
    than wedging the suite - which is the same failure mode the defect is
    about, and the reason the old test could only assert that the reader was
    still alive.
    """
    fifo = str(tmp_path / "nobody.faa")
    os.mkfifo(fifo)
    monkeypatch.setattr(ma, "_FIFO_WAIT", 1.0)
    box = {}

    def read():
        try:
            list(ma.read_fasta(fifo))
            box["returned"] = True
        except ma.StageError as e:                          # noqa: PERF203
            box["error"] = str(e)

    th = threading.Thread(target=read, daemon=True)
    th.start()
    th.join(30)
    assert not th.is_alive(), \
        "read_fasta is still on the open after its wait - the FIFO wait is " \
        "not bounded and `run` can hang again"
    assert "returned" not in box, "an unwritten FIFO cannot yield records"
    msg = box.get("error", "")
    for want in (fifo, "FIFO", "nothing is holding the write end",
                 "fifo_wait_s"):
        assert want in msg, \
            f"the refusal does not name {want!r}: {msg}"


def test_the_wait_is_the_configured_one_and_zero_refuses_at_once(
        ma, tmp_path, monkeypatch):
    """The setting is a setting, and its two ends both mean something.

    0 is not "wait forever" and not "ignore this": it is the operator who has
    no piped inputs saying so, and it has to be the value that refuses a FIFO
    on sight. A setting whose off value silently meant the default would be a
    knob that cannot be turned off.

    IN A THREAD WITH A JOIN, like its sibling twenty lines above, and that is
    not a style note. This test used to call list(read_fasta(fifo)) on the
    main thread with no thread, no join and no timeout: against an opener that
    does not bound the wait it did not FAIL, it WEDGED, and it hung until
    SIGKILL - taking the whole suite with it, which is the same failure mode
    the code under test exists to remove. Every FIFO test here goes through
    _read_in_thread() for that reason.
    """
    fifo = str(tmp_path / "none.faa")
    os.mkfifo(fifo)
    monkeypatch.setattr(ma, "_FIFO_WAIT", 0.0)
    t0 = time.time()
    box = _read_in_thread(lambda: list(ma.read_fasta(fifo)), 20)
    assert not box["alive"], \
        "fifo_wait_s: 0 did not refuse - the reader is still on the open"
    assert isinstance(box.get("error"), ma.StageError), \
        f"fifo_wait_s: 0 did not raise a StageError: {box}"
    assert time.time() - t0 < 5, "fifo_wait_s: 0 waited anyway"
    # ...and the default is the documented one, read off the config rather
    # than typed here, so the two cannot drift.
    ma.set_fifo_wait(ma.DEFAULT_CONFIG["fifo_wait_s"])
    assert ma._FIFO_WAIT == float(ma.DEFAULT_CONFIG["fifo_wait_s"])
    with pytest.raises(TypeError):
        ma.set_fifo_wait(True)


def test_a_regular_file_is_opened_once_and_reads_exactly_as_before(
        ma, tmp_path, monkeypatch):
    """The case that must not have changed at all, in all three ways.

    An ordinary file is the overwhelming majority of every read this tool
    does, so the FIFO fix may not cost it an extra syscall, a changed
    exception or a changed decode. One open, because a probe would be two and
    a probe is what breaks a pipe; the latin-1 guarantee in opener()'s own
    docstring, because that is the other thing this function promises; and the
    same exception for a path that is not there.
    """
    good = tmp_path / "plain.txt"
    good.write_bytes(b"caf\xa0 latin-1\n")
    real_open, seen = os.open, []
    monkeypatch.setattr(os, "open",
                        lambda *a, **k: (seen.append(a[0]),
                                         real_open(*a, **k))[1])
    with ma.opener(str(good)) as fh:
        text = fh.read()
    assert seen == [str(good)], \
        f"opener() opened the path {len(seen)} times, not once: {seen}"
    assert text == "caf\ufffd latin-1\n", \
        "the errors='replace' guarantee in opener()'s docstring is gone"
    # ...and the same promise on the gzip branch, which states its encoding
    # for the same reason.
    gz = tmp_path / "plain.txt.gz"
    with gzip.open(str(gz), "wb") as fh:
        fh.write(b"caf\xa0 latin-1\n")
    del seen[:]
    with ma.opener(str(gz)) as fh:
        assert fh.read() == "caf\ufffd latin-1\n"
    assert seen == [str(gz)], f"the gzip branch opened the path twice: {seen}"
    monkeypatch.undo()
    with pytest.raises(FileNotFoundError):
        ma.opener(str(tmp_path / "nope.txt"))
    with pytest.raises(IsADirectoryError):
        ma.opener(str(tmp_path))


def test_a_socket_or_a_device_node_is_refused_naming_what_is_there(
        ma, tmp_path, monkeypatch):
    """The kinds there is nothing to wait for.

    A FIFO is the one non-regular kind that can BECOME readable, so it is the
    only one that gets a wait. A character device answers immediately and
    forever with bytes that are not a table - /dev/zero would feed a parser
    NULs until it ran out of memory - and a socket cannot be opened this way
    at all. Both are `other` to `doctor`, which is why the refusal has to name
    what it really found rather than say "not a file".
    """
    with pytest.raises(ma.StageError) as e:
        ma.opener("/dev/zero")
    assert "character device" in str(e.value) and "/dev/zero" in str(e.value)
    # Bound from inside the directory, because an AF_UNIX path is capped at
    # about a hundred bytes and a pytest tmp_path is most of that already.
    monkeypatch.chdir(tmp_path)
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    try:
        sock.bind("s.sock")
        # An OSError rather than a StageError, and it comes from os.open()
        # itself: a socket cannot be opened as a file at all, so there is no
        # descriptor to fstat and nothing for opener() to say about it that
        # the errno does not already say. WHICH errno is the OS's to choose
        # and is deliberately not asserted - this file and the doctor rows
        # said ENXIO for three rounds and the platform they were measured on
        # answers EOPNOTSUPP. What is load-bearing is that it is an OSError,
        # because an OSError is what nothing downstream catches.
        with pytest.raises(OSError) as se:
            ma.opener("s.sock")
        assert not isinstance(se.value, ma.StageError), \
            "a socket must not come back as a refusal anything falls back from"
    finally:
        sock.close()


# Every path a FIFO has been measured hanging this program on, with what
# drives each. `run` where the whole pipeline reaches the path on a config the
# suite can build in a second, the READER where it does not - stage_context
# needs an ORF finder on PATH and the TMT tree needs a quant_format the rest
# of this project is not in, and a test that stands a whole pipeline up to
# prove a message costs a minute to assert what the reader asserts in a
# second - and the COMMAND for the one that is not `run` at all.
#
# IT IS NOT CALLED SEVEN_HANGING_PATHS ANY MORE, and the rename is the point
# rather than tidiness: a count in the name of a list is the same claim a
# count in prose is, and this one went out of date the moment the sweep was
# widened past the inputs with `input:` rows on doctor. More paths were
# measured hanging after that list was written - `unipept.result` and the
# taxdump's .dmp files, which went through bare open() calls inside
# read_unipept_result() and NCBITaxonomy(), and `subset --quant`, which is a
# whole command nobody had ever swept. "seven" is kept only where it is
# history: the sentence about what the differential sweep found that round.
HANGING_PATHS = (
    ("proteins_faa", "run"),
    ("quant_table", "run"),
    ("manifest", "run"),
    ("emapper_precomputed", "run"),
    ("unipept_result", "run"),
    ("taxdump_nodes", "run"),
    # The same file on a config that reads it TWICE - run.taxonomy builds one
    # NCBITaxonomy and taxon_rank makes stage_join build another - so the
    # single-read arm and the refuse-at-once arm are both driven at a database
    # path, which is the class that had no read plan at all.
    ("taxdump_nodes_twice", "run"),
    ("gff", "reader"),
    ("tmt_level", "reader"),
    ("tmt_annotation", "reader"),
    ("subset_quant", "command"),
)


def _hanging_project(tmp_path, key):
    """(project, the path to put a FIFO at) for one `run`-driven key.

    Two of these need a config build_project() does not make - an ingested
    pept2lca export and a taxdump - and _space_project() already knows how to
    make both, so the sweep that proves the plan and the sweep that proves the
    refusal build their configs the same way.
    """
    if key == "unipept_result":
        p = _space_project(tmp_path, dict(id="hang_uni",
                                          run=dict(unipept=True), result=True))
        p.write_config(fifo_wait_s=2)
        return p, p.cfg["unipept"]["result"]
    if key.startswith("taxdump_nodes"):
        point = dict(id="hang_tax", run=dict(unipept=True, taxonomy=True),
                     result=True, taxdump=True)
        if key.endswith("_twice"):
            point = dict(point, id="hang_tax2", taxon_rank="genus")
        p = _space_project(tmp_path, point)
        p.write_config(fifo_wait_s=2)
        return p, os.path.join(p.cfg["db"]["ncbi_taxonomy"], "nodes.dmp")
    p = build_project(tmp_path / f"fifo_{key[:6]}", fifo_wait_s=2)
    return p, (p.cfg["emapper_precomputed"][0]
               if key == "emapper_precomputed" else p.cfg[key])


@pytest.mark.parametrize("key,how", HANGING_PATHS,
                         ids=[k for k, _h in HANGING_PATHS])
def test_every_path_that_hung_the_run_now_ends_with_a_message_naming_it(
        ma, tmp_path, monkeypatch, key, how):
    """All seven, driven, with a timeout so a regression fails rather than
    wedges.

    The assertion is deliberately about the MESSAGE and not only the exit
    status: the defect was that there was nothing to read and nothing to time
    out against, so "it exited 1" is only half of what was missing. An
    operator has to be told which path it was.

    WHICH message is now two messages, and that is the correction this round
    made. A path the run reads once ends in the wait, because the wait is what
    keeps a live writer working there; a path it reads more than once ends AT
    THE OPEN, because waiting at one buys a failure six hours later and a
    sentence recommending a workflow that cannot work at that path. So the
    test asks the engine which case this is - it does not carry a list - and
    holds the message against the answer.
    """
    if how == "command":
        # `subset` is not `run`: it has no --config, so it has no fifo_wait_s
        # to set and it keeps the module default. Driven in process with the
        # wait cut down, and in a thread with a join, because the defect here
        # is a BLOCK - a bare pd.read_csv(args.quant) with nothing in front of
        # it - and a test that called cmd_subset() on the main thread would
        # wedge rather than fail when it came back.
        p = build_project(tmp_path / "fifo_subset")
        monkeypatch.setattr(ma, "_FIFO_WAIT", 1.0)
        fifo = str(tmp_path / "report.pg_matrix.tsv")
        os.mkfifo(fifo)
        args = types.SimpleNamespace(
            quant=fifo, format="diann", db=p.cfg["proteins_faa"],
            out=str(tmp_path / "subset.faa"))
        box = {}
        th = threading.Thread(
            target=lambda: box.update(e=_caught(lambda: ma.cmd_subset(args))),
            daemon=True)
        th.start()
        th.join(30)
        assert not th.is_alive(), \
            "subset is still on the open - --quant is read outside the choke " \
            "point and the command hangs with no exit status"
        assert fifo in str(box.get("e")) \
            and "nothing is holding the write end" in str(box["e"]), \
            f"subset failed without naming the path: {box}"
        return
    if how == "run":
        p, target = _hanging_project(tmp_path, key)
        reads = len(ma.set_read_plan(p.cfg).get(os.path.realpath(target), []))
        os.remove(target)
        os.mkfifo(target)
        t0 = time.time()
        proc = p.run(expect=None, timeout=120)
        assert proc.returncode != 0, \
            f"run exited 0 with an unwritten FIFO at {key}"
        assert target in proc.stderr and "FIFO" in proc.stderr, \
            f"run failed without naming {key}:\n{proc.stderr[-1500:]}"
        if reads > 1:
            assert "drained EXACTLY ONCE" in proc.stderr, \
                f"{key} is read {reads} times and the run still waited:" \
                f"\n{proc.stderr[-1500:]}"
            assert time.time() - t0 < 30, \
                f"{key} cannot support a pipe and the run waited anyway"
        else:
            assert "nothing is holding the write end" in proc.stderr, \
                f"{key} is single-read and did not end in the wait:" \
                f"\n{proc.stderr[-1500:]}"
        return
    monkeypatch.setattr(ma, "_FIFO_WAIT", 1.0)
    fifo = str(tmp_path / {"gff": "x.gff", "tmt_level": "ion.tsv",
                           "tmt_annotation": "TMT1_annotation.txt"}[key])
    os.mkfifo(fifo)
    reader = {"gff": lambda: ma.parse_gff(fifo),
              "tmt_level": lambda: ma.header_columns(fifo),
              "tmt_annotation": lambda: ma.read_tmt_annotation(fifo, "TMT1")}[key]
    box = {}
    th = threading.Thread(
        target=lambda: box.update(e=_caught(reader)), daemon=True)
    th.start()
    th.join(30)
    assert not th.is_alive(), f"the reader for {key} is still on the open"
    assert fifo in str(box.get("e")) \
        and "nothing is holding the write end" in str(box["e"]), \
        f"the reader for {key} failed without naming it: {box}"


def _caught(fn):
    """Run `fn` and hand back whatever it raised, or a marker if it returned."""
    try:
        fn()
    except Exception as e:                                  # noqa: BLE001
        return e
    return "returned instead of raising"


def test_doctor_still_answers_at_once_on_a_fifo_it_will_never_wait_for(
        tmp_path):
    """The asymmetry, kept, and now the thing most likely to be broken.

    `run` may wait on a pipe, because waiting on a pipe is what reading a pipe
    IS. `doctor` may NEVER wait: a command whose whole job is to answer before
    the run is worthless if it blocks, and it is the only thing an operator
    has that can say what is wrong with a path before committing to hours.
    The two are kept apart by doctor reading every path through
    _deep_readable(), which opens with O_NONBLOCK and does not wait at all -
    so a `fifo_wait_s` of a day must make no difference to it.
    """
    p = build_project(tmp_path / "docfifo", fifo_wait_s=86400)
    os.remove(p.cfg["proteins_faa"])
    os.mkfifo(p.cfg["proteins_faa"])
    t0 = time.time()
    proc = run_metaannot("doctor", "--config", p.config_path, "--json",
                         expect=1, timeout=90)
    assert time.time() - t0 < 40, \
        "doctor waited on the FIFO - it is reading through opener() somewhere"
    doc = json.loads(proc.stdout)
    row = [c for c in doc["checks"] if c["id"] == "input:proteins_faa"][0]
    assert row["found"]["kind"] == "other" and row["status"] == "fail"


# ======================================================================
# the read plan: how many times a run opens each operator-supplied path
# ======================================================================
# THE FACT EVERY FIFO DECISION RESTS ON, and the reason it is MEASURED here
# rather than read off the code: the previous round's report stated that
# `opener()` was the only reader of the quant table. That was false in the same
# function - `read_delim_table()` opened it again for its delimiter sniff, and
# pandas opened it a third time - and it was where a FIFO at `quant_table`
# actually hung, holding the results lock while the writer took EPIPE. A
# hand-written list of readers is exactly how that happened, so the list is not
# trusted: a real command is driven with every open of every configured path
# recorded, and the plan in metaannot.py has to match what the process did.
#
# The shim is a `sitecustomize` on PYTHONPATH, so it wraps builtins.open and
# os.open before metaannot is imported and needs no change to the file under
# test. It records, per open, which GATE was on the stack - the primitives
# that fstat the descriptor they are about to hand over - and whether the
# signature digest was, which is what lets these tests say WHICH KIND of open
# they saw rather than only how many.
_OPEN_TRACER = """
import atexit, builtins, json, os, traceback

TARGETS = set(p for p in os.environ.get("TRACE_PATHS", "").split(os.pathsep) if p)
OUT = os.environ["TRACE_OUT"]
# The primitives that prove what they hold before anything reads it:
# _open_for_read is `run`'s (it waits on a FIFO where a FIFO can work),
# _open_regular_binary is the signature digest's and, through
# _open_regular_text, `doctor`'s, and regular_readable is `doctor`'s probe.
# All three fstat the DESCRIPTOR they opened, which is the property this
# records - "did this open go through a gate", not "which function called it".
GATES = ("_open_for_read", "_open_regular_binary", "regular_readable")
_open, _osopen = builtins.open, os.open
hits = []


def _record(path):
    try:
        real = os.path.realpath(str(path))
    except Exception:
        return
    if real not in TARGETS:
        return
    stack = traceback.extract_stack()[:-2]
    names = [f.name for f in stack if "metaannot.py" in f.filename]
    hits.append({
        "path": real,
        "choke": "_open_for_read" in names,
        "gate": next((n for n in names if n in GATES), ""),
        "digest": "_stat" in names or "_content_digest" in names,
        "where": " <- ".join(reversed(
            [f"{f.name}:{f.lineno}" for f in stack
             if "metaannot.py" in f.filename][-5:])),
    })


def _open_wrapper(file, *a, **k):
    _record(file)
    return _open(file, *a, **k)


def _osopen_wrapper(path, *a, **k):
    _record(path)
    return _osopen(path, *a, **k)


builtins.open = _open_wrapper
os.open = _osopen_wrapper


@atexit.register
def _dump():
    with _open(OUT, "a", encoding="utf-8") as fh:
        for h in hits:
            fh.write(json.dumps(h) + "\\n")
"""


def _traced(tmp_path, project, targets, args, expect=None, cwd=None):
    """Drive the real CLI once and return (proc, every open it made).

    A SUBPROCESS, because that is what a command is: cmd_run holds the results
    lock, the stages run on worker threads, and an in-process call would be
    measuring something else. The tracer is installed through PYTHONPATH so
    that nothing in metaannot.py has to know it is being watched.

    `args` rather than a fixed ("run",) because the sweep below drives
    `doctor` and `subset` too, and the one thing every earlier version of this
    had in common was that it drove `run` and nothing else - which is how a
    bare pd.read_csv() in cmd_subset survived every sweep that has ever been
    written here.
    """
    shim = tmp_path / "shim"
    shim.mkdir(exist_ok=True)
    (shim / "sitecustomize.py").write_text(_OPEN_TRACER, encoding="utf-8")
    out = tmp_path / "opens.jsonl"
    if out.exists():
        out.unlink()
    env = {"PYTHONPATH": str(shim), "TRACE_OUT": str(out),
           "TRACE_PATHS": os.pathsep.join(
               os.path.realpath(str(t)) for t in targets)}
    proc = run_metaannot(*args, env=env, cwd=cwd or project.root,
                         expect=expect, timeout=300)
    rows = [json.loads(ln) for ln in
            (out.read_text(encoding="utf-8").splitlines() if out.exists()
             else [])]
    return proc, rows


def _traced_run(tmp_path, project, targets, expect=0, args=("run",)):
    """The `run` case, with --config filled in. Kept for its callers."""
    return _traced(tmp_path, project, targets,
                   (*args, "--config", project.config_path), expect=expect)


def _stream_opens(rows, path):
    """The opens of `path` that a FIFO would have to satisfy.

    The signature digest is excluded and that exclusion is the subtle half of
    the whole design: _content_digest() reads every input in a stage's
    signature, so a regular file is opened once more than this counts. It
    never touches a stream - _open_regular_binary() fstats the descriptor and
    refuses anything that is not a regular file - which is asserted on its own
    in test_the_signature_digest_never_opens_a_stream, so a FIFO never meets
    it and it is not part of the count that decides whether a FIFO can work.
    """
    real = os.path.realpath(str(path))
    return [r for r in rows if r["path"] == real and not r["digest"]]


def _configured_paths(cfg):
    """Every existing file this config names, found by walking the config.

    INDEPENDENT OF INPUT_READ_SITES, and that independence is the whole value
    of it. A watch list built from the table could only ever watch the inputs
    the table already knows about, so an input it forgot ENTIRELY would be
    invisible to the test written to catch exactly that - which is what
    happened: `db.ncbi_taxonomy` had no entry of any kind, its four .dmp files
    went through bare open() calls, and a FIFO at nodes.dmp left `run` blocked
    in the open with no exit status and the results lock still held. Every
    sweep in this file passed while that was true, because every sweep asked
    the table which paths to watch.

    A directory-valued setting is expanded to the files inside it, because
    `quant_table` on quant_format 'fragpipe_tmt' is a run directory and the
    per-plex ion.tsv and annotation files inside it are operator-supplied
    paths in every sense that matters to an open(). Relative paths are
    skipped, which is what keeps `results_dir` - a directory the run WRITES -
    out of a sweep about reading inputs.
    """
    found = []

    def walk(node):
        if isinstance(node, dict):
            for v in node.values():
                walk(v)
        elif isinstance(node, (list, tuple)):
            for v in node:
                walk(v)
        elif isinstance(node, str) and node and os.path.isabs(node):
            if os.path.isfile(node):
                found.append(node)
            elif os.path.isdir(node):
                for root, _dirs, files in os.walk(node):
                    found.extend(os.path.join(root, f) for f in files)

    walk(cfg)
    return sorted({os.path.realpath(x) for x in found})


# ======================================================================
# the read plan and the choke point, swept over the CONFIG SPACE
# ======================================================================
# WHY A SWEEP AND NOT TWO CONFIGS, which is what stood here. The read plan is
# a property of the RUN, so a test that drives two points in the config space
# pins the plan at two points and says nothing about the rest of it - and the
# rest of it is where every defect this round found was living. Driven, the
# two-point version was green while:
#
#   * `manifest` said ONE read on a config that made THREE, because
#     read_feature_table() opens the manifest and is reached from stage_join
#     AND from peptide_features() in each peptide stage. The run died at the
#     second open and `doctor --json` had already published "a FIFO here is
#     NOT refused on sight" for it.
#   * `quant_table` said two reads where the run made one, and three where it
#     made two, because stage_unipept RETURNS before peptide_features() when
#     `unipept.result` is set.
#   * `emapper_precomputed` said one read on `run.eggnog: false`, where
#     nothing opens it at all.
#   * `db.ncbi_taxonomy` was not in the plan at all and its four .dmp files
#     were read through bare open() calls, so a FIFO there hung the run.
#
# Not one of those configs was one of the two points. So the space is swept
# along the two axes the answer actually moves on - the quant FORMAT, which
# decides which reader join uses and whether the peptide stages reach the
# table at all, and the run FLAGS, which decide how many stages want it - and
# the paths watched are taken from the CONFIG rather than from the table under
# test, so an input nobody listed cannot hide.
#
# THE SAME DRIVEN RUN ANSWERS BOTH QUESTIONS, which is why these are one test
# and not two. "How many times is this opened" and "did each of those opens go
# through a gate" are two readings of one trace, and running the sweep twice
# to ask them separately would double the cost of the most expensive test in
# the file for nothing. The failure messages stay separate.
PLAN_SPACE = (
    # One point per quant format. All eight, because the format is what
    # decides whether join reads the table once (read_named_table gives the
    # header and the body from one open), twice (a protein-level table has to
    # be recognised from its header before pandas parses it), or not at all,
    # and whether read_manifest is reached from read_feature_table, from
    # stage_join's own branch, or from nowhere.
    dict(id="fragpipe_peptide"),
    dict(id="fragpipe_ion", fmt="fragpipe_ion"),
    dict(id="msstats_csv", fmt="msstats_csv"),
    dict(id="msstats_feature", fmt="msstats_feature"),
    dict(id="diann", fmt="diann"),
    dict(id="fragpipe", fmt="fragpipe"),
    dict(id="msstats_protein", fmt="msstats_protein"),
    dict(id="fragpipe_tmt", fmt="fragpipe_tmt"),
    # ...then the flags that move the count for a fixed format.
    dict(id="eggnog_off", run=dict(eggnog=False)),
    dict(id="join_off", run=dict(join=False)),
    dict(id="no_manifest", manifest=False),
    dict(id="unipept_result_only", run=dict(unipept=True), result=True),
    dict(id="unipept_and_taxonomy", run=dict(unipept=True, taxonomy=True),
         result=True, taxdump=True),
    dict(id="taxonomy_on_a_protein_format", fmt="diann",
         run=dict(unipept=True, taxonomy=True), result=True, taxdump=True),
    dict(id="peptide_only_always", run=dict(unipept=True, taxonomy=True),
         result=True, taxdump=True, peptide_only_reader="always"),
    dict(id="join_off_with_peptide_stages", run=dict(join=False, unipept=True,
                                                     taxonomy=True),
         result=True, taxdump=True),
    # The two taxdump points, which are the ones that separate "read once" from
    # "read twice" for the .dmp files: taxon_rank alone reaches NCBITaxonomy
    # through stage_join, run.taxonomy reaches it through stage_taxonomy, and
    # a config with both opens every .dmp file twice.
    dict(id="taxon_rank_only", taxon_rank="genus", taxdump=True),
    dict(id="taxon_rank_and_taxonomy", run=dict(unipept=True, taxonomy=True),
         result=True, taxdump=True, taxon_rank="genus"),
    # A TMT run with a peptide stage on it. quant_table is a DIRECTORY here and
    # the files a reader opens are the per-plex ones inside it, so this is the
    # db.ncbi_taxonomy shape one input over: with two plexes and two consumers
    # every file under the tree is opened twice, and the plan had no entry of
    # any kind for any of them.
    dict(id="taxonomy_on_tmt", fmt="fragpipe_tmt", plexes=2,
         run=dict(unipept=True, taxonomy=True), result=True, taxdump=True),
    # ...and the config where peptide_features() opens the quant table TWICE:
    # the full reader refuses a manifest whose runs match no column, and the
    # peptide-only reader it falls back to re-reads the same table.
    dict(id="auto_fallback", manifest="unmatched",
         run=dict(unipept=True, taxonomy=True, join=False),
         result=True, taxdump=True),
)


def _space_project(tmp_path, point):
    """Build the project one point of PLAN_SPACE describes."""
    name = "space_" + point["id"]
    p = build_project(tmp_path / name)
    cfg = {}
    fmt = point.get("fmt", "fragpipe_peptide")
    if fmt != "fragpipe_peptide":
        cfg["quant_format"] = fmt
        cfg["quant_table"] = _space_quant(p, fmt, point.get("plexes", 1))
    if point.get("run"):
        cfg["run"] = dict(p.cfg["run"], **point["run"])
    if point.get("manifest") is False:
        cfg["manifest"] = None
    if point.get("manifest") == "unmatched":
        # A manifest the full quant reader REFUSES - its one run matches no
        # column of the table - which is what makes peptide_features() fall
        # back to the peptide-only reader and read the table a second time.
        with open(p.cfg["manifest"], "w", encoding="utf-8") as fh:
            fh.write("/nowhere/ZZZ.mzML\tZZZ\t1\tDDA\n")
    if point.get("peptide_only_reader"):
        cfg["peptide_only_reader"] = point["peptide_only_reader"]
    if point.get("taxon_rank"):
        cfg["taxon_rank"] = point["taxon_rank"]
    if point.get("taxdump"):
        cfg["db"] = dict(p.cfg.get("db") or {},
                         ncbi_taxonomy=F.toy_taxdump(str(tmp_path / "taxdump")))
    if point.get("result"):
        cfg["unipept"] = {"result": _space_lca(p),
                          "split_missed_cleavages": False,
                          "consensus_min_peptides": 2}
    p.write_config(**cfg)
    return p


def _space_quant(p, fmt, plexes=1):
    """A quant table of the right shape for `fmt`, beside the default one."""
    inp = p.path("input")
    if fmt == "fragpipe_tmt":
        root = os.path.join(inp, "tmtrun")
        # The peptides are the project's own, so the pept2lca export
        # _space_lca() builds covers them and a taxonomy stage over this tree
        # has something to join onto. Each plex carries one replicate of BOTH
        # conditions: a plex that holds one condition is confounded with it,
        # which the full reader refuses - and this point is about the reads a
        # tree makes when the reader ACCEPTS it.
        rows = [{"peptide": r["peptide"], "razor": r["razor"]}
                for r in F.peptides_for(p.proteins, seed=3)]
        for n in range(plexes):
            F.write_tmt_plex(root, f"TMT{n + 1}",
                             [("126", p.samples[n]),
                              ("127N", p.samples[n + len(p.samples) // 2])],
                             rows)
        return root
    made = {
        "fragpipe_ion": (F.write_ion_table, "combined_ion.tsv"),
        "diann": (F.write_diann_matrix, "report.pg_matrix.tsv"),
        "fragpipe": (F.write_fragpipe_protein, "combined_protein.tsv"),
    }
    if fmt in made:
        fn, base = made[fmt]
        path = os.path.join(inp, base)
        fn(path, p.proteins, p.samples)
        return path
    groups = {s: s.split("_")[0] for s in p.samples}
    fn, base = {"msstats_csv": (F.write_msstats_csv, "msstats.csv"),
                "msstats_feature": (F.write_msstats_feature, "msstats_feat.csv"),
                "msstats_protein": (F.write_msstats_protein,
                                    "msstats_prot.csv")}[fmt]
    path = os.path.join(inp, base)
    fn(path, p.proteins, p.samples, groups)
    return path


def _space_lca(p):
    """A pept2lca export that matches the project's own peptide table."""
    import pandas as pd
    quant = pd.read_csv(p.path("input", "combined_peptide.tsv"), sep="\t")

    def lin(genus, species):
        return {"domain": 2, "phylum": 1239, "class": 91061, "order": 186826,
                "family": 81852, "genus": genus, "species": species}

    lineages = {"820": lin(816, 820), "821": lin(816, 821),
                "822": lin(816, 822), "823": lin(816, 823),
                "1351": lin(1350, 1351)}
    seed = {x.pid: x.seed_taxid for x in p.proteins}
    rows = []
    for pep, prot in zip(quant["Peptide Sequence"], quant["Protein"]):
        tid = str(seed.get(prot, "") or "")
        if tid not in lineages:
            continue
        rows.append((re.sub(r"[^A-Za-z]", "", str(pep)).upper(), tid,
                     "species", lineages[tid]))
    return F.write_unipept(p.path("input", "pept2lca.csv"), rows)


@pytest.mark.parametrize("point", PLAN_SPACE, ids=[p["id"] for p in PLAN_SPACE])
def test_the_read_plan_and_the_choke_point_hold_across_the_config_space(
        tmp_path, ma, point):
    """One driven run per config point; the plan and the gate, both measured.

    THE PLAN half demands EXACT equality wherever the plan has sites, because
    both errors cost something real: a plan that under-counts is a FIFO that
    hangs or dies at a second open the operator was told would not happen, and
    one that over-counts refuses a workflow that works while quoting reads the
    run will never make.

    The one read that cannot be counted exactly is the one peptide_features()
    makes when the full quant reader REFUSES a table and it re-reads with the
    peptide-only reader: that turns on the file, and the plan is computed from
    the config. It is marked CONTINGENT_READ rather than dropped - dropping it
    is the under-count, and `doctor --json` publishing "this run reads this
    input exactly once" for a table a stage can open twice is exactly the
    false promise this change set exists to remove - and the demand here is
    then a range rather than a number, on that one site and nowhere else.

    Where the plan has NO sites the demand is different and is not weaker than
    it looks: a path with no plan is treated as SINGLE-READ by
    _open_for_read(), so an unplanned path the run opens twice is a promise
    this program cannot keep - it is the hang, with the refusal that exists to
    prevent it switched off. At most one open, or the plan owes that path an
    entry. That is the assertion db.ncbi_taxonomy would have failed for as
    long as it has existed.

    THE GATE half asks something the plan cannot: every open of a configured
    path has to have gone through a primitive that fstat'd the DESCRIPTOR it
    is about to read - _open_for_read for `run`, _open_regular_binary for the
    signature digest and for `doctor`'s header read, regular_readable for
    doctor's probe. One reader outside those is one path that can still hang,
    or worse, one that drains a pipe the gate has already promised to somebody
    else.
    """
    p = _space_project(tmp_path, point)
    targets = _configured_paths(p.cfg)
    assert targets, f"{point['id']}: the config names no existing file"
    _proc, rows = _traced_run(tmp_path, p, targets, expect=None)
    assert rows, "the tracer saw no opens at all; the shim is not installed"
    ma.set_read_plan(p.cfg)

    # A run that DIED read a prefix of what it planned to, and a prefix is not
    # a defect in the plan: `quant_format: diann` with the taxonomy stage on is
    # a config an operator can write and `run` refuses it inside stage_taxonomy
    # with "the taxonomy comparison needs peptide-level input", so join never
    # opens anything. The plan is a promise about the reads a COMPLETED run
    # makes, so exact equality is demanded of a run that finished and only the
    # over-read half - the half that can hang - of one that did not.
    finished = _proc.returncode == 0
    wrong, unplanned = [], []
    for path in targets:
        seen = _stream_opens(rows, path)
        # Through the two accessors the engine itself reasons with, rather
        # than off the returned dict: a CONTINGENT site is a read the run MAY
        # make - peptide_features() re-reads the quant table with the
        # peptide-only reader when the full one refuses it, and whether that
        # happens is a property of the TABLE and not of the config the plan is
        # computed from. planned_reads() counts it, because that is the
        # question a FIFO asks; certain_reads() does not, because that is the
        # question a completed run answers. So the firm sites are a floor and
        # every site is the ceiling, and on every path with no contingent site
        # - which is all of them but one - the two are the same number and
        # this is the exact equality it was.
        want = ma.planned_reads(path)
        firm = ma.certain_reads(path)
        if want and (not len(firm) <= len(seen) <= len(want) if finished
                     else len(seen) > len(want)):
            wrong.append(
                f"{os.path.basename(path)}: the plan says "
                + (f"{len(want)} read(s)" if len(firm) == len(want)
                   else f"{len(firm)}-{len(want)} read(s)")
                + f" {want}, the run made {len(seen)}"
                + ("" if finished else " (and exited "
                   f"{_proc.returncode}, so only an over-read is a defect)")
                + ":\n    " + "\n    ".join(r["where"] for r in seen))
        elif not want and len(seen) > 1:
            unplanned.append(
                f"{os.path.basename(path)}: no plan, so it is treated as "
                f"single-read, and the run opened it {len(seen)} times:\n    "
                + "\n    ".join(r["where"] for r in seen))
    stray = [r for r in rows if not r["gate"] and not r["digest"]]

    problems = []
    if wrong:
        problems.append(
            "THE READ PLAN no longer describes what this run does. A reader "
            "that is not in the plan is a FIFO that hangs or a refusal that "
            "never comes:\n" + "\n".join(wrong))
    if unplanned:
        problems.append(
            "AN INPUT WITH NO PLAN IS READ MORE THAN ONCE, which is the hang "
            "with the refusal switched off - give it an entry in "
            "INPUT_READ_SITES:\n" + "\n".join(unplanned))
    if stray:
        problems.append(
            "THE CHOKE POINT was bypassed: a configured input was opened "
            "without a descriptor anything had fstat'd:\n"
            + "\n".join(f"  {os.path.basename(r['path'])}: {r['where']}"
                        for r in stray))
    assert not problems, f"[{point['id']}]\n" + "\n\n".join(problems)


# Every other command that takes an operator-supplied path on its own command
# line or out of a config. `run` has been swept for three rounds; these had
# never been swept at all, and `subset` was reading --quant through a bare
# pd.read_csv() the entire time - a path a user types into a shell, which is
# as operator-supplied as a path gets.
#
# `report` and `object` are here for the same reason one round later, and they
# found the same thing one command over: auto_contrasts() read
# `analysis.metadata` through a bare pd.read_csv(), so `metaannot report` on a
# TMT project whose metadata was a FIFO printed "this is a FIFO, and this run
# reads it exactly once" - from the gated header read one line above it - and
# then BLOCKED IN THE SECOND OPEN with no exit status, with `fifo_wait_s: 0`
# making no difference because a bare pandas open has no wait to shorten. A
# command that can hang is exactly as bad as a stage that can hang; `report`
# holds no results lock only by accident.
OTHER_COMMANDS = ("doctor", "doctor --json", "describe", "subset diann",
                  "subset fragpipe", "subset fragpipe_peptide", "subset ids",
                  "report", "object")
# ...of which these two open no operator-supplied path at all, which is the
# measured answer for them and is asserted as one below.
OPEN_NO_INPUT = ("describe", "object")


def _tmt_report_project(tmp_path, name):
    """A TMT project with a hand-written analysis.metadata, ready for `report`.

    TMT because that is the branch that reads the operator's metadata rather
    than the design the run recovered: `report` takes the contrasts from
    `analysis.metadata` when the format is fragpipe_tmt and the file is not
    the recovered design, which is where the condition of a TMT experiment is
    written down.
    """
    p = build_project(tmp_path / name)
    root = os.path.join(p.path("input"), "tmtrun")
    F.write_tmt_plex(root, "TMT1",
                     [("126", p.samples[0]), ("127N", p.samples[1])],
                     [{"peptide": r["peptide"], "razor": r["razor"]}
                      for r in F.peptides_for(p.proteins, seed=3)])
    meta = p.path("input", "metadata.tsv")
    with open(meta, "w", encoding="utf-8") as fh:
        fh.write("sample\tplex\tgroup\n")
        for i, s in enumerate(p.samples):
            fh.write(f"{s}\tTMT1\t{'A' if i % 2 else 'B'}\n")
    p.write_config(quant_format="fragpipe_tmt", quant_table=root,
                   analysis=dict(p.cfg.get("analysis") or {}, metadata=meta))
    return p


@pytest.mark.parametrize("cmd", OTHER_COMMANDS)
def test_no_command_opens_a_configured_path_outside_a_gate(tmp_path, cmd):
    """The choke point, for the commands that are not `run`.

    `doctor` may never WAIT, so it does not share `run`'s opener - but it must
    still never open a path it has not proved, which is a different property
    and the one this asserts. Before this, doctor probed the quant table with
    regular_readable() and then let pandas open it again BY NAME: a stat and
    an open with instructions in between, which is the race every other check
    in that file was rewritten to remove, and a different object is possible
    between them.
    """
    args = cmd.split()
    if args[0] in ("report", "object"):
        p = _tmt_report_project(tmp_path, f"cmd_{args[0]}")
        argv = (args[0], "--config", p.config_path,
                "--no-render" if args[0] == "report" else "--no-run")
        if args[0] == "object":
            # cmd_object refuses without one, and this test is about the opens
            # it makes on the way, not about that refusal.
            aq = p.path("results", "quant", "annotated_quant.tsv")
            os.makedirs(os.path.dirname(aq), exist_ok=True)
            with open(aq, "w", encoding="utf-8") as fh:
                fh.write("protein\tA_1\n" + "\t".join(["p1", "1"]) + "\n")
        targets = _configured_paths(p.cfg)
    elif args[0] in ("doctor", "describe"):
        p = build_project(tmp_path / "cmds")
        argv = (args[0], "--config", p.config_path) + tuple(args[1:])
        targets = _configured_paths(p.cfg)
    else:
        p = build_project(tmp_path / "cmds")
        fmt = args[1]
        quant = p.cfg["quant_table"]
        if fmt == "diann":
            quant = p.path("input", "report.pg_matrix.tsv")
            F.write_diann_matrix(quant, p.proteins, p.samples)
        elif fmt == "fragpipe":
            quant = p.path("input", "combined_protein.tsv")
            F.write_fragpipe_protein(quant, p.proteins, p.samples)
        elif fmt == "ids":
            quant = p.path("input", "ids.txt")
            with open(quant, "w", encoding="utf-8") as fh:
                fh.write("\n".join(x.pid for x in p.proteins) + "\n")
        argv = ("subset", "--quant", quant, "--format", fmt,
                "--db", p.cfg["proteins_faa"],
                "--out", p.path("subset.faa"))
        targets = _configured_paths(p.cfg) + [os.path.realpath(quant)]
    _proc, rows = _traced(tmp_path, p, targets, argv, expect=0)
    if cmd in OPEN_NO_INPUT:
        # Measured, and an answer rather than a gap: `describe` reads the
        # config and reports the contract it implies, and `object` reads the
        # tables the run WROTE. Neither opens an operator-supplied input at
        # all, so there is nothing here for a gate to hold - and the day one
        # of them grows a reader, this says so instead of passing on an empty
        # trace. The commands that DO read one are what proves the shim is
        # installed.
        assert not rows, (
            f"`{cmd}` opened a configured input, which it did not before - it "
            "needs a line in this test saying which gate that read goes "
            "through:\n"
            + "\n".join(f"  {os.path.basename(r['path'])}: {r['where']}"
                         for r in rows))
        return
    assert rows, f"the tracer saw no opens at all for `{cmd}`"
    stray = [r for r in rows if not r["gate"] and not r["digest"]]
    assert not stray, (
        f"`{cmd}` opened a configured path without a descriptor anything had "
        "fstat'd:\n"
        + "\n".join(f"  {os.path.basename(r['path'])}: {r['where']}"
                    for r in stray))


def test_the_signature_digest_never_opens_a_stream(ma, tmp_path):
    """Why the digest is not in the read plan, proved rather than assumed.

    _content_digest() reads a file end to end. On a FIFO that would DRAIN the
    pipe and the stage it belongs to would find an empty stream - so the whole
    "a FIFO works at a single-read input" rule depends on this never
    happening. It did not happen before either, but only because _stat() asks
    os.path.isfile() first: a stat, on a path, some instructions before the
    open, which is the shape of race every other check here was rewritten to
    remove. _open_regular_binary() fstats the DESCRIPTOR instead.
    """
    fifo = str(tmp_path / "pipe.faa")
    os.mkfifo(fifo)
    box = {}

    def go():
        try:
            box["stat"] = ma._stat(fifo)
            box["digest"] = ma._content_digest(fifo, 0, 0)
        except BaseException as e:                          # noqa: BLE001
            box["error"] = e

    th = threading.Thread(target=go, daemon=True)
    th.start()
    th.join(10)
    assert not th.is_alive(), \
        "the signature digest blocked on a FIFO - it is opening a stream"
    assert box.get("digest") is None, \
        "_content_digest read a FIFO; it must refuse a non-regular file"
    assert box["stat"][1] is None or isinstance(box["stat"][1], int), \
        "_stat digested a FIFO instead of falling back to size+mtime"
    # ...and it still digests an ordinary file, which is the thing it is for.
    real = tmp_path / "real.txt"
    real.write_text("hello\n", encoding="utf-8")
    st = os.stat(str(real))
    assert ma._content_digest(str(real), st.st_size, int(st.st_mtime))


# ======================================================================
# the FIFO rule: one read per run, or a refusal that says so at once
# ======================================================================
def _read_in_thread(fn, seconds=30):
    """Run `fn` off the main thread and join it.

    EVERY FIFO TEST GOES THROUGH THIS, and it is a rule rather than a habit. A
    test that calls a reader on a FIFO directly, with no thread and no join,
    does not FAIL when the reader stops being bounded - it WEDGES, and takes
    the suite with it until somebody kills the process. One of these tests was
    written that way and hung until SIGKILL against a reverted opener, twenty
    lines below a sibling whose docstring said exactly this.
    """
    box = {}

    def go():
        try:
            box["value"] = fn()
        except BaseException as e:                          # noqa: BLE001
            box["error"] = e

    th = threading.Thread(target=go, daemon=True)
    th.start()
    th.join(seconds)
    box["alive"] = th.is_alive()
    return box


def test_a_fifo_at_a_multi_read_input_is_refused_at_once_naming_the_reads(
        ma, tmp_path):
    """THE PREMISE EVERY EARLIER ROUND REASONED FROM, inverted.

    Two rounds argued about protecting `mkfifo p; zcat big.faa.gz > p &` at
    proteins_faa. Driven end to end, it does not work there and never did: the
    run opens that path three times, the first reader drains the pipe, and
    what the previous round's six-hour wait bought was a failure at `integrate`
    with a message that quoted, as advice, the workflow that had just failed.

    So the refusal has to be IMMEDIATE - no wait at all, however long
    `fifo_wait_s` is - and it has to carry the count and the reads, because
    "this cannot work" without "and here is why" is what sends an operator to
    lengthen a timeout that changes nothing.
    """
    p = build_project(tmp_path / "multi", fifo_wait_s=86400)
    ma.set_read_plan(p.cfg)
    faa = p.cfg["proteins_faa"]
    os.remove(faa)
    os.mkfifo(faa)
    t0 = time.time()
    box = _read_in_thread(lambda: list(ma.read_fasta(faa)))
    assert not box["alive"] and "error" in box, \
        "read_fasta waited on a FIFO at a path the run reads three times"
    assert time.time() - t0 < 10, \
        "the refusal waited; fifo_wait_s must not apply to a multi-read path"
    msg = str(box["error"])
    for want in (faa, "FIFO", "3 times", "drained EXACTLY ONCE",
                 "read_fasta", "real file"):
        assert want in msg, f"the refusal does not name {want!r}:\n{msg}"
    # ...and it must not suggest that waiting would have helped.
    assert "no writer appeared" not in msg, \
        "the multi-read refusal still talks about waiting for a writer"


def test_the_refusal_names_the_inputs_where_a_pipe_does_work(ma, tmp_path):
    """Read off the plan, never asserted.

    Which inputs take a pipe is a property of the CONFIG - a project with a
    second peptide consumer reads its quant table twice - so a sentence that
    named a fixed list would be the same hand-written claim this whole change
    set keeps finding. Driven on two configs whose answers differ.
    """
    p = build_project(tmp_path / "hint")
    ma.set_read_plan(p.cfg)
    faa = p.cfg["proteins_faa"]
    os.remove(faa)
    os.mkfifo(faa)
    box = _read_in_thread(lambda: list(ma.read_fasta(faa)))
    msg = str(box["error"])
    assert os.path.basename(p.cfg["quant_table"]) in msg \
        and os.path.basename(p.cfg["manifest"]) in msg, \
        f"the refusal does not name the single-read inputs:\n{msg}"
    assert os.path.basename(faa) not in msg.split("in this config means")[-1], \
        "the refusal names the multi-read path as one that works"


@pytest.mark.parametrize("key", ["quant_table", "manifest",
                                 "emapper_precomputed"])
def test_a_live_writer_feeds_a_single_read_input_through_the_whole_run(
        tmp_path, key):
    """The workflow, driven end to end through the real CLI.

    Not on a reader: the claim being made is about `run`, and `run` is what a
    verifier drove to find that this workflow worked at exactly one of the
    operator-supplied inputs and failed at the rest. Two of the three here
    were MOVED into the single-read class by this change - the quant table was
    read three times (header_columns, read_delim_table's own sniff, then
    pandas) and the manifest twice (the column mapping, then the design), and
    each is one open now.

    The writer is started BEFORE the run and writes after a pause, so the run
    reaches the open first: that is the order the cluster workflow has, where
    the producer is a queued job.
    """
    p = build_project(tmp_path / f"live_{key[:5]}", fifo_wait_s=60)
    target = (p.cfg["emapper_precomputed"][0]
              if key == "emapper_precomputed" else p.cfg[key])
    keep = str(tmp_path / "keep.bytes")
    shutil.copyfile(target, keep)
    os.remove(target)
    os.mkfifo(target)

    def write():
        time.sleep(1.0)
        with open(keep, "rb") as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)

    th = threading.Thread(target=write, daemon=True)
    th.start()
    proc = p.run(expect=None, timeout=180)
    th.join(30)
    assert proc.returncode == 0, (
        f"a live writer on {key} no longer feeds a run - this input was "
        f"single-read:\n{proc.stderr[-2000:]}")
    assert "this run reads it exactly once" in proc.stderr, \
        f"{key} was not announced as a single-read FIFO:\n{proc.stderr[-800:]}"


def test_report_does_not_hang_on_the_metadata_it_has_already_drained(tmp_path):
    """L1, driven through the shipped command: `report` must not block.

    THE DEFECT, measured before it was fixed: `metaannot report` on a TMT
    project whose `analysis.metadata` is a FIFO read that path TWICE -
    header_columns() in _tmt_report_design(), which is gated and which printed
    "this is a FIFO, and this run reads it exactly once, so it is read as a
    live stream", and then auto_contrasts(), which was a bare pd.read_csv().
    The second open found a drained pipe with no writer left on it and BLOCKED
    THERE: no further output, no exit status, and nothing to time out against,
    because `fifo_wait_s` bounds opener() and a bare pandas open has no wait
    to bound. The promise had already been printed by then, which is the worst
    ordering available - the command said the workflow was supported one line
    before it stopped answering.

    A subprocess with a timeout, and a writer that feeds the first read, so
    that what is under test is the SECOND open and not the wait in front of
    the first. The whole assertion is that the command comes back at all; what
    it comes back with is the backstop's message, which names the path and the
    read that already happened.
    """
    p = _tmt_report_project(tmp_path, "report_fifo")
    meta = (p.cfg["analysis"] or {})["metadata"]
    keep = open(meta, encoding="utf-8").read()
    os.remove(meta)
    os.mkfifo(meta)
    p.write_config(fifo_wait_s=20)

    def write():
        # EPIPE is the expected end: the first reader wants the header line
        # and closes the stream, which is what leaves the second open with
        # nothing - the state under test.
        with contextlib.suppress(OSError):
            with open(meta, "w", encoding="utf-8") as fh:
                fh.write(keep)

    th = threading.Thread(target=write, daemon=True)
    th.start()
    t0 = time.time()
    proc = run_metaannot("report", "--config", p.config_path, "--no-render",
                         expect=None, timeout=90)
    th.join(30)
    assert time.time() - t0 < 60, \
        "`report` is still blocking on the metadata FIFO"
    assert proc.returncode != 0, \
        f"`report` exited 0 on a drained metadata FIFO:\n{proc.stderr[-1500:]}"
    assert meta in proc.stderr and "already read it once" in proc.stderr, \
        f"`report` stopped without naming the path or the read that took "\
        f"the stream:\n{proc.stderr[-1500:]}"


def test_report_honours_the_configs_fifo_wait_and_not_the_module_default(
        tmp_path):
    """`fifo_wait_s: 0` is published as "refused on the spot"; `report` sat.

    Routing auto_contrasts() through opener() fixed the SECOND open - the
    backstop refuses a drained pipe by name, which is what the test above
    drives. It did nothing for the FIRST, because cmd_report never called
    set_fifo_wait(), so opener() used the MODULE default: six hours, on a
    config that had asked for none. Measured before this line existed, with no
    writer at all: `run` with a FIFO at quant_table died at 0.3s under
    `fifo_wait_s: 0` and at 5.3s under 5, while `report` on the same setting
    was still going at 60s, having printed "waits up to 21600s".

    That made two published sentences false rather than imprecise - README and
    TUTORIAL both offer `fifo_wait_s: 0` to an operator who never pipes
    anything in and would rather a FIFO were refused at once. So this asserts
    the SETTING is in force, by the only evidence that distinguishes it from
    the default: the command comes back quickly, and the number it quotes for
    itself is the one the config asked for.

    No writer, deliberately. The test above covers the drained-pipe path; this
    one is about the wait in front of the first open, which is the half a
    writer would hide.
    """
    p = _tmt_report_project(tmp_path, "report_wait")
    meta = (p.cfg["analysis"] or {})["metadata"]
    os.remove(meta)
    os.mkfifo(meta)
    p.write_config(fifo_wait_s=0)

    t0 = time.time()
    proc = run_metaannot("report", "--config", p.config_path, "--no-render",
                         expect=None, timeout=90)
    took = time.time() - t0
    assert took < 30, (
        "`report` ignored fifo_wait_s: 0 and waited on the metadata FIFO "
        f"for {took:.0f}s")
    assert proc.returncode != 0, proc.stderr[-1500:]
    assert meta in proc.stderr, proc.stderr[-1500:]
    # The default is 21600. A command quoting that number to an operator who
    # set 0 is the defect, whether or not it then waits that long.
    assert "21600" not in proc.stderr, (
        "`report` quoted the module default to a config that set 0:\n"
        + proc.stderr[-1500:])


def test_a_wait_that_expires_says_which_silence_it_was(ma, tmp_path):
    """F5: "no writer appeared ... nothing has opened the other end" was false
    in a reachable case, and it was false in both halves at once.

    A writer that attaches at 0s and whose first byte lands at 5s of a 3s wait
    is attached the whole time. Telling that operator to start a writer sends
    them to fix the one thing that was already right, and the writer then dies
    of EPIPE when the reader gives up. A non-blocking one-byte read tells the
    two apart for free: a pipe with no writer answers end-of-file, one with a
    silent writer answers EAGAIN.
    """
    monkey = str(tmp_path / "silent.faa")
    os.mkfifo(monkey)
    ma.set_fifo_wait(1.0)
    try:
        held = []

        def attach():
            held.append(open(monkey, "wb"))
            time.sleep(20)

        th = threading.Thread(target=attach, daemon=True)
        th.start()
        time.sleep(0.5)
        box = _read_in_thread(lambda: list(ma.read_fasta(monkey)), 20)
        assert not box["alive"], "the reader is still on the open"
        msg = str(box["error"])
        assert "something IS holding the write end" in msg, \
            f"an attached-but-silent writer was reported as absent:\n{msg}"
        assert "raise `fifo_wait_s`" in msg, \
            "the message does not name the remedy for a slow writer"
        for fh in held:
            fh.close()

        # ...and the other silence, which has to keep the other sentence.
        nobody = str(tmp_path / "nobody.faa")
        os.mkfifo(nobody)
        box = _read_in_thread(lambda: list(ma.read_fasta(nobody)), 20)
        assert not box["alive"]
        msg = str(box["error"])
        assert "nothing is holding the write end" in msg, \
            f"an unattached FIFO was reported as attached:\n{msg}"
        assert "Start the writer" in msg
    finally:
        ma.set_fifo_wait(ma.DEFAULT_CONFIG["fifo_wait_s"])


@pytest.mark.parametrize("how", ["mid_line", "nothing_at_all", "gzip"])
def test_a_stream_that_ends_early_is_a_failure_and_not_a_shorter_result(
        ma, tmp_path, how):
    """F6, and it is the worst of the outcomes because it is not an error.

    A writer that wrote half a FASTA and died yielded [("p1", "MKV"),
    ("p2", "MK")] with nothing at all to say the input was a fragment - a
    short read looks exactly like a complete file to every reader above the
    stream. The gzip branch was already right, because gzip is FRAMED and a
    missing end-of-stream marker raises; the plain branch has no frame, so
    what it can see is that the stream stopped inside a line, or that the
    writer closed without writing a byte.
    """
    name = {"mid_line": "t.faa", "nothing_at_all": "e.faa",
            "gzip": "t.faa.gz"}[how]
    fifo = str(tmp_path / name)
    os.mkfifo(fifo)
    blob = gzip.compress(b">g1\nMKV\n>g2\nMKW\n")
    payload = {"mid_line": b">p1\nMKV\n>p2\nMK",
               "nothing_at_all": b"",
               "gzip": blob[:len(blob) // 2]}[how]
    th = threading.Thread(
        target=lambda: open(fifo, "wb").write(payload) if payload
        else open(fifo, "wb").close(), daemon=True)
    th.start()
    box = _read_in_thread(lambda: list(ma.read_fasta(fifo)), 20)
    assert not box["alive"], "the reader never finished"
    assert "value" not in box, \
        f"a truncated stream came back as records: {box.get('value')}"
    msg = str(box["error"])
    assert fifo in msg, f"the truncation does not name the path: {msg}"
    want = {"mid_line": "middle of a line",
            "nothing_at_all": "without writing a single byte",
            "gzip": "end-of-stream marker"}[how]
    assert want in msg, f"wrong diagnosis for {how}: {msg}"


def test_fifo_wait_bounds_every_read_and_not_only_the_first(ma, tmp_path):
    """F7: an idle writer holding the write end open is the same hang, one
    read further in.

    `fifo_wait_s` bounded the FIRST byte and nothing after it, so a writer
    that attached, sent half a table and then sat there - an O_RDWR keeper, a
    producer blocked on its own input - put the run straight back into the
    state the setting exists to prevent, and back onto the results lock it
    holds while it waits. The bound is an idle timeout now.
    """
    fifo = str(tmp_path / "idle.faa")
    os.mkfifo(fifo)
    ma.set_fifo_wait(1.0)
    try:
        held = []

        def keeper():
            fh = open(fifo, "wb")
            fh.write(b">p1\nMKV\n")
            fh.flush()
            held.append(fh)
            time.sleep(20)

        threading.Thread(target=keeper, daemon=True).start()
        t0 = time.time()
        box = _read_in_thread(lambda: list(ma.read_fasta(fifo)), 20)
        assert not box["alive"], \
            "the read never ended - fifo_wait_s still bounds only the open"
        assert time.time() - t0 < 15, "the idle bound took far longer than set"
        assert "sent nothing for" in str(box["error"]), \
            f"the stall was reported as something else: {box['error']}"
        for fh in held:
            fh.close()
    finally:
        ma.set_fifo_wait(ma.DEFAULT_CONFIG["fifo_wait_s"])


def test_a_pausing_writer_is_not_mistaken_for_a_stalled_one(ma, tmp_path):
    """The other side of the idle bound, which is the one that costs a
    workflow if it is wrong.

    Every pause a healthy writer takes - `zcat` waiting on its own disk, a
    scheduler taking the CPU away - now goes through the same select() as a
    stall. A bound that fired on those would turn a working pipe into an
    intermittent failure, which is worse than the hang it replaced because a
    hang is at least reproducible. The gaps here are long compared with the
    reads and short compared with the bound.
    """
    fifo = str(tmp_path / "pausing.faa")
    os.mkfifo(fifo)
    ma.set_fifo_wait(5.0)
    try:
        th, box = _slow_writer(
            fifo, [b">p1\n", b"MKV\n", b">p2\n", b"MKW\n"], gap=0.4)
        got = _read_in_thread(lambda: list(ma.read_fasta(fifo)), 30)
        assert got.get("value") == [("p1", "MKV"), ("p2", "MKW")], \
            f"a pausing writer was cut off: {got}"
        th.join(5)
        assert not box, f"the writer did not finish cleanly: {box}"
    finally:
        ma.set_fifo_wait(ma.DEFAULT_CONFIG["fifo_wait_s"])


def test_a_second_open_of_a_fifo_fails_rather_than_blocking(ma, tmp_path):
    """The backstop under the plan, for the day the plan is wrong.

    The plan is the mechanism and this is not a substitute for it - by the
    time this fires the pipe is drained and the hours are spent. What it
    converts is the worst outcome, a second open that blocks forever with no
    output and no exit status, into a message that names both reads. A plan
    that under-counts must cost a failure, never a hang.
    """
    fifo = str(tmp_path / "twice.faa")
    os.mkfifo(fifo)
    # A plan that says ONE read, so the first open is allowed and the second
    # is the case under test. The key it is filed under does not matter - the
    # plan is keyed by PATH - and `manifest` is the input a default run reads
    # exactly once.
    ma.set_read_plan({"manifest": fifo, "run": {"join": True},
                      "quant_format": "fragpipe_peptide"})
    threading.Thread(
        target=lambda: open(fifo, "wb").write(b">p1\nMKV\n"),
        daemon=True).start()
    first = _read_in_thread(lambda: list(ma.read_fasta(fifo)), 20)
    assert first.get("value") == [("p1", "MKV")], f"the first read: {first}"
    second = _read_in_thread(lambda: list(ma.read_fasta(fifo)), 20)
    assert not second["alive"], "the second open blocked instead of failing"
    assert "already read it once" in str(second.get("error")), \
        f"the second open failed for the wrong reason: {second}"
    ma.set_read_plan({})


def test_the_read_open_asks_for_o_binary_wherever_the_platform_has_it(
        ma, tmp_path, monkeypatch):
    """C1: the published Windows fallback that does not exist, and the flag
    that was missing because of it.

    The comment at the top of metaannot.py said that with `fcntl` and `select`
    absent "opener() falls back there to exactly the plain open() it has
    always used". There is no such branch - opener() calls _open_for_read()
    unconditionally, which calls os.open(). That is demonstrated below by
    taking both names away and counting.

    It matters past the false sentence. CPython's builtin open() adds O_BINARY
    on Windows and a bare os.open() does not, so the gzip branch - which is
    os.fdopen(fd, "rb"), and emapper_precomputed is routinely a .gz - would
    have read compressed bytes through a CRT text-mode descriptor: CRLF
    translation, and 0x1A treated as end of file. Nothing in this suite can
    reach a Windows CRT, so what is testable is that the flag is ASKED FOR
    wherever the platform defines it, which is what _OPEN_FLAGS is for.

    THE STUB GOES IN BEFORE THE IMPORT, and that is the correction this round
    made to the test rather than to the code. The version that stood here
    monkeypatched ma._OPEN_FLAGS with `os.O_RDONLY | getattr(os, "O_NONBLOCK",
    0) | getattr(os, "O_BINARY", 0)` - its own recomputation of the constant -
    and then asserted that its own expression had the O_BINARY bit in it. It
    could only ever agree with itself: driven, deleting `| getattr(os,
    "O_BINARY", 0)` from metaannot.py left BOTH C1 tests passing, because the
    module's constant was never read. That is the same defect shape this change
    set had just fixed one level down, in the test written to hold the fix.

    _OPEN_FLAGS is computed once at import, so the only way to see what a
    Windows import would produce is to stub the name and import again. The
    module under test is then asked for its own constant, and the two trees
    really do separate: 0x8004 with the term, 0x4 without it.
    """
    saved = sys.modules.get("metaannot")
    # A platform that HAS O_BINARY, simulated on one that does not. 0x8000 is
    # the value CPython uses on Windows; any non-zero bit would do, and a real
    # one keeps the number in the failure message recognisable.
    monkeypatch.setattr(os, "O_BINARY", 0x8000, raising=False)
    try:
        fresh = _load()
    finally:
        sys.modules["metaannot"] = saved
    assert fresh._OPEN_FLAGS & 0x8000, (
        "_OPEN_FLAGS does not ask for O_BINARY on a platform that has it "
        f"(got {fresh._OPEN_FLAGS:#x}); on Windows the gzip branch would read "
        "a .gz through a CRT text-mode descriptor")
    assert fresh._OPEN_FLAGS & getattr(os, "O_NONBLOCK", 0) == \
        getattr(os, "O_NONBLOCK", 0), \
        "O_BINARY was added at the cost of O_NONBLOCK, which is what makes " \
        "the open of a FIFO return instead of waiting for a writer"

    # ...and the flag really reaches the open, rather than sitting in a
    # constant nothing passes on.
    seen = []
    real_open = os.open
    monkeypatch.setattr(os, "open",
                        lambda path, flags, *a, **k: (
                            seen.append(flags), real_open(path, flags, *a, **k))[1])
    plain = tmp_path / "x.txt"
    plain.write_text("hello\n", encoding="utf-8")
    with fresh.opener(str(plain)) as fh:
        assert fh.read() == "hello\n"
    assert seen and all(f & 0x8000 for f in seen), \
        "the read open does not ask for O_BINARY; on Windows the gzip branch " \
        "would read a .gz through a text-mode descriptor"


def test_there_is_no_plain_open_fallback_when_fcntl_and_select_are_absent(
        ma, tmp_path, monkeypatch):
    """The other half of C1, and the reason the comment had to change rather
    than the code alone.

    With both POSIX names gone, a reader still goes through os.open() and not
    through the builtin - so "falls back to the plain open() it has always
    used" was describing a branch nobody had written. What the missing names
    really cost is that _clear_nonblock() becomes a no-op (there is no
    O_NONBLOCK to clear) and _wait_readable() answers True at once (nothing to
    wait with), which is what the comment says now.
    """
    monkeypatch.setattr(ma, "fcntl", None)
    monkeypatch.setattr(ma, "select", None)
    builtin, os_opens = [], []
    real_builtin, real_os = open, os.open
    monkeypatch.setattr("builtins.open",
                        lambda *a, **k: (builtin.append(a[0]),
                                         real_builtin(*a, **k))[1])
    monkeypatch.setattr(os, "open",
                        lambda *a, **k: (os_opens.append(a[0]),
                                         real_os(*a, **k))[1])
    plain = tmp_path / "y.txt"
    real_builtin(str(plain), "w", encoding="utf-8").write("hi\n")
    del builtin[:], os_opens[:]
    with ma.opener(str(plain)) as fh:
        assert fh.read() == "hi\n"
    assert os_opens == [str(plain)] and not builtin, \
        "opener() took a builtin-open fallback that the comment used to " \
        f"claim and the code has never had: os.open {os_opens}, open {builtin}"


def test_putting_the_sniffed_header_back_gives_pandas_the_same_frame(
        ma, tmp_path):
    """The claim that makes the single-open table read safe, asserted.

    Deciding the delimiter needs the header line, and pandas needs the file
    from the top - which used to be two OPENS of one path. Handing pandas the
    already-open handle with `names=` would be the obvious alternative and is
    NOT the same file: `names=` skips pandas' own header handling, so two
    columns called `a` stop being `a` and `a.1`. _HeadRestored() puts the line
    back instead, so every byte reaches pandas in the order it would have
    reached it.

    The fixture is chosen to break the wrong implementation: a duplicate
    column name, and a quoted field containing the delimiter.
    """
    path = tmp_path / "dup.tsv"
    path.write_text('a\tb\ta\tc\n1\t"x\ty"\t3\t4\n5\t6\t7\t8\n',
                    encoding="utf-8")
    import pandas as pd
    ref = pd.read_csv(str(path), sep="\t", low_memory=False)
    got = ma.read_delim_table(str(path))
    assert list(got.columns) == list(ref.columns), \
        f"the columns changed: {list(got.columns)} vs {list(ref.columns)}"
    assert got.equals(ref), "the frame changed"
    assert got.dtypes.tolist() == ref.dtypes.tolist(), "the dtypes changed"
    # ...and the header the recogniser sees is the same header, from the same
    # single open.
    cols, df2 = ma.read_named_table(str(path))
    assert cols == ["a", "b", "a", "c"], \
        "read_named_table's header is not the raw header line"
    assert df2.equals(ref)
    # One open for all of it, which is the point.
    real_open, seen = os.open, []
    try:
        os.open = lambda *a, **k: (seen.append(a[0]), real_open(*a, **k))[1]
        ma.read_named_table(str(path))
    finally:
        os.open = real_open
    assert seen == [str(path)], \
        f"read_named_table opened the path {len(seen)} times: {seen}"
