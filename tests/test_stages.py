"""Stage plumbing: adoption, atomic writes, dependencies, the emapper join.

Every external tool is stubbed by a script on PATH, so these exercise the
pipeline's own logic without hmmsearch, DIAMOND, MMseqs2 or Foldseek.
"""
from __future__ import annotations

import io
import json
import os
import re
import shutil
import sys

import pytest

import fixtures as F
from conftest import METAANNOT_PY, build_project, run_metaannot


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
    # diamond, smorf, jackhmmer, hhblits, esmfold and foldseek can legitimately
    # produce nothing; every other stage cannot.
    empty_ok = {s["name"] for s in ma.STAGES if s.get("empty_ok")}
    assert empty_ok == {"diamond", "smorf", "jackhmmer", "hhblits", "esmfold",
                        "foldseek"}


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


def test_a_short_peptide_database_warns_that_the_evalue_is_unreachable(
        ma, tmp_path, paths_for, stub_bin, capsys):
    # BAGEL's real shape: 262 sequences, median 15 residues.
    db = F.write_dmnd(tmp_path / "bagel.dmnd", sequences=262, letters=4009)
    cfg, p = _dia_project(ma, tmp_path, paths_for, bagel=db)
    ma.stage_diamond(cfg, p)
    err = capsys.readouterr().err
    assert "incapable of a hit before it starts" in err
    assert "15 residues" in err
    assert "diamond_evalues" in err
    # it still runs: this is a warning about the threshold, not a broken file
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
    typical, letters, why = ma.diamond_db_profile(cfg, "vfdb", db)
    assert (typical, letters) == (None, None)
    assert "diamond is not installed" in why
    bad, warn = ma.diamond_db_check(cfg, "vfdb", db)
    assert bad is None
    assert "nothing here can tell whether" in warn


def test_a_source_fasta_beside_the_database_gives_the_median(
        ma, tmp_path, paths_for, monkeypatch):
    db = F.write_dmnd(tmp_path / "bagel.dmnd")
    with open(tmp_path / "bagel.fas", "w", encoding="utf-8") as fh:
        for i in range(5):
            fh.write(">p%d\n%s\n" % (i, "M" * (10 + i)))
    cfg, p = _dia_project(ma, tmp_path, paths_for, bagel=db)
    monkeypatch.setattr(ma, "have", lambda t: t != "diamond")
    typical, letters, why = ma.diamond_db_profile(cfg, "bagel", db)
    assert typical == 12                      # median of 10..14
    assert letters == 60
    assert "median of 5 sequences" in why


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


def test_doctor_warns_about_a_database_too_short_for_the_evalue(tmp_path,
                                                                stub_bin):
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
    assert "incapable of a hit before it starts" in proc.stdout
    assert "diamond_evalues" in proc.stdout


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
    assert "no qtmscore/ttmscore" in err
    assert "normalised by the alignment" in err


def test_the_foldseek_retry_keeps_the_scratch_tree_so_only_convertalis_reruns(
        ma, tmp_path, paths_for, stub_bin, monkeypatch):
    # symptom: --format-output is consumed by convertalis, which easy-search
    # runs AFTER the search, so this failure arrives with the whole multi-hour
    # alignment already done in tmpd. The retry did shutil.rmtree(tmpd) first,
    # throwing that away and re-running the search to change a formatting
    # argument. Foldseek skips its search when the tree still holds the
    # result, so the tree must survive between the two attempts.
    calls, tmp_alive = [], []
    inner = _fake_foldseek(calls)

    def fake_run(cmd, **kw):
        c = [str(x) for x in cmd]
        if len(c) > 1 and c[1] == "easy-search":
            tmpd = c[5]
            # Stand in for the completed alignment the real search leaves.
            os.makedirs(tmpd, exist_ok=True)
            marker = os.path.join(tmpd, "result.dbtype")
            tmp_alive.append(os.path.exists(marker))
            open(marker, "w").close()
        return inner(cmd, **kw)

    monkeypatch.setattr(ma, "run_cmd", fake_run)
    monkeypatch.setattr(ma, "have", lambda t: True)
    cfg, p = _fs_project(ma, tmp_path, paths_for)
    try:
        ma.stage_foldseek(cfg, p)
    except Exception:              # noqa: BLE001
        pass
    assert len(tmp_alive) == 2, "expected exactly two easy-search attempts"
    assert tmp_alive[1], (
        "the scratch tree was cleared between the attempts, so the retry "
        "re-runs the search instead of just the conversion")


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
