"""Stage plumbing: adoption, atomic writes, dependencies, the emapper join.

Every external tool is stubbed by a script on PATH, so these exercise the
pipeline's own logic without hmmsearch, DIAMOND, MMseqs2 or Foldseek.
"""
from __future__ import annotations

import json
import os
import shutil

import pytest

import fixtures as F
from conftest import build_project, run_metaannot


def _searchable(tmp_path, root, **over):
    """A project with the search stages on and fake database files present."""
    db = tmp_path / "db"
    db.mkdir(exist_ok=True)
    for name in ("Pfam-A.hmm", "dbCAN.txt", "hmm_PGAP.LIB"):
        (db / name).write_text("HMMER3/f\n", encoding="utf-8")
    (db / "vfdb.dmnd").write_text("fake", encoding="utf-8")
    (db / "merops.dmnd").write_text("fake", encoding="utf-8")
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


def test_adding_a_diamond_database_invalidates_integrate(ma, tmp_path,
                                                         stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    sig_before = proj.state()["integrate"]["signature"]
    db = tmp_path / "db" / "card.dmnd"
    db.write_text("fake", encoding="utf-8")
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
