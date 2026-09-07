"""Invariants that must hold on a real run, and the reproducibility claims
the README makes but nothing verified."""
from __future__ import annotations

import hashlib
import os
import random
import shutil

import numpy as np
import pandas as pd
import pytest

import fixtures as F
from conftest import build_project, run_metaannot
from test_stages import _searchable


RESULT_FILES = ("annotation_final.tsv", "bin_summary.tsv",
                "quant/annotated_quant.tsv", "quant/peptide_evidence.tsv",
                "quant/feature_quant.tsv", "quant/taxon_intensity.tsv",
                "quant/taxon_size_factors.tsv",
                "quant/design_from_input.tsv", "quant/group_conflicts.tsv")


def _digest(proj, files=RESULT_FILES):
    h = hashlib.sha1()
    for rel in files:
        path = proj.rpath(*rel.split("/"))
        h.update(rel.encode())
        h.update(open(path, "rb").read() if os.path.exists(path) else b"<absent>")
    return h.hexdigest()


# --- finding 45 -------------------------------------------------------
@pytest.fixture(scope="module")
def finished(tmp_path_factory):
    """One complete run, reused by every invariant below."""
    root = tmp_path_factory.mktemp("invariants")
    proj = build_project(root / "p", n_extra=20)
    proj.run()
    return proj


def test_every_quantified_group_is_annotated(finished):
    aq = pd.read_csv(finished.rpath("quant", "annotated_quant.tsv"), sep="\t")
    ann = pd.read_csv(finished.rpath("annotation_final.tsv"), sep="\t")
    assert set(aq["group_id"]) <= set(ann["protein_id"].astype(str))
    assert aq["bin"].notna().all()


def test_every_assigned_protein_is_quantified(finished):
    fq = pd.read_csv(finished.rpath("quant", "feature_quant.tsv"), sep="\t")
    aq = pd.read_csv(finished.rpath("quant", "annotated_quant.tsv"), sep="\t")
    assigned = set(fq["assigned_protein"].dropna().astype(str))
    assert assigned <= set(aq["group_id"].astype(str))


def test_the_rollup_equals_the_sum_of_its_assigned_features(finished):
    fq = pd.read_csv(finished.rpath("quant", "feature_quant.tsv"), sep="\t")
    aq = pd.read_csv(finished.rpath("quant", "annotated_quant.tsv"), sep="\t")
    samples = finished.samples
    want = fq.dropna(subset=["assigned_protein"]).groupby(
        "assigned_protein")[samples].sum(min_count=1)
    got = aq.set_index("group_id")[samples]
    common = want.index.intersection(got.index)
    assert len(common) == len(want)
    np.testing.assert_allclose(got.loc[common].values, want.loc[common].values,
                               rtol=1e-9)


def test_n_features_used_matches_the_feature_table(finished):
    fq = pd.read_csv(finished.rpath("quant", "feature_quant.tsv"), sep="\t")
    ev = pd.read_csv(finished.rpath("quant", "peptide_evidence.tsv"), sep="\t")
    counted = fq.dropna(subset=["assigned_protein"]).groupby(
        "assigned_protein").size()
    ev = ev.set_index("protein_id")
    for pid, n in counted.items():
        assert int(ev.loc[pid, "n_features_used"]) == int(n)


def test_taxon_intensity_equals_the_sum_over_its_members(finished):
    aq = pd.read_csv(finished.rpath("quant", "annotated_quant.tsv"), sep="\t",
                     dtype={"effective_taxid": str})
    ti = pd.read_csv(finished.rpath("quant", "taxon_intensity.tsv"), sep="\t",
                     dtype={"effective_taxid": str}).set_index("effective_taxid")
    samples = finished.samples
    have = aq[aq["effective_taxid"].fillna("").ne("")]
    want = have.groupby("effective_taxid")[samples].sum()
    assert set(want.index) == set(ti.index)
    np.testing.assert_allclose(ti.loc[want.index, samples].values,
                               want.values, rtol=1e-9)


def test_every_taxon_with_members_has_a_finite_positive_size_factor(finished):
    sf = pd.read_csv(finished.rpath("quant", "taxon_size_factors.tsv"),
                     sep="\t", dtype={"effective_taxid": str})
    ti = pd.read_csv(finished.rpath("quant", "taxon_intensity.tsv"), sep="\t",
                     dtype={"effective_taxid": str})
    assert set(sf["effective_taxid"]) == set(ti["effective_taxid"])
    vals = sf[finished.samples].to_numpy(dtype=float)
    finite = vals[~np.isnan(vals)]
    assert finite.size, "no taxon got a size factor at all"
    assert (finite > 0).all()
    assert np.isfinite(finite).all()


def test_every_bin_is_one_of_the_seven(finished, ma):
    ann = pd.read_csv(finished.rpath("annotation_final.tsv"), sep="\t")
    assert set(ann["bin"]) <= set(ma.BIN_ORDER)
    aq = pd.read_csv(finished.rpath("quant", "annotated_quant.tsv"), sep="\t")
    assert set(aq["bin"].dropna()) <= set(ma.BIN_ORDER)


def test_no_duplicate_identifiers_anywhere(finished):
    for rel, col in (("annotation_final.tsv", "protein_id"),
                     ("quant/annotated_quant.tsv", "group_id"),
                     ("quant/peptide_evidence.tsv", "protein_id"),
                     ("quant/feature_quant.tsv", "feature_id")):
        path = finished.rpath(*rel.split("/"))
        df = pd.read_csv(path, sep="\t")
        c = col if col in df.columns else df.columns[0]
        dup = df[c][df[c].duplicated()]
        assert dup.empty, f"{rel} has duplicate {c}: {list(dup)[:5]}"


def test_no_negative_intensities(finished):
    for rel in ("quant/annotated_quant.tsv", "quant/feature_quant.tsv",
                "quant/taxon_intensity.tsv"):
        df = pd.read_csv(finished.rpath(*rel.split("/")), sep="\t")
        vals = df[[c for c in finished.samples if c in df.columns]]
        assert not (vals.fillna(0) < 0).any().any(), rel


def test_the_bin_summary_totals_match_the_annotation_table(finished):
    ann = pd.read_csv(finished.rpath("annotation_final.tsv"), sep="\t")
    s = pd.read_csv(finished.rpath("bin_summary.tsv"), sep="\t")
    assert int(s.loc[s["bin"] == "TOTAL", "n"].iloc[0]) == len(ann)
    per = s[s["bin"] != "TOTAL"].set_index("bin")["n"]
    assert per.sum() == len(ann)
    for b, n in ann["bin"].value_counts().items():
        assert int(per[b]) == int(n)


# --- finding 46 -------------------------------------------------------
FORMAT_BUILDERS = {}


def _project_for(root, fmt, **over):
    """A project whose quant_table is in `fmt`."""
    proj = build_project(root, **over)
    ps, samples = proj.proteins, proj.samples
    groups = {s: s.split("_")[0] for s in samples}
    inp = proj.path("input")
    if fmt == "fragpipe_peptide":
        return proj
    if fmt == "fragpipe_ion":
        q = F.write_ion_table(os.path.join(inp, "combined_ion.tsv"), ps,
                              samples)
    elif fmt == "fragpipe":
        q = F.write_fragpipe_protein(os.path.join(inp, "combined_protein.tsv"),
                                     ps, samples)
    elif fmt == "diann":
        q = F.write_diann_matrix(os.path.join(inp, "report.pg_matrix.tsv"), ps,
                                 samples)
    elif fmt == "msstats_csv":
        q = F.write_msstats_csv(os.path.join(inp, "MSstats.csv"), ps, samples,
                                groups)
    elif fmt == "msstats_feature":
        q = F.write_msstats_feature(os.path.join(inp, "feature.tsv"), ps,
                                    samples, groups)
    elif fmt == "msstats_protein":
        q = F.write_msstats_protein(os.path.join(inp, "protein.tsv"), ps,
                                    samples, groups)
    else:
        raise AssertionError(fmt)
    cfg = {"quant_table": q, "quant_format": fmt}
    if fmt.startswith("msstats"):
        cfg["manifest"] = ""          # the design comes from the long table
    proj.write_config(**cfg)
    return proj


ALL_FORMATS = ["diann", "fragpipe", "fragpipe_peptide", "fragpipe_ion",
               "msstats_csv", "msstats_feature", "msstats_protein"]


@pytest.mark.parametrize("fmt", ALL_FORMATS)
def test_parallel_and_serial_produce_identical_output(tmp_path, fmt):
    # symptom: the README says parallel and serial execution are INTENDED to
    # produce identical output and that no shipped test verifies it.
    par = _project_for(tmp_path / f"par_{fmt}", fmt)
    par.run()
    ser = _project_for(tmp_path / f"ser_{fmt}", fmt)
    ser.run("--serial")
    assert _digest(par) == _digest(ser)


@pytest.mark.parametrize("mode", ["protein_unique", "taxon_unique", "razor"])
def test_parallel_and_serial_agree_for_every_assignment_mode(tmp_path, mode):
    par = _project_for(tmp_path / f"par_{mode}", "fragpipe_peptide",
                       peptide_assignment=mode)
    par.run()
    ser = _project_for(tmp_path / f"ser_{mode}", "fragpipe_peptide",
                       peptide_assignment=mode)
    ser.run("--serial")
    assert _digest(par) == _digest(ser)


@pytest.mark.slow
@pytest.mark.parametrize("fmt", ALL_FORMATS)
@pytest.mark.parametrize("mode", ["protein_unique", "taxon_unique", "razor"])
def test_parallel_equals_serial_across_every_format_and_mode(tmp_path, fmt,
                                                             mode):
    par = _project_for(tmp_path / f"p_{fmt}_{mode}", fmt,
                       peptide_assignment=mode)
    par.run()
    ser = _project_for(tmp_path / f"s_{fmt}_{mode}", fmt,
                       peptide_assignment=mode)
    ser.run("--serial")
    assert _digest(par) == _digest(ser)


# --- finding 47 -------------------------------------------------------
def _shuffle_fasta(path, seed=7):
    recs = []
    cur = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">") and cur:
            recs.append(cur)
            cur = []
        cur.append(line)
    if cur:
        recs.append(cur)
    random.Random(seed).shuffle(recs)
    with open(path, "w", encoding="utf-8") as fh:
        for r in recs:
            fh.writelines(r)


def _shuffle_rows(path, seed=7, sep="\t"):
    lines = open(path, encoding="utf-8").read().splitlines()
    head, body = lines[0], lines[1:]
    random.Random(seed).shuffle(body)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join([head] + body) + "\n")


def _sorted_digest(proj):
    """Content digest with every table sorted, so a legitimate change of row
    ORDER is not counted as a change of result."""
    h = hashlib.sha1()
    for rel in RESULT_FILES:
        path = proj.rpath(*rel.split("/"))
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        df = df.sort_values(list(df.columns)).reset_index(drop=True)
        h.update(rel.encode())
        h.update(df.to_csv(index=False).encode())
    return h.hexdigest()


def test_shuffling_the_fasta_record_order_changes_nothing(tmp_path):
    ref = build_project(tmp_path / "ref")
    ref.run()
    alt = build_project(tmp_path / "alt")
    _shuffle_fasta(alt.path("input", "proteins.faa"))
    alt.run()
    assert _sorted_digest(alt) == _sorted_digest(ref)


def test_shuffling_the_peptide_row_order_changes_nothing(tmp_path):
    ref = build_project(tmp_path / "ref")
    ref.run()
    alt = build_project(tmp_path / "alt")
    _shuffle_rows(alt.path("input", "combined_peptide.tsv"))
    alt.run()
    assert _sorted_digest(alt) == _sorted_digest(ref)


def test_shuffling_the_manifest_row_order_changes_nothing(tmp_path):
    ref = build_project(tmp_path / "ref")
    ref.run()
    alt = build_project(tmp_path / "alt")
    path = alt.path("input", "experiment.fp-manifest")
    lines = open(path, encoding="utf-8").read().splitlines()
    random.Random(3).shuffle(lines)
    open(path, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    alt.run()
    assert _sorted_digest(alt) == _sorted_digest(ref)


# --- finding 48 -------------------------------------------------------
def test_five_repeated_runs_give_one_hash(tmp_path):
    seen = set()
    for i in range(5):
        proj = build_project(tmp_path / f"run{i}")
        proj.run()
        seen.add(_digest(proj))
    assert len(seen) == 1, f"{len(seen)} distinct results from 5 runs"


def test_repeated_runs_agree_across_hash_seeds(tmp_path):
    a = build_project(tmp_path / "h0")
    a.run(env={"PYTHONHASHSEED": "0"})
    b = build_project(tmp_path / "h1")
    b.run(env={"PYTHONHASHSEED": "12345"})
    assert _digest(a) == _digest(b)


# --- finding 49 -------------------------------------------------------
def test_a_second_run_changes_no_output_and_runs_no_stage(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    before = _digest(proj)
    mtimes = {rel: os.path.getmtime(proj.rpath(*rel.split("/")))
              for rel in RESULT_FILES
              if os.path.exists(proj.rpath(*rel.split("/")))}
    proc = proj.run()
    assert "done: 0 run, 0 adopted" in proc.stderr
    assert _digest(proj) == before
    for rel, t in mtimes.items():
        assert os.path.getmtime(proj.rpath(*rel.split("/"))) == t, \
            f"{rel} was rewritten by a run that reported doing nothing"


def test_a_third_run_is_still_a_no_op(tmp_path, stub_bin):
    proj = _searchable(tmp_path, tmp_path / "p")
    proj.run()
    proj.run()
    proc = proj.run()
    assert "done: 0 run" in proc.stderr


# --- text encoding ----------------------------------------------------
def test_every_text_output_is_utf8(tmp_path):
    # symptom: text files were written with the platform codec, so on Windows
    # the report died on a Unicode character and left an empty file.
    ps = F.protein_set()
    ps[0].description = "β-glucosidase — Bacteroides sp."
    proj = build_project(tmp_path / "utf8", proteins=ps)
    proj.run()
    for root, _dirs, files in os.walk(proj.results):
        for f in files:
            if not f.endswith((".tsv", ".faa", ".txt", ".log")):
                continue
            open(os.path.join(root, f), encoding="utf-8").read()
