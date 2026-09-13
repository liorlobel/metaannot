"""Quantification input: manifests, identifiers, zeros, and the shared-peptide
rule that decides which peptide quantifies which protein."""
from __future__ import annotations

import os
import re
import shutil

import numpy as np
import pandas as pd
import pytest

import fixtures as F
from conftest import build_project, run_metaannot


def _cfg(ma, **over):
    import json
    c = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    c.update(over)
    return c


# --- finding 30 -------------------------------------------------------
def test_fragpipe_zeros_are_read_as_missing_by_default(ma, tmp_path, capsys):
    # symptom: FragPipe writes 0 for "not quantified", not for "measured as
    # zero". Summed as a real zero it turns missingness into fold change.
    ps = F.protein_set()[:2]
    path = str(tmp_path / "combined_peptide.tsv")
    F.write_peptide_table(path, ps, ["S1", "S2"],
                          zeros=[(0, "S1"), (1, "S2"), (2, "S1")])
    feats, int_cols, _ = ma.read_feature_table(path, "fragpipe_peptide",
                                               _cfg(ma))
    assert feats["S1 Intensity"].isna().sum() == 2
    assert feats["S2 Intensity"].isna().sum() == 1
    err = capsys.readouterr().err
    assert "3 intensity cell(s) are 0" in err
    assert "missingness into fold change" in err


def test_zero_intensity_is_missing_false_restores_summing(ma, tmp_path):
    ps = F.protein_set()[:2]
    path = str(tmp_path / "combined_peptide.tsv")
    F.write_peptide_table(path, ps, ["S1", "S2"], zeros=[(0, "S1")])
    feats, _, _ = ma.read_feature_table(path, "fragpipe_peptide",
                                        _cfg(ma, zero_intensity_is_missing=False))
    assert feats["S1 Intensity"].isna().sum() == 0
    assert float(feats["S1 Intensity"].iloc[0]) == 0.0


def test_a_missing_feature_is_not_summed_as_a_zero(ma, tmp_path):
    # the consequence: a protein whose only feature is missing in a sample
    # must be NA there, not 0.
    ps = [F.Protein("P1", "MKV" * 40, seed_taxid="820")]
    rows = [{"peptide": "PEPTIDEK", "razor": "P1", "candidates": ["P1"]}]
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, ps, ["S1", "S2"], rows=rows,
                          zeros=[(0, "S1")])
    feats, int_cols, _ = ma.read_feature_table(path, "fragpipe_peptide",
                                               _cfg(ma))
    q, ev, _ = ma.rollup_features(feats, int_cols, {"P1": "820"},
                                  "taxon_unique", 1)
    assert pd.isna(q.loc[0, "S1 Intensity"])
    assert q.loc[0, "S2 Intensity"] > 0


# --- finding 31 -------------------------------------------------------
def test_a_fractionated_manifest_collapses_to_one_row_per_sample(ma, tmp_path,
                                                                 capsys):
    # symptom: a manifest is one row per raw FILE. A fractionated acquisition
    # repeats the sample name across its fraction rows, which is exactly how
    # FragPipe denotes fractions — and it used to be refused as a duplicate.
    path = F.simple_manifest(str(tmp_path / "e.fp-manifest"),
                             ["A_1", "A_2", "B_1", "B_2"], fractions=4)
    m = ma.read_manifest(path)
    assert len(m) == 4
    assert sorted(m["sample"]) == ["A_1", "A_2", "B_1", "B_2"]
    assert set(m["n_fractions"]) == {4}
    err = capsys.readouterr().err
    assert "4 sample(s) are split across multiple fraction files" in err


def test_fraction_rows_that_disagree_about_the_design_are_refused(ma, tmp_path):
    # collapsing rows that disagree would invent a sample that never existed.
    path = F.write_manifest(str(tmp_path / "e.fp-manifest"), [
        ("/raw/a_f1.mzML", "A", "1", "DDA"),
        ("/raw/a_f2.mzML", "A", "1", "DIA"),
    ])
    with pytest.raises(ma.StageError) as e:
        ma.read_manifest(path)
    assert "different data_type" in str(e.value)


def test_two_different_samples_resolving_to_one_quant_column_are_refused(ma):
    # symptom: last one wins used to be silent, producing a design with fewer
    # samples than the manifest and the wrong replicate labels.
    m = pd.DataFrame([
        {"file": "/raw/x.mzML", "experiment": "A", "bioreplicate": "1",
         "basename": "x", "sample": "A_1"},
        {"file": "/raw/y.mzML", "experiment": "A", "bioreplicate": "",
         "basename": "A_1", "sample": "A"},
    ])
    with pytest.raises(ma.StageError) as e:
        ma.map_manifest_to_columns(m, ["A_1 Intensity"], " Intensity")
    assert "several samples map to the same quant column" in str(e.value)


def test_a_windows_manifest_path_is_split_on_backslashes(ma):
    # symptom: os.path.basename on Linux does not split 'D:\runs\EX1.mzML', so
    # the whole path became the sample name on the analysis machine.
    assert ma.manifest_basename(r"D:\runs\EX1.mzML") == "EX1"
    assert ma.manifest_basename("/raw/EX1.mzML") == "EX1"


def test_a_single_space_is_not_a_manifest_separator(ma, tmp_path):
    # symptom: raw-file paths routinely contain spaces; guessing there
    # silently produced a design with one nameless group.
    path = tmp_path / "e.fp-manifest"
    path.write_text("/raw/my file.mzML\n", encoding="utf-8")
    with pytest.raises(ma.StageError) as e:
        ma.read_manifest(str(path))
    assert "only one field" in str(e.value)


def test_a_fractionated_run_produces_the_same_quant_as_an_unfractionated_one(
        tmp_path):
    # end to end: fractions must change the design, not the numbers.
    plain = build_project(tmp_path / "plain", fractions=1)
    frac = build_project(tmp_path / "frac", fractions=3)
    plain.run()
    frac.run()
    a = pd.read_csv(plain.rpath("quant", "annotated_quant.tsv"), sep="\t")
    b = pd.read_csv(frac.rpath("quant", "annotated_quant.tsv"), sep="\t")
    pd.testing.assert_frame_equal(a, b)


# --- the two invisibility rates, and the subtraction between them --------
# symptom: the pipeline has printed the database rate (from finalise) and the
# quantified rate (from the join) for as long as it has had bins, on separate
# pages of a 9,000-line log, and never once said that they differ. On the
# first full real run it was 43.6% against 29.4% -- the KO-less fraction is
# over-represented among the proteins that were actually expressed, which is
# the claim this tool exists to make.


def test_the_join_states_the_gap_between_the_two_invisibility_rates(tmp_path):
    """Quantify only the KO-less half, so the ratio cannot be read backwards.

    The knitted report's own version of this assertion runs on a fixture where
    every database protein is also quantified, so the two rates are equal and
    the ratio is 1.00x -- which passes whichever way round it is divided. Here
    the populations are made to differ on purpose.
    """
    import fixtures as F
    proj = build_project(tmp_path / "gap")
    koless = [x for x in proj.proteins if x.expect_bin != "1_ko_pathway"]
    assert 0 < len(koless) < len(proj.proteins), "the fixture stopped skewing"
    F.write_peptide_table(proj.path("input", "combined_peptide.tsv"),
                          koless, proj.samples)
    err = proj.run().stderr

    ann = pd.read_csv(proj.rpath("annotation_final.tsv"), sep="\t", dtype=str)
    db = ann["bin"].dropna().astype(str)
    db_pct = 100 * (db != "1_ko_pathway").mean()
    line = [l for l in err.splitlines() if "invisible to KEGG" in l and "join:" in l]
    assert len(line) == 1, line
    line = line[0]
    assert "100.0% of the" in line, \
        f"every quantified protein here is KO-less: {line}"
    assert f"against {db_pct:.1f}% of the {len(db):,} protein(s) in the " \
           f"search database" in line, line
    assert f"{100.0 / db_pct:.2f}x as likely" in line, \
        f"the ratio is quantified-over-database, not the other way: {line}"
    assert db_pct < 100.0, "the two rates must differ or this proves nothing"


# --- the funnel: 455,571 -> 1,282, in one file rather than a dozen lines ---
# symptom: the first full real run narrowed by two and a half orders of
# magnitude and said why across a dozen lines of an 8,990-line log, in three
# stages, some of them counting features and some counting proteins.


def _funnel(proj):
    return pd.read_csv(proj.rpath("quant", "quant_funnel.tsv"), sep="\t",
                       dtype=str).fillna("")


def test_the_quant_funnel_is_written_beside_the_table_it_explains(
        ma, tmp_path):
    proj = build_project(tmp_path / "funnel", fractions=1)
    proj.run()
    f = _funnel(proj)
    assert list(f.columns) == list(ma._QuantFunnel.COLUMNS)
    rows = len(pd.read_csv(proj.rpath("quant", "annotated_quant.tsv"),
                           sep="\t"))
    assert f["step"].iloc[-1] == "rows written to annotated_quant.tsv"
    assert int(f["after"].iloc[-1]) == rows, \
        "the funnel's last count is the file it sits next to, or it is fiction"
    # What this cannot discriminate, said rather than implied: the last row
    # counts `out` and not `q` because the merges above are left joins and a
    # duplicate key on a right-hand side ADDS rows. On this fixture nothing
    # duplicates, so the two are equal and swapping them passes here. The
    # assertion is still the right one — it fails on any run where they do
    # differ — but it is not evidence that the right count was chosen.


def test_the_quant_funnel_reconciles_within_a_unit_and_never_across_one(
        tmp_path):
    """The arithmetic, and the place the arithmetic is not allowed to run.

    A funnel that chained a feature count straight into a protein count would
    read as one number shrinking while being a different claim at every step.
    So a row's `before` is the last `after` recorded FOR ITS OWN UNIT, not the
    row above it — this stage interrupts the protein chain with the feature
    rows in the middle of it — and the first row of each unit has no `before`
    at all.
    """
    proj = build_project(tmp_path / "reconcile", fractions=1)
    proj.run()
    f = _funnel(proj)
    assert len(f) >= 4, f.to_string()
    last, starts = {}, 0
    for _, r in f.iterrows():
        unit = r["unit"]
        if r["before"] == "":
            starts += 1
            assert unit not in last, \
                f"a second chain for {unit!r} silently breaks the first: {dict(r)}"
            assert r["dropped"] == "", \
                f"a row with no `before` cannot have dropped anything: {dict(r)}"
        else:
            assert int(r["before"]) == last[unit], \
                f"{dict(r)} does not start where this unit was left ({last[unit]})"
            assert (int(r["dropped"])
                    == int(r["before"]) - int(r["after"])), dict(r)
        last[unit] = int(r["after"])
    assert starts == len(set(f["unit"])) == 2, \
        "one chain start per unit, and this fixture counts proteins and features"
    # the interruption is the point: the protein chain has to step OVER the
    # feature rows rather than restart after them
    units = list(f["unit"])
    assert units.index("feature") > 0 and units[-1] == "protein"
    assert units.count("protein") >= 4 and "feature" in units[1:-1]


def test_the_funnel_prices_a_filter_that_really_removes_something(tmp_path):
    # The fixture drops nothing anywhere, so every count above reconciles at
    # zero and would reconcile at zero against a funnel that had subtracted
    # the wrong pair. min_features_per_protein is the one filter a config can
    # make bite without a second fixture.
    proj = build_project(tmp_path / "bites", fractions=1)
    proj.write_config(min_features_per_protein=99)
    proj.run()
    f = _funnel(proj)
    row = f[f["step"] == "min_features_per_protein=99"].iloc[0]
    assert int(row["dropped"]) == int(row["before"]) - int(row["after"])
    assert int(row["dropped"]) > 0, \
        "a threshold of 99 features per protein removed nothing at all"
    assert int(row["after"]) == 0
    assert "absent from" in row["why"]
    # and the row after it still starts where this one left the chain
    nxt = f.iloc[f.index.get_loc(row.name) + 1]
    assert int(nxt["before"]) == int(row["after"])


def test_the_funnel_prices_min_features_per_protein_even_at_the_default(
        tmp_path):
    # A funnel silent about a filter is read as a funnel with no such filter,
    # and this is the one a reader most wants priced before they raise it.
    proj = build_project(tmp_path / "priced", fractions=1)
    proj.run()
    f = _funnel(proj)
    hit = f[f["step"].str.startswith("min_features_per_protein=")]
    assert len(hit) == 1, f["step"].tolist()
    assert "removes nothing" in hit["why"].iloc[0]
    assert int(hit["dropped"].iloc[0]) == 0


def test_the_funnel_says_that_the_report_narrows_again_after_it(tmp_path):
    # It ends at annotated_quant.tsv, which is the middle of the narrowing and
    # not the end of it: min_valid_per_group and analysis.min_plexes run in R
    # over that file. A funnel that stopped without saying so would be read as
    # the whole story.
    proj = build_project(tmp_path / "downstream", fractions=1)
    err = proj.run().stderr
    assert "min_valid_per_group" in err and "quant_funnel.tsv" in err


def test_the_written_quant_table_carries_the_dominance_column_the_report_reads(
        ma, tmp_path):
    """The column the report's dominance line is computed from, on disk.

    Everything else about that line is pinned on frames built inside the
    tests, and the report's own half is pinned against an R literal - so a
    rename on the Python side would leave the report reading a column that is
    no longer written, which is the exact defect this change was made to fix,
    with every test still green. This is the one assertion that couples the
    two: the file a run really writes has to carry the name the report really
    asks for.

    Read off the R source rather than repeated here, so that the coupling is
    to the report and not to a string in this file.
    """
    proj = build_project(tmp_path / "written", fractions=1)
    proj.run()
    cols = pd.read_csv(proj.rpath("quant", "annotated_quant.tsv"),
                       sep="\t", nrows=0).columns
    assert "taxon_unique_dominated" in cols, sorted(cols)
    assert "taxon_unique_dominated" in ma.RMD_TEMPLATE, \
        "the report no longer reads the column this test exists to couple"


def test_a_manifest_run_matching_no_column_names_the_columns_present(ma,
                                                                     tmp_path):
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, F.protein_set()[:2], ["S1", "S2"])
    man = F.simple_manifest(str(tmp_path / "e.fp-manifest"), ["Z_9"])
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(path, "fragpipe_peptide",
                              _cfg(ma, manifest=man))
    assert "match no column" in str(e.value)
    assert "columns present" in str(e.value)


# --- finding 32 -------------------------------------------------------
def test_a_protein_id_carrying_a_description_is_reduced_to_its_first_token(
        ma, tmp_path, capsys):
    # symptom: FragPipe fills 'Protein ID' with '<id> <description>' for a
    # metagenome database, a value that can never match the fasta.
    ps = F.protein_set()[:3]
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, ps, ["S1"], protein_id_with_desc=True)
    feats, _, _ = ma.read_feature_table(path, "fragpipe_peptide", _cfg(ma))
    assert set(feats["razor_protein"]) <= {p.pid for p in ps}
    assert not any(" " in x for x in feats["razor_protein"])


def test_the_run_reports_how_many_values_carried_a_description(ma, tmp_path,
                                                               capsys):
    ps = F.protein_set()[:2]
    path = str(tmp_path / "q.tsv")
    # drop the bare 'Protein' column so the reader falls back to 'Protein ID'
    df = pd.read_csv(F.write_peptide_table(path, ps, ["S1"]), sep="\t")
    df = df.drop(columns=["Protein"])
    df.to_csv(path, sep="\t", index=False)
    ma.read_feature_table(path, "fragpipe_peptide", _cfg(ma))
    err = capsys.readouterr().err
    assert "contain a space" in err
    assert "using the first token" in err


def test_the_protein_column_is_preferred_over_protein_id(ma, tmp_path):
    # 'Protein' is the fasta header; 'Protein ID' is FragPipe's decorated one.
    ps = F.protein_set()[:2]
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, ps, ["S1"])
    feats, _, _ = ma.read_feature_table(path, "fragpipe_peptide", _cfg(ma))
    assert feats["razor_protein"].iloc[0] == ps[0].pid


def test_first_token_leaves_a_bare_identifier_alone(ma):
    assert ma.first_token("CDPNAMPK_339076 hypothetical protein") == \
        "CDPNAMPK_339076"
    assert ma.first_token("P12345") == "P12345"
    assert ma.first_token("") == ""


def test_decoy_and_contaminant_rows_are_dropped_and_counted(ma, tmp_path,
                                                            capsys):
    ps = F.protein_set()[:2] + [
        F.Protein("rev_P1", "MKV" * 30), F.Protein("CON__TRYP", "MKV" * 30)]
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, ps, ["S1"])
    feats, _, _ = ma.read_feature_table(path, "fragpipe_peptide", _cfg(ma))
    assert not any(x.startswith(("rev_", "CON__"))
                   for x in feats["razor_protein"])
    assert "decoy/contaminant row(s) dropped" in capsys.readouterr().err


def test_a_decoy_candidate_does_not_veto_a_shared_peptide(ma, tmp_path):
    # symptom: a decoy candidate can never carry eggNOG taxonomy, so leaving
    # it in the candidate list made every peptide it touches
    # shared_unknown_taxon — vetoing real quantification.
    ps = [F.Protein("P1", "MKV" * 40, seed_taxid="820"),
          F.Protein("P2", "MKV" * 40, seed_taxid="820")]
    rows = [{"peptide": "PEPTIDEK", "razor": "P1", "candidates": ["P1"]}]
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, ps, ["S1"], rows=rows,
                          shared=[(0, ["P2", "rev_P9"])])
    feats, _, _ = ma.read_feature_table(path, "fragpipe_peptide", _cfg(ma))
    assert "rev_P9" not in feats["candidates"].iloc[0]
    assert set(feats["candidates"].iloc[0]) == {"P1", "P2"}


# --- finding 37 -------------------------------------------------------
def _shared_features(ma, tmp_path, name="q.tsv"):
    """Four features: unique, taxon-unique, cross-taxon shared, and shared
    with a candidate that has no taxonomy at all."""
    ps = [F.Protein(f"P{i}", "MKV" * 40) for i in range(1, 6)]
    rows = [
        {"peptide": "UNIQUEK", "razor": "P1", "candidates": ["P1"]},
        {"peptide": "SAMETAXK", "razor": "P1", "candidates": ["P1"]},
        {"peptide": "CROSSTAXK", "razor": "P1", "candidates": ["P1"]},
        {"peptide": "NOTAXK", "razor": "P1", "candidates": ["P1"]},
    ]
    path = str(tmp_path / name)
    F.write_peptide_table(path, ps, ["S1", "S2"], rows=rows,
                          shared=[(1, ["P2"]), (2, ["P3"]), (3, ["P4"])])
    feats, int_cols, _ = ma.read_feature_table(path, "fragpipe_peptide",
                                               _cfg(ma))
    taxon_of = {"P1": "820", "P2": "820", "P3": "999", "P4": ""}
    return feats, int_cols, taxon_of


def test_retention_is_monotonic_across_the_assignment_modes(ma, tmp_path):
    # razor >= taxon_unique >= protein_unique, always.
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    kept = {}
    for mode in ("protein_unique", "taxon_unique", "razor"):
        _, ev, fc = ma.rollup_features(feats, int_cols, taxon_of, mode, 1)
        kept[mode] = int(fc["_assigned"].notna().sum())
    assert kept["razor"] >= kept["taxon_unique"] >= kept["protein_unique"]
    assert kept["protein_unique"] == 1
    assert kept["taxon_unique"] == 2
    assert kept["razor"] == 4


def test_razor_conserves_total_intensity(ma, tmp_path):
    # razor keeps every feature, so nothing may be lost in the roll-up.
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    q, _, _ = ma.rollup_features(feats, int_cols, taxon_of, "razor", 1)
    for c in int_cols:
        assert q[c].sum() == pytest.approx(feats[c].sum())


def test_a_candidate_with_no_taxonomy_is_never_taxon_unique(ma, tmp_path):
    # symptom: treating it as taxon-unique would readmit exactly the
    # unannotated proteins this pipeline exists to scrutinise.
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    _, _, fc = ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    klass = dict(zip(fc["peptide"], fc["_class"]))
    assert klass["UNIQUEK"] == "unique"
    assert klass["SAMETAXK"] == "taxon_unique"
    assert klass["CROSSTAXK"] == "shared"
    assert klass["NOTAXK"] == "shared_unknown_taxon"
    assert "taxon_unique" not in klass["NOTAXK"]


def test_the_family_fallback_is_opt_in_and_named_separately(ma, tmp_path):
    # a family is a SEQUENCE CLUSTER, not an organism, so the class it
    # produces must stay visible rather than being called taxon_unique.
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    family_of = {"P1": "fam1", "P2": "fam1", "P3": "fam2", "P4": "fam1"}
    _, ev, fc = ma.rollup_features(feats, int_cols, taxon_of,
                                   "taxon_or_family_unique", 1,
                                   family_of=family_of)
    klass = dict(zip(fc["peptide"], fc["_class"]))
    assert klass["NOTAXK"] == "family_unique"
    assert klass["CROSSTAXK"] == "shared"


def test_peptide_evidence_records_the_counts_per_protein(ma, tmp_path):
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    _, ev, _ = ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    row = ev.set_index("protein_id").loc["P1"]
    assert row["n_unique"] == 1
    assert row["n_taxon_unique"] == 1
    assert row["n_features_dropped"] == 2
    assert bool(row["taxon_unique_dominated"]) is False


def _dominance_frame(ma, tmp_path, n_dominated, n_clean, n_dropped,
                     name="dom.tsv"):
    """A frame whose two candidate denominators differ, on purpose.

    `n_dominated` proteins carry one unique feature and two taxon-unique ones
    (flagged); `n_clean` carry one unique feature (quantified, not flagged);
    `n_dropped` carry nothing but cross-taxon shared features, so every one of
    their features is dropped and they arrive in peptide_evidence.tsv as a row
    with n_features_used == 0 - quantified nowhere, and unflaggable whatever
    the assignment rule decided.
    """
    ps, rows, shared, taxon_of = [], [], [], {}
    i = 0
    for k in range(n_dominated + n_clean):
        pid = f"Q{k}"
        ps.append(F.Protein(pid, "MKV" * 40)); taxon_of[pid] = str(800 + k)
        rows.append({"peptide": f"UQ{k}AAAAK", "razor": pid,
                     "candidates": [pid]})
        i += 1
        if k < n_dominated:
            for j in range(2):
                nb = f"N{k}x{j}"
                ps.append(F.Protein(nb, "MKV" * 40)); taxon_of[nb] = str(800 + k)
                rows.append({"peptide": f"TX{k}x{j}AAAK", "razor": pid,
                             "candidates": [pid]})
                shared.append((i, [nb])); i += 1
    for k in range(n_dropped):
        pid, other = f"Z{k}", f"W{k}"
        ps.append(F.Protein(pid, "MKV" * 40)); taxon_of[pid] = "700"
        ps.append(F.Protein(other, "MKV" * 40)); taxon_of[other] = "701"
        rows.append({"peptide": f"ZX{k}AAAAK", "razor": pid,
                     "candidates": [pid]})
        shared.append((i, [other])); i += 1
    path = str(tmp_path / name)
    F.write_peptide_table(path, ps, ["S1", "S2"], rows=rows, shared=shared)
    feats, int_cols, _ = ma.read_feature_table(path, "fragpipe_peptide",
                                               _cfg(ma))
    return feats, int_cols, taxon_of


def _dom_line(err):
    """The one logged line that states the dominance rate."""
    hits = [ln for ln in err.splitlines() if "rest more on shared-but-taxon-" in ln]
    assert len(hits) <= 1, f"more than one dominance line: {hits}"
    return hits[0] if hits else ""


def test_the_dominance_rate_is_over_the_proteins_that_could_carry_the_flag(
        ma, tmp_path, capsys):
    # symptom: the denominator was len(ev), and ev comes from an OUTER join
    # with the dropped-feature counts, so every protein whose features were
    # ALL dropped sat in it as a row of zeros. (0 + 0) > 0 is False, so such a
    # row could never be flagged while still enlarging the denominator - on
    # the first full real run 5,039 of 8,238 rows, which reported a 54.7%
    # finding as 21.2%.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 6)
    _, ev, _ = ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    assert len(ev) == 18                      # 12 assessable + 6 with nothing
    line = _dom_line(capsys.readouterr().err)
    assert "7/12 (58.3%)" in line
    assert "/18" not in line.split("peptide_evidence.tsv has")[0]
    assert "with at least one assigned feature" in line


def test_a_protein_with_every_feature_dropped_is_not_in_the_dominance_denominator(
        ma, tmp_path, capsys):
    # the mechanism, stated on its own: adding rows that hold no measurement
    # must not move the rate. It is a positivity rate over patients who were
    # never tested.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 0,
                                                 name="a.tsv")
    ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    none_dropped = _dom_line(capsys.readouterr().err)
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 40,
                                                 name="b.tsv")
    ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    many_dropped = _dom_line(capsys.readouterr().err)
    assert "7/12 (58.3%)" in none_dropped
    assert "7/12 (58.3%)" in many_dropped


def test_the_dominance_line_names_the_rows_it_left_out_of_its_denominator(
        ma, tmp_path, capsys):
    # the duty owed for narrowing a denominator: say how many rows were taken
    # out of it, in the sentence that used it, so nobody recomputes the rate
    # off peptide_evidence.tsv and gets a different number.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 6)
    ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    line = _dom_line(capsys.readouterr().err)
    assert "peptide_evidence.tsv has 18 rows" in line
    assert "the other 6 had every feature dropped" in line
    assert "could never be flagged" in line


def test_the_proteins_with_no_assigned_feature_are_counted_in_their_own_line(
        ma, tmp_path, capsys):
    # nothing in this pipeline used to state this count - the features line
    # counts FEATURES, and the min_features line has already excluded these
    # proteins - so a reader filtering peptide_evidence.tsv had no way to know
    # a third of it held no number at all.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 6)
    ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    err = capsys.readouterr().err
    hits = [ln for ln in err.splitlines()
            if "protein(s) had every feature dropped under" in ln]
    assert len(hits) == 1
    assert "6/18 protein(s)" in hits[0]
    assert "under 'taxon_unique'" in hits[0]
    assert "no number in annotated_quant.tsv" in hits[0]
    # INFO and never WARN: in a strain-redundant database this is the rule the
    # user chose doing what it says, and a WARN on every real run is the line
    # a reader learns to skip.
    assert hits[0].split("]")[1].strip().startswith("INFO")


def test_a_run_that_drops_nothing_keeps_the_old_dominance_denominator(
        ma, tmp_path, capsys):
    # the property most worth having: where the defect does not exist the
    # number does not move. With nothing dropped, len(ev) IS the assessable
    # count, and the trailing sentence is not emitted at all.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 0)
    _, ev, _ = ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    line = _dom_line(capsys.readouterr().err)
    assert f"7/{len(ev)} (58.3%)" in line
    assert "peptide_evidence.tsv has" not in line


def test_the_dominance_denominator_is_the_min_features_line_s_own_denominator(
        ma, tmp_path, capsys):
    # the argument that settled which population to report over: `before` in
    # the min_features WARN is len(quant), i.e. the proteins with at least one
    # assigned feature. Two adjacent lines of one stage must not print
    # different denominators with nothing saying they are different
    # populations.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 6)
    ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 2)
    err = capsys.readouterr().err
    retained = [ln for ln in err.splitlines() if "proteins retained with" in ln]
    assert len(retained) == 1
    # Parsed out of both lines rather than written twice here: the invariant
    # is that they ARE one denominator, printed one way - "3,199" beside
    # "3199" in one funnel reads as two populations.
    mine = re.search(r"(\d[\d,]*)/(\d[\d,]*) \(", _dom_line(err))
    theirs = re.search(r"(\d[\d,]*)/(\d[\d,]*) proteins retained", retained[0])
    assert mine and theirs
    assert mine.group(2) == theirs.group(2) == "12"
    assert mine.group(1) == "7" and theirs.group(1) == "7"


def test_a_handful_of_assessable_proteins_states_counts_without_a_percentage(
        ma, tmp_path, capsys):
    # below ASSESSABLE_MIN_N the counts are still printed; the quotable
    # fraction is not, because a percentage over a handful of proteins is what
    # gets pasted into a methods section as a claim about a population.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 3, 4, 6)
    ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    line = _dom_line(capsys.readouterr().err)
    assert "3/7 protein(s)" in line
    assert "%" not in line
    assert "too few to state as a rate" in line
    assert line.split("]")[1].strip().startswith("INFO")


def test_enough_assessable_proteins_raises_the_dominance_rate_to_a_warning(
        ma, tmp_path, capsys):
    # the other side of the same floor, so the tier boundary is pinned from
    # both directions rather than assumed.
    assert ma.ASSESSABLE_MIN_N == 10
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 5, 5, 6)
    ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    line = _dom_line(capsys.readouterr().err)
    assert "5/10 (50.0%)" in line
    assert line.split("]")[1].strip().startswith("WARN")


def test_no_taxonomy_keeps_the_dominance_line_silent_rather_than_reporting_zero(
        ma, tmp_path, capsys):
    # under protein_unique - which is also where stage_join falls back when no
    # seed_taxid exists - both taxon columns are structurally 0, so the flag
    # can never be True. A line reading "0/12 (0.0%)" on every such run is the
    # line that fires on everything and takes the real one with it.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 6)
    _, ev, _ = ma.rollup_features(feats, int_cols, taxon_of, "protein_unique",
                                  1)
    assert not bool(ev["taxon_unique_dominated"].any())
    assert _dom_line(capsys.readouterr().err) == ""


def test_a_run_where_nothing_is_assigned_states_no_rate_and_does_not_divide(
        ma, tmp_path, capsys):
    # the degenerate case the guard already covers, pinned rather than
    # assumed: with every feature dropped the assessable count is 0, and a
    # percentage computed before the guard would turn a reporting line into a
    # ZeroDivisionError that takes the whole join stage with it.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 0, 0, 6)
    q, ev, _ = ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    assert len(q) == 0
    assert int((ev["n_features_used"] == 0).sum()) == len(ev) == 6
    err = capsys.readouterr().err
    assert _dom_line(err) == ""
    assert "6/6 protein(s) had every feature dropped" in err


def test_every_row_flagged_taxon_unique_dominated_rests_more_on_shared_features(
        ma, tmp_path):
    # the part that was already right, pinned so a fix to the SUMMARY cannot
    # quietly widen the per-protein classification. The flag is written into
    # peptide_evidence.tsv and annotated_quant.tsv and is true of exactly the
    # rows whose taxon- plus family-unique features outnumber their own unique
    # ones - a tie is False, because "rests MORE on" is strict.
    feats, int_cols, taxon_of = _dominance_frame(ma, tmp_path, 7, 5, 6)
    _, ev, _ = ma.rollup_features(feats, int_cols, taxon_of, "taxon_unique", 1)
    flagged = ev["n_taxon_unique"] + ev["n_family_unique"] > ev["n_unique"]
    assert list(ev["taxon_unique_dominated"]) == list(flagged)
    assert int(flagged.sum()) == 7
    dom = ev[ev["taxon_unique_dominated"]]
    assert (dom["n_taxon_unique"] == 2).all() and (dom["n_unique"] == 1).all()
    # and the rows that hold no measurement are flagged by neither.
    empty = ev[ev["n_features_used"] == 0]
    assert len(empty) == 6 and not bool(empty["taxon_unique_dominated"].any())


def test_peptide_evidence_survives_a_run_where_nothing_is_assigned(ma,
                                                                   tmp_path):
    # symptom: with every feature shared across taxa — which a
    # strain-redundant metagenome database really does produce — the evidence
    # table's index came from the dropped side of the join and was named
    # 'razor_protein', so building it raised KeyError('protein_id') instead
    # of writing an empty table and letting the run say so.
    ps = [F.Protein(f"P{i}", "MKV" * 40) for i in range(1, 4)]
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, ps, ["S1", "S2"],
                          rows=[{"peptide": "CROSSTAXK", "razor": "P1"}],
                          shared=[(0, ["P2"])])
    feats, int_cols, _ = ma.read_feature_table(path, "fragpipe_peptide",
                                               _cfg(ma))
    q, ev, fc = ma.rollup_features(feats, int_cols,
                                   {"P1": "820", "P2": "999"},
                                   "taxon_unique", 1)
    assert len(q) == 0
    assert list(ev["protein_id"]) == ["P1"]
    assert int(ev.iloc[0]["n_features_used"]) == 0
    assert int(ev.iloc[0]["n_features_dropped"]) == 1


def test_an_unknown_assignment_mode_is_refused(ma, tmp_path):
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    with pytest.raises(ma.StageError) as e:
        ma.rollup_features(feats, int_cols, taxon_of, "razr", 1)
    assert "peptide_assignment must be one of" in str(e.value)


def test_the_feature_count_filter_reports_what_it_dropped(ma, tmp_path,
                                                          capsys):
    # symptom: proteins removed by min_features_per_protein were absent from
    # every report table with nothing said.
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    q, _, _ = ma.rollup_features(feats, int_cols, taxon_of, "protein_unique", 3)
    assert len(q) == 0
    err = capsys.readouterr().err
    assert "assigned features; the" in err and "dropped are absent" in err


def test_median_polish_preserves_the_summed_total_but_not_the_ratios(ma,
                                                                     tmp_path):
    # opt-in alternative: magnitudes stay comparable with "sum", ratios do not.
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    s, _, _ = ma.rollup_features(feats, int_cols, taxon_of, "razor", 1,
                                 rollup_method="sum")
    m, _, _ = ma.rollup_features(feats, int_cols, taxon_of, "razor", 1,
                                 rollup_method="median_polish")
    assert m[int_cols].sum().sum() == pytest.approx(
        s[int_cols].sum().sum(), rel=1e-6)


def test_an_unknown_rollup_method_is_refused(ma, tmp_path):
    feats, int_cols, taxon_of = _shared_features(ma, tmp_path)
    with pytest.raises(ma.StageError) as e:
        ma.rollup_features(feats, int_cols, taxon_of, "razor", 1,
                           rollup_method="medianpolish")
    assert "rollup_method must be one of" in str(e.value)


# --- finding 38 -------------------------------------------------------
def test_a_fragpipe_tmt_table_is_refused_naming_the_channels(ma, tmp_path):
    # symptom: reporter channels are named 'Intensity <sample>', which the
    # suffix rule cannot see, so the only column it matched was the bare MS1
    # 'Intensity' — one precursor value pooled over every channel. The run
    # exited 0 and reported ONE sample.
    path = str(tmp_path / "peptide.tsv")
    F.write_tmt_peptide_table(path, F.protein_set()[:3],
                              ["126", "127N", "127C", "128N"])
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(path, "fragpipe_peptide", _cfg(ma))
    msg = str(e.value)
    assert "isobaric (TMT/iTRAQ) output" in msg
    assert "Intensity 126" in msg
    assert "does not read reporter-ion" in msg


def test_the_tmt_refusal_reaches_the_cli_as_a_stage_failure(tmp_path):
    proj = build_project(tmp_path / "tmt")
    F.write_tmt_peptide_table(proj.path("input", "combined_peptide.tsv"),
                              proj.proteins, ["126", "127N"])
    proc = proj.run(expect=1)
    assert "isobaric (TMT/iTRAQ) output" in proc.stderr
    assert "Traceback" not in proc.stderr


def test_a_label_free_table_with_a_bare_intensity_column_still_works(ma,
                                                                     tmp_path):
    # the refusal must not fire on label-free input that happens to carry a
    # bare 'Intensity' column alongside real sample columns.
    ps = F.protein_set()[:2]
    path = str(tmp_path / "q.tsv")
    F.write_peptide_table(path, ps, ["S1", "S2"],
                          extra_cols={"Intensity": [1, 2, 3, 4, 5, 6]})
    feats, int_cols, _ = ma.read_feature_table(path, "fragpipe_peptide",
                                               _cfg(ma))
    assert "S1 Intensity" in int_cols and "S2 Intensity" in int_cols




# --- FragPipe TMT: the per-plex reader --------------------------------
def _tmt_run(tmp_path, **kw):
    """Two plexes, four samples, one feature in both and one in TMT1 only.

    The pool sits at a DIFFERENT channel in each plex, which is what a real
    8-plex run does and what no single reference_channel can describe.
    """
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("131C", "Pool01")],
                     [{"peptide": "SHAREDPEPK", "razor": "P_ko_path",
                       "values": {"A1": 100, "A2": 200, "Pool01": 50}},
                      {"peptide": "ONLYONEK", "razor": "P_dark1",
                       "values": {"A1": 0, "A2": 400, "Pool01": 50}}], **kw)
    F.write_tmt_plex(root, "TMT2",
                     [("126", "B1"), ("127N", "B2"), ("131N", "Pool02")],
                     [{"peptide": "SHAREDPEPK", "razor": "P_ko_path",
                       "mapped": ["P_dark1"],
                       "values": {"B1": 300, "B2": 600, "Pool02": 100}}], **kw)
    return root


def _tmt_cfg(ma, root, **tmt):
    cfg = _cfg(ma, quant_table=root, quant_format="fragpipe_tmt")
    # within_plex_normalise defaults to "median", which rescales every cell.
    # The tests below are about which cells exist and where they came from, so
    # they turn it off and read FragPipe's own numbers; the normalisation has
    # its own tests further down.
    cfg["tmt"]["within_plex_normalise"] = "none"
    cfg["tmt"].update(tmt)
    return cfg


def test_tmt_reporter_columns_are_mapped_through_the_plex_annotation(ma,
                                                                     tmp_path):
    # symptom: the reporter columns are named 'Intensity <sample>', the PREFIX
    # form, so the label-free suffix rule saw only the bare MS1 'Intensity'
    # and a TMT run reported one sample. The annotation is what names them.
    root = _tmt_run(tmp_path)
    feats, int_cols, design = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root))
    assert int_cols == ["A1", "A2", "Pool01", "B1", "B2", "Pool02"]
    assert list(design["plex"]) == ["TMT1"] * 3 + ["TMT2"] * 3
    assert list(design["channel"])[:3] == ["126", "127N", "131C"]
    assert set(feats.columns) >= {"feature_id", "peptide", "razor_protein",
                                  "candidates"}


def test_a_tmt_feature_missing_from_a_plex_is_na_not_zero(ma, tmp_path):
    # symptom: only ~45% of ion keys are shared between two plexes, so an
    # outer join that filled 0 for "not identified here" would turn
    # plex-shaped missingness into fold change. FragPipe's own 0 means "not
    # quantified" and is missing too.
    root = _tmt_run(tmp_path)
    feats, int_cols, _ = ma.read_feature_table(root, "fragpipe_tmt",
                                               _tmt_cfg(ma, root))
    row = feats.set_index("feature_id").loc["ONLYONEK_n[230]ONLYONEK_2"]
    assert row["A2"] == 400
    assert pd.isna(row["A1"])                    # a literal 0 in the table
    assert pd.isna(row["B1"]) and pd.isna(row["B2"])   # absent from TMT2


def test_the_tmt_reference_is_resolved_per_plex_by_sample_name(ma, tmp_path):
    # symptom: the pool sits at 131C in six plexes of a real run and at 131N
    # in the other two, so a single reference_channel cannot express the
    # dataset. The stable signal is the sample NAME.
    root = _tmt_run(tmp_path)
    cfg = _tmt_cfg(ma, root, reference_name="Pool*",
                   use_reference_ratios=True)
    feats, int_cols, design = ma.read_feature_table(root, "fragpipe_tmt", cfg)
    assert int_cols == ["A1", "A2", "B1", "B2"]     # the references are gone
    row = feats.set_index("feature_id").loc["SHAREDPEPK_n[230]SHAREDPEPK_2"]
    assert row["A1"] == 2.0 and row["B2"] == 6.0


def test_a_fixed_tmt_reference_channel_that_misses_a_plex_is_refused(ma,
                                                                     tmp_path):
    root = _tmt_run(tmp_path)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt",
                              _tmt_cfg(ma, root, reference_channel="131C"))
    msg = str(e.value)
    assert "TMT2" in msg and "matches 0 of its channels" in msg
    assert "reference_name" in msg


def test_the_tmt_annotation_file_is_named_after_its_plex(ma, tmp_path):
    # symptom: FragPipe writes TMT1/TMT1_annotation.txt, never a plain
    # annotation.txt, so a reader looking for the plain name maps nothing.
    root = _tmt_run(tmp_path)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt",
                              _tmt_cfg(ma, root, annotation="annotation.txt"))
    assert "has no annotation file" in str(e.value)
    assert "TMT1_annotation.txt" in str(e.value)   # what is actually there


def test_a_relative_tmt_annotation_is_read_from_the_plex_not_the_cwd(
        ma, tmp_path, monkeypatch):
    # symptom: a relative tmt.annotation was resolved against the process's
    # working directory, so a same-named file beside the shell won over the
    # plex's own annotation and the reporter columns were mapped through it -
    # silently, because the CHANNELS still line up and only the sample names
    # are wrong.
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1", [("126", "A1"), ("127N", "A2")],
                     [{"peptide": "PEPTIDEK", "razor": "P_ko_path"}],
                     columns_named="channel")
    pdir = os.path.join(root, "TMT1")
    os.makedirs(os.path.join(pdir, "ann"))
    shutil.move(os.path.join(pdir, "TMT1_annotation.txt"),
                os.path.join(pdir, "ann", "TMT1.txt"))
    cwd = tmp_path / "cwd"
    (cwd / "ann").mkdir(parents=True)
    for decoy in (cwd / "TMT1_annotation.txt", cwd / "ann" / "TMT1.txt"):
        decoy.write_text("126 WRONG1\n127N WRONG2\n", encoding="utf-8")
    monkeypatch.chdir(cwd)
    # the pattern form, with a directory in it, and the {plex: path} map form
    for ann in ("ann/{plex}.txt", {"TMT1": "ann/TMT1.txt"}):
        _f, int_cols, _d = ma.read_feature_table(
            root, "fragpipe_tmt", _tmt_cfg(ma, root, annotation=ann))
        assert int_cols == ["A1", "A2"], ann
    # and the decoy is named when the plex really has no such file
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt",
                              _tmt_cfg(ma, root,
                                       annotation="TMT1_annotation.txt"))
    assert "resolved against the plex directory" in str(e.value)


def test_two_tmt_plexes_cannot_claim_the_same_sample_name(ma, tmp_path):
    # symptom: sample names become the columns of the joined matrix, so a
    # collision would silently merge two channels of two plexes.
    root = str(tmp_path / "run")
    for plex in ("TMT1", "TMT2"):
        F.write_tmt_plex(root, plex, [("126", "A1"), ("127N", "A2")],
                         [{"peptide": "PEPTIDEK", "razor": "P_ko_path"}])
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt", _tmt_cfg(ma, root))
    assert "both claim the sample name 'A1'" in str(e.value)


def test_tmt_reporter_columns_that_match_no_annotation_are_refused(ma,
                                                                   tmp_path):
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1", [("126", "A1"), ("127N", "A2")],
                     [{"peptide": "PEPTIDEK", "razor": "P_ko_path"}])
    ann = os.path.join(root, "TMT1", "TMT1_annotation.txt")
    with open(ann, "w", encoding="utf-8") as fh:
        fh.write("126 OTHER1\n127N OTHER2\n")
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt", _tmt_cfg(ma, root))
    msg = str(e.value)
    assert "cannot be mapped to samples" in msg and "Intensity A1" in msg


def test_a_tmt_channel_with_no_sample_is_dropped(ma, tmp_path):
    # symptom: FragPipe names an unassigned channel '<PLEX>_<CHANNEL>'. Its
    # signal is isotope carry-over from its neighbours, not a sample.
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "TMT1_127N")],
                     [{"peptide": "PEPTIDEK", "razor": "P_ko_path"}])
    feats, int_cols, _ = ma.read_feature_table(root, "fragpipe_tmt",
                                               _tmt_cfg(ma, root))
    assert int_cols == ["A1"]
    feats, int_cols, _ = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, drop_empty_channels=False))
    assert int_cols == ["A1", "TMT1_127N"]


def test_min_plexes_filters_on_how_many_plexes_saw_the_feature(ma, tmp_path):
    root = _tmt_run(tmp_path)
    feats, _, _ = ma.read_feature_table(root, "fragpipe_tmt",
                                        _tmt_cfg(ma, root, min_plexes=2))
    assert list(feats["feature_id"]) == ["SHAREDPEPK_n[230]SHAREDPEPK_2"]


def test_the_tmt_reader_never_falls_back_to_the_tmt_report_matrices(ma,
                                                                    tmp_path):
    # symptom: a plex whose per-plex table is missing is a broken run, and
    # quantifying tmt-report/ instead would be undetectable in the output.
    root = _tmt_run(tmp_path)
    os.remove(os.path.join(root, "TMT2", "ion.tsv"))
    F.write_tmt_report_matrix(os.path.join(root, "tmt-report",
                                           "abundance_protein_MD.tsv"),
                              ["A1", "A2", "B1", "B2"])
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt", _tmt_cfg(ma, root))
    assert "has no ion.tsv" in str(e.value)


def test_a_tmt_report_matrix_is_refused_by_name(ma, tmp_path):
    # symptom: it reads perfectly as a wide table, and its values are already
    # log2 and median-centred, so nothing downstream would ever notice.
    path = str(tmp_path / "abundance_protein_MD.tsv")
    F.write_tmt_report_matrix(path, ["A1", "A2"])
    for fmt in ("fragpipe_peptide", "fragpipe_tmt"):
        with pytest.raises(ma.StageError) as e:
            ma.read_feature_table(path, fmt, _cfg(ma))
        assert "tmt-report matrix" in str(e.value)
        assert "ReferenceIntensity" in str(e.value)


def test_the_tmt_flavour_of_msstats_csv_is_refused_by_name(ma, tmp_path):
    # symptom: unquoted commas in Protein.Description killed the C parser
    # with a tokenising error that named neither TMT nor the file, so the
    # recogniser has to run on the header before pandas parses the body.
    path = str(tmp_path / "msstats.csv")
    F.write_tmt_msstats_csv(path)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(path, "msstats_csv", _cfg(ma))
    msg = str(e.value)
    assert "TMT flavour of msstats.csv" in msg and "Channel 126" in msg


def test_the_tmt_peptide_only_reader_covers_every_plex(ma, tmp_path):
    # the taxonomy stages need peptides and candidates, never an intensity,
    # so they must survive a TMT project whose channels the quant reader
    # would refuse.
    root = _tmt_run(tmp_path)
    cfg = _tmt_cfg(ma, root)
    cfg["peptide_only_reader"] = "always"
    feats = ma.peptide_features(cfg, "unipept")
    assert set(feats["peptide"]) == {"SHAREDPEPK", "ONLYONEK"}


# --- TMT: the design, the condition and the reference -----------------
def _tmt_named_run(tmp_path, plexes, seed=5):
    """plexes: {plex: [(channel, sample), ...]}. Every sample sees the shared
    peptide, so the design is what the test is about, not the missingness."""
    root = str(tmp_path / "run")
    for plex, chans in plexes.items():
        F.write_tmt_plex(root, plex, chans,
                         [{"peptide": "SHAREDPEPK", "razor": "P_ko_path"},
                          {"peptide": f"ONLY{plex}K", "razor": "P_dark1"}],
                         seed=seed)
    return root


def test_the_tmt_reference_is_not_a_biological_sample(ma, tmp_path, capsys):
    # symptom: a pooled bridge left in the sample columns acquires a condition
    # in the design, joins a group's mean, and drags the size factors towards
    # a pool that is in every plex by construction.
    root = _tmt_run(tmp_path)
    feats, int_cols, design = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, reference_name="Pool*"))
    assert int_cols == ["A1", "A2", "B1", "B2"]
    assert "Pool01" not in set(design["sample"])
    assert "is_reference" not in design.columns   # nothing left to flag
    err = capsys.readouterr().err
    assert "covariate treatment" in err
    notes = "\n".join(design.attrs["design_notes"])
    assert "reference:" in notes and "covariate" in notes
    assert "TMT1=131C/Pool01" in notes and "TMT2=131N/Pool02" in notes


def test_the_ratio_treatment_reports_the_values_it_deleted(ma, tmp_path,
                                                           capsys):
    # symptom: a feature with no reference value in a plex becomes NA for
    # every channel of that plex. Whether that costs 0.5% of the matrix or
    # most of it is a property of the data, so it has to be counted, not
    # assumed.
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("131C", "Pool01")],
                     [{"peptide": "HASREFK", "razor": "P_ko_path",
                       "values": {"A1": 100, "A2": 200, "Pool01": 50}},
                      # the pool itself was not quantified here: FragPipe's 0
                      {"peptide": "NOREFK", "razor": "P_dark1",
                       "values": {"A1": 300, "A2": 400, "Pool01": 0}}])
    cfg = _tmt_cfg(ma, root, reference_name="Pool*", use_reference_ratios=True)
    feats, int_cols, design = ma.read_feature_table(root, "fragpipe_tmt", cfg)
    row = feats.set_index("feature_id").loc["NOREFK_n[230]NOREFK_2"]
    assert pd.isna(row["A1"]) and pd.isna(row["A2"])
    err = capsys.readouterr().err
    assert "costing 2 of 4 value(s)" in err
    assert "cost 2 of 4 non-reference value(s) (50.00%)" in err
    assert "substantially sparser" in err
    notes = "\n".join(design.attrs["design_notes"])
    assert "ratios" in notes and "2/4 value(s) lost" in notes


def test_a_zero_reference_is_missing_not_a_divide_into_infinity(ma, tmp_path,
                                                                capsys):
    # symptom: with zero_intensity_is_missing false a reference FragPipe wrote
    # as 0 stays a number, so x/0 put +inf into every other channel of the
    # plex while notna() counted each inf as an observed value - the cost line
    # reported nothing lost over a matrix whose infinities survive the roll-up
    # and log2 and turn the taxon size factors into NaN.
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("131C", "Pool01")],
                     [{"peptide": "HASREFK", "razor": "P_ko_path",
                       "values": {"A1": 100, "A2": 200, "Pool01": 50}},
                      {"peptide": "NOREFK", "razor": "P_dark1",
                       "values": {"A1": 300, "A2": 400, "Pool01": 0}}])
    cfg = _tmt_cfg(ma, root, reference_name="Pool*", use_reference_ratios=True)
    cfg["zero_intensity_is_missing"] = False     # the 0 reaches the divide
    feats, int_cols, design = ma.read_feature_table(root, "fragpipe_tmt", cfg)
    assert not np.isinf(feats[int_cols].to_numpy(dtype=float)).any()
    row = feats.set_index("feature_id").loc["NOREFK_n[230]NOREFK_2"]
    assert pd.isna(row["A1"]) and pd.isna(row["A2"])
    err = capsys.readouterr().err
    assert "costing 2 of 4 value(s)" in err      # counted, not silently kept
    assert "cost 2 of 4 non-reference value(s) (50.00%)" in err
    assert "2/4 value(s) lost" in "\n".join(design.attrs["design_notes"])


def test_the_condition_is_derived_from_unambiguous_sample_names(ma, tmp_path,
                                                                capsys):
    # the condition can be in the names. When it is, deriving it saves the
    # user a metadata file - but it is a guess about their experiment, so it
    # is logged as one and recorded as one.
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "resp_1"), ("127N", "nonresp_1")],
        "TMT2": [("126", "resp_2"), ("127N", "nonresp_2")]})
    _, _, design = ma.read_feature_table(root, "fragpipe_tmt",
                                         _tmt_cfg(ma, root))
    assert list(design["group"]) == ["resp", "nonresp", "resp", "nonresp"]
    assert list(design["plex"]) == ["TMT1", "TMT1", "TMT2", "TMT2"]
    err = capsys.readouterr().err
    assert "This is a GUESS" in err
    assert "derived from the sample name before the first '_'" in \
        "\n".join(design.attrs["design_notes"])


def test_sample_names_that_carry_no_condition_leave_the_design_without_one(
        ma, tmp_path, capsys):
    # the real dataset: MF#### codes. Nothing in them is a condition, and the
    # plex is a batch, so there is nothing to fall back on.
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "MF0030"), ("127N", "MF0071")],
        "TMT2": [("126", "MF0084"), ("127N", "MF0202")]})
    _, _, design = ma.read_feature_table(root, "fragpipe_tmt",
                                         _tmt_cfg(ma, root))
    assert "group" not in design.columns
    err = capsys.readouterr().err
    assert "could not be derived from the sample names" in err
    assert "NOT taken from the plex" in err
    assert "not derived" in "\n".join(design.attrs["design_notes"])


def test_the_advice_says_metadata_is_missing_only_when_it_actually_is(
        ma, tmp_path, capsys):
    """The warning used to say "supply it in analysis.metadata" whether or
    not analysis.metadata was set and already named every sample.

    On the real Pittsburgh run it was set and did, and the line read as a
    defect when it was only a division of labour: design_from_input.tsv
    records what the INPUT knew, and the metadata is merged at report time.
    """
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "MF0030"), ("127N", "MF0071")],
        "TMT2": [("126", "MF0084"), ("127N", "MF0202")]})
    md = tmp_path / "meta.tsv"
    md.write_text("sample\tgroup\n"
                  "MF0030\tresponder\n"
                  "MF0071\tnon_responder\n"
                  "MF0084\tresponder\n"
                  "MF0202\tnon_responder\n", encoding="utf-8")
    cfg = _tmt_cfg(ma, root)
    cfg["analysis"] = dict(cfg.get("analysis") or {},
                           metadata=str(md), sample_col="sample")
    ma.read_feature_table(root, "fragpipe_tmt", cfg)
    err = capsys.readouterr().err
    assert "does name every sample" in err
    assert "supply it in analysis.metadata" not in err


def test_the_advice_names_the_samples_the_metadata_leaves_out(ma, tmp_path,
                                                             capsys):
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "MF0030"), ("127N", "MF0071")],
        "TMT2": [("126", "MF0084"), ("127N", "MF0202")]})
    md = tmp_path / "meta.tsv"
    md.write_text("sample\tgroup\n"
                  "MF0030\tresponder\n", encoding="utf-8")
    cfg = _tmt_cfg(ma, root)
    cfg["analysis"] = dict(cfg.get("analysis") or {},
                           metadata=str(md), sample_col="sample")
    ma.read_feature_table(root, "fragpipe_tmt", cfg)
    err = capsys.readouterr().err
    assert "does not name 3 of these samples" in err


def test_the_advice_says_so_when_the_metadata_path_does_not_exist(ma, tmp_path,
                                                                 capsys):
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "MF0030"), ("127N", "MF0071")],
        "TMT2": [("126", "MF0084"), ("127N", "MF0202")]})
    cfg = _tmt_cfg(ma, root)
    cfg["analysis"] = dict(cfg.get("analysis") or {},
                           metadata=str(tmp_path / "absent.tsv"),
                           sample_col="sample")
    ma.read_feature_table(root, "fragpipe_tmt", cfg)
    err = capsys.readouterr().err
    assert "which does not exist" in err


def test_an_ambiguous_split_of_the_sample_names_is_refused_not_picked(ma,
                                                                     tmp_path):
    # 'a-x_1' groups as {a-x, a-y, b-z} on '_' and as {a, b} on '-'. Both are
    # legal splits of every name and they disagree, so which one is the
    # condition is not knowable from the names.
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "a-x_1"), ("127N", "a-x_2"), ("128N", "a-y_1")],
        "TMT2": [("126", "a-y_2"), ("127N", "b-z_1"), ("128N", "b-z_2")]})
    _, _, design = ma.read_feature_table(root, "fragpipe_tmt",
                                         _tmt_cfg(ma, root))
    assert "group" not in design.columns
    assert "group differently" in "\n".join(design.attrs["design_notes"])


def test_condition_from_name_takes_a_regex_with_one_capture_group(ma,
                                                                  tmp_path):
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "R01T"), ("127N", "N02T")],
        "TMT2": [("126", "R03T"), ("127N", "N04T")]})
    _, _, design = ma.read_feature_table(
        root, "fragpipe_tmt",
        _tmt_cfg(ma, root, condition_from_name=r"^([RN])\d+"))
    assert list(design["group"]) == ["R", "N", "R", "N"]
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt",
                              _tmt_cfg(ma, root,
                                       condition_from_name=r"^(X)(\d+)"))
    assert "exactly one" in str(e.value)


def test_condition_from_name_off_never_guesses(ma, tmp_path, capsys):
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "resp_1"), ("127N", "nonresp_1")],
        "TMT2": [("126", "resp_2"), ("127N", "nonresp_2")]})
    _, _, design = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, condition_from_name=""))
    assert "group" not in design.columns
    assert "no condition is derived" in capsys.readouterr().err


def test_a_condition_confounded_with_the_plex_stops_the_run(ma, tmp_path):
    # symptom: with one condition per plex the batch effect and the biology
    # are the same vector. Every fold change would be both, and the failure
    # otherwise surfaces as an unestimable coefficient deep in limma.
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "resp_1"), ("127N", "resp_2")],
        "TMT2": [("126", "nonresp_1"), ("127N", "nonresp_2")]})
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt", _tmt_cfg(ma, root))
    msg = str(e.value)
    assert "perfectly confounded with the plex" in msg
    assert "each of the 2 plexes contains exactly one condition" in msg


def test_the_join_stage_records_how_the_tmt_design_was_made(ma, tmp_path):
    # design_from_input.tsv is a table of samples; it cannot say where its
    # condition came from or what happened to a reference channel that is no
    # longer in it. The report copies these lines into design_record.txt.
    root = _tmt_named_run(tmp_path, {
        "TMT1": [("126", "resp_1"), ("127N", "nonresp_1"), ("131C", "Pool01")],
        "TMT2": [("126", "resp_2"), ("127N", "nonresp_2"), ("131N", "Pool02")]})
    proj = build_project(tmp_path / "p", quant_table=root,
                         quant_format="fragpipe_tmt", manifest="",
                         tmt={"reference_name": "Pool*"})
    proj.run()
    d = pd.read_csv(proj.rpath("quant", "design_from_input.tsv"), sep="\t")
    assert list(d.columns) == ["sample", "plex", "channel", "group"]
    assert "Pool01" not in set(d["sample"])
    notes = open(proj.rpath("quant", "design_notes.txt"),
                 encoding="utf-8").read()
    assert "condition source: derived from the sample name" in notes
    assert "reference:        covariate" in notes


# --- column detection -------------------------------------------------

def test_boolean_metadata_columns_are_not_detected_as_samples(ma, tmp_path):
    # symptom: pandas calls a bool column numeric, so 'Is Decoy' and
    # 'Is Contaminant' were auto-detected as sample channels.
    path = tmp_path / "combined_protein.tsv"
    pd.DataFrame({"Protein": ["P_ko_path", "P_dark1"],
                  "Is Decoy": [False, False],
                  "Is Contaminant": [False, False],
                  "S1": [10.0, 20.0], "S2": [30.0, 40.0]}).to_csv(
        path, sep="\t", index=False)
    df = ma.read_delim_table(str(path))
    meta = [c for c in ma.FRAGPIPE_META if c in df.columns]
    int_cols = [c for c in df.columns
                if c not in meta and pd.api.types.is_numeric_dtype(df[c])
                and not pd.api.types.is_bool_dtype(df[c])]
    assert int_cols == ["S1", "S2"]


def test_the_delimiter_is_taken_from_the_header_not_sniffed(ma, tmp_path):
    # symptom: sep=None forces the pure-python parser, ~18x slower and 4-6x
    # the memory of the C parser on a 400k-row table, for no benefit.
    path = tmp_path / "x.csv"
    path.write_text("a,b,c\n1,2,3\n", encoding="utf-8")
    assert list(ma.read_delim_table(str(path)).columns) == ["a", "b", "c"]
    path2 = tmp_path / "x.tsv"
    path2.write_text("a\tb\tc\n1\t2\t3\n", encoding="utf-8")
    assert list(ma.read_delim_table(str(path2)).columns) == ["a", "b", "c"]


def test_the_peptide_only_reader_survives_a_table_the_quant_reader_refuses(
        ma, tmp_path, capsys):
    # symptom: a table whose QUANTIFICATION this tool cannot use aborted the
    # taxonomy work too, although its peptide column was perfectly readable.
    path = str(tmp_path / "peptide.tsv")
    F.write_tmt_peptide_table(path, F.protein_set()[:3], ["126", "127N"])
    cfg = _cfg(ma, quant_table=path, quant_format="fragpipe_peptide",
               peptide_only_reader="auto")
    feats = ma.peptide_features(cfg, "unipept")
    assert len(feats) == 3
    assert "falling back to the peptide-only reader" in capsys.readouterr().err


def test_peptide_only_reader_never_hides_the_refusal(ma, tmp_path):
    path = str(tmp_path / "peptide.tsv")
    F.write_tmt_peptide_table(path, F.protein_set()[:3], ["126", "127N"])
    cfg = _cfg(ma, quant_table=path, quant_format="fragpipe_peptide",
               peptide_only_reader="never")
    with pytest.raises(ma.StageError):
        ma.peptide_features(cfg, "unipept")


# --- TMT: within-plex normalisation, purity, and the size factors -----
def _tmt_loaded_run(tmp_path, name, loads, missing=False, seed=13):
    """One plex per entry of `loads`, four channels plus a pool, thirty
    proteins over three taxa, and NO biology: every sample has the same true
    abundance, so anything that survives into a per-sample number is an
    artefact.

    `loads` scales every channel of a plex, which is what a loading difference
    between plexes is. With missing=True, 7 of taxon T2's 10 proteins are seen
    in the LAST plex only, which is the plex-shaped missingness the whole
    design has to be read against; the last plex is the one whose loading is
    furthest from the rest, because that is where mixing references inside one
    plex with references across all of them costs the most.
    """
    import math
    import random
    rng = random.Random(seed)
    root = str(tmp_path / name)
    prots = {f"T{t}_P{i}": f"T{t}" for t in (0, 1, 2) for i in range(10)}
    for plex, load in loads.items():
        chans = [(c, f"{plex}_S{j}") for j, c in
                 enumerate(("126", "127N", "128N", "129N"))]
        chans.append(("131C", f"Pool_{plex}"))
        rows = []
        for prot, tax in prots.items():
            if (missing and tax == "T2" and int(prot.rsplit("P", 1)[1]) >= 3
                    and plex != list(loads)[-1]):  # 7 of T2's 10
                continue
            base = 20000 * (1.0 + 0.4 * rng.random())
            rows.append({"peptide": f"{prot}PEPK", "razor": prot,
                         "values": {s: round(base * load
                                             * math.exp(rng.gauss(0, 0.06)))
                                    for _, s in chans}})
        F.write_tmt_plex(root, plex, chans, rows, seed=seed)
    return root, prots


def _size_factor_parts(ma, root, prots, norm="none"):
    """-> (common per-sample log2 factor, taxon-specific residual, plex map,
    protein table, size factors).

    The taxon size factor is one number per taxon per sample. Split it into
    the part every taxon shares (the sample's loading, which a size factor is
    FOR) and the part that differs between taxa: only the second can carry a
    plex effect into a taxon's ratio model.
    """
    # condition_from_name off: these sample names carry the plex, and the
    # derivation is right to refuse them. This fixture is about the size
    # factors, not the design.
    cfg = _tmt_cfg(ma, root, reference_name="Pool*", condition_from_name="",
                   within_plex_normalise=norm)
    feats, int_cols, design = ma.read_fragpipe_tmt(root, cfg)
    prot, _ev, _fc = ma.rollup_features(feats, int_cols, {}, "razor", 0)
    prot["tax"] = prot["group_id"].map(prots)
    sf = ma.taxon_size_factors(prot, "tax", int_cols, 4)
    lg = np.log2(sf.set_index("tax")[int_cols].astype(float))
    # the shared part is read off the taxa that HAVE a factor everywhere: a
    # taxon that is NA in most samples cannot contribute to a per-sample
    # median without dragging it wherever it happens to exist.
    whole = lg.dropna(axis=0, how="any")
    common = (whole if len(whole) else lg).median(axis=0)
    plex = dict(zip(design["sample"], design["plex"]))
    return common, lg.sub(common, axis=1), plex, prot, sf


def test_a_plex_effect_does_not_leak_into_the_taxon_size_factor(ma, tmp_path):
    # The size factor is a median of ratios computed across ALL samples, so an
    # isobaric run has to be asked whether the batch got into it. Measured
    # here by running the SAME random draw twice, once with equal plex loading
    # and once with a 7.5x spread between plexes: the plex belongs in the
    # per-sample part (that is what a size factor is for), and what must not
    # move is the taxon-by-taxon part, which is what the ratio model divides
    # by.
    flat, prots = _tmt_loaded_run(tmp_path, "flat",
                                  {"TMT1": 1.0, "TMT2": 1.0, "TMT3": 1.0})
    skew, _ = _tmt_loaded_run(tmp_path, "skew",
                              {"TMT1": 1.0, "TMT2": 3.0, "TMT3": 0.4})
    c0, r0, _p0, _t0, _s0 = _size_factor_parts(ma, flat, prots)
    c1, r1, plex, _t1, sf1 = _size_factor_parts(ma, skew, prots)
    # every taxon could use the complete-case reference, so nothing here is
    # the poscounts fallback the next test is about
    assert set(sf1["method"]) == {"median_of_ratios"}

    def by_plex(c):
        return {p: float(c[[s for s in c.index if plex[s] == p]].mean())
                for p in sorted(set(plex.values()))}

    got = by_plex(c1)
    assert abs((got["TMT2"] - got["TMT1"]) - np.log2(3.0)) < 0.15
    assert abs((got["TMT3"] - got["TMT1"]) - np.log2(0.4)) < 0.15
    assert max(abs(v) for v in by_plex(c0).values()) < 0.15
    # and it left NO trace in the taxon-specific part: the same numbers with
    # and without the plex effect, to within the rounding of the integer
    # intensities the fixture writes - 1e-3 log2 against a plex effect of
    # 1.58 log2, three orders of magnitude smaller than what it would have to
    # be to reach a taxon's ratio model.
    assert np.nanmax(np.abs(r0.to_numpy() - r1.to_numpy())) < 1e-3


def test_a_plex_confined_taxon_is_where_the_size_factor_is_exposed(ma,
                                                                   tmp_path):
    # The exposure that IS real, and it is missingness rather than the plex
    # effect: with too few members observed in every plex a taxon loses the
    # complete-case reference, falls to the poscounts variant, and its median
    # then mixes proteins referenced inside one plex with proteins referenced
    # across all of them.
    root, prots = _tmt_loaded_run(tmp_path, "conf",
                                  {"TMT1": 1.0, "TMT2": 3.0, "TMT3": 0.4},
                                  missing=True)
    _c, resid, plex, prot, sf = _size_factor_parts(ma, root, prots)
    assert "median_of_ratios_poscounts" in set(sf["method"])
    exposed, total = ma.tmt_size_factor_plex_exposure(
        prot, "tax", list(plex), plex, 4)
    assert (exposed, total) == (1, 3)      # T2 only, and it is detected
    # T0 and T1 are unaffected; T2's own factor moves by far more than their
    # noise, and only in the plex its members are confined to
    spread = (resid.max(axis=1) - resid.min(axis=1)).abs()
    assert spread["T0"] < 0.2 and spread["T1"] < 0.2
    assert abs(resid.loc["T2"].dropna().mean()) > 0.4


def test_the_plex_exposure_docstring_quoted_other_numbers_than_this_fixture(
        ma, tmp_path):
    # symptom: tmt_size_factor_plex_exposure's docstring said the taxon-
    # specific part of the factor is "plex-free to 0.08 log2" and is "pulled
    # 0.77 log2" for a taxon with "half" its proteins confined to one plex.
    # This fixture confines 7 of 10 and measures 0.03 and 1.40 — which is what
    # README.md and the v0.3.0 CHANGELOG entry quote. Three sources, two
    # answers, one fixture, and no way for a reader to tell which was measured.
    root, prots = _tmt_loaded_run(tmp_path, "docstr",
                                  {"TMT1": 1.0, "TMT2": 3.0, "TMT3": 0.4},
                                  missing=True)
    _c, resid, plex, _prot, _sf = _size_factor_parts(ma, root, prots)
    # per plex, because that is the comparison the prose makes: what the
    # batch did to one taxon's factor inside the plex its members live in.
    per_plex = {(t, p): float(resid.loc[t, [s for s in resid.columns
                                            if plex[s] == p]].mean())
                for t in resid.index for p in sorted(set(plex.values()))}
    clean = max(abs(v) for (t, _p), v in per_plex.items()
                if t != "T2" and not np.isnan(v))
    confined = max(abs(v) for (t, _p), v in per_plex.items()
                   if t == "T2" and not np.isnan(v))
    doc = " ".join((ma.tmt_size_factor_plex_exposure.__doc__ or "").split())
    assert f"{clean:.2f} log2" in doc, \
        f"the docstring does not quote the measured {clean:.2f} log2"
    assert f"{confined:.2f} log2" in doc, \
        f"the docstring does not quote the measured {confined:.2f} log2"
    assert "7 of its 10" in doc, "the docstring misstates how many are confined"


def test_the_exposure_check_is_silent_when_there_is_one_plex(ma, tmp_path):
    root, prots = _tmt_loaded_run(tmp_path, "conf2",
                                  {"TMT1": 1.0, "TMT2": 2.0}, missing=True)
    _c, _r, plex, prot, _sf = _size_factor_parts(ma, root, prots)
    assert ma.tmt_size_factor_plex_exposure(
        prot, "tax", list(plex), plex, 4) == (1, 3)
    # one plex is not a batch, and there is nothing to be exposed to
    assert ma.tmt_size_factor_plex_exposure(
        prot, "tax", list(plex), {s: "TMT1" for s in plex}, 4) == (0, 0)


def test_within_plex_normalisation_brings_the_channels_to_one_scale(ma,
                                                                    tmp_path,
                                                                    capsys):
    # symptom: channels of one plex differ by how much peptide was loaded and
    # how completely it was labelled - a per-channel constant with no biology
    # in it, which the roll-up would sum straight into the protein.
    root = str(tmp_path / "run")
    rows = [{"peptide": f"PEP{i}K", "razor": "P_ko_path",
             "values": {"A1": 100 * (i + 1), "A2": 400 * (i + 1),
                        "A3": 200 * (i + 1)}} for i in range(5)]
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("128N", "A3")], rows)
    feats, int_cols, design = ma.read_feature_table(
        root, "fragpipe_tmt",
        _tmt_cfg(ma, root, within_plex_normalise="median"))
    med = feats[int_cols].median()
    assert abs(med["A1"] - med["A2"]) < 1e-6
    assert abs(med["A2"] - med["A3"]) < 1e-6
    # centred on the plex's own median channel, so the values stay on the
    # linear intensity scale the size factors need rather than becoming ~1
    assert 100 < float(med["A1"]) < 100000
    err = capsys.readouterr().err
    assert "within-plex median centring" in err
    notes = "\n".join(design.attrs["design_notes"])
    assert "within-plex norm: median centring per channel" in notes


def test_within_plex_normalisation_leaves_the_between_plex_difference(
        ma, tmp_path):
    # Deliberate: the between-plex difference is the batch, and the plex term
    # in the model (or the report's own median normalisation) is what removes
    # it. Removing it here would hide it from both.
    root, _prots = _tmt_loaded_run(tmp_path, "keep",
                                   {"TMT1": 1.0, "TMT2": 4.0})
    cfg = _tmt_cfg(ma, root, within_plex_normalise="median")
    feats, int_cols, design = ma.read_fragpipe_tmt(root, cfg)
    plex = dict(zip(design["sample"], design["plex"]))
    med = feats[int_cols].median()
    m1 = np.median([med[s] for s in int_cols if plex[s] == "TMT1"])
    m2 = np.median([med[s] for s in int_cols if plex[s] == "TMT2"])
    assert 3.0 < m2 / m1 < 5.0


def test_within_plex_normalisation_can_be_switched_off_and_says_so(
        ma, tmp_path, capsys):
    root = _tmt_run(tmp_path)
    feats, _int_cols, design = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, within_plex_normalise="none"))
    row = feats.loc[feats["peptide"].eq("SHAREDPEPK")].iloc[0]
    assert row["A1"] == 100 and row["A2"] == 200      # FragPipe's own numbers
    err = capsys.readouterr().err
    assert "within_plex_normalise is 'none'" in err
    assert "none (tmt.within_plex_normalise: none)" in \
        "\n".join(design.attrs["design_notes"])


def test_an_unknown_within_plex_normalisation_is_refused(ma, tmp_path):
    root = _tmt_run(tmp_path)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(
            root, "fragpipe_tmt",
            _tmt_cfg(ma, root, within_plex_normalise="quantile"))
    assert "within_plex_normalise" in str(e.value)
    assert "quantile" in str(e.value)


# --- TMT: min_purity, which only psm.tsv can answer -------------------
def _tmt_purity_run(tmp_path, purity, **plexkw):
    root = str(tmp_path / "run")
    rows = [{"peptide": "CLEANPEPK", "razor": "P_ko_path"},
            {"peptide": "DIRTYPEPK", "razor": "P_dark1"},
            {"peptide": "NOPSMPEPK", "razor": "P_dark2"}]
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("131C", "Pool01")],
                     rows, psm=purity, **plexkw)
    return root


def test_min_purity_filters_on_the_median_of_a_features_psms(ma, tmp_path,
                                                             capsys):
    # symptom: purity exists ONLY in psm.tsv (SCHEMA.md), so a min_purity key
    # honoured at ion level without that join would be a filter the user
    # believes in and the numbers never saw.
    root = _tmt_purity_run(tmp_path, {"CLEANPEPK": [0.99, 0.95, 0.20],
                                      "DIRTYPEPK": [0.30, 0.40, 0.99]})
    feats, _ic, _d = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, min_purity=0.5))
    kept = set(feats["peptide"])
    # median 0.95 survives, median 0.40 does not: one bad spectrum does not
    # condemn a feature and one good one does not save it
    assert "CLEANPEPK" in kept and "DIRTYPEPK" not in kept
    assert "min_purity=0.5 drops 1 of 3" in capsys.readouterr().err


def test_a_feature_with_no_psm_row_is_kept_and_counted(ma, tmp_path, capsys):
    # symptom: an unmatched key is a JOIN failure (FragPipe leaves 'Modified
    # Peptide' empty on rows ion.tsv writes a modified sequence for), and
    # dropping those would look exactly like a purity filter working.
    root = _tmt_purity_run(tmp_path, {"CLEANPEPK": 0.99, "DIRTYPEPK": 0.10})
    feats, _ic, design = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, min_purity=0.5))
    assert "NOPSMPEPK" in set(feats["peptide"])
    err = capsys.readouterr().err
    assert "match no row of psm.tsv" in err and "they are KEPT" in err
    assert "1 unjudged and kept" in "\n".join(design.attrs["design_notes"])


def test_min_purity_without_a_psm_table_is_refused_naming_it(ma, tmp_path):
    root = _tmt_run(tmp_path)                 # written without psm.tsv
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt",
                              _tmt_cfg(ma, root, min_purity=0.5))
    assert "psm.tsv" in str(e.value) and "TMT1" in str(e.value)


def test_min_purity_outside_zero_to_one_is_refused(ma, tmp_path):
    root = _tmt_run(tmp_path)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt",
                              _tmt_cfg(ma, root, min_purity=50))
    assert "between 0 and 1" in str(e.value)


def test_psm_tables_are_stage_inputs_only_when_min_purity_reads_them(
        ma, tmp_path):
    # symptom: listing psm.tsv unconditionally would invalidate every cached
    # join the moment FragPipe rewrote a file the run never opened.
    root = _tmt_purity_run(tmp_path, {"CLEANPEPK": 0.9})
    off = ma.quant_inputs(_tmt_cfg(ma, root))
    on = ma.quant_inputs(_tmt_cfg(ma, root, min_purity=0.5))
    assert not any(f.endswith("psm.tsv") for f in off)
    assert any(f.endswith("psm.tsv") for f in on)


def test_a_non_numeric_min_purity_dies_in_the_reader_not_in_the_signature(
        ma, tmp_path):
    # symptom: the stage signature float()ed tmt.min_purity before the reader
    # could validate anything, so a value like '90%' surfaced as a bare
    # ValueError traceback instead of the message that names the key, the
    # range and the value.
    root = _tmt_run(tmp_path)
    cfg = _tmt_cfg(ma, root, min_purity="90%")
    assert ma.quant_inputs(cfg) == [root]        # tolerant, not a traceback
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt", cfg)
    assert "tmt.min_purity must be a number between 0 and 1" in str(e.value)
    assert "'90%'" in str(e.value)


def test_a_weak_and_empty_channel_is_named_as_under_corrected(ma, tmp_path,
                                                              capsys):
    # symptom: the median is taken over OBSERVED values, so a channel whose
    # low end went missing has a median sitting above its true centre and is
    # scaled up too little. An unequal but COMPLETE load is what median
    # centring is for and must not raise the same alarm.
    root = str(tmp_path / "run")
    rows = []
    for i in range(20):
        v = {"A1": 10000 + i, "A2": 10000 + i,
             # A3 is a tenth of the others and half of its values are gone
             "A3": (1000 + i) if i % 2 else 0}
        rows.append({"peptide": f"PEP{i}K", "razor": "P_ko_path",
                     "values": v})
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("128N", "A3")], rows)
    ma.read_feature_table(root, "fragpipe_tmt",
                          _tmt_cfg(ma, root, within_plex_normalise="median"))
    err = capsys.readouterr().err
    assert "UNDER-corrected" in err and "A3 (128N)" in err
    assert "50% missing" in err


def test_an_unequal_but_complete_load_is_normalised_without_an_alarm(ma,
                                                                     tmp_path,
                                                                     capsys):
    root = str(tmp_path / "run")
    rows = [{"peptide": f"PEP{i}K", "razor": "P_ko_path",
             "values": {"A1": 10000 + i, "A2": 1000 + i, "A3": 5000 + i}}
            for i in range(20)]
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("128N", "A3")], rows)
    ma.read_feature_table(root, "fragpipe_tmt",
                          _tmt_cfg(ma, root, within_plex_normalise="median"))
    err = capsys.readouterr().err
    assert "within-plex median centring" in err
    assert "UNDER-corrected" not in err


# --- TMT: the effect that was planted, and the batch that was not ------
def _tmt_planted_cfg(ma, root, **tmt):
    """The planted run's config: a real reference, and every other default
    left alone — this is the path a user gets, median centring included."""
    cfg = _cfg(ma, quant_table=root, quant_format="fragpipe_tmt")
    cfg["tmt"].update(dict(reference_name="Pool*"))
    cfg["tmt"].update(tmt)
    return cfg


def _fit_group_plex(prot, int_cols, group_of, plex_of):
    """-> {protein: (with the plex in the model, without it)}, in log2.

    The first number is what '~ 0 + group + plex' estimates for
    group_b - group_a: an ordinary least-squares fit of an additive
    condition-plus-batch model, which is the model the report writes for an
    isobaric run. The second is the same contrast with no plex term, i.e. the
    difference of the two group means. In a CROSSED design the two agree and
    neither says anything about the covariate; the fixture is deliberately
    unbalanced so that they do not.
    """
    g = np.array([group_of[s] for s in int_cols])
    p = np.array([plex_of[s] for s in int_cols])
    levels, plexes = sorted(set(g)), sorted(set(p))
    X = np.column_stack([(g == x).astype(float) for x in levels]
                        + [(p == x).astype(float) for x in plexes[1:]])
    out = {}
    for _i, row in prot.iterrows():
        y = np.log2(row[int_cols].astype(float).to_numpy(dtype=float))
        ok = np.isfinite(y)
        beta, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
        with_plex = beta[levels.index("b")] - beta[levels.index("a")]
        blind = float(np.nanmean(y[(g == "b") & ok])
                      - np.nanmean(y[(g == "a") & ok]))
        out[row["group_id"]] = (float(with_plex), blind)
    return out


def _planted_estimates(ma, root, truth, **tmt):
    """Read the planted run, roll it up to proteins and fit both models."""
    cfg = _tmt_planted_cfg(ma, root, **tmt)
    feats, int_cols, design = ma.read_feature_table(root, "fragpipe_tmt", cfg)
    prot, _ev, _fc = ma.rollup_features(feats, int_cols, {}, "razor", 0)
    group_of = dict(zip(design["sample"], design["group"]))
    plex_of = dict(zip(design["sample"], design["plex"]))
    return feats, int_cols, design, prot, _fit_group_plex(
        prot, int_cols, group_of, plex_of)


def test_the_planted_condition_effect_is_recovered_and_the_plex_is_absorbed(
        ma, tmp_path):
    # symptom: this is the whole reason the plex goes into the model as a
    # covariate. Two plexes with a 3x loading difference and an UNBALANCED
    # condition split — 3 'a' + 1 'b' in one plex, the reverse in the other —
    # so a planted 2x up and a planted 2x down have to come back at their
    # planted size, and the batch must not appear as a fold change of its own.
    root, truth = F.tmt_planted_run(str(tmp_path / "run"))
    feats, int_cols, design, _prot, est = _planted_estimates(ma, root, truth)
    assert int_cols == truth.samples                 # the pools are not here
    assert list(design["group"]) == [truth.group_of[s] for s in int_cols]

    # the plex effect is really in the matrix: 3x, as planted
    med = feats[int_cols].median()
    per_plex = {p: np.median([med[s] for s in int_cols
                              if truth.plex_of[s] == p])
                for p in sorted(set(truth.plex_of.values()))}
    assert np.log2(per_plex["TMT2"] / per_plex["TMT1"]) == pytest.approx(
        truth.plex_log2["TMT2"], abs=0.15)

    up = [est[p][0] for p in truth.regulated if truth.effect_of[p] > 0]
    down = [est[p][0] for p in truth.regulated if truth.effect_of[p] < 0]
    null = [est[p][0] for p in truth.null]
    assert up and down and len(null) == len(truth.null)
    # recovered, not merely non-zero: +1 and -1 log2, to within 0.12
    assert np.mean(up) == pytest.approx(+truth.effect_log2, abs=0.12)
    assert np.mean(down) == pytest.approx(-truth.effect_log2, abs=0.12)
    assert max(abs(v) for v in null) < 0.12
    assert min(up) > 0 > max(down)

    # and the covariate is what did it. Without the plex term the SAME data
    # gives every protein the same spurious lift, because the plex and the
    # condition are correlated: (n_b/n - n_a/n) x the plex effect, which is
    # half of it here.
    bias = 0.5 * truth.plex_log2["TMT2"]
    assert np.mean([est[p][1] for p in truth.null]) == pytest.approx(
        bias, abs=0.12)
    # a protein that truly moved DOWN 2x is then reported as barely moving
    assert np.mean([est[p][1] for p in truth.regulated
                    if truth.effect_of[p] < 0]) > -0.4


def test_the_ratio_and_the_covariate_route_recover_the_same_planted_effect(
        ma, tmp_path):
    # the two reference treatments are different arithmetic — one divides
    # every channel by its plex pool, the other keeps the pool out of the
    # matrix and puts the plex in the model — so they have to be shown to
    # agree on the answer rather than assumed to.
    root, truth = F.tmt_planted_run(str(tmp_path / "run"))
    _f, cov_cols, _d, _p, cov = _planted_estimates(ma, root, truth)
    _f2, rat_cols, _d2, _p2, rat = _planted_estimates(
        ma, root, truth, use_reference_ratios=True)
    assert cov_cols == rat_cols                      # the same samples
    for pid in sorted(truth.regulated):
        c, r = cov[pid][0], rat[pid][0]
        assert np.sign(c) == np.sign(truth.effect_of[pid])
        assert np.sign(r) == np.sign(c), f"{pid}: {c:+.3f} vs {r:+.3f}"
        assert abs(r - c) < 0.10, f"{pid}: covariate {c:+.3f}, ratio {r:+.3f}"
        assert abs(r - truth.effect_of[pid]) < 0.12
    # what the ratios buy, and why both are offered: dividing by a master
    # pool takes the plex out of the DATA, so even the plex-blind estimate is
    # right — under the covariate route only the model can fix it.
    assert np.mean([rat[p][1] for p in truth.null]) == pytest.approx(
        0, abs=0.12)
    assert np.mean([cov[p][1] for p in truth.null]) > 0.5


def test_a_protein_seen_in_one_plex_only_is_dropped_by_min_plexes(ma,
                                                                  tmp_path):
    # symptom: min_valid_per_group counts SAMPLES and cannot see which batch
    # they were in, so a protein quantified in every channel of one plex
    # satisfies it completely — while every one of its numbers carries that
    # plex's loading and none of them is comparable with anything else.
    root, truth = F.tmt_planted_run(
        str(tmp_path / "run"), n_confined=3,
        layout=(("a", "a", "b", "b"), ("a", "a", "b", "b")))
    feats, int_cols, design = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_planted_cfg(ma, root))
    prot = ma.rollup_features(feats, int_cols, {}, "razor", 0)[0] \
        .set_index("group_id")
    group_of = dict(zip(design["sample"], design["group"]))

    # the count min_valid_per_group would see, computed here rather than
    # assumed: every confined protein passes it at 2 in BOTH groups, on one
    # plex alone.
    for pid in sorted(truth.confined):
        row = prot.loc[pid]
        counts = {g: int(sum(pd.notna(row[s]) for s in int_cols
                             if group_of[s] == g))
                  for g in sorted(set(group_of.values()))}
        assert counts == {"a": 2, "b": 2}
        assert {truth.plex_of[s] for s in int_cols
                if pd.notna(row[s])} == {"TMT1"}

    kept = set(ma.read_feature_table(
        root, "fragpipe_tmt",
        _tmt_planted_cfg(ma, root, min_plexes=2))[0]["razor_protein"])
    assert not (kept & truth.confined)
    assert kept == set(truth.null | truth.regulated)


def test_a_condition_confounded_with_the_plex_stops_the_whole_run(tmp_path):
    # the guard has its own unit test; this is what the user sees. One
    # condition per plex is not a design this pipeline can analyse, and it
    # has to fail at the join rather than deep inside limma.
    proteins = F.protein_set()
    root, _truth = F.tmt_planted_run(
        str(tmp_path / "run"), proteins=proteins, n_proteins=6,
        layout=(("a", "a", "a", "a"), ("b", "b", "b", "b")))
    proj = build_project(tmp_path / "p", proteins=proteins, quant_table=root,
                         quant_format="fragpipe_tmt", manifest="",
                         tmt={"reference_name": "Pool*"})
    proc = proj.run(expect=1)
    assert "perfectly confounded with the plex" in proc.stderr
    assert "each of the 2 plexes contains exactly one condition" in proc.stderr
    assert "TMT1" in proc.stderr and "TMT2" in proc.stderr
    assert "Traceback" not in proc.stderr


# --- TMT: a reporter zero is not a measurement ------------------------
def _tmt_zero_run(tmp_path):
    """One plex, three channels. P_one has a single peptide and loses one
    channel to a zero; P_many keeps two of its three peptides there."""
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1",
                     [("126", "A1"), ("127N", "A2"), ("128N", "A3")],
                     [{"peptide": "ONLYPEPK", "razor": "P_one",
                       "values": {"A1": 100, "A2": 0, "A3": 300}},
                      {"peptide": "MANYPEP1K", "razor": "P_many",
                       "values": {"A1": 100, "A2": 0, "A3": 300}},
                      {"peptide": "MANYPEP2K", "razor": "P_many",
                       "values": {"A1": 100, "A2": 40, "A3": 300}},
                      {"peptide": "MANYPEP3K", "razor": "P_many",
                       "values": {"A1": 100, "A2": 60, "A3": 300}}])
    return root


def test_a_reporter_zero_is_read_as_missing_and_never_summed_as_a_number(
        ma, tmp_path, capsys):
    # symptom: FragPipe writes 0 into 5.9%-16% of the reporter cells of a
    # real plex for "not quantified". Kept as a number it becomes a real
    # measurement of zero: a protein seen in one channel only is reported as
    # ABSENT rather than as unmeasured, which is a fold change.
    root = _tmt_zero_run(tmp_path)
    feats, int_cols, _d = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root))
    f = feats.set_index("peptide")
    assert pd.isna(f.loc["ONLYPEPK", "A2"])
    assert f.loc["ONLYPEPK", "A1"] == 100 and f.loc["ONLYPEPK", "A3"] == 300
    prot = ma.rollup_features(feats, int_cols, {}, "razor", 0)[0] \
        .set_index("group_id")
    assert pd.isna(prot.loc["P_one", "A2"])
    # and the protein with three peptides is the sum of the two measured ones
    assert prot.loc["P_many", "A2"] == 100.0
    assert prot.loc["P_many", "A1"] == 300.0
    assert "2 of 12 reporter cell(s) (16.7%) are 0" in capsys.readouterr().err


def test_keeping_the_reporter_zeros_changes_what_it_is_meant_to_change(
        ma, tmp_path):
    # the opt-out has to be worth having: with the zeros kept, the same
    # protein reads as measured-at-zero, which log2 cannot represent.
    root = _tmt_zero_run(tmp_path)
    cfg = _tmt_cfg(ma, root)
    cfg["zero_intensity_is_missing"] = False
    feats, int_cols, _d = ma.read_feature_table(root, "fragpipe_tmt", cfg)
    prot = ma.rollup_features(feats, int_cols, {}, "razor", 0)[0] \
        .set_index("group_id")
    assert prot.loc["P_one", "A2"] == 0.0
    assert prot.loc["P_many", "A2"] == 100.0       # 0 adds nothing to a sum


def test_a_reporter_zero_never_enters_the_within_plex_median(ma, tmp_path):
    # the centring divides each channel by its OWN median, so a zero counted
    # as a measurement drags that median down and the channel is then scaled
    # up by far more than its loading deserves — the zeros would be paid for
    # a second time, by every real value in the channel.
    root = _tmt_zero_run(tmp_path)

    def a2_of(zeros_are_missing):
        cfg = _tmt_cfg(ma, root, within_plex_normalise="median")
        cfg["zero_intensity_is_missing"] = zeros_are_missing
        feats, int_cols, _d = ma.read_feature_table(root, "fragpipe_tmt", cfg)
        med = feats[int_cols].median()
        # whatever the scale, centring puts every channel on the same one
        assert float(med["A1"]) == pytest.approx(float(med["A3"]))
        return float(feats.set_index("peptide").loc["MANYPEP2K", "A2"])

    # A2 holds 0, 0, 40, 60. Observed, its median is 50 and the plex's median
    # channel is A1 at 100, so A2 is scaled 2x and its 40 becomes 80. With the
    # zeros counted the median is 20, the scale factor is 5, and the same
    # measurement is reported as 200.
    assert a2_of(True) == pytest.approx(80.0)
    assert a2_of(False) == pytest.approx(200.0)


# --- TMT: the files that are still refused ----------------------------
TMT_REPORT_FILES = [(kind, level)
                    for kind in ("abundance", "ratio")
                    for level in ("gene", "protein", "peptide",
                                  "modified-peptide")]


@pytest.mark.parametrize("kind,level", TMT_REPORT_FILES)
def test_every_tmt_report_matrix_is_refused_whatever_its_level(ma, tmp_path,
                                                               kind, level):
    # symptom: a real run writes eight of these, and the peptide-level ones
    # carry 'Peptide' and 'Mapped Proteins', so they look more like readable
    # feature input than the protein ones do. Every one of them is already
    # log2 and median-centred and already carries TMT-Integrator's own
    # inference, so every one has to be refused by name — not just the one
    # that happened to be tested.
    path = str(tmp_path / f"{kind}_{level}_MD.tsv")
    F.write_tmt_report_matrix(path, ["A1", "A2"], level=level)
    for fmt in ("fragpipe_tmt", "fragpipe_peptide"):
        with pytest.raises(ma.StageError) as e:
            ma.read_feature_table(path, fmt, _cfg(ma))
        msg = str(e.value)
        assert "tmt-report matrix" in msg
        assert "ReferenceIntensity" in msg           # what it found
        assert "quant_format: fragpipe_tmt" in msg   # and what to do instead
    # and the taxonomy stages' reader, which never touches an intensity, must
    # not read it either: its rows are TMT-Integrator's protein groups.
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_peptides(path, "fragpipe_peptide", _cfg(ma))
    assert "tmt-report matrix" in str(e.value)


def test_a_tmt_report_matrix_is_refused_at_protein_level_too(tmp_path):
    # the protein-level formats do not go through read_feature_table at all,
    # so the refusal has to sit in the join stage as well — and this is the
    # format a user reaches for when the file is called 'abundance_protein'.
    proj = build_project(tmp_path / "p")
    path = proj.path("input", "abundance_protein_MD.tsv")
    F.write_tmt_report_matrix(path, proj.samples, level="protein")
    proj.write_config(quant_table=path, quant_format="fragpipe", manifest="")
    proc = proj.run(expect=1)
    assert "tmt-report matrix" in proc.stderr
    assert "ReferenceIntensity" in proc.stderr
    assert "Traceback" not in proc.stderr


def test_a_per_plex_protein_table_is_refused_not_read_as_one_plex(ma,
                                                                  tmp_path):
    # symptom: TMTn/protein.tsv has neither of the markers the other isobaric
    # files are caught by — no ReferenceIntensity, no 'Channel <mass>' — so
    # the protein-level column detector took every numeric column and
    # quantified ONE plex as the whole experiment, with 'Length',
    # 'Protein Qvalue' and 'Razor Intensity' sitting in the matrix as if they
    # were samples.
    path = str(tmp_path / "protein.tsv")
    F.write_tmt_protein_table(path, F.protein_set()[:3],
                              [("126", "A1"), ("127N", "A2"),
                               ("131C", "Pool01")])
    with pytest.raises(ma.StageError) as e:
        ma.refuse_per_plex_reporter_table(path, ma.header_columns(path))
    msg = str(e.value)
    assert "per-plex FragPipe TMT table" in msg
    assert "Intensity A1" in msg                     # what it found
    assert "SINGLE plex" in msg and "Razor Intensity" in msg
    assert "quant_format: fragpipe_tmt" in msg


def test_the_per_plex_protein_refusal_reaches_the_cli(tmp_path):
    proj = build_project(tmp_path / "p")
    path = proj.path("input", "protein.tsv")
    F.write_tmt_protein_table(path, proj.proteins,
                              [("126", "A1"), ("131C", "Pool01")])
    proj.write_config(quant_table=path, quant_format="fragpipe", manifest="")
    proc = proj.run(expect=1)
    assert "per-plex FragPipe TMT table" in proc.stderr
    assert "Traceback" not in proc.stderr


def test_a_label_free_protein_table_is_not_caught_by_that_refusal(ma,
                                                                  tmp_path):
    # the discriminator is the PREFIX form: FragPipe writes 'Intensity
    # <sample>' for a reporter channel and '<sample> Intensity' for a
    # label-free run, so a combined_protein.tsv can never match.
    path = str(tmp_path / "combined_protein.tsv")
    F.write_fragpipe_protein(path, F.protein_set()[:3], ["A_1", "B_1"])
    assert ma.refuse_per_plex_reporter_table(path,
                                             ma.header_columns(path)) is None


def test_a_per_plex_psm_table_is_refused_naming_its_channels(ma, tmp_path):
    # psm.tsv carries the reporter columns too, and it is the table the
    # purity filter reads — so it is exactly the file a user might point
    # quant_table at. It has one row per spectrum, not per feature.
    root = str(tmp_path / "run")
    d = F.write_tmt_plex(root, "TMT1",
                         [("126", "A1"), ("127N", "A2")],
                         [{"peptide": "PEPTIDEK", "razor": "P_ko_path"}],
                         psm={"PEPTIDEK": 0.9})
    path = os.path.join(d, "psm.tsv")
    for fmt in ("fragpipe_peptide", "fragpipe_ion"):
        with pytest.raises(ma.StageError) as e:
            ma.read_feature_table(path, fmt, _cfg(ma))
        msg = str(e.value)
        assert "isobaric (TMT/iTRAQ) output" in msg
        assert "Intensity A1" in msg                 # what it found
        assert "quant_format: fragpipe_tmt" in msg


@pytest.mark.parametrize("what", ["ion.tsv", "peptide.tsv", "protein.tsv",
                                  "psm.tsv"])
def test_fragpipe_tmt_pointed_at_one_file_names_the_run_directory(ma,
                                                                  tmp_path,
                                                                  what):
    # quant_table is a DIRECTORY for this format, and the natural mistake is
    # to point it at the table inside one plex — which would be one plex read
    # as the experiment if it were allowed through.
    root = _tmt_run(tmp_path)
    path = os.path.join(root, "TMT1", "ion.tsv")
    target = os.path.join(root, "TMT1", what)
    if not os.path.exists(target):
        shutil.copyfile(path, target)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(target, "fragpipe_tmt", _tmt_cfg(ma, root))
    msg = str(e.value)
    assert "must be the run directory" in msg
    assert what in msg                               # names what it was given


def test_fragpipe_tmt_pointed_at_one_plex_directory_says_so(ma, tmp_path):
    # the other half of the same mistake: the plex directory itself holds no
    # plex directories, and falling through to "no features" would be worse
    # than saying which glob matched nothing.
    root = _tmt_run(tmp_path)
    plex = os.path.join(root, "TMT1")
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(plex, "fragpipe_tmt", _tmt_cfg(ma, plex))
    msg = str(e.value)
    assert "no plex directory matches tmt.plex_glob 'TMT*'" in msg
    assert "ion.tsv" in msg                          # what is there instead


# --- TMT: the reference channel is in no output that names samples -----
def _tmt_planted_project(tmp_path, name="p", **over):
    """A whole project over a planted two-plex run, ready to run()."""
    proteins = F.protein_set()
    root, truth = F.tmt_planted_run(str(tmp_path / "run"), proteins=proteins,
                                    n_proteins=8,
                                    layout=(("a", "a", "b", "b"),
                                            ("a", "a", "b", "b")))
    cfg = dict(quant_table=root, quant_format="fragpipe_tmt", manifest="",
               tmt={"reference_name": "Pool*"})
    cfg.update(over)
    proj = build_project(tmp_path / name, proteins=proteins, **cfg)
    proj.truth = truth
    return proj


def test_the_reference_channel_reaches_no_table_that_names_samples(tmp_path):
    # symptom: a pooled bridge left among the samples acquires a condition in
    # the design, joins a group's mean and drags the size factors towards a
    # channel that is in every plex by construction. It has to be gone from
    # every file the report and the R object read samples from — those are
    # what become colData.
    proj = _tmt_planted_project(tmp_path)
    proj.run()
    pools = set(proj.truth.reference_of.values())
    assert pools == {"Pool01", "Pool02"}

    design = pd.read_csv(proj.rpath("quant", "design_from_input.tsv"), sep="\t")
    assert sorted(design["sample"]) == sorted(proj.truth.samples)
    assert not (pools & set(design["sample"]))
    assert "is_reference" not in design.columns      # nothing left to flag

    recorded = [l.strip() for l in open(
        proj.rpath("quant", "sample_columns.txt"), encoding="utf-8")
        if l.strip()]
    assert recorded == proj.truth.samples
    quant = pd.read_csv(proj.rpath("quant", "annotated_quant.tsv"), sep="\t",
                        nrows=1)
    assert not (pools & set(quant.columns))
    feat = pd.read_csv(proj.rpath("quant", "feature_quant.tsv"), sep="\t",
                       nrows=1)
    assert not (pools & set(feat.columns))
    # and it is recorded that a channel WAS dropped, and which one
    notes = open(proj.rpath("quant", "design_notes.txt"),
                 encoding="utf-8").read()
    assert "TMT1=131C/Pool01" in notes and "TMT2=131N/Pool02" in notes


# --- TMT: the peptide layer behaves exactly as it does label-free ------
_SHARED_ROWS = [
    {"peptide": "UNIQUEK", "razor": "P1"},
    {"peptide": "SAMETAXK", "razor": "P1", "mapped": ["P2"]},
    {"peptide": "CROSSTAXK", "razor": "P1", "mapped": ["P3"]},
    {"peptide": "NOTAXK", "razor": "P1", "mapped": ["P4"]},
]
_SHARED_TAXA = {"P1": "820", "P2": "820", "P3": "999", "P4": ""}


def _tmt_shared_features(ma, tmp_path):
    """The four features of _shared_features, in two TMT plexes.

    tmt.level 'peptide', so the feature id is the peptide sequence and is
    literally the same string the label-free reader produces — the two can
    then be compared row for row.
    """
    root = str(tmp_path / "run")
    for plex, chans in (("TMT1", [("126", "S1"), ("127N", "S2")]),
                        ("TMT2", [("126", "S3"), ("127N", "S4")])):
        F.write_tmt_plex(root, plex, chans, _SHARED_ROWS, level="peptide")
    feats, int_cols, _d = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, level="peptide"))
    return feats, int_cols, _SHARED_TAXA


def test_the_shared_peptide_rule_classifies_tmt_features_as_it_does_others(
        ma, tmp_path):
    # the point of reading the per-plex tables rather than tmt-report/: the
    # shared-peptide rule is the layer this tool exists for, and it has to
    # behave identically on isobaric input, not merely run.
    lf_feats, lf_cols, taxon_of = _shared_features(ma, tmp_path, "lf.tsv")
    tmt_feats, tmt_cols, _t = _tmt_shared_features(ma, tmp_path)
    assert list(tmt_feats["feature_id"]) == list(lf_feats["feature_id"])
    assert [sorted(c) for c in tmt_feats["candidates"]] == \
           [sorted(c) for c in lf_feats["candidates"]]

    kept = {}
    for mode in ("protein_unique", "taxon_unique", "razor"):
        lf = ma.rollup_features(lf_feats, lf_cols, taxon_of, mode, 1)[2]
        tm = ma.rollup_features(tmt_feats, tmt_cols, taxon_of, mode, 1)[2]
        assert dict(zip(tm["peptide"], tm["_class"])) == \
               dict(zip(lf["peptide"], lf["_class"]))
        kept[mode] = int(tm["_assigned"].notna().sum())
    # the same numbers the label-free tests pin, on isobaric input
    assert (kept["protein_unique"], kept["taxon_unique"], kept["razor"]) == \
           (1, 2, 4)


def test_peptide_evidence_counts_the_same_things_on_tmt_input(ma, tmp_path):
    lf_feats, lf_cols, taxon_of = _shared_features(ma, tmp_path, "lf.tsv")
    tmt_feats, tmt_cols, _t = _tmt_shared_features(ma, tmp_path)
    lf_ev = ma.rollup_features(lf_feats, lf_cols, taxon_of, "taxon_unique",
                               1)[1].set_index("protein_id")
    ev = ma.rollup_features(tmt_feats, tmt_cols, taxon_of, "taxon_unique",
                            1)[1].set_index("protein_id")
    assert list(ev.columns) == list(lf_ev.columns)
    row = ev.loc["P1"]
    assert row["n_unique"] == 1
    assert row["n_taxon_unique"] == 1
    assert row["n_features_dropped"] == 2
    assert row["n_features_used"] == 2
    assert bool(row["taxon_unique_dominated"]) is False
    assert row.to_dict() == lf_ev.loc["P1"].to_dict()


def test_a_feature_shared_in_one_plex_only_is_shared_in_all_of_them(ma,
                                                                    tmp_path):
    # TMT-specific: FragPipe maps a peptide per plex, so a peptide can be
    # razor-only in one plex and mapped to a second protein in another.
    # Taking the first plex's word for it would make the rule depend on which
    # plex happened to be read first, and quantify a cross-taxon peptide as
    # if it were unique.
    root = str(tmp_path / "run")
    F.write_tmt_plex(root, "TMT1", [("126", "S1"), ("127N", "S2")],
                     [{"peptide": "CROSSTAXK", "razor": "P1"}],
                     level="peptide")
    F.write_tmt_plex(root, "TMT2", [("126", "S3"), ("127N", "S4")],
                     [{"peptide": "CROSSTAXK", "razor": "P1",
                       "mapped": ["P3"]}], level="peptide")
    feats, int_cols, _d = ma.read_feature_table(
        root, "fragpipe_tmt", _tmt_cfg(ma, root, level="peptide"))
    assert sorted(feats["candidates"].iloc[0]) == ["P1", "P3"]
    fc = ma.rollup_features(feats, int_cols, _SHARED_TAXA, "taxon_unique",
                            1)[2]
    assert list(fc["_class"]) == ["shared"]
    assert fc["_assigned"].isna().all()


def test_a_tmt_run_writes_the_same_peptide_evidence_file_as_a_label_free_one(
        tmp_path):
    # end to end, because the file on disk is what the report and the R
    # object read: the columns, not just the in-memory frame, have to match.
    lf = build_project(tmp_path / "lf")
    lf.run()
    tmt = _tmt_planted_project(tmp_path, "tmt")
    tmt.run()
    a = pd.read_csv(lf.rpath("quant", "peptide_evidence.tsv"), sep="\t")
    b = pd.read_csv(tmt.rpath("quant", "peptide_evidence.tsv"), sep="\t")
    assert list(a.columns) == list(b.columns)
    assert set(b["rollup_method"]) == {"sum"}
    assert set(b["peptide_assignment"]) == {"taxon_unique"}
    # two peptides per protein in the planted fixture, both unique, and every
    # protein of the fixture is there
    assert len(b) == 8
    assert set(b["n_features_used"]) == {2}


# --- the log level has to follow the content -------------------------
# symptom: on a 60-hour TMT run the design recovery logged
#   WARN tmt: the condition could not be derived from the sample names ...
#        analysis.metadata (...) does name every sample, so the report
#        supplies the condition; this affects only design_from_input.tsv
# A warning that ends by saying nothing is affected is a false alarm by its own
# admission, and a WARN on a run that long is a thing you stop and
# investigate. _metadata_note already got the PROSE right -- its own docstring
# calls this case a false alarm -- but every caller logged WARN regardless.
def _no_separator_design(ma):
    """Sample names carrying no condition, as the real cohort's do not:
    MF0030, MF0071, MF001A3A817."""
    return pd.DataFrame({"sample": [f"MF{i:04d}" for i in range(6)],
                         "plex": ["TMT1"] * 3 + ["TMT2"] * 3})


def _levels_for(ma, capsys, cfg):
    capsys.readouterr()
    design, note = ma._tmt_add_condition(_no_separator_design(ma), cfg)
    lines = [l for l in capsys.readouterr().err.splitlines()
             if "condition could not be derived" in l
             or "no condition is" in l]
    assert len(lines) == 1, lines
    return lines[0], note


def test_a_condition_the_metadata_supplies_is_information_not_a_warning(
        ma, tmp_path, capsys):
    md = tmp_path / "metadata.tsv"
    md.write_text("sample\tgroup\n"
                  + "".join(f"MF{i:04d}\tresponder\n" for i in range(6)),
                  encoding="utf-8")
    line, note = _levels_for(ma, capsys, {
        "tmt": {"condition_from_name": "auto"},
        "analysis": {"metadata": str(md)}})
    assert "INFO" in line, line
    assert "WARN" not in line
    assert "does name every sample" in line
    # the note still records that it was NOT derived, which is the claim
    # design_record.txt has to carry
    assert "not derived" in note


def test_no_metadata_at_all_is_still_a_warning(ma, capsys):
    line, _ = _levels_for(ma, capsys, {"tmt": {"condition_from_name": "auto"}})
    assert "WARN" in line, line
    assert "supply it in analysis.metadata" in line


def test_metadata_that_misses_samples_is_still_a_warning(ma, tmp_path,
                                                        capsys):
    md = tmp_path / "metadata.tsv"
    md.write_text("sample\tgroup\nMF0000\tresponder\nMF0001\tresponder\n",
                  encoding="utf-8")
    line, _ = _levels_for(ma, capsys, {
        "tmt": {"condition_from_name": "auto"},
        "analysis": {"metadata": str(md)}})
    assert "WARN" in line, line
    assert "does not name 4 of these samples" in line


def test_metadata_pointing_at_nothing_is_still_a_warning(ma, tmp_path,
                                                         capsys):
    line, _ = _levels_for(ma, capsys, {
        "tmt": {"condition_from_name": "auto"},
        "analysis": {"metadata": str(tmp_path / "absent.tsv")}})
    assert "WARN" in line, line
    assert "which does not exist" in line


def test_metadata_without_the_sample_column_is_still_a_warning(ma, tmp_path,
                                                               capsys):
    md = tmp_path / "metadata.tsv"
    md.write_text("subject\tgroup\nMF0000\tresponder\n", encoding="utf-8")
    line, _ = _levels_for(ma, capsys, {
        "tmt": {"condition_from_name": "auto"},
        "analysis": {"metadata": str(md)}})
    assert "WARN" in line, line
    assert "has no 'sample' column" in line


def test_an_empty_condition_spec_follows_the_same_rule(ma, tmp_path, capsys):
    # the other caller of _metadata_note: tmt.condition_from_name left empty
    md = tmp_path / "metadata.tsv"
    md.write_text("sample\tgroup\n"
                  + "".join(f"MF{i:04d}\tresponder\n" for i in range(6)),
                  encoding="utf-8")
    line, _ = _levels_for(ma, capsys, {
        "tmt": {"condition_from_name": ""},
        "analysis": {"metadata": str(md)}})
    assert "INFO" in line, line
    line, _ = _levels_for(ma, capsys, {"tmt": {"condition_from_name": ""}})
    assert "WARN" in line, line


def test_the_note_is_a_pair_so_the_level_cannot_drift_from_the_prose(ma):
    # the defect was that the prose and the level were decided separately.
    # _metadata_note returns both now, so a new branch cannot add a message
    # without also saying whether it is recoverable.
    import inspect
    src = inspect.getsource(ma._tmt_add_condition)
    body = src[src.index("def _metadata_note"):src.index("if not spec:")]
    for ret in [l for l in body.splitlines() if "return" in l]:
        pass                      # returns span lines; check the pairing below
    assert body.count("return") == 6, body.count("return")
    assert "False)" in body and "True)" in body
    assert '"WARN"' not in body, "the note must not choose its own level"
