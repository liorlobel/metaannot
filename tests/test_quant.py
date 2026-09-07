"""Quantification input: manifests, identifiers, zeros, and the shared-peptide
rule that decides which peptide quantifies which protein."""
from __future__ import annotations

import os

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
