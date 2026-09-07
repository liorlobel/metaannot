"""Every parser, against the malformed input a real tool actually produces.

The pipeline's parsers are deliberately tolerant, which is what makes them
dangerous: a file in the wrong format becomes zero records rather than an
error, and every protein silently loses that evidence.
"""
from __future__ import annotations

import contextlib
import os

import numpy as np
import pandas as pd
import pytest

import fixtures as F


# --- finding 14 -------------------------------------------------------
def test_malformed_tblout_line_leaves_no_phantom_entry(ma, tmp_path):
    # symptom: a defaultdict entry was created before the fields were parsed,
    # so a protein with NO Pfam hit acquired an empty hit list — and all([])
    # is True, which binned it 3d_duf_only.
    path = tmp_path / "pfam.tblout"
    path.write_text(
        F.tblout_line("good_prot", "Peptidase_S8", "PF00082.1") + "\n"
        # 18 columns but a non-numeric E-value
        + " ".join(["bad_prot", "-", "Fam", "PF1", "not-a-number", "x", "0.0",
                    "1e-5", "10", "0", "1.0", "1", "0", "0", "1", "1", "1",
                    "1", "-"]) + "\n"
        # truncated line
        + "short_prot -\n", encoding="utf-8")
    hits = ma.parse_hmm_tblout(str(path))
    assert "good_prot" in hits
    assert "bad_prot" not in hits, "a malformed row created an empty entry"
    assert "short_prot" not in hits
    assert all(v for v in hits.values()), "no entry may be an empty list"


def test_a_protein_with_no_pfam_hit_is_not_duf_only(ma, tmp_path, paths_for):
    # the consequence of finding 14, asserted where it bites: the DUF test is
    # `bool(hits) and all(...)`, so an empty list can never mean "DUF only".
    cfg, p = paths_for()
    proteins = F.protein_set()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), proteins)
    F.write_emapper(p.emapper, proteins)
    F.write_tblout(p.pfam, [("P_ko_path", "Peptidase_S8", "PF00082.1")])
    df = ma.build_annotation(cfg, p)
    assert not bool(df.loc["P_dark1", "duf_only"])
    assert df.loc["P_dark1", "bin"] == "4_dark"


# --- finding 15 -------------------------------------------------------
@pytest.mark.parametrize("writer,name,kind", [
    (lambda p: open(p, "w").write("this is not a tblout at all\n"),
     "pfam.tblout", "HMM hits"),
    (lambda p: open(p, "w").write("garbage\nmore garbage\n"),
     "dbcan.domtblout", "domain hits"),
    (lambda p: open(p, "w").write("nonsense\n"), "kofam.tsv", "KO assignments"),
    (lambda p: open(p, "w").write("nonsense\n"), "interpro.tsv",
     "InterPro matches"),
])
def test_zero_records_from_a_non_empty_file_warns(ma, tmp_path, capsys,
                                                  writer, name, kind):
    # symptom: a corrupted pfam.tblout parsed to zero hits in silence,
    # stripping every Pfam annotation and shifting every bin.
    path = str(tmp_path / name)
    writer(path)
    parse = {"HMM hits": ma.parse_hmm_tblout,
             "domain hits": lambda x: ma.parse_hmm_domtblout(x, 0.35, 1e-15),
             "KO assignments": ma.parse_kofam,
             "InterPro matches": ma.parse_interproscan}[kind]
    assert parse(path) == {} or len(parse(path)) == 0
    err = capsys.readouterr().err
    assert f"parsed 0 {kind} from a non-empty file" in err


def test_an_empty_file_does_not_warn(ma, tmp_path, capsys):
    # no records from an EMPTY file is not a format failure, it is no hits.
    path = str(tmp_path / "pfam.tblout")
    open(path, "w").close()
    ma.parse_hmm_tblout(path)
    assert "parsed 0" not in capsys.readouterr().err


# --- finding 16 -------------------------------------------------------
def test_a_bare_greater_than_in_the_fasta_is_refused(ma, tmp_path):
    # symptom: an empty identifier was accepted, and every such record then
    # collided on the same empty id.
    path = tmp_path / "p.faa"
    path.write_text(">good\nMKV\n>\nMKA\n", encoding="utf-8")
    with pytest.raises(ma.StageError) as e:
        list(ma.read_fasta(str(path)))
    assert "empty identifier" in str(e.value)


def test_fasta_id_is_the_first_whitespace_token(ma, tmp_path):
    path = tmp_path / "p.faa"
    path.write_text(">P1 hypothetical protein [Bacteroides]\nMKV\n",
                    encoding="utf-8")
    assert [pid for pid, _ in ma.read_fasta(str(path))] == ["P1"]


def test_duplicate_fasta_ids_are_reported_and_the_first_is_kept(ma, tmp_path,
                                                                paths_for,
                                                                capsys):
    cfg, p = paths_for()
    path = tmp_path / "p.faa"
    path.write_text(">P1\nMKVAA\n>P1\nMKAAA\n>P2\nMKW\n", encoding="utf-8")
    cfg["proteins_faa"] = str(path)
    df = ma.build_annotation(cfg, p)
    assert list(df.index.sort_values()) == ["P1", "P2"]
    assert "1 duplicate FASTA ids" in capsys.readouterr().err


# --- finding 17 -------------------------------------------------------
MALFORMED = {
    "empty": "",
    "header_only": "# only a header line\n",
    "truncated": "P1\t-\tFam\n",
    "ragged": "P1\t-\tFam\tPF1\t1e-5\t10\nP2\n\nP3\t\t\t\t\t\n",
    "non_numeric": "P1 - Fam PF1 x y 0 0 0 0 0 0 0 0 0 0 0 0 -\n",
    "crlf": "P1 - Fam PF1 1e-5 10.0 0.0 1e-5 10.0 0.0 1.0 1 0 0 1 1 1 1 -\r\n",
    "binary_ish": "\x00\x01\x02 not text\n",
}


@pytest.mark.parametrize("case", sorted(MALFORMED))
@pytest.mark.parametrize("fn", [
    "parse_hmm_tblout", "parse_kofam", "parse_interproscan", "parse_cluster",
    "parse_signalp6", "parse_tmbed", "parse_hmm_lib_desc",
])
def test_every_parser_returns_cleanly_on_malformed_input(ma, tmp_path, fn,
                                                         case):
    # symptom: a parser that raises turns a bad file into a traceback three
    # stages later instead of a warning at the point it was read.
    path = tmp_path / f"{fn}_{case}.txt"
    path.write_bytes(MALFORMED[case].encode("utf-8", "replace"))
    out = getattr(ma, fn)(str(path))
    assert out is not None


@pytest.mark.parametrize("case", ["truncated", "ragged", "non_numeric",
                                  "binary_ish"])
def test_emapper_without_a_query_header_says_so(ma, tmp_path, case):
    path = tmp_path / f"e_{case}.annotations"
    path.write_bytes(MALFORMED[case].encode("utf-8", "replace"))
    with pytest.raises(ValueError) as e:
        ma.parse_emapper(str(path))
    assert "#query" in str(e.value)


@pytest.mark.xfail(reason="LIVE DEFECT: parse_emapper raises a bare "
                          "KeyError(\"None of ['protein_id'] are in the "
                          "columns\") on an emapper.annotations file with no "
                          "data rows, instead of the clear 'no #query header' "
                          "message it gives every other malformed file. "
                          "Reachable by adopting a zero-row eggnog table.",
                   raises=KeyError, strict=True)
@pytest.mark.parametrize("case", ["empty", "header_only"])
def test_emapper_with_no_rows_returns_cleanly(ma, tmp_path, case):
    path = tmp_path / f"e_{case}.annotations"
    path.write_text(MALFORMED[case], encoding="utf-8")
    assert ma.parse_emapper(str(path)) is not None


@pytest.mark.parametrize("case", sorted(MALFORMED))
def test_thresholded_parsers_return_cleanly_on_malformed_input(ma, tmp_path,
                                                               case):
    path = tmp_path / f"x_{case}.txt"
    path.write_bytes(MALFORMED[case].encode("utf-8", "replace"))
    assert ma.parse_hmm_domtblout(str(path), 0.35, 1e-15) is not None
    assert ma.parse_diamond(str(path), 1e-10, 50, 30) is not None
    assert ma.parse_foldseek(str(path), 1e-3, 0.9, 0.5) is not None
    assert ma.parse_hhr_dir(str(tmp_path), 90.0) is not None


def test_crlf_tblout_still_parses(ma, tmp_path):
    path = tmp_path / "pfam.tblout"
    body = F.tblout_line("P1", "Peptidase_S8", "PF00082.1")
    path.write_bytes((body + "\r\n").encode("utf-8"))
    hits = ma.parse_hmm_tblout(str(path))
    assert "P1" in hits
    assert hits["P1"][0][0] == "Peptidase_S8"


def test_missing_files_are_not_an_error_for_the_optional_parsers(ma, tmp_path):
    absent = str(tmp_path / "nope.tsv")
    assert ma.parse_diamond(absent, 1e-10, 50, 30) == {}
    assert ma.parse_foldseek(absent, 1e-3, 0.9, 0.5) == {}
    assert ma.parse_context(absent) == {}


# --- specific parser behaviours the bins depend on --------------------
def test_hmm_tblout_description_is_only_read_when_all_fixed_columns_exist(ma,
                                                                          tmp_path):
    # a truncated line must not have its 7th field read as free text.
    path = tmp_path / "jack.tblout"
    path.write_text(
        F.tblout_line("t1", "q1", "-", desc="Uncharacterized protein") + "\n"
        + "t2 - q2 - 1e-5 10.0 0.0\n", encoding="utf-8")
    hits = ma.parse_hmm_tblout(str(path))
    assert hits["t1"][0][4] == "Uncharacterized protein"
    assert hits["t2"][0][4] == ""


def test_hmm_tblout_dash_description_becomes_empty(ma, tmp_path):
    path = tmp_path / "p.tblout"
    path.write_text(F.tblout_line("t1", "q1", "-", desc="-") + "\n",
                    encoding="utf-8")
    assert ma.parse_hmm_tblout(str(path))["t1"][0][4] == ""


def test_dbcan_hmm_suffix_is_stripped(ma, tmp_path):
    # symptom: dbCAN NAME fields end in '.hmm', which made dbcan_hits
    # string-incomparable to eggNOG's CAZy column.
    path = str(tmp_path / "d.domtblout")
    F.write_domtblout(path, [("P1", "GH13.hmm", 100, 1, 90, 1e-20)])
    hits = ma.parse_hmm_domtblout(path, 0.35, 1e-15)
    assert hits["P1"][0][0] == "GH13"


def test_diamond_rows_in_the_default_12_column_format_are_refused(ma, tmp_path,
                                                                  capsys):
    # symptom: a hand-run DIAMOND with the default --outfmt 6 has 12 numeric
    # columns, so every field would be read from the wrong position.
    path = tmp_path / "vfdb.tsv"
    path.write_text("\t".join(["P1", "VFG1", "88.0", "150", "2", "0", "1",
                               "150", "1", "150", "1e-40", "300"]) + "\n",
                    encoding="utf-8")
    assert ma.parse_diamond(str(path), 1e-10, 50, 30) == {}
    assert "not in the column format" in capsys.readouterr().err


def test_kofam_rows_below_the_family_threshold_are_not_ko_calls(ma, tmp_path):
    path = tmp_path / "kofam.tsv"
    path.write_text("*\tP1\tK00001\t100\t150.0\t1e-40\t\"alcohol dehydrogenase\"\n"
                    " \tP2\tK00002\t100\t 10.0\t1e-01\t\"below threshold\"\n",
                    encoding="utf-8")
    out = ma.parse_kofam(str(path))
    assert "P1" in out and "P2" not in out
    assert out["P1"][0][3] == "alcohol dehydrogenase", "quotes must be stripped"


def test_tmbed_counts_abutting_helices_as_two_segments(ma, tmp_path):
    # symptom: TMbed encodes orientation in the CASE of the label, so a
    # hairpin's two helices abut; counting runs of the class merged them.
    path = str(tmp_path / "t.pred")
    F.write_tmbed(path, {"P1": "iiHHHHhhhhiii", "P2": "iiBBBBbbbbiii"})
    out = ma.parse_tmbed(path)
    assert out["P1"] == (2, 0)
    assert out["P2"] == (0, 2)


def test_interproscan_rows_without_iprlookup_still_parse(ma, tmp_path):
    path = tmp_path / "ip.tsv"
    path.write_text("P1\tmd5\t300\tPfam\tPF00082\tPeptidase S8\t10\t120\t"
                    "1.0E-20\tT\t01-01-2026\n", encoding="utf-8")
    rec = ma.parse_interproscan(str(path))["P1"][0]
    assert rec["analysis"] == "Pfam" and rec["ipr"] == ""


def test_hhr_best_hit_is_found_even_when_the_name_fills_the_column(ma,
                                                                   tmp_path):
    # symptom: the hit table is printed as '%-30.30s', so a name that fills
    # all 30 columns left one space before Prob and a mandatory description
    # group could not match — the BEST hit was silently skipped.
    d = tmp_path / "hh"
    F.write_hhr(str(d), "P1", [
        ("A" * 30, "", 100.0, 1e-30),
        ("PF00082.20", "Peptidase S8", 99.0, 1e-20),
    ])
    out = ma.parse_hhr_dir(str(d), 90.0)
    assert out["P1"][0] == "A" * 30
    assert out["P1"][1] == 100.0


def test_foldseek_target_priority_beats_raw_bitscore(ma, tmp_path, capsys):
    # symptom: an unannotated AFDB50 model routinely outscores a described
    # Swiss-Prot hit, and the description is what the report reads.
    path = str(tmp_path / "hits.tsv")
    F.write_foldseek(path, [
        ("P1", "AF-Q9X0-F1", 1e-20, 900, 0.99, 0.8, "AF-Q9X0-F1", "afdb50"),
        ("P1", "sp|P0A0|TOX", 1e-15, 400, 0.95, 0.7, "aerolysin toxin",
         "swissprot"),
    ])
    best = ma.parse_foldseek(path, 1e-3, 0.9, 0.5,
                             target_priority=["swissprot", "afdb50"])
    assert best["P1"][0] == "sp|P0A0|TOX"
    # without a preference the higher bitscore wins, as before
    best2 = ma.parse_foldseek(path, 1e-3, 0.9, 0.5)
    assert best2["P1"][0] == "AF-Q9X0-F1"


def test_foldseek_query_chain_suffix_is_stripped(ma, tmp_path):
    path = str(tmp_path / "hits.tsv")
    F.write_foldseek(path, [("P1.pdb_A", "T1", 1e-20, 900, 0.99, 0.8, "hdr",
                             "pdb")])
    assert "P1" in ma.parse_foldseek(path, 1e-3, 0.9, 0.5)


def test_legacy_foldseek_rows_warn_that_the_tm_gate_is_weaker(ma, tmp_path,
                                                              capsys):
    path = str(tmp_path / "hits.tsv")
    F.write_foldseek(path, [("P1", "T1", 1e-20, 900, 0.99, 0.8, "hdr", "pdb")])
    ma.parse_foldseek(path, 1e-3, 0.9, 0.5)
    assert "10-column Foldseek format" in capsys.readouterr().err


# --- finding 18 -------------------------------------------------------
@contextlib.contextmanager
def string_dtype(enabled):
    """Run a block under pandas' new string dtype, or under the old one.

    The defect this section covers is a difference BETWEEN pandas versions, so
    testing it on whichever pandas happens to be installed tests half of it.
    `future.infer_string` turns on the pandas-3 behaviour today, which is what
    lets one machine cover both.
    """
    fut = getattr(pd.options, "future", None)
    if enabled and not hasattr(fut, "infer_string"):
        pytest.skip("this pandas has no future.infer_string switch")
    old = getattr(fut, "infer_string", None) if fut is not None else None
    try:
        if fut is not None and hasattr(fut, "infer_string"):
            pd.options.future.infer_string = enabled
        yield
    finally:
        if old is not None:
            pd.options.future.infer_string = old


@pytest.mark.parametrize("new_string_dtype", [False, True])
def test_all_na_columns_do_not_break_the_feature_key(ma, new_string_dtype):
    # symptom: pandas >= 3 gives astype(str) a StringDtype whose NA is a real
    # float nan, so "_".join(row) raised TypeError. FragmentIon and
    # ProductCharge in an MSstats export are usually all-NA.
    with string_dtype(new_string_dtype):
        _assert_feature_key_survives_all_na(ma)


def _assert_feature_key_survives_all_na(ma):
    df = pd.DataFrame({"PeptideSequence": ["PEPK", "TIDEK", "PEPK"],
                       "PrecursorCharge": [2, 3, 2],
                       "FragmentIon": [np.nan, np.nan, np.nan],
                       "ProductCharge": [None, None, None]})
    cols = ["PeptideSequence", "PrecursorCharge", "FragmentIon",
            "ProductCharge"]
    out = ma.join_cols(df, cols)
    # What is asserted is that a key comes back at all, and that it still has
    # one slot per column. How a MISSING value is spelled inside it depends on
    # the pandas version — older ones render it "nan"/"None", newer ones keep
    # it NA so fillna("") empties the slot — and nothing in the pipeline reads
    # the filler back, so pinning the spelling only breaks the test on a
    # different pandas.
    assert all(isinstance(x, str) for x in out)
    assert all(x.count("_") == len(cols) - 1 for x in out), \
        "a missing column must still occupy its slot in the feature key"
    assert out.iloc[0] == out.iloc[2], "identical rows must give one key"
    assert out.iloc[0] != out.iloc[1]


@pytest.mark.parametrize("new_string_dtype", [False, True])
def test_a_partly_missing_column_still_separates_two_features(
        ma, new_string_dtype):
    # the slot matters: two features differing only in a column that is
    # missing for one of them must not collapse into one.
    with string_dtype(new_string_dtype):
        df = pd.DataFrame({"PeptideSequence": ["PEPK", "PEPK"],
                           "FragmentIon": ["y1", None]})
        out = ma.join_cols(df, ["PeptideSequence", "FragmentIon"])
        assert out.iloc[0] != out.iloc[1]
        assert all(isinstance(x, str) for x in out)


def test_msstats_csv_with_all_na_fragment_columns_reads(ma, tmp_path):
    # the same defect where it bites: the whole MSstats path.
    path = tmp_path / "MSstats.csv"
    rows = ["ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,"
            "ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,Intensity"]
    for i, run in enumerate(["r1", "r2"]):
        rows.append(f"P1,PEPTIDEK,2,NA,NA,L,A,{i+1},{run},{1000*(i+1)}")
        rows.append(f"P2,TIDEPEPK,2,NA,NA,L,B,{i+1},{run},{2000*(i+1)}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    cfg = dict(ma.DEFAULT_CONFIG)
    feats, int_cols, design = ma.read_feature_table(str(path), "msstats_csv",
                                                    cfg)
    assert sorted(int_cols) == ["r1", "r2"]
    assert set(feats["razor_protein"]) == {"P1", "P2"}
    assert design is not None and "group" in design.columns


# --- finding 19 -------------------------------------------------------
@pytest.mark.parametrize("modern", [True, False])
@pytest.mark.parametrize("sep", [",", "\t"])
def test_unipept_reads_both_the_cli_and_the_legacy_gem_columns(ma, tmp_path,
                                                               modern, sep):
    # symptom: the current CLI writes domain_id where the old Ruby gem wrote
    # superkingdom_id, which is what the parser read.
    path = str(tmp_path / "pept2lca.csv")
    F.write_unipept(path, [("PEPTIDEK", 820, "species",
                            {"domain": 2, "genus": 816, "species": 820})],
                    modern=modern, sep=sep)
    df = ma.read_unipept_result(path)
    assert "domain_id" in df.columns
    assert int(df["domain_id"].iloc[0]) == 2
    assert int(df["species_id"].iloc[0]) == 820


def test_unipept_export_without_taxid_columns_is_refused(ma, tmp_path):
    # symptom: the web export is name-based; it used to pass silently and give
    # every protein an empty lineage, which the verdict logic called a
    # disagreement.
    path = tmp_path / "web_export.csv"
    path.write_text("peptide,taxon_name,genus,species\n"
                    "PEPTIDEK,Bacteroides uniformis,Bacteroides,B. uniformis\n",
                    encoding="utf-8")
    with pytest.raises(ma.StageError) as e:
        ma.read_unipept_result(str(path))
    assert "not a pept2lca result" in str(e.value)


def test_a_bare_rank_column_is_never_taken_as_a_taxid(ma, tmp_path):
    # symptom: a bare `genus` column holds the NAME; accepting it would put
    # species names into the taxid comparison.
    path = tmp_path / "p.csv"
    path.write_text("peptide,taxon_id,genus,species\n"
                    "PEPTIDEK,820,Bacteroides,Bacteroides uniformis\n",
                    encoding="utf-8")
    with pytest.raises(ma.StageError) as e:
        ma.read_unipept_result(str(path))
    assert "no rank taxid columns" in str(e.value)


def test_pept2lca_without_the_all_flag_is_refused(ma, tmp_path):
    path = tmp_path / "p.csv"
    path.write_text("peptide,taxon_id,taxon_name,taxon_rank\n"
                    "PEPTIDEK,820,B. uniformis,species\n", encoding="utf-8")
    with pytest.raises(ma.StageError) as e:
        ma.read_unipept_result(str(path))
    assert "--all" in str(e.value)


# --- finding 20 -------------------------------------------------------
def test_positive_values_without_a_column_to_test_is_refused(ma, tmp_path):
    # symptom: positive_values was tested against a None column when score_col
    # was omitted, so it matched nothing and the predictor was silently zero.
    path = tmp_path / "bastion.tsv"
    path.write_text("id\tlabel\nP1\teffector\nP2\tnon-effector\n",
                    encoding="utf-8")
    spec = {"bastion6": {"file": str(path), "id_col": "id",
                         "positive_values": ["effector"]}}
    with pytest.raises(ma.StageError) as e:
        ma.parse_external_predictions(spec, pd.Index(["P1", "P2"]))
    assert "positive_values needs score_col" in str(e.value)


def test_positive_values_with_a_column_selects_only_those_rows(ma, tmp_path):
    path = tmp_path / "bastion.tsv"
    path.write_text("id\tlabel\nP1\teffector\nP2\tnon-effector\n",
                    encoding="utf-8")
    spec = {"bastion6": {"file": str(path), "id_col": "id",
                         "score_col": "label",
                         "positive_values": ["effector"]}}
    cols = ma.parse_external_predictions(spec, pd.Index(["P1", "P2"]))
    assert list(cols["pred_bastion6"]) == [True, False]


def test_a_misspelled_score_column_is_refused_not_treated_as_a_positives_list(
        ma, tmp_path):
    # symptom: the fallthrough marked every listed protein an effector.
    path = tmp_path / "b.tsv"
    path.write_text("id\tscore\nP1\t0.9\n", encoding="utf-8")
    spec = {"b": {"file": str(path), "id_col": "id", "score_col": "scoer"}}
    with pytest.raises(ma.StageError) as e:
        ma.parse_external_predictions(spec, pd.Index(["P1"]))
    assert "is not in" in str(e.value)


def test_a_spec_with_no_usable_column_is_refused_unless_all_positive_is_set(
        ma, tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text("id\tscore\nP1\t0.9\n", encoding="utf-8")
    spec = {"b": {"file": str(path), "id_col": "id"}}
    with pytest.raises(ma.StageError):
        ma.parse_external_predictions(spec, pd.Index(["P1"]))


def test_a_headerless_id_list_is_not_split_on_a_character_of_its_own_ids(
        ma, tmp_path, capsys):
    # symptom: csv.Sniffer on a one-id-per-line file picked a character out of
    # the ids themselves, so every protein became a negative.
    path = tmp_path / "ids.txt"
    path.write_text("OIDECCNN_00001\nOIDECCNN_00002\n", encoding="utf-8")
    df = ma._read_prediction_table(str(path), "t4sepp")
    assert list(df["id"]) == ["OIDECCNN_00001", "OIDECCNN_00002"]


def test_a_predictor_keyed_on_foreign_ids_is_reported(ma, tmp_path, capsys):
    path = tmp_path / "b.tsv"
    path.write_text("id\tscore\nOTHER_1\t0.9\n", encoding="utf-8")
    spec = {"b": {"file": str(path), "id_col": "id", "score_col": "score",
                  "threshold": 0.5}}
    ma.parse_external_predictions(spec, pd.Index(["P1", "P2"]))
    assert "zero id overlap" in capsys.readouterr().err


# --- emapper header handling -----------------------------------------
def test_emapper_2_0_column_names_are_canonicalised(ma, tmp_path):
    # symptom: eggnog-mapper 2.0.x names its first column '#query_name' and
    # spells three annotation columns differently; set_index raised KeyError.
    path = tmp_path / "old.emapper.annotations"
    path.write_text("#query_name\tseed_eggNOG_ortholog\tCOG Functional cat."
                    "\teggNOG free text desc.\n"
                    "P1\t820.ABC\tS\tsomething\n", encoding="utf-8")
    df = ma.parse_emapper(str(path))
    assert df.index.name == "protein_id"
    assert "COG_category" in df.columns and "Description" in df.columns
    assert df.loc["P1", "seed_ortholog"] == "820.ABC"


def test_emapper_dashes_become_empty_strings(ma, tmp_path):
    path = str(tmp_path / "e.annotations")
    F.write_emapper(path, F.protein_set()[:3])
    df = ma.parse_emapper(path)
    assert df.loc["P_ko_orphan", "KEGG_Pathway"] == ""


# --- a defect found by the first full run on real data ----------------
# 0xa0 is a latin-1 non-breaking space and is not valid UTF-8. VFDB subject
# titles, InterPro signature descriptions, HMM DESC lines and FASTA headers
# all carry latin-1 in the wild.
NON_UTF8 = b"\xa0"


@pytest.mark.parametrize("name,payload", [
    ("vfdb.tsv",
     b"P1\tVFG0001\t88.0\t150\t1e-40\t300.0\t90\t80\themolysin" + NON_UTF8 + b"BL\n"),
    ("interproscan.tsv",
     b"P1\tmd5\t300\tPfam\tPF00082\tPeptidase" + NON_UTF8
     + b"S8\t10\t120\t1e-20\tT\t01-01-2026\n"),
    ("kofam.tsv",
     b"*\tP1\tK00001\t100\t150.0\t1e-40\t\"alcohol" + NON_UTF8 + b"dh\"\n"),
    ("pfam.tblout",
     b"P1 - Fam PF1 1e-5 10.0 0.0 1e-5 10.0 0.0 1.0 1 0 0 1 1 1 1 d"
     + NON_UTF8 + b"esc\n"),
    ("ncbifam.lib", b"NAME  TIGR1\nDESC  hypothetical" + NON_UTF8 + b"protein\n//\n"),
    ("proteins.faa", b">P1 desc" + NON_UTF8 + b"here\nMKVAA\n"),
])
def test_a_non_utf8_byte_in_a_tool_output_is_not_fatal(ma, tmp_path, name,
                                                       payload):
    # symptom: one 0xa0 in a search result killed the integrate stage of a real
    # run — after InterProScan had already spent three hours — with a message
    # that named neither the file nor the stage:
    #   FATAL stage 'integrate' failed: 'utf-8' codec can't decode byte 0xa0
    path = tmp_path / name
    path.write_bytes(payload)
    parse = {
        "vfdb.tsv": lambda p: ma.parse_diamond(p, 1e-10, 50, 30),
        "interproscan.tsv": ma.parse_interproscan,
        "kofam.tsv": ma.parse_kofam,
        "pfam.tblout": ma.parse_hmm_tblout,
        "ncbifam.lib": ma.parse_hmm_lib_desc,
        "proteins.faa": lambda p: list(ma.read_fasta(p)),
    }[name]
    out = parse(str(path))
    assert out, f"{name} parsed to nothing"


def test_a_bad_byte_costs_one_character_not_the_record(ma, tmp_path):
    # errors="replace", so the damage is bounded: the record still parses and
    # only the offending character of the free text is lost.
    path = tmp_path / "vfdb.tsv"
    path.write_bytes(b"P1\tVFG0001\t88.0\t150\t1e-40\t300.0\t90\t80\t"
                     b"hemolysin" + NON_UTF8 + b"BL\n")
    hit = ma.parse_diamond(str(path), 1e-10, 50, 30)["P1"]
    assert hit[0] == "VFG0001"
    assert hit[3].startswith("hemolysin") and hit[3].endswith("BL")


def test_an_identifier_is_not_silently_altered_by_the_replacement(ma,
                                                                   tmp_path):
    # the one place a replacement could do real harm is an id, so check that
    # a CLEAN id beside a dirty description comes through untouched.
    path = tmp_path / "p.faa"
    path.write_bytes(b">CDPNAMPK_339076 hypothetical" + NON_UTF8
                     + b" protein\nMKVAA\n")
    assert [pid for pid, _ in ma.read_fasta(str(path))] == ["CDPNAMPK_339076"]


def test_a_gzipped_input_is_read_as_utf8_whatever_the_locale(ma, tmp_path):
    # symptom: the .gz branch passed no encoding at all, so TextIOWrapper used
    # the LOCALE's codec — never UTF-8 by construction. emapper_precomputed is
    # routinely a .gz, making it the input most likely to be read differently
    # on the server than on the laptop.
    import gzip
    path = tmp_path / "cat.emapper.annotations.gz"
    with gzip.open(path, "wb") as fh:
        fh.write("#query\tseed_ortholog\tDescription\n"
                 "P1\t820.SEED\tβ-glucosidase\n".encode("utf-8"))
    df = ma.parse_emapper(str(path))
    assert df.loc["P1", "Description"] == "β-glucosidase"


def test_a_gzipped_input_with_a_bad_byte_is_also_survivable(ma, tmp_path):
    import gzip
    path = tmp_path / "cat.emapper.annotations.gz"
    with gzip.open(path, "wb") as fh:
        fh.write(b"#query\tseed_ortholog\tDescription\n"
                 b"P1\t820.SEED\tprotein" + NON_UTF8 + b"alpha\n")
    assert ma.parse_emapper(str(path)).loc["P1", "seed_ortholog"] == "820.SEED"
