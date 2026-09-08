"""Do the sources that reach the same protein agree about it?

Coverage overlap and concordance are different questions, and the tool used to
answer only the first. These tests pin the second, and in particular pin the
namespace trap that made two searches of ONE library read as 97% conflict.
"""
import pandas as pd
import pytest


def _df(rows):
    """A frame shaped like annotation_final.tsv, with only the columns the
    agreement check reads."""
    cols = ["pfam_accs", "pfam_hits", "pfams_emapper", "ncbifam_accs",
            "ncbifam_hits", "interpro_sigs", "ko", "kofam_ko"]
    return pd.DataFrame([{c: r.get(c, "") for c in cols} for r in rows])


def _row(tab, needle):
    m = tab[tab["comparison"].str.contains(needle)]
    assert len(m) == 1, f"expected one {needle!r} row, got {len(m)}"
    return m.iloc[0]


# ----------------------------------------------------------------------
def test_two_searches_that_agree_are_reported_as_identical(ma):
    tab = ma.source_agreement(_df([
        {"pfam_accs": "PF00005;PF00664",
         "interpro_sigs": "Gene3D:G3DSA:3.40.50.300;Pfam:PF00005;Pfam:PF00664"},
        {"pfam_accs": "PF01234", "interpro_sigs": "Pfam:PF01234"},
    ]))
    r = _row(tab, "pfam: hmmsearch")
    assert r["both"] == 2 and r["identical"] == 2
    assert r["disjoint"] == 0 and r["pct_agree"] == 100.0


def test_one_side_finding_more_is_overlap_not_conflict(ma):
    tab = ma.source_agreement(_df([
        {"pfam_accs": "PF00005;PF00664", "interpro_sigs": "Pfam:PF00005"},
    ]))
    r = _row(tab, "pfam: hmmsearch")
    assert r["overlapping"] == 1 and r["disjoint"] == 0
    assert r["a_superset"] == 1, "the direction has to be recorded, not just the fact"
    assert r["b_superset"] == 0


def test_each_side_holding_something_the_other_lacks_is_counted_apart(ma):
    tab = ma.source_agreement(_df([
        {"pfam_accs": "PF00005;PF11111", "interpro_sigs": "Pfam:PF00005;Pfam:PF22222"},
    ]))
    r = _row(tab, "pfam: hmmsearch")
    assert r["overlapping"] == 1
    assert r["mutually_exclusive"] == 1
    assert r["a_superset"] == 0 and r["b_superset"] == 0


def test_a_real_conflict_is_reported_as_disjoint(ma):
    tab = ma.source_agreement(_df([
        {"pfam_accs": "PF00005", "interpro_sigs": "Pfam:PF99999"},
    ]))
    r = _row(tab, "pfam: hmmsearch")
    assert r["disjoint"] == 1 and r["identical"] == 0
    assert r["pct_disjoint"] == 100.0


def test_a_version_suffix_is_not_a_disagreement(ma):
    """PF00005.29 and PF00005 are one family; counting them apart would make
    the whole table an artefact of formatting."""
    tab = ma.source_agreement(_df([
        {"pfam_accs": "PF00005.29", "interpro_sigs": "Pfam:PF00005"},
    ]))
    assert _row(tab, "pfam: hmmsearch")["identical"] == 1


def test_the_eggnog_ko_prefix_is_not_a_disagreement(ma):
    """eggNOG writes 'ko:K02014'; KOfamScan writes 'K02014'."""
    tab = ma.source_agreement(_df([
        {"ko": "ko:K02014", "kofam_ko": "K02014"},
        {"ko": "ko:K01347,ko:K12684", "kofam_ko": "K12684,K01347"},
    ]))
    r = _row(tab, "ko: eggnog")
    assert r["identical"] == 2 and r["disjoint"] == 0


def test_coverage_by_one_source_alone_is_not_counted_as_agreement(ma):
    tab = ma.source_agreement(_df([
        {"ko": "ko:K00001", "kofam_ko": ""},
        {"ko": "", "kofam_ko": "K00002"},
        {"ko": "ko:K00003", "kofam_ko": "K00003"},
    ]))
    r = _row(tab, "ko: eggnog")
    assert r["a_only"] == 1 and r["b_only"] == 1
    assert r["both"] == 1, "only proteins BOTH sources called can be compared"


# ----------------------------------------------------------------------
# the namespace trap
# ----------------------------------------------------------------------
def test_ncbifam_accessions_make_the_comparison_possible_at_all(ma):
    """Before ncbifam_accs existed there was nothing to compare.

    ncbifam_hits holds family NAMES ('PorV_fam'); InterProScan holds
    ACCESSIONS ('NF033709'). Comparing those reads as total conflict between
    two searches of one library that in fact agree.
    """
    tab = ma.source_agreement(_df([
        {"ncbifam_hits": "PorV_fam;T9SS_OM_PorV",
         "ncbifam_accs": "NF033709;NF033710",
         "interpro_sigs": "NCBIfam:NF033709;NCBIfam:NF033710"},
    ]))
    r = _row(tab, "ncbifam: hmmsearch")
    assert r["identical"] == 1 and r["disjoint"] == 0, \
        "the accession column is what makes these comparable"


def test_a_near_total_conflict_is_called_out_as_a_namespace_problem(ma, capsys,
                                                                   tmp_path):
    """A handful of disjoint calls is ordinary; ~100% is not real disagreement.

    The log has to say so, or the number invites exactly the wrong conclusion.
    """
    rows = [{"ncbifam_accs": f"NAME{i}", "interpro_sigs": f"NCBIfam:NF{i:06d}"}
            for i in range(60)]
    ma.write_source_agreement(_df(rows), str(tmp_path / "agreement.tsv"))
    err = capsys.readouterr().err
    assert "not a credible rate of real disagreement" in err
    assert "names against accessions" in err


def test_an_ordinary_conflict_rate_is_reported_without_the_namespace_warning(
        ma, capsys, tmp_path):
    rows = ([{"ko": f"ko:K{i:05d}", "kofam_ko": f"K{i:05d}"} for i in range(95)]
            + [{"ko": "ko:K00001", "kofam_ko": "K99999"} for _ in range(5)])
    ma.write_source_agreement(_df(rows), str(tmp_path / "agreement.tsv"))
    err = capsys.readouterr().err
    assert "not a credible rate" not in err
    assert "95.0% agree" in err


def test_the_table_is_written_even_when_nothing_is_comparable(ma, tmp_path):
    """A declared stage output that is only sometimes created makes the stage
    look unfinished and rerun forever."""
    out = tmp_path / "agreement.tsv"
    ma.write_source_agreement(_df([{"pfam_accs": "", "interpro_sigs": ""}]),
                              str(out))
    assert out.exists()
    tab = pd.read_csv(out, sep="\t")
    assert list(tab.columns)[:3] == ["comparison", "tests", "a"]
    assert len(tab) == 0


def test_a_missing_column_is_skipped_rather_than_crashing(ma):
    """Not every run has every source; kofam may never have run."""
    df = pd.DataFrame([{"pfam_accs": "PF00005", "interpro_sigs": "Pfam:PF00005"}])
    tab = ma.source_agreement(df)
    assert len(tab) == 1 and "pfam" in tab.iloc[0]["comparison"]
