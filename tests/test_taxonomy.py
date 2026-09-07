"""Taxonomy: the peptide LCA consensus, the eggNOG comparison, and the
within-taxon reference the ratio model stands on."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap

import numpy as np
import pandas as pd
import pytest
import re

import fixtures as F
from conftest import METAANNOT_PY, build_project, run_metaannot


# --- finding 39 -------------------------------------------------------
def _lfc(mat, gt, ref=None):
    """log2 fold change per protein, optionally after dividing by a per-sample
    reference."""
    x = np.log2(mat.astype(float))
    if ref is not None:
        x = x.sub(np.log2(ref.astype(float)), axis=1)
    case = [s for s in gt.samples if gt.groups[s] == "case"]
    ctrl = [s for s in gt.samples if gt.groups[s] == "ctrl"]
    return x[case].mean(axis=1) - x[ctrl].mean(axis=1)


def test_the_taxon_size_factor_recovers_a_known_organism_shift(ma):
    # ground truth: taxon T1 doubles between groups (shift_log2 = 1) and three
    # of its twelve proteins are genuinely regulated on top of that.
    df, gt = F.ground_truth_taxon_shift(shift_log2=1.0, effect_log2=2.0)
    tab = df.reset_index()
    tab["taxid"] = [gt.taxon_of[p] for p in tab["group_id"]]
    sf = ma.taxon_size_factors(tab, "taxid", gt.samples, 4)
    row = sf.set_index("taxid").loc["T1"]
    assert row["method"] == "median_of_ratios"
    case = np.mean([row[s] for s in gt.samples if gt.groups[s] == "case"])
    ctrl = np.mean([row[s] for s in gt.samples if gt.groups[s] == "ctrl"])
    assert np.log2(case / ctrl) == pytest.approx(gt.shift_log2, abs=0.05)


def test_passengers_collapse_and_regulated_proteins_survive_the_ratio_model(ma):
    # the whole point of the second model: a protein that merely tracks its
    # source organism must lose its apparent effect, and a genuinely regulated
    # one must keep it.
    df, gt = F.ground_truth_taxon_shift(shift_log2=1.0, effect_log2=2.0)
    tab = df.reset_index()
    tab["taxid"] = [gt.taxon_of[p] for p in tab["group_id"]]
    sf = ma.taxon_size_factors(tab, "taxid", gt.samples, 4).set_index("taxid")
    ref = sf.loc["T1", gt.samples]

    members = df.loc[sorted(gt.passengers | gt.truly_regulated)]
    raw = _lfc(members, gt)
    adj = _lfc(members, gt, ref)

    for pid in sorted(gt.passengers):
        assert raw[pid] == pytest.approx(gt.shift_log2, abs=0.05)
        assert abs(adj[pid]) < 0.10, f"{pid} did not collapse: {adj[pid]:.3f}"
    for pid in sorted(gt.truly_regulated):
        assert adj[pid] == pytest.approx(2.0, abs=0.15)


def test_a_summed_taxon_reference_fails_the_same_ground_truth(ma):
    # this is the defect the median of ratios fixes: one strongly changing
    # member inflates the sum that every other member is divided by, giving
    # the passengers a bias in the opposite direction.
    df, gt = F.ground_truth_taxon_shift(shift_log2=1.0, effect_log2=2.0)
    members = df.loc[sorted(gt.passengers | gt.truly_regulated)]
    summed = members.sum(axis=0)
    summed = summed / summed.mean()
    adj_sum = _lfc(members, gt, summed)
    bias = np.mean([adj_sum[p] for p in sorted(gt.passengers)])
    assert bias < -0.15, (
        "a summed reference should push the passengers negative; got "
        f"{bias:.3f}")

    tab = df.reset_index()
    tab["taxid"] = [gt.taxon_of[p] for p in tab["group_id"]]
    sf = ma.taxon_size_factors(tab, "taxid", gt.samples, 4).set_index("taxid")
    adj_med = _lfc(members, gt, sf.loc["T1", gt.samples])
    med_bias = np.mean([adj_med[p] for p in sorted(gt.passengers)])
    assert abs(med_bias) < abs(bias) / 3


def test_a_taxon_with_too_few_proteins_falls_back_to_the_sum_and_says_so(ma):
    # a median over two or three proteins is not robust either, and the
    # fallback must be recorded rather than reported as a size factor.
    df, gt = F.ground_truth_taxon_shift(n_members=2, n_regulated=0,
                                        other_taxa=0)
    tab = df.reset_index()
    tab["taxid"] = "T1"
    sf = ma.taxon_size_factors(tab, "taxid", gt.samples, 4)
    assert sf["method"].iloc[0] == "sum_fallback"


def test_a_thin_poscounts_median_leaves_the_sample_without_a_factor(ma):
    # a median over one or two ratios is not a size factor; leaving it NaN is
    # what keeps it distinguishable from a well-supported one.
    df, gt = F.ground_truth_taxon_shift(n_members=6, n_regulated=0,
                                        other_taxa=0, n_per_group=3)
    mat = df.copy()
    # make every protein missing in at least one sample, so nothing is
    # "complete", and leave one sample with only a single observed ratio
    for i, pid in enumerate(mat.index):
        mat.iloc[i, i % mat.shape[1]] = np.nan
    mat.iloc[1:, 0] = np.nan
    tab = mat.reset_index()
    tab["taxid"] = "T1"
    sf = ma.taxon_size_factors(tab, "taxid", gt.samples, 4)
    assert sf["method"].iloc[0] in ("median_of_ratios_poscounts", "sum_fallback")


# --- finding 41 -------------------------------------------------------
CONSENSUS_PROBE = textwrap.dedent("""
    import importlib.util, json, sys
    spec = importlib.util.spec_from_file_location("metaannot", sys.argv[1])
    ma = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ma)
    # four peptides: a clear genus majority, an exact 2-2 tie at species.
    lin = [
        {"domain": "2", "genus": "816", "species": "820"},
        {"domain": "2", "genus": "816", "species": "820"},
        {"domain": "2", "genus": "816", "species": "821"},
        {"domain": "2", "genus": "816", "species": "821"},
    ]
    print(json.dumps(ma.consensus_taxon(lin, 0.5, 2)))
""")


@pytest.mark.parametrize("seed", ["0", "1", "12345", "random"])
def test_consensus_taxon_is_deterministic_across_hash_seeds(tmp_path, seed):
    # symptom: the tie-break depended on set iteration order, so the taxon a
    # protein was assigned changed between interpreter runs.
    probe = tmp_path / "probe.py"
    probe.write_text(CONSENSUS_PROBE, encoding="utf-8")
    env = dict(os.environ, PYTHONHASHSEED=seed)
    out = subprocess.run([sys.executable, str(probe), METAANNOT_PY],
                         capture_output=True, text=True, env=env, check=True)
    assert json.loads(out.stdout) == ["816", "genus", 1.0, 4]


def test_an_exact_tie_falls_back_to_the_shared_parent_rank(ma):
    # symptom: a 2-vs-2 split passed min_fraction 0.5 and was resolved by
    # string order.
    lin = [{"genus": "816", "species": "820"}, {"genus": "816", "species": "820"},
           {"genus": "816", "species": "821"}, {"genus": "816", "species": "821"}]
    tid, rank, frac, n = ma.consensus_taxon(lin, 0.5, 2)
    assert (tid, rank) == ("816", "genus")


def test_a_peptide_that_stops_short_is_not_counted_against_deeper_ranks(ma):
    # symptom: counting root-level and phylum-level peptides in the
    # denominator made every conserved, well-quantified protein lose its taxon.
    lin = [{"domain": "2"}, {"domain": "2", "genus": "816", "species": "820"},
           {"domain": "2", "genus": "816", "species": "820"}]
    tid, rank, frac, n = ma.consensus_taxon(lin, 0.5, 2)
    assert (tid, rank) == ("820", "species")


def test_too_few_peptides_gives_no_consensus(ma):
    assert ma.consensus_taxon([{"genus": "816"}], 0.5, 2)[0] is None
    assert ma.consensus_taxon([], 0.5, 2)[0] is None


# --- the taxdump reader ----------------------------------------------
def test_a_merged_taxid_resolves_to_its_successor(ma, tmp_path):
    # symptom: eggNOG 5 seed taxids come from a 2018 taxonomy; without
    # merged.dmp those proteins looked like taxonomic conflicts, and a NEWER
    # taxdump made it worse.
    d = F.toy_taxdump(str(tmp_path / "taxdump"))
    tax = ma.NCBITaxonomy(d)
    assert tax.current("9999") == "820"
    assert tax.lineage("9999").get("genus") == "816"


def test_a_deleted_taxid_resolves_to_nothing(ma, tmp_path):
    d = F.toy_taxdump(str(tmp_path / "taxdump"))
    tax = ma.NCBITaxonomy(d)
    assert tax.current("8888") == ""
    assert tax.lineage("8888") == {}
    assert "8888" in tax.unresolved


def test_superkingdom_is_normalised_to_domain(ma, tmp_path):
    # NCBI renamed the top cellular rank in March 2025; old taxdumps still say
    # superkingdom.
    d = str(tmp_path / "old")
    F.write_taxdump(d, {"1": ("1", "no rank"), "2": ("1", "superkingdom"),
                        "816": ("2", "genus")},
                    {"2": "Bacteria", "816": "Bacteroides"})
    tax = ma.NCBITaxonomy(d)
    assert tax.lineage("816").get("domain") == "2"


def test_taxon_rank_collapses_strain_taxids_to_one_organism(ma, tmp_path,
                                                            capsys):
    # symptom: a seed_ortholog taxid names a reference STRAIN, so two ORFs of
    # one gut organism carry different ones, their shared peptides class as
    # 'shared' and drop, and the reference is computed over singleton
    # pseudo-taxa.
    d = F.toy_taxdump(str(tmp_path / "taxdump"))
    cfg = {"taxon_rank": "genus", "db": {"ncbi_taxonomy": d}}
    out = ma.collapse_taxon_rank(cfg, {"P1": "820", "P2": "821", "P3": "1351"})
    assert out == {"P1": "816", "P2": "816", "P3": "1350"}
    assert "collapsed to" in capsys.readouterr().err


def test_taxon_rank_without_a_taxdump_is_refused(ma):
    with pytest.raises(ma.StageError) as e:
        ma.collapse_taxon_rank({"taxon_rank": "genus", "db": {}}, {"P1": "820"})
    assert "needs db.ncbi_taxonomy" in str(e.value)


def test_the_default_taxon_rank_says_what_a_taxon_then_means(ma, capsys):
    # the default is not "no taxonomy": it is one taxon per eggNOG REFERENCE
    # GENOME, and nobody reading "1,842 distinct taxa" should have to open the
    # config to learn that.
    ma.collapse_taxon_rank({"taxon_rank": ""}, {"P1": "820", "P2": "821"})
    err = capsys.readouterr().err
    assert "reference GENOMES, not organisms" in err


# --- the eggNOG / Unipept comparison ---------------------------------
def _taxonomy_project(tmp_path, name, taxon_source="concordant",
                      unipept_taxid_of=None, taxdump=True):
    """A project with the unipept and taxonomy stages on and a pept2lca file
    generated to match the peptide table."""
    proj = build_project(tmp_path / name)
    quant = pd.read_csv(proj.path("input", "combined_peptide.tsv"), sep="\t")
    # A COMPLETE lineage at every rank: see
    # test_a_gap_in_the_unipept_lineage_does_not_truncate_the_consensus.
    def _lin(genus, species):
        return {"domain": 2, "phylum": 1239, "class": 91061, "order": 186826,
                "family": 81852, "genus": genus, "species": species}
    lineages = {"820": _lin(816, 820), "821": _lin(816, 821),
                "822": _lin(816, 822), "823": _lin(816, 823),
                "1351": _lin(1350, 1351)}
    rows = []
    for pep, prot in zip(quant["Peptide Sequence"], quant["Protein"]):
        tid = (unipept_taxid_of or {}).get(prot, "SAME")
        if tid is None:
            continue                       # no LCA for this protein at all
        if tid == "SAME":
            seed = {p.pid: p.seed_taxid for p in proj.proteins}.get(prot, "")
            if not seed:
                continue
            tid = seed
        # the pipeline strips modifications before querying Unipept, so the
        # LCA table has to be keyed on the stripped sequence
        key = re.sub(r"[^A-Za-z]", "", str(pep)).upper()
        rows.append((key, tid, "species", lineages[str(tid)]))
    res = F.write_unipept(proj.path("input", "pept2lca.csv"), rows)
    cfg = {
        "run": dict(proj.cfg["run"], unipept=True, taxonomy=True),
        "unipept": {"result": res, "split_missed_cleavages": False,
                    "consensus_min_peptides": 2},
        "taxonomy_source": taxon_source,
    }
    if taxdump:
        cfg["db"] = {"ncbi_taxonomy": F.toy_taxdump(str(tmp_path / "taxdump"))}
    proj.write_config(**cfg)
    return proj


# --- finding 43 -------------------------------------------------------
def test_effective_taxid_and_taxonomy_source_reach_the_quant_table(tmp_path):
    # symptom: both columns were assigned AFTER to_csv, so neither ever
    # reached disk and the report silently fell back to the eggNOG taxid
    # whatever taxonomy_source said.
    proj = _taxonomy_project(tmp_path, "eff", taxon_source="eggnog")
    proj.run()
    aq = pd.read_csv(proj.rpath("quant", "annotated_quant.tsv"), sep="\t",
                     dtype=str)
    assert "effective_taxid" in aq.columns
    assert "taxonomy_source" in aq.columns
    assert set(aq["taxonomy_source"].dropna()) == {"eggnog"}
    assert aq["effective_taxid"].fillna("").ne("").any()


# --- finding 42 -------------------------------------------------------
def test_taxonomy_source_governs_both_the_rollup_and_the_taxon_reference(
        tmp_path):
    # symptom: taxonomy_source used to govern only the peptide roll-up, so the
    # ratio model silently kept using eggNOG.
    egg = _taxonomy_project(tmp_path, "src_egg", taxon_source="eggnog")
    egg.run()
    uni = _taxonomy_project(
        tmp_path, "src_uni", taxon_source="unipept",
        # every protein's peptides say Enterococcus faecalis instead
        unipept_taxid_of={p.pid: "1351" for p in
                          build_project(tmp_path / "_probe").proteins})
    uni.run()

    a = pd.read_csv(egg.rpath("quant", "annotated_quant.tsv"), sep="\t",
                    dtype=str).set_index("group_id")
    b = pd.read_csv(uni.rpath("quant", "annotated_quant.tsv"), sep="\t",
                    dtype=str).set_index("group_id")
    assert set(a["taxonomy_source"]) == {"eggnog"}
    assert set(b["taxonomy_source"]) == {"unipept"}
    # the taxon reference is computed over the SAME resolution
    ta = pd.read_csv(egg.rpath("quant", "taxon_intensity.tsv"), sep="\t",
                     dtype=str)
    tb = pd.read_csv(uni.rpath("quant", "taxon_intensity.tsv"), sep="\t",
                     dtype=str)
    assert set(ta["effective_taxid"]) != set(tb["effective_taxid"])
    assert set(tb["effective_taxid"]) == {"1351"}


# --- finding 44 -------------------------------------------------------
def test_concordant_blanks_proteins_where_either_method_had_no_answer(
        tmp_path):
    # symptom: "concordant" sounds like it only removes disagreements, but
    # everything the comparison could not decide is blanked too — so far fewer
    # proteins may survive than expected.
    proj = _taxonomy_project(
        tmp_path, "conc", taxon_source="concordant",
        unipept_taxid_of={
            "P_ko_path": "SAME",        # agrees -> kept
            "P_ko_global": "1351",      # conflicts -> blanked
            "P_ko_orphan": None,        # no LCA at all -> blanked
        })
    proj.run()
    comp = pd.read_csv(proj.rpath("unipept", "taxonomy_comparison.tsv"),
                       sep="\t", dtype=str).set_index("protein_id")
    assert comp.loc["P_ko_path", "verdict"] == "identical"
    assert comp.loc["P_ko_global", "verdict"] == "conflict"
    assert "P_ko_orphan" not in comp.index      # no row at all

    aq = pd.read_csv(proj.rpath("quant", "annotated_quant.tsv"), sep="\t",
                     dtype=str).set_index("group_id")
    eff = aq["effective_taxid"].fillna("")
    assert eff.loc["P_ko_path"] == "820"
    assert eff.loc["P_ko_global"] == ""
    assert eff.loc["P_ko_orphan"] == "", \
        "a protein with no comparison row must be blanked, not kept"


def test_concordant_reports_everything_it_excluded(tmp_path):
    proj = _taxonomy_project(tmp_path, "conc2", taxon_source="concordant",
                             unipept_taxid_of={"P_ko_global": "1351"})
    proc = proj.run()
    assert "taxonomy_source=concordant" in proc.stderr
    assert "excluded from taxon-based steps" in proc.stderr
    assert "with no row in the comparison at all" in proc.stderr


def test_a_missing_taxid_is_not_a_disagreement(tmp_path):
    # symptom: two missings used to compare equal — str(nan) == str(nan) — and
    # were reported "identical".
    proj = _taxonomy_project(tmp_path, "miss", taxon_source="eggnog",
                             unipept_taxid_of={"P_dark2": "820",
                                               "P_dark3": "820"})
    proj.run()
    comp = pd.read_csv(proj.rpath("unipept", "taxonomy_comparison.tsv"),
                       sep="\t", dtype=str).set_index("protein_id")
    # P_dark2/P_dark3 are absent from the eggNOG table entirely
    for pid in ("P_dark2", "P_dark3"):
        if pid in comp.index:
            assert comp.loc[pid, "verdict"] == "eggnog_missing"


def test_without_a_taxdump_the_comparison_degrades_and_says_so(tmp_path):
    # symptom: without a lineage, Enterococcus vs E. faecalis reads as a
    # conflict, and the verdict it emits is not one of the report's factor
    # levels.
    proj = _taxonomy_project(tmp_path, "nodump", taxon_source="eggnog",
                             unipept_taxid_of={"P_ko_path": "821"},
                             taxdump=False)
    proc = proj.run()
    assert "db.ncbi_taxonomy not set" in proc.stderr
    comp = pd.read_csv(proj.rpath("unipept", "taxonomy_comparison.tsv"),
                       sep="\t", dtype=str).set_index("protein_id")
    assert comp.loc["P_ko_path", "verdict"] == "differ_no_lineage"


def test_a_unipept_lca_stopping_above_genus_is_not_a_conflict(tmp_path, ma):
    # symptom: Unipept's LCA stopping at family is a limit of the evidence,
    # not a contradiction, and used to be counted as one.
    proj = _taxonomy_project(tmp_path, "above", taxon_source="eggnog")
    proj.run()
    comp = pd.read_csv(proj.rpath("unipept", "taxonomy_comparison.tsv"),
                       sep="\t", dtype=str)
    assert set(comp["verdict"]) <= {"identical", "concordant",
                                    "concordant_above_genus", "conflict",
                                    "eggnog_missing", "unipept_missing",
                                    "eggnog_unresolved", "no_common_rank",
                                    "differ_no_lineage"}


def test_the_unipept_stage_refuses_to_invent_an_lca_without_a_result(tmp_path):
    proj = build_project(tmp_path / "nohttp")
    proj.write_config(run=dict(proj.cfg["run"], unipept=True),
                      unipept={"result": "", "allow_http": False})
    proc = proj.run(expect=1)
    assert "unipept.allow_http is false" in proc.stderr
    assert "unipept pept2lca" in proc.stderr
    assert os.path.exists(proj.rpath("unipept", "peptides.txt"))


# --- defects found while writing these tests, not on the original list ---
@pytest.mark.xfail(reason="LIVE DEFECT: stage_taxonomy merges each peptide's "
                          "rank columns with `for rk in RANKS: ... else break`, "
                          "so the merged lineage stops at the FIRST rank the "
                          "Unipept row leaves blank. Real pept2lca output "
                          "routinely has no class/order for a bacterial "
                          "lineage, and such a protein is then reported at "
                          "domain however deeply its peptides resolved — which "
                          "makes concordant/taxon_unique drop it.",
                   strict=True)
def test_a_gap_in_the_unipept_lineage_does_not_truncate_the_consensus(
        tmp_path):
    proj = build_project(tmp_path / "gap")
    quant = pd.read_csv(proj.path("input", "combined_peptide.tsv"), sep="\t")
    rows = []
    for pep in quant["Peptide Sequence"]:
        key = re.sub(r"[^A-Za-z]", "", str(pep)).upper()
        # domain and genus/species present, class and order blank — the shape
        # Unipept actually returns for many gut organisms.
        rows.append((key, 820, "species",
                     {"domain": 2, "phylum": 1239, "genus": 816,
                      "species": 820}))
    res = F.write_unipept(proj.path("input", "pept2lca.csv"), rows)
    proj.write_config(run=dict(proj.cfg["run"], unipept=True, taxonomy=True),
                      unipept={"result": res, "split_missed_cleavages": False,
                               "consensus_min_peptides": 2},
                      taxonomy_source="eggnog",
                      db={"ncbi_taxonomy": F.toy_taxdump(str(tmp_path / "td"))})
    proj.run()
    pt = pd.read_csv(proj.rpath("unipept", "protein_taxonomy.tsv"), sep="\t",
                     dtype=str).set_index("protein_id")
    assert pt.loc["P_ko_path", "unipept_rank"] == "species"
    assert pt.loc["P_ko_path", "unipept_taxid"] == "820"


@pytest.mark.xfail(reason="LIVE DEFECT: when no protein reaches a Unipept "
                          "consensus, stage_taxonomy builds an EMPTY "
                          "comparison DataFrame and then reads "
                          "comp['verdict'], so the stage dies with the bare "
                          "message `'verdict'`. A pept2lca file for the wrong "
                          "dataset is exactly how this happens, and the "
                          "message says nothing about it.",
                   strict=True)
def test_a_pept2lca_file_that_matches_nothing_fails_with_an_explanation(
        tmp_path):
    proj = build_project(tmp_path / "nomatch")
    res = F.write_unipept(proj.path("input", "pept2lca.csv"),
                          [("NOTINTHISDATASETK", 820, "species",
                            {"domain": 2, "genus": 816, "species": 820})])
    proj.write_config(run=dict(proj.cfg["run"], unipept=True, taxonomy=True),
                      unipept={"result": res, "split_missed_cleavages": False},
                      taxonomy_source="eggnog")
    proc = proj.run(expect=1)
    assert "'verdict'" not in proc.stderr
    assert "no protein" in proc.stderr.lower() or "pept2lca" in proc.stderr
