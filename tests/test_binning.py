"""The evidence bins: the one thing this tool exists to compute."""
from __future__ import annotations

import itertools
import os

import pytest

import fixtures as F


def _annotate(ma, tmp_path, paths_for, proteins, **evidence):
    """Seed a results directory with tool outputs and run build_annotation."""
    cfg, p = paths_for(f"r{abs(hash(tuple(sorted(evidence)))) % 100000}")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), proteins)
    F.write_emapper(p.emapper, proteins)
    for key, val in evidence.items():
        {"pfam": lambda v: F.write_tblout(p.pfam, v),
         "ncbifam": lambda v: F.write_tblout(p.ncbifam, v),
         "dbcan": lambda v: F.write_domtblout(p.dbcan, v),
         "foldseek": lambda v: F.write_foldseek(p.foldseek, v),
         "jackhmmer": lambda v: F.write_tblout(p.jackhmmer, v),
         "kofam": lambda v: F.write_kofam(p.kofam, v),
         "interpro": lambda v: F.write_interpro(p.interpro, v),
         "hhblits": lambda v: [F.write_hhr(p.hhr_dir, q, h) for q, h in v],
         }[key](val)
    return cfg, p, ma.build_annotation(cfg, p)


# --- finding 34 -------------------------------------------------------
def test_an_eggnog_pfam_alone_is_domain_evidence(ma, tmp_path, paths_for):
    # symptom: binning on the pfam stage's hmmsearch hits alone meant that
    # with run.pfam false, a protein eggNOG had already assigned a domain to
    # was reported as having no evidence at all — 36% of a real dark bin.
    _, _, df = _annotate(ma, tmp_path, paths_for, F.protein_set())
    assert df.loc["P_eggpfam", "bin"] == "3_annotated_no_ko"
    assert bool(df.loc["P_eggpfam", "has_seq_annotation"])


@pytest.mark.parametrize("pid", ["P_eggduf", "P_eggupf"])
def test_an_eggnog_pfam_that_is_only_a_duf_lands_in_the_duf_bin(ma, tmp_path,
                                                                paths_for, pid):
    # the other half: DUF1234 and UPF0102 name a family, not a function, so
    # they must not promote a protein to 3_annotated_no_ko.
    _, _, df = _annotate(ma, tmp_path, paths_for, F.protein_set())
    assert df.loc[pid, "bin"] == "3d_duf_only"
    assert not bool(df.loc[pid, "has_seq_annotation"])


def test_an_hmmsearch_upf_family_is_not_counted_as_annotation(ma, tmp_path,
                                                              paths_for):
    # symptom: DUF_RE never matched Pfam's UPF0xxx models, so a UPF-only
    # protein was reported as annotated.
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps,
                         pfam=[("P_dark1", "UPF0102", "PF01894.1")])
    assert df.loc["P_dark1", "bin"] == "3d_duf_only"


def test_an_uninformative_ncbifam_family_does_not_leave_the_dark_bin(
        ma, tmp_path, paths_for):
    # symptom: NCBIfam carries families whose whole DESC is "hypothetical
    # protein"; counting those as annotation promotes a protein on a family
    # name that says nothing. The DESC is read from the library, because
    # --tblout carries the TARGET's description, never the query HMM's.
    ps = F.protein_set()
    cfg, p = paths_for("ncbifam")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    F.write_tblout(p.ncbifam, [("P_dark1", "TIGR00001", "TIGR00001"),
                               ("P_dark2", "TIGR00002", "TIGR00002")])
    cfg["db"]["ncbifam_hmm"] = F.write_hmm_library(
        str(tmp_path / "ncbifam.lib"),
        [("TIGR00001", "TIGR00001", "hypothetical protein"),
         ("TIGR00002", "TIGR00002", "ribosomal protein S12")])
    df = ma.build_annotation(cfg, p)
    assert df.loc["P_dark1", "bin"] == "3d_duf_only"
    assert df.loc["P_dark2", "bin"] == "3_annotated_no_ko"


def test_ncbifam_uninformative_test_can_be_switched_off(ma, tmp_path,
                                                        paths_for, capsys):
    ps = F.protein_set()
    cfg, p = paths_for("ncbifam_off")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    cfg["ncbifam_uninformative_test"] = False
    F.write_emapper(p.emapper, ps)
    F.write_tblout(p.ncbifam, [("P_dark1", "TIGR00001", "TIGR00001")])
    df = ma.build_annotation(cfg, p)
    assert df.loc["P_dark1", "bin"] == "3_annotated_no_ko"
    assert "every NCBIfam hit counts" in capsys.readouterr().err


def test_a_structural_interpro_analysis_cannot_promote_out_of_the_dark_bin(
        ma, tmp_path, paths_for):
    # symptom: InterProScan reports SignalP as SignalP_GRAM_POSITIVE etc.,
    # none of which EQUALS the configured "SignalP", so a signal-peptide-only
    # protein was promoted by the very analysis meant to be ignored.
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps, interpro=[
        ("P_dark1", "SignalP_GRAM_POSITIVE", "SignalP-noTM", "SignalP-noTM",
         "-", "-", "-"),
        ("P_dark2", "Gene3D", "G3DSA:3.40.50", "Rossmann fold", "IPR001",
         "NAD(P)-binding domain", "-"),
    ])
    assert df.loc["P_dark1", "bin"] == "4_dark"
    assert df.loc["P_dark2", "bin"] == "3_annotated_no_ko"


def test_an_interpro_row_with_no_description_is_not_informative(ma, tmp_path,
                                                                paths_for):
    # symptom: about half a real InterProScan TSV carries an empty signature
    # description, and an empty string matches no DUF pattern, so every such
    # row counted as informative on its own.
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps, interpro=[
        ("P_dark1", "SUPERFAMILY", "SSF52540", "-", "-", "-", "-"),
        ("P_dark2", "PANTHER", "PTHR1", "FAMILY NOT NAMED", "-", "-", "-"),
    ])
    assert df.loc["P_dark1", "bin"] == "4_dark"
    assert df.loc["P_dark2", "bin"] == "4_dark"


# --- finding 36 -------------------------------------------------------
def test_a_ko_whose_only_map_is_global_is_an_orphan(ma, tmp_path, paths_for):
    # symptom: "Metabolic pathways" (map01100) is not pathway annotation in
    # any useful sense, so a protein whose only map is a global one must not
    # count as KEGG-visible.
    _, _, df = _annotate(ma, tmp_path, paths_for, F.protein_set())
    assert df.loc["P_ko_global", "bin"] == "2_ko_orphan"
    assert not bool(df.loc["P_ko_global", "kegg_enrichment_visible"])
    assert df.loc["P_ko_path", "bin"] == "1_ko_pathway"
    assert bool(df.loc["P_ko_path", "kegg_enrichment_visible"])


def test_every_global_map_in_the_list_is_excluded(ma):
    assert "map01100" in ma.GLOBAL_MAPS
    # map01130 was retired by KEGG but eggnog-mapper 2.1.x still emits it
    assert "map01130" in ma.GLOBAL_MAPS


# --- finding 35 -------------------------------------------------------
PREDICATES = ("ko", "path", "seq", "duf", "struct", "prof")


def _combination_proteins():
    """One protein per combination of the six deciding predicates."""
    out = []
    for bits in itertools.product([0, 1], repeat=len(PREDICATES)):
        d = dict(zip(PREDICATES, bits))
        pid = "C_" + "".join(str(b) for b in bits)
        out.append((pid, d))
    return out


def _expected_bin(d):
    """BIN_ORDER, applied by hand: first match wins."""
    if d["ko"] and d["path"]:
        return "1_ko_pathway"
    if d["ko"]:
        return "2_ko_orphan"
    if d["seq"]:
        return "3_annotated_no_ko"
    if d["duf"]:
        return "3d_duf_only"
    if d["struct"]:
        return "3s_structure_only"
    if d["prof"]:
        return "3p_profile_only"
    return "4_dark"


def test_bin_assignment_is_total_and_mutually_exclusive(ma, tmp_path,
                                                        paths_for):
    # symptom: the bins are the tool's whole output vocabulary. Every
    # combination of the deciding predicates must land in exactly one of the
    # seven, and no combination carrying evidence may land in 4_dark.
    combos = _combination_proteins()
    proteins, pfam, dbcan, fs, jh = [], [], [], [], []
    for pid, d in combos:
        proteins.append(F.Protein(
            pid, "M" + "A" * 120,
            ko="ko:K01234" if d["ko"] else "",
            pathway=("ko00010,map00010" if d["path"] else ""),
            seed_taxid="820"))
        if d["seq"]:
            dbcan.append((pid, "GH13.hmm", 100, 1, 90, 1e-20))
        if d["duf"]:
            pfam.append((pid, "DUF1234", "PF01234.1"))
        if d["struct"]:
            fs.append((pid, "sp|P0A0|TOX", 1e-20, 400, 0.99, 0.8,
                       "aerolysin toxin", "swissprot"))
        if d["prof"]:
            jh.append((f"UniRef50_{pid}", pid, "-"))

    cfg, p = paths_for("combos")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "combo.faa"), proteins)
    F.write_emapper(p.emapper, proteins)
    F.write_tblout(p.pfam, pfam)
    F.write_domtblout(p.dbcan, dbcan)
    F.write_foldseek(p.foldseek, fs)
    F.write_tblout(p.jackhmmer, jh)
    df = ma.build_annotation(cfg, p)

    assert set(df["bin"]) <= set(ma.BIN_ORDER)
    wrong = []
    for pid, d in combos:
        got, want = df.loc[pid, "bin"], _expected_bin(d)
        if got != want:
            wrong.append((pid, d, got, want))
    assert wrong == [], f"{len(wrong)} combinations bin wrongly: {wrong[:5]}"

    # no evidence may land in 4_dark. `path` alone is not evidence: a
    # KEGG_Pathway cell with no KEGG_ko is contradictory input, so the two
    # combinations that carry only it are dark, and nothing else is.
    dark = set(df.index[df["bin"] == "4_dark"])
    assert dark == {pid for pid, d in combos
                    if not any(d[k] for k in PREDICATES if k != "path")}


def test_the_module_docstring_the_report_and_bin_order_use_one_vocabulary(ma):
    # symptom: a bin missing from the report's BIN_LEVELS becomes NA in every
    # figure, and one missing from BIN_COLS drops out of the legend. Both are
    # silent, so the copies are checked at import.
    levels, cols = ma._rmd_bin_vocabulary(ma.RMD_TEMPLATE)
    assert levels == ma.BIN_ORDER
    assert cols == ma.BIN_ORDER
    for b in ma.BIN_ORDER:
        assert b in ma.__doc__, f"{b} is not in the module docstring"


def test_bin_summary_has_a_row_for_every_bin_in_classifier_order(ma, tmp_path,
                                                                 paths_for):
    # symptom: groupby drops a bin with no proteins, so an empty
    # 3s_structure_only was indistinguishable from a Foldseek stage that never
    # ran — and the rows came out alphabetically, 3p before 3s.
    import pandas as pd
    cfg, p, df = _annotate(ma, tmp_path, paths_for, F.protein_set())
    ma.write_summary(df, p.summary)
    s = pd.read_csv(p.summary, sep="\t")
    assert list(s["bin"])[:-1] == list(ma.BIN_ORDER)
    assert list(s["bin"])[-1] == "TOTAL"
    assert s["n"].dtype.kind == "i", "counts must not be upcast to float"
    assert int(s.loc[s["bin"] == "TOTAL", "n"].iloc[0]) == len(df)


def test_a_bin_name_outside_bin_order_is_fatal(ma, tmp_path, paths_for):
    cfg, p, df = _annotate(ma, tmp_path, paths_for, F.protein_set())
    df = df.copy()
    df.loc[df.index[0], "bin"] = "5_invented"
    with pytest.raises(ma.StageError) as e:
        ma.write_summary(df, p.summary)
    assert "not in BIN_ORDER" in str(e.value)


# --- the rescue bins --------------------------------------------------
def test_an_uncharacterised_fold_or_profile_hit_is_not_a_rescue(ma, tmp_path,
                                                                paths_for,
                                                                capsys):
    # symptom: the default hhblits DB is Pfam and AFDB50 is mostly TrEMBL
    # "Uncharacterized protein", so a dark protein's best hit is often a DUF —
    # and it used to leave 4_dark and be counted as a structure rescue.
    ps = F.protein_set()
    _, _, df = _annotate(
        ma, tmp_path, paths_for, ps,
        foldseek=[("P_dark1", "AF-Q9X0-F1", 1e-20, 400, 0.99, 0.8,
                   "Uncharacterized protein", "afdb50"),
                  ("P_dark2", "sp|P0A0|X", 1e-20, 400, 0.99, 0.8,
                   "aerolysin toxin", "swissprot")])
    assert df.loc["P_dark1", "bin"] == "4_dark"
    assert df.loc["P_dark2", "bin"] == "3s_structure_only"
    assert "do not count as a rescue" in capsys.readouterr().err


def test_a_duf_only_protein_keeps_its_bin_despite_a_fold(ma, tmp_path,
                                                         paths_for, capsys):
    # deliberate: a DUF still names a family, so duf_only is tested before the
    # structural condition. The run counts how many, so the 3s/3p rescue
    # numbers are read with that in mind.
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps,
                         foldseek=[("P_eggduf", "sp|P0A0|X", 1e-20, 400, 0.99,
                                    0.8, "aerolysin toxin", "swissprot")])
    assert df.loc["P_eggduf", "bin"] == "3d_duf_only"
    assert df.loc["P_eggduf", "rescued_by"] == "structure"
    assert "have structure or profile evidence" in capsys.readouterr().err


def test_jackhmmer_hits_are_keyed_on_the_query_not_the_target(ma, tmp_path,
                                                              paths_for):
    # symptom: parse_hmm_tblout keys on column 1, which is right for hmmsearch
    # and inverted for jackhmmer; the dict was full of UniRef ids and every
    # protein got "".
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps,
                         jackhmmer=[("UniRef50_A0A1", "P_dark1", "-"),
                                    ("P_dark2", "P_dark2", "-")])
    assert df.loc["P_dark1", "jackhmmer_hit"] == "UniRef50_A0A1"
    assert df.loc["P_dark1", "bin"] == "3p_profile_only"
    assert df.loc["P_dark2", "jackhmmer_hit"] == "", "a self-hit is not a homolog"


def test_kofam_can_rescue_a_protein_from_the_ko_less_bins(ma, tmp_path,
                                                          paths_for, capsys):
    # symptom, and the reason the kofam stage is a CONTROL: if KOfam rescues a
    # large slice, the KO-less bins were partly an artefact of eggNOG's search.
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps, kofam=[
        ("P_dark1", "K01234", 200.0, 1e-50, "phosphoglucomutase")])
    assert bool(df.loc["P_dark1", "has_ko"])
    assert df.loc["P_dark1", "ko_source"] == "kofam"
    # the run's own KO->map pairs promote it to 1_ko_pathway, because K01234
    # is mapped in this eggNOG table
    assert df.loc["P_dark1", "bin"] == "1_ko_pathway"
    err = capsys.readouterr().err
    assert "leave the KO-less bins on KOfam evidence alone" in err


def test_disagreeing_ko_calls_are_recorded_per_protein(ma, tmp_path,
                                                       paths_for):
    # symptom: ko_source "both" said nothing about agreement, so a protein
    # eggNOG called K00001 and KOfam called K99999 looked confirmed.
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps, kofam=[
        ("P_ko_path", "K99999", 200.0, 1e-50, "something else"),
        ("P_ko_orphan", "K03456", 200.0, 1e-50, "orphan KO")])
    assert bool(df.loc["P_ko_path", "ko_conflict"])
    assert not bool(df.loc["P_ko_orphan", "ko_conflict"])
    assert df.loc["P_ko_orphan", "ko_source"] == "both"


# --- the dark work-list ----------------------------------------------
def test_the_structure_worklist_and_the_profile_query_set_are_separate_files(
        ma, tmp_path, paths_for, capsys):
    # symptom: a GPU budget (max_dark_structures) silently truncated the
    # unrelated profile searches, because both read dark.faa.
    ps = F.protein_set(n_extra=0)
    cfg, p = paths_for("dark")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    cfg["max_dark_structures"] = 1
    F.write_emapper(p.emapper, ps)
    ma.build_annotation(cfg, p, emit_dark=p.dark, emit_dark_all=p.dark_all)
    n_dark = sum(1 for _ in ma.read_fasta(p.dark))
    n_all = sum(1 for _ in ma.read_fasta(p.dark_all))
    assert n_dark == 1
    assert n_all > n_dark
    assert "hhblits and jackhmmer query" in capsys.readouterr().err


def test_decoys_and_contaminants_are_kept_off_the_dark_worklist(ma, tmp_path,
                                                                paths_for,
                                                                capsys):
    # symptom: entrapment and contaminant sequences are unannotated by
    # construction, score highly (no_ko + small_protein) and floated to the
    # top of the fold budget, then reappeared in the report as novel dark
    # proteins.
    ps = F.protein_set() + [
        F.Protein("rev_P_decoy", "M" + "A" * 60, in_emapper=False),
        F.Protein("HUMANHOST_P1", "M" + "A" * 60, in_emapper=False)]
    cfg, p = paths_for("decoy")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    ma.build_annotation(cfg, p, emit_dark=p.dark, emit_dark_all=p.dark_all)
    ids = {pid for pid, _ in ma.read_fasta(p.dark)}
    assert not any(i.startswith(("rev_", "HUMANHOST_")) for i in ids)
    assert "excluded from" in capsys.readouterr().err


def test_a_requested_structure_that_was_never_folded_is_distinguishable(
        ma, tmp_path, paths_for, capsys):
    # symptom: reporting the request as "attempted" made "we folded it and
    # found nothing" indistinguishable from "we never folded it".
    ps = F.protein_set()
    cfg, p = paths_for("struct")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    ma.build_annotation(cfg, p, emit_dark=p.dark, emit_dark_all=p.dark_all)
    df = ma.build_annotation(cfg, p)
    assert bool(df["structure_requested"].any())
    assert not bool(df["structure_attempted"].any())
    assert "requested structures exist" in capsys.readouterr().err


# --- the effector score ----------------------------------------------
def test_a_weak_diamond_hit_scores_half_the_weight_of_a_strong_one(
        ma, tmp_path, paths_for):
    # symptom: a 31%-identity VFDB hit scored the same as a 90% one.
    ps = F.protein_set()
    cfg, p = paths_for("dia")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    F.write_diamond(os.path.join(p.diamond_dir, "vfdb.tsv"), [
        ("P_dark1", "VFG1", 95.0, 1e-40, 300.0, 90, "hemolysin"),
        ("P_dark2", "VFG2", 31.0, 1e-40, 300.0, 90, "hemolysin")])
    df = ma.build_annotation(cfg, p)
    w = cfg["diamond_weights"]["vfdb"]
    assert df.loc["P_dark1", "effector_score"] - \
        df.loc["P_dark2", "effector_score"] == w - w // 2


def test_a_bare_lpxtg_motif_alone_is_not_a_sortase_substrate(ma, tmp_path,
                                                             paths_for):
    # symptom: a bare LP.TG 4-mer occurs by chance in ~0.03% of proteins and
    # on its own used to set surface_or_secreted, which gates the shortlist.
    real = "M" + "A" * 100 + "LPKTG" + "AVILMFWCGPAVILMFWCGP" + "KRK"
    fake = "M" + "A" * 100 + "LPKTG" + "DEDEDEDEDEDEDEDEDEDE" + "DED"
    ps = [F.Protein("P_real", real, in_emapper=False),
          F.Protein("P_fake", fake, in_emapper=False)]
    cfg, p = paths_for("lpxtg")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    df = ma.build_annotation(cfg, p)
    assert bool(df.loc["P_real", "lpxtg"])
    assert not bool(df.loc["P_fake", "lpxtg"])


def test_a_toxin_pattern_only_matches_whole_words(ma, tmp_path, paths_for):
    # symptom: unanchored, "RTX" matched inside "MRTXase" and "Rhs" inside
    # longer names, each adding the largest single effector weight.
    ps = F.protein_set()
    _, _, df = _annotate(ma, tmp_path, paths_for, ps, foldseek=[
        ("P_dark1", "T1", 1e-20, 400, 0.99, 0.8, "MRTXase family protein",
         "pdb"),
        ("P_dark2", "T2", 1e-20, 400, 0.99, 0.8, "RTX toxin", "pdb")])
    assert not bool(df.loc["P_dark1", "toxin_fold"])
    assert bool(df.loc["P_dark2", "toxin_fold"])


def test_a_pilin_signal_peptide_is_scored(ma, tmp_path, paths_for):
    # symptom: SignalP 6 has six classes; PILIN was parsed and then ignored
    # everywhere, so type IV pilins scored below an OTHER-class protein.
    ps = F.protein_set()
    cfg, p = paths_for("pilin")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    F.write_signalp(p.signalp, {"P_dark1": "PILIN", "P_dark2": "OTHER"})
    df = ma.build_annotation(cfg, p)
    assert bool(df.loc["P_dark1", "surface_or_secreted"])
    assert not bool(df.loc["P_dark2", "surface_or_secreted"])
    assert df.loc["P_dark1", "effector_score"] > \
        df.loc["P_dark2", "effector_score"]


def test_an_integer_context_column_is_not_reported_as_a_detected_flag(
        ma, tmp_path, paths_for):
    # symptom: `v == 1` was also true for the integer n_cazymes_in_window, so
    # the literal string "n_cazymes_in_window" appeared in context_flags for
    # every protein with exactly one CAZyme neighbour.
    import pandas as pd
    ps = F.protein_set()
    cfg, p = paths_for("ctx")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    pd.DataFrame([{"protein_id": "P_dark1", "n_cazymes_in_window": 1,
                   "prophage": True, "pul": False}]).to_csv(
        p.context, sep="\t", index=False)
    df = ma.build_annotation(cfg, p)
    flags = df.loc["P_dark1", "context_flags"]
    assert "n_cazymes_in_window" not in flags
    assert "prophage" in flags and "pul" not in flags


def test_the_shortlist_gate_survives_without_signalp_or_tmbed(ma, tmp_path,
                                                              paths_for):
    # symptom (a doc claim, not a code one): docs/signalp-6.md and README said
    # surface_or_secreted was "False for everything" without topology, so an
    # empty effector shortlist was "empty by construction". It is not — the
    # gate is an OR of four terms and two of them, the LPxTG motif and an
    # anchor domain, need no topology tool at all. On the real UC run 604
    # proteins passed it with topology off.
    anchored = "M" + "A" * 150
    sortase = ("M" + "A" * 100 + "LPKTG" + "AVILMFWCGPAVILMFWCGP" + "KRK")
    ps = [F.Protein("P_anchor", anchored, in_emapper=False),
          F.Protein("P_lpxtg", sortase, in_emapper=False),
          F.Protein("P_plain", "M" + "D" * 150, in_emapper=False)]
    cfg, p = paths_for("gate")
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    # an anchor Pfam, from the pfam stage — nothing to do with topology
    F.write_tblout(p.pfam, [("P_anchor", "Gram_pos_anchor",
                             cfg["anchor_pfams"][0])])
    assert not os.path.exists(p.signalp) and not os.path.exists(p.tmbed)

    df = ma.build_annotation(cfg, p)
    assert (df["sp_class"].fillna("") == "").all(), "no SignalP ran"
    assert (df["n_tmb"].fillna(0) == 0).all(), "no tmbed ran"
    assert bool(df.loc["P_anchor", "surface_or_secreted"])
    assert bool(df.loc["P_lpxtg", "surface_or_secreted"])
    assert not bool(df.loc["P_plain", "surface_or_secreted"])
    assert int(df["surface_or_secreted"].sum()) == 2
