"""The effectors stage is gone, and the things its removal exposed.

It ingested predictions from Bastion/EffectiveDB/T4SEpp - web services this
tool cannot invoke, and which are now unreachable - so it asked the user to go
and do work elsewhere. Removal was a provable no-op: the stage contributed
nothing unless configured, and the score never depended on it.
"""
import pandas as pd
import pytest


# ----------------------------------------------------------------------
# the stage is gone, and its absence is not an error
# ----------------------------------------------------------------------
def test_there_is_no_effectors_stage(ma):
    names = [s["name"] for s in ma.STAGES]
    assert "effectors" not in names
    assert not hasattr(ma, "stage_effectors")
    assert not hasattr(ma, "parse_external_predictions")


def test_no_stage_still_depends_on_the_effectors_output(ma):
    for s in ma.STAGES:
        assert "effectors" not in (s.get("deps") or []), \
            f"stage {s['name']} still lists effectors as a dependency"


def test_the_retired_config_keys_are_gone_from_the_defaults(ma):
    assert "effector_predictions" not in ma.DEFAULT_CONFIG
    assert "effector_prediction_weight" not in ma.DEFAULT_CONFIG
    assert "effectors" not in ma.DEFAULT_CONFIG["run"]


# ----------------------------------------------------------------------
# an old config gets an explanation, not a spelling suggestion
# ----------------------------------------------------------------------
def test_a_retired_key_is_named_as_removed_not_as_a_typo(ma, capsys):
    """Without this the user goes looking for their own mistake.

    `run` is a fixed key set, so `effectors: false` in a config written last
    month falls through to the unknown-key path and nearest_config_key helpfully
    suggests the closest surviving name.
    """
    ma.report_unknown_keys({"run": {"effectors": False}}, "config.yaml")
    err = capsys.readouterr().err
    assert "no longer a setting" in err
    assert "did you mean" not in err.lower(), \
        "there is no key to mean instead; a suggestion would send the user hunting"


def test_a_real_typo_still_gets_its_suggestion(ma, capsys):
    ma.report_unknown_keys({"run": {"unipep": True}}, "config.yaml")
    err = capsys.readouterr().err
    assert "unrecognised key" in err
    assert "no longer a setting" not in err


def test_every_retired_key_names_where_to_read_about_it(ma):
    for key, msg in ma.RETIRED_KEYS.items():
        assert "CHANGELOG" in msg, f"{key} does not say where the reason is"


# ----------------------------------------------------------------------
# VFDB by category, not by a flat weight
# ----------------------------------------------------------------------
def _vfdb_run(ma, tmp_path, paths_for, name, hits, **over):
    """build_annotation over a project whose only extra evidence is a VFDB
    DIAMOND table, so what is asserted is the number the real scoring path
    produced and not a second implementation of it.

    An earlier version of these tests scored the row itself, in six lines
    beside the assertions, and so covered nothing: build_annotation could stop
    reading vfdb_category_weights altogether and they would all still pass.
    That is how the missing `vfdb_category_weights` entry in the integrate
    stage's key list survived - re-weighting VFDB was a silent no-op on a
    re-run and nothing here could see it.

    `hits` is {protein_id: stitle}. Every hit is written at 95% identity, well
    above diamond_strong_pident, so the half-weight rule for a weak hit is not
    in play and a difference between two scores is a difference between two
    category weights. P_dark1, P_dark2 and P_dark3 all score the same without
    a DIAMOND hit, which is what makes that difference readable.
    """
    import os
    import fixtures as F
    ps = F.protein_set()
    cfg, p = paths_for(name)
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / f"{name}.faa"), ps)
    cfg.update(over)
    F.write_emapper(p.emapper, ps)
    F.write_diamond(os.path.join(p.diamond_dir, "vfdb.tsv"),
                    [(pid, f"VFG{i:04d}", 95.0, 1e-40, 300.0, 90, title)
                     for i, (pid, title) in enumerate(hits.items())])
    return ma.build_annotation(cfg, p)


def test_an_exported_effector_outweighs_a_housekeeping_virulence_hit(
        ma, tmp_path, paths_for):
    """The defect this fixes: on a real gut metaproteome the two largest VFDB
    categories were immune modulation (965 hits) and nutritional/metabolic
    (903 - GroEL, ClpP, GuaA), each collecting the largest DIAMOND weight in
    the config, while effector delivery and exotoxin together were 11%."""
    cw = ma.DEFAULT_CONFIG["vfdb_category_weights"]
    df = _vfdb_run(ma, tmp_path, paths_for, "vfcat", {
        "P_dark1": "(gb|AAL20571) (sopB) inositol phosphate phosphatase SopB "
                   "(VFC0086) [Effector delivery system (VF0367)]",
        "P_dark2": "(gb|AAA57975) (groEL) chaperonin GroEL (VFC0272) "
                   "[Nutritional/Metabolic factor (VF0483)]"})
    eff = int(df.loc["P_dark1", "export_score"])
    hk = int(df.loc["P_dark2", "export_score"])
    assert eff > hk, "housekeeping in a virulence coat must not outrank an effector"
    assert eff - hk == cw["VFC0086"] - cw["VFC0272"]


def test_a_hit_with_no_category_code_falls_back_to_the_flat_weight(
        ma, tmp_path, paths_for):
    """Both halves of the fallback: a title VFDB wrote without a code at all,
    and a code this config has no opinion about. Neither may score zero - the
    hit is still evidence - and neither may score a category weight."""
    flat = ma.DEFAULT_CONFIG["diamond_weights"]["vfdb"]
    df = _vfdb_run(ma, tmp_path, paths_for, "vfflat", {
        "P_dark1": "(gb|AAB12345) (hlyA) haemolysin, category not stated",
        "P_dark2": "(gb|AAB67890) (xyzA) a category this config does not "
                   "list (VFC9999) [Something new]"})
    base = int(df.loc["P_dark3", "export_score"])      # no VFDB hit at all
    assert int(df.loc["P_dark1", "export_score"]) - base == flat
    assert int(df.loc["P_dark2", "export_score"]) - base == flat


def test_the_vfdb_category_parse_rate_is_logged_once(ma, tmp_path, paths_for,
                                                     capsys):
    """Said once per run, because the failure it is there to expose is silent:
    if VFDB rewords its titles so the code no longer sits in parentheses,
    every hit takes the fallback weight and the only visible sign is this
    line reporting 0/N."""
    _vfdb_run(ma, tmp_path, paths_for, "vflog", {
        "P_dark1": "(gb|AAL20571) (sopB) effector SopB (VFC0086) [T3SS]",
        "P_dark2": "(gb|AAA57975) (groEL) chaperonin GroEL (VFC0272) [NMF]",
        "P_dark3": "(gb|AAB12345) (hlyA) haemolysin, category not stated"})
    flat = ma.DEFAULT_CONFIG["diamond_weights"]["vfdb"]
    err = capsys.readouterr().err
    assert "vfdb: 2/3 hit(s) carry a VFC category code" in err
    assert f"diamond_weights.vfdb={flat}" in err
    assert err.count("carry a VFC category code") == 1, \
        "one line per run, not one per hit"


def test_the_categories_are_keyed_on_the_numeric_code_not_the_prose(ma):
    """VFDB can reword a category name; it will not renumber the code."""
    for k in ma.DEFAULT_CONFIG["vfdb_category_weights"]:
        assert k.startswith("VFC") and k[3:].isdigit(), \
            f"{k} is not a VFC code, so a reworded category name would break it"


def test_clearing_the_category_map_restores_the_flat_weight(ma, tmp_path,
                                                            paths_for):
    """The escape hatch has to exist: an empty map is the old behaviour."""
    flat = ma.DEFAULT_CONFIG["diamond_weights"]["vfdb"]
    df = _vfdb_run(ma, tmp_path, paths_for, "vfoff", {
        "P_dark1": "(gb|AAA57975) (groEL) chaperonin GroEL (VFC0272) [NMF]"},
        vfdb_category_weights={})
    base = int(df.loc["P_dark3", "export_score"])      # no VFDB hit at all
    assert int(df.loc["P_dark1", "export_score"]) - base == flat, \
        "with the map cleared a VFC0272 hit must score the flat weight again"
    assert "vfdb_category_weights" in ma.FREEFORM, \
        "an unlisted VFDB category must not read as a typo"


# ----------------------------------------------------------------------
# toxin_fold, which carried weight 3 and never fired
# ----------------------------------------------------------------------
def _tox_re(ma):
    import re
    pats = ma.DEFAULT_CONFIG["toxin_fold_patterns"]
    return re.compile(r"\b(?:" + "|".join(pats) + r")\b", re.I)


def test_a_holotoxin_is_matched(ma):
    """The real miss: PDB 2vse is a Tc-family holotoxin, and the whole-word
    rule means the shipped 'Tc toxin' pattern cannot match inside 'holotoxin'.
    On 38,204 real proteins the old list fired zero times."""
    desc = "2vse-assembly1_A Structure and mode of action of a mosquitocidal holotoxin"
    assert _tox_re(ma).search(desc)


# Each case names the ONE pattern it exists to exercise, and the test below
# removes that pattern to prove the match came from it. The first version of
# this test asserted only "something matched" on
# "Ntox47 nuclease toxin domain" — which matches on `nuclease toxin`, so it
# passed while `Ntox` matched nothing at all. A description that a neighbouring
# pattern also covers cannot tell you whether the pattern you meant to add
# works.
@pytest.mark.parametrize("pattern,desc", [
    ("CdiA", "CdiA-CT toxin domain"),
    ("LXG", "LXG domain-containing protein"),
    (r"Ntox\d*", "Ntox47 domain-containing protein"),
    (r"Ntox\d*", "Bacterial toxin 28 domain (Ntox28)"),
    ("zeta toxin", "zeta toxin family protein"),
    ("nuclease toxin", "HNH nuclease toxin"),
    ("contact-dependent", "contact-dependent growth inhibition system"),
])
def test_the_families_gut_commensals_actually_carry_are_matched(ma, pattern,
                                                                desc):
    """Every real Ntox family is Ntox followed by a number, and a digit is a
    word character — so the trailing \\b of the whole-word rule meant a bare
    `Ntox` matched none of Ntox15/Ntox28/Ntox47. Same failure as 'Tc toxin'
    inside 'holotoxin', in a pattern added to fix that one."""
    import re
    pats = ma.DEFAULT_CONFIG["toxin_fold_patterns"]
    assert pattern in pats, f"{pattern!r} is no longer in toxin_fold_patterns"
    assert _tox_re(ma).search(desc), f"{desc!r} scores no toxin fold"

    # The whole point: without THIS pattern the description must stop matching.
    # Otherwise the case is being carried by one of its neighbours and says
    # nothing about the pattern it was written for.
    without = re.compile(
        r"\b(?:" + "|".join(p for p in pats if p != pattern) + r")\b", re.I)
    assert not without.search(desc), (
        f"{desc!r} still matches with {pattern!r} removed, so this case does "
        f"not actually exercise {pattern!r}")


@pytest.mark.parametrize("desc", [
    "Cytidine deaminase",
    "Patatin-like phospholipase",
    "Lysine--tRNA ligase",
    "Hemolysin III family protein",
])
def test_housekeeping_lookalikes_are_still_not_matched(ma, desc):
    """The word-boundary comment in the config exists because bare
    'deaminase' and 'hemolysin' already burned this once."""
    assert not _tox_re(ma).search(desc), f"{desc} must not score a toxin fold"


# ----------------------------------------------------------------------
# the rename
# ----------------------------------------------------------------------
def test_export_score_and_effector_score_are_the_same_numbers(ma, tmp_path,
                                                              paths_for):
    """The old name is carried one release so existing scripts keep working.
    It is an alias, not a second opinion."""
    import fixtures as F
    cfg, p = paths_for()
    ps = F.protein_set()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    df = ma.build_annotation(cfg, p)
    assert "export_score" in df.columns
    assert "effector_score" in df.columns
    assert (df["export_score"] == df["effector_score"]).all()


def test_the_bin_summary_reports_the_new_name(ma, tmp_path, paths_for):
    import fixtures as F
    cfg, p = paths_for()
    ps = F.protein_set()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    df = ma.build_annotation(cfg, p)
    ma.write_summary(df, p.summary)
    cols = pd.read_csv(p.summary, sep="\t").columns
    assert "median_export_score" in cols


# ----------------------------------------------------------------------
# upgrading a results directory written before the rename
# ----------------------------------------------------------------------
def _v02_results(ma, tmp_path, paths_for):
    """A finished results directory whose pass1 is v0.2-shaped.

    Built by running this version and then removing from annotation_pass1.tsv
    exactly the two columns v0.2 did not have: export_score, which had not yet
    been split out of effector_score, and ncbifam_accs.
    """
    import fixtures as F
    cfg, p = paths_for()
    ps = F.protein_set()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    ma.stage_integrate_pass1(cfg, p)
    ma.stage_integrate_final(cfg, p)
    old = pd.read_csv(p.pass1, sep="\t", dtype=str).drop(
        columns=["export_score", "ncbifam_accs"])
    old.to_csv(p.pass1, sep="\t", index=False)
    return cfg, p


def test_a_pass1_from_before_the_rename_is_refused_by_name(ma, tmp_path,
                                                           paths_for):
    # symptom: `run --only finalise` on a v0.2 results directory died with
    # KeyError('export_score') out of write_summary, naming neither the file
    # it could not use nor the stage that would rebuild it.
    cfg, p = _v02_results(ma, tmp_path, paths_for)
    with pytest.raises(ma.StageError) as e:
        ma.stage_integrate_final(cfg, p)
    msg = str(e.value)
    assert "export_score" in msg, "the missing column has to be named"
    assert p.pass1 in msg, "the file that cannot be reused has to be named"
    assert "integrate" in msg, "the stage that rebuilds it has to be named"


def test_a_refused_pass1_leaves_the_previous_final_table_alone(ma, tmp_path,
                                                               paths_for):
    # symptom: the reuse branch wrote annotation_final.tsv BEFORE
    # write_summary raised, so the upgrade half-applied - the final table came
    # out v0.2-shaped with no export_score at all, over a good one.
    cfg, p = _v02_results(ma, tmp_path, paths_for)
    before = open(p.final, encoding="utf-8").read()
    try:
        ma.stage_integrate_final(cfg, p)
    except Exception:
        pass
    assert open(p.final, encoding="utf-8").read() == before, \
        "a pass1 the stage goes on to refuse must not reach annotation_final.tsv"
    assert "export_score" in before.split("\n", 1)[0]


def test_ann_core_cols_lists_every_column_build_annotation_writes(ma, tmp_path,
                                                                  paths_for):
    # symptom: the check above is only as general as this list. v0.2 was
    # missing ncbifam_accs as well as export_score, and a list that drifts
    # from build_annotation lets the next such column through in silence.
    import fixtures as F
    cfg, p = paths_for()
    ps = F.protein_set()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"), ps)
    F.write_emapper(p.emapper, ps)
    df = ma.build_annotation(cfg, p)
    assert list(df.reset_index().columns) == list(ma.ANN_CORE_COLS)
