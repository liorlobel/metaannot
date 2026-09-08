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
def _score_row(ma, desc, pident=90.0):
    cfg = ma.DEFAULT_CONFIG
    cw = cfg["vfdb_category_weights"]
    import re
    m = re.search(r"\((VFC\d+)\)", desc)
    return cw.get(m.group(1), cfg["diamond_weights"]["vfdb"]) if m else \
        cfg["diamond_weights"]["vfdb"]


def test_an_exported_effector_outweighs_a_housekeeping_virulence_hit(ma):
    """The defect this fixes: on a real gut metaproteome the two largest VFDB
    categories were immune modulation (965 hits) and nutritional/metabolic
    (903 - GroEL, ClpP, GuaA), each collecting the largest DIAMOND weight in
    the config, while effector delivery and exotoxin together were 11%."""
    eff = _score_row(ma, "Type III secretion effector SopB (VFC0086)")
    hk = _score_row(ma, "GroEL chaperonin (VFC0272)")
    assert eff > hk, "housekeeping in a virulence coat must not outrank an effector"


def test_a_hit_with_no_category_code_falls_back_to_the_flat_weight(ma):
    fallback = _score_row(ma, "some description with no code")
    assert fallback == ma.DEFAULT_CONFIG["diamond_weights"]["vfdb"]


def test_the_categories_are_keyed_on_the_numeric_code_not_the_prose(ma):
    """VFDB can reword a category name; it will not renumber the code."""
    for k in ma.DEFAULT_CONFIG["vfdb_category_weights"]:
        assert k.startswith("VFC") and k[3:].isdigit(), \
            f"{k} is not a VFC code, so a reworded category name would break it"


def test_clearing_the_category_map_restores_the_flat_weight(ma):
    """The escape hatch has to exist: an empty map is the old behaviour."""
    assert isinstance(ma.DEFAULT_CONFIG["vfdb_category_weights"], dict)
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


@pytest.mark.parametrize("desc", [
    "contact-dependent growth inhibition toxin CdiA",
    "LXG domain-containing toxin",
    "Ntox47 nuclease toxin domain",
    "zeta toxin family protein",
])
def test_the_families_gut_commensals_actually_carry_are_matched(ma, desc):
    assert _tox_re(ma).search(desc)


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
