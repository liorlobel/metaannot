"""Documentation against the code.

Every number and every flag in README.md and TUTORIAL.md is an instruction
someone will follow on a shared server, so a doc that has drifted from the
code is a defect, not a typo.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys

import pytest

from conftest import METAANNOT_PY, ROOT

README = os.path.join(ROOT, "README.md")
TUTORIAL = os.path.join(ROOT, "TUTORIAL.md")
CHANGELOG = os.path.join(ROOT, "CHANGELOG.md")
CLAUDE = os.path.join(ROOT, "CLAUDE.md")
DOCS = [README, TUTORIAL, CLAUDE]


def _text(path):
    return open(path, encoding="utf-8").read()


def _norm(s):
    """Whitespace-collapsed, so a phrase that the markdown happens to wrap
    across two lines still matches."""
    return re.sub(r"\s+", " ", s)


def _code_lines(text):
    """Lines inside fenced code blocks."""
    out, inside = [], False
    for line in text.splitlines():
        if line.strip().startswith("```"):
            inside = not inside
            continue
        if inside:
            out.append(line)
    return out


@pytest.fixture(scope="module")
def cli_flags():
    """{subcommand: set of flags}, read from the real --help output."""
    top = subprocess.run([sys.executable, METAANNOT_PY, "--help"],
                         capture_output=True, text=True).stdout
    subs = set(re.findall(r"^\s{4}(\w+)\s", top, re.M))
    out = {"": set(re.findall(r"(--[\w-]+)", top))}
    for s in sorted(subs):
        h = subprocess.run([sys.executable, METAANNOT_PY, s, "--help"],
                           capture_output=True, text=True).stdout
        out[s] = set(re.findall(r"(--?[\w-]+)", h))
    return out


# --- finding 59 -------------------------------------------------------
def test_the_readme_stage_chain_matches_stage_names(ma):
    # symptom: the chain is what a reader uses to decide what --from means.
    txt = _text(README)
    block = re.search(r"```\n(emapper.*?)\n```", txt, re.S)
    assert block, "the stage chain block is gone from README.md"
    named = re.findall(r"[a-z_]+", block.group(1))
    assert named == list(ma.STAGE_NAMES)


def test_every_bin_is_documented(ma):
    txt = _text(README)
    for b in ma.BIN_ORDER:
        assert b in txt, f"{b} is not documented in README.md"
    documented = set(re.findall(r"`(\d[a-z]?_[a-z_]+)`", txt))
    unknown = documented - set(ma.BIN_ORDER)
    assert not unknown, f"README documents bins that do not exist: {unknown}"


def test_every_subcommand_is_documented(cli_flags):
    txt = _text(README) + _text(TUTORIAL)
    for sub in cli_flags:
        if not sub:
            continue
        assert re.search(rf"metaannot\.py {sub}\b", txt), \
            f"the '{sub}' subcommand appears in no doc"


def test_every_run_toggle_is_documented(ma):
    txt = _text(README) + _text(TUTORIAL)
    for flag in ma.DEFAULT_CONFIG["run"]:
        assert re.search(rf"\b{flag}\b", txt), \
            f"run.{flag} is documented nowhere"


def test_stages_without_a_run_flag_are_documented_as_always_running(ma):
    always = [s["name"] for s in ma.STAGES if s["enabled"] is None]
    assert sorted(always) == ["finalise", "integrate"]
    assert "`integrate` and `finalise` have no flag and always run" in \
        _norm(_text(README))


def test_the_documented_default_stage_set_matches_default_config(ma):
    # README: "A fresh init config turns on six of them: eggnog, pfam, dbcan,
    # diamond, cluster and join."
    on = {k for k, v in ma.DEFAULT_CONFIG["run"].items() if v}
    assert on == {"eggnog", "pfam", "dbcan", "diamond", "cluster", "join"}
    assert "turns on **six** of them" in _norm(_text(README))


@pytest.mark.parametrize("claim,key,want", [
    ("aborts below 50% overlap", "emapper_min_coverage", 0.50),
    ("warn coverage", "emapper_warn_coverage", 0.90),
])
def test_the_documented_coverage_thresholds_match_the_defaults(ma, claim, key,
                                                               want):
    assert ma.DEFAULT_CONFIG[key] == want
    txt = _text(README) + _text(TUTORIAL) + _text(CLAUDE)
    pct = f"{int(want * 100)}%"
    assert pct in txt, f"{claim}: {pct} is not stated anywhere"


def test_documented_numeric_defaults_match_default_config(ma):
    d = ma.DEFAULT_CONFIG
    checks = {
        "min_features_per_protein": 1,
        "taxon_min_proteins_for_factor": 4,
        "stage_workers": 4,
        "effector_prediction_weight": 3,
        "emapper_dbmem_min_gb": 64,
        "max_dark_structures": 2000,
        "max_len_structure": 700,
    }
    for k, v in checks.items():
        assert d[k] == v, f"{k} changed to {d[k]}; update the docs too"
    a = d["analysis"]
    assert a["fdr"] == 0.05 and a["min_valid_per_group"] == 3
    assert a["taxon_min_proteins"] == d["taxon_min_proteins_for_factor"], \
        "the README documents these two as deliberately equal"
    assert d["thresholds"]["smorf_max_len"] == 100
    assert "hard floor of 90 nt" in _text(README)


def test_every_peptide_assignment_mode_is_documented_somewhere(ma):
    txt = _text(README) + _text(TUTORIAL) + _text(CHANGELOG) + _text(CLAUDE)
    for mode in ma.ASSIGNMENT_MODES:
        assert mode in txt, f"{mode} is documented nowhere"
    assert set(ma.ASSIGNMENT_MODES) == {"protein_unique", "taxon_unique",
                                        "taxon_or_family_unique", "razor"}


@pytest.mark.xfail(reason="LIVE DOC GAP: README's 'The shared-peptide rule' "
                          "section lists protein_unique, taxon_unique and "
                          "razor but not taxon_or_family_unique, which is a "
                          "real accepted value of peptide_assignment and the "
                          "only one that can silently attribute a peptide to "
                          "one member of a SEQUENCE CLUSTER. It appears only "
                          "in CHANGELOG.md.",
                   strict=True)
def test_the_readme_shared_peptide_section_lists_every_mode(ma):
    txt = _text(README)
    for mode in ma.ASSIGNMENT_MODES:
        assert f"`{mode}`" in txt, f"{mode} is not in README.md"


def test_the_documented_quant_formats_are_the_real_ones(ma):
    txt = _text(README)
    documented = set(re.findall(r"`(diann|fragpipe|fragpipe_peptide|"
                                r"fragpipe_ion|msstats_csv|msstats_feature|"
                                r"msstats_protein)`", txt))
    assert documented == ma.ALL_FORMATS


def test_the_documented_free_form_blocks_are_the_real_ones(ma):
    txt = _text(README)
    for block in ma.FREEFORM:
        assert block in txt, f"{block} is not named as free-form in README.md"


def test_the_documented_default_diamond_weights_match(ma):
    # the README tells the user a database without a weight scores 0
    w = ma.DEFAULT_CONFIG["diamond_weights"]
    assert set(w) == set(ma.DEFAULT_CONFIG["db"]["diamond"])
    assert "scores 0 and the run warns" in _norm(_text(README))


def test_the_signature_version_is_documented_as_separate_from_the_tool_version(
        ma):
    assert "SIGNATURE_VERSION" in _text(README)
    assert ma.SIGNATURE_VERSION != ma.__version__


def test_the_changelog_names_the_released_version(ma):
    assert f"## v{ma.__version__}" in _text(CHANGELOG)


# --- finding 60 -------------------------------------------------------
def _documented_invocations(text):
    """(line, subcommand, [flags]) for every metaannot.py COMMAND in a doc.

    Only lines inside fenced code blocks count: prose mentions the file name
    in backticks, and `python metaannot.py --help`. inside a sentence carries
    the full stop into the token.
    """
    out = []
    for line in _code_lines(text):
        if "metaannot.py" not in line or line.strip().startswith("#"):
            continue
        tail = line.split("metaannot.py", 1)[1]
        # stop at a shell pipe or redirect: what follows belongs to another
        # program (tail, grep, tee).
        tail = re.split(r"[|>]|\d>&\d", tail)[0]
        words = [w.strip("`.,") for w in tail.split()]
        sub = words[0] if words and not words[0].startswith("-") else ""
        flags = [w for w in words if w.startswith("--")]
        out.append((line.strip(), sub, flags))
    return out


@pytest.mark.parametrize("doc", ["README.md", "TUTORIAL.md", "CLAUDE.md",
                                 "examples/server-run-plan/README.md"])
def test_every_documented_command_and_flag_exists_in_the_cli(cli_flags, doc):
    # symptom: the docs are a runbook. A flag that no longer exists sends the
    # reader into an argparse error six hours into a session.
    text = _text(os.path.join(ROOT, doc))
    bad = []
    for line, sub, flags in _documented_invocations(text):
        if sub and sub not in cli_flags:
            bad.append((line, f"unknown subcommand '{sub}'"))
            continue
        known = cli_flags.get(sub, set()) | cli_flags[""]
        for f in flags:
            name = f.split("=")[0]
            if name not in known:
                bad.append((line, f"unknown flag '{name}' for '{sub or 'top'}'"))
    assert bad == [], "\n".join(f"{w}: {l}" for l, w in bad)


def test_the_troubleshooting_table_quotes_messages_the_code_can_emit():
    # a symptom column that no longer matches any string in the tool is a
    # runbook entry nobody can find.
    code = _text(METAANNOT_PY)
    quoted = [
        "another metaannot is already running here",
        "is an HTML page, not the database",
        "not adopting",
        "parsed 0",
        "unrecognised key",
        "this looks like isobaric (TMT/iTRAQ) output",
        "cannot run:",
        "deadlock:",
        # split across two source lines in the die() call, so only the first
        # half is a contiguous literal
        "have rows with different ",
        "matches no fasta id",
    ]
    tut = _norm(_text(TUTORIAL))
    for q in quoted:
        assert _norm(q) in tut, f"the troubleshooting table lost '{q}'"
        assert q in code, f"TUTORIAL quotes '{q}', which the code cannot emit"


def test_the_readme_does_not_claim_untested_things_are_verified():
    # the honesty of the "what has actually been run" section is load-bearing:
    # it is what stops a reader treating a search-stage run as production.
    txt = _text(README)
    assert "Not on real data" in txt
    assert "Report and R object: synthetic data only" in txt


def test_the_docs_and_the_tool_agree_on_the_signalp_command_line():
    # docs/signalp-6.md tells the user which SignalP version to install.
    doc = _text(os.path.join(ROOT, "docs", "signalp-6.md"))
    code = _text(METAANNOT_PY)
    for flag in ("--fastafile", "--organism", "--output_dir", "--format",
                 "--mode"):
        assert flag in doc and flag in code
    assert "prediction_results.txt" in doc and "prediction_results.txt" in code


def test_the_tutorial_bin_reference_numbers_are_internally_consistent():
    # the tutorial quotes a real run's bin table; the counts must add up to
    # the protein total it also quotes.
    txt = _text(TUTORIAL)
    m = re.search(r"^1_ko_pathway\s+(\d+).*?^4_dark\s+(\d+)", txt,
                  re.S | re.M)
    assert m, "the reference bin table is gone from TUTORIAL.md"
    counts = [int(x) for x in re.findall(r"^\S+\s+(\d+)\s+\d+\.\d%", txt, re.M)]
    assert sum(counts) == 38204, \
        f"the quoted bin counts sum to {sum(counts)}, not the 38,204 proteins"


@pytest.mark.xfail(reason="LIVE DOC DRIFT: README and TUTORIAL still describe "
                          "a duplicate top-level `run:` or `db:` key as being "
                          "SILENTLY collapsed by yaml.safe_load, with the "
                          "first block discarded. load_config now installs a "
                          "no-duplicate loader and exits naming the key and "
                          "both line numbers, so the documented failure mode "
                          "can no longer happen and the advice reads as if it "
                          "still can.",
                   strict=True)
def test_the_docs_describe_duplicate_keys_as_refused_not_silently_collapsed():
    txt = _norm(_text(README) + " " + _text(TUTORIAL))
    assert "keeps only the last" not in txt
    assert "duplicate key" in txt


def test_a_duplicate_key_really_is_refused():
    # whatever the docs say, this is the behaviour to protect.
    src = _text(METAANNOT_PY)
    assert "duplicate key" in src
    assert "_NoDupLoader" in src


def test_the_docs_do_not_claim_the_shortlist_is_empty_by_construction():
    # symptom: docs/signalp-6.md and README said surface_or_secreted is "False
    # for everything" without SignalP/tmbed, so an empty effector shortlist was
    # a missing-tool artefact. On the first real run 604 proteins passed the
    # gate with topology off, and the shortlist was empty because nothing was
    # significant — a different conclusion entirely.
    signalp_doc = _text(os.path.join(ROOT, "docs", "signalp-6.md"))
    txt = _norm(_text(README) + " " + signalp_doc)
    assert "empty by construction" not in txt
    assert "surface_or_secreted is False for everything" not in txt
    # and the two terms that keep working must be named where the claim was
    assert "LPxTG" in signalp_doc and "anchor domain" in signalp_doc


def test_the_changelog_does_not_assert_the_claim_either():
    # the same wrong sentence survived in CHANGELOG's v0.1.0 entry, which the
    # check above did not cover. A released entry is a historical record, so it
    # carries a marked correction rather than a silent rewrite — which is why
    # this looks for the ASSERTING form only: v0.2.0's entry and the correction
    # both quote the phrase in order to disown it.
    txt = _norm(_text(CHANGELOG))
    assert "is empty by construction" not in txt
    assert "shortlist is empty" not in txt
    # the correction must still be visible rather than the claim just deleted
    assert "Corrected after v0.2.0" in txt


def test_the_documented_gate_terms_match_the_code():
    # the doc describes surface_or_secreted as an OR of four things; if the
    # code's definition changes, the description has to change with it.
    code = _text(METAANNOT_PY)
    start = code.find('df["surface_or_secreted"] = (')
    assert start != -1, "surface_or_secreted is no longer assigned in one place"
    # to the matching paren, not to the first one: the expression spans lines
    # and contains nested calls
    i = code.index("(", start)
    depth, j = 0, i
    while j < len(code):
        if code[j] == "(":
            depth += 1
        elif code[j] == ")":
            depth -= 1
            if depth == 0:
                break
        j += 1
    expr = code[i + 1:j]
    for term in ("sp_class", "lpxtg", "anchor_domain", "n_tmb"):
        assert term in expr, f"{term} is no longer part of the gate"


# --- the worked example under examples/ --------------------------------
EXAMPLE = os.path.join(ROOT, "examples", "server-run-plan")


def test_the_example_runbook_scripts_are_valid_shell():
    # they are meant to be copied to a server and run; a syntax error there is
    # found at hour six, not at review time.
    import glob
    scripts = sorted(glob.glob(os.path.join(EXAMPLE, "*.sh")))
    assert scripts, "the example lost its scripts"
    for s in scripts:
        r = subprocess.run(["bash", "-n", s], capture_output=True, text=True)
        assert r.returncode == 0, f"{os.path.basename(s)}: {r.stderr}"


def test_the_example_configs_have_no_unrecognised_keys(ma):
    # a typo in a shipped example config is a typo every reader copies. The
    # same check `doctor` runs, applied at review time.
    import glob
    import yaml as _yaml
    cfgs = sorted(glob.glob(os.path.join(EXAMPLE, "*", "config.yaml")))
    assert len(cfgs) >= 8, f"expected the eight example configs, found {len(cfgs)}"
    for c in cfgs:
        body = _yaml.safe_load(open(c, encoding="utf-8"))
        bad = ma.unknown_keys(body, ma.DEFAULT_CONFIG)
        assert bad == [], f"{os.path.relpath(c, ROOT)}: {bad}"


def test_the_example_configs_enable_the_stages_the_readme_claims(ma):
    # the README says nine stages are on and names why each of the others is
    # off. If a config drifts from that, the prose is wrong.
    import glob
    import yaml as _yaml
    want = {"eggnog", "pfam", "dbcan", "diamond", "cluster", "ncbifam",
            "kofam", "interpro", "join"}
    for c in sorted(glob.glob(os.path.join(EXAMPLE, "*", "config.yaml"))):
        run = _yaml.safe_load(open(c, encoding="utf-8"))["run"]
        on = {k for k, v in run.items() if v}
        assert on == want, f"{os.path.relpath(c, ROOT)} enables {sorted(on)}"
        # every stage that exists must be stated one way or the other
        missing = set(ma.DEFAULT_CONFIG["run"]) - set(run)
        assert missing == set(), f"{os.path.relpath(c, ROOT)} is silent about {missing}"
