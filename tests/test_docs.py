"""Documentation against the code.

Every number and every flag in README.md and TUTORIAL.md is an instruction
someone will follow on a shared server, so a doc that has drifted from the
code is a defect, not a typo.
"""
from __future__ import annotations

import ast
import io
import json
import os
import re
import subprocess
import sys
import tokenize

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


def _code_only(path):
    """`path`'s source with every comment and string literal removed.

    For asserting that a property holds of the CODE when the same words are
    also, legitimately, all over the prose. This file's job is to keep the docs
    and the code in step, and the docs live in the code too: a rule this
    codebase writes down at the point it is enforced will be written down in a
    comment or a docstring, in the same words the README uses. A raw `in src`
    scan cannot tell "the console binds 127.0.0.1" from "a 127.0.0.1 listener
    was rejected because...", and it is the second one that this project keeps
    writing.

    Tokenized rather than regexed, so a `#` inside a string and a quote inside
    a comment both land on the right side of the line. FSTRING_* is named
    because 3.12 stopped emitting f-strings as one STRING token; on 3.9, where
    this suite runs, that branch is simply never taken.
    """
    out = []
    with io.open(path, encoding="utf-8") as fh:
        for tok in tokenize.generate_tokens(fh.readline):
            name = tokenize.tok_name.get(tok.type, "")
            if name in ("COMMENT", "STRING") or name.startswith("FSTRING"):
                continue
            out.append(tok.string)
    return "\n".join(out)


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
    # README: "A fresh init config turns on six run: flags: eggnog (which runs
    # the emapper stage), pfam, dbcan, diamond, cluster and join."
    on = {k for k, v in ma.DEFAULT_CONFIG["run"].items() if v}
    assert on == {"eggnog", "pfam", "dbcan", "diamond", "cluster", "join"}
    txt = _norm(_text(README))
    assert "turns on **six**" in txt
    # The count and the list have to agree, or the prose drifts from the config
    # one name at a time.
    for name in on:
        assert f"`{name}`" in txt, f"{name} is on by default but unlisted"


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


def test_the_tmbed_length_cap_is_documented(ma):
    # symptom: tmbed_max_len drops every protein over 3,000 residues from the
    # topology stage — no n_tmh, no n_tmb, a changed stage output and a forced
    # re-run of it — and appeared in no document, so the first anybody heard
    # of the exclusion was a WARN in the middle of a multi-hour run.
    txt = _norm(_text(README))
    assert "tmbed_max_len" in txt, "tmbed_max_len is documented nowhere"
    assert ma.DEFAULT_CONFIG["tmbed_max_len"] == 3000
    assert f"tmbed_max_len: {ma.DEFAULT_CONFIG['tmbed_max_len']}" in txt, \
        "the documented default has drifted from DEFAULT_CONFIG"
    assert "tmbed_excluded.tsv" in txt, \
        "a reader has to be told where the proteins that were cut are listed"
    tmbed = next(s for s in ma.STAGES if s["name"] == "tmbed")
    assert "tmbed_max_len" in tmbed["keys"], \
        "the README says changing the cap re-runs the stage"


def test_every_peptide_assignment_mode_is_documented_somewhere(ma):
    txt = _text(README) + _text(TUTORIAL) + _text(CHANGELOG) + _text(CLAUDE)
    for mode in ma.ASSIGNMENT_MODES:
        assert mode in txt, f"{mode} is documented nowhere"
    assert set(ma.ASSIGNMENT_MODES) == {"protein_unique", "taxon_unique",
                                        "taxon_or_family_unique", "razor"}


def test_the_readme_shared_peptide_section_lists_every_mode(ma):
    txt = _text(README)
    for mode in ma.ASSIGNMENT_MODES:
        assert f"`{mode}`" in txt, f"{mode} is not in README.md"


def test_the_documented_quant_formats_are_the_real_ones(ma):
    txt = _text(README)
    documented = set(re.findall(r"`(diann|fragpipe|fragpipe_peptide|"
                                r"fragpipe_ion|fragpipe_tmt|msstats_csv|"
                                r"msstats_feature|msstats_protein)`", txt))
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
    """The NEWEST section is this tree's version — not merely a section
    somewhere in the file.

    `f"## v{__version__}" in text` was satisfied by ANY section, so a CHANGELOG
    carrying a written-up `## v0.5.0` above a `metaannot.py` still saying
    `__version__ = "0.4.0"` passed it: the old section matched the substring and
    the new one went unread. That is exactly the state this guard is for. It is
    the state a release is cut FROM — notes written, version not bumped — and
    the tag then ships a tool that reports itself as the previous release, in
    `--version`, in `describe --json`, and in the `_run` record every results
    directory keeps for ever.

    So the FIRST heading is the one that has to match, and the rest have to
    descend from it, which is the same defect one line further down: a section
    filed in the wrong place is a section nobody will find. Ordering is
    compared on the parsed numbers, because `v0.10.0` sorts before `v0.9.0` as
    a string and this project will get there.
    """
    heads = re.findall(r"^## v(\d+)\.(\d+)\.(\d+)\b", _text(CHANGELOG), re.M)
    assert heads, "CHANGELOG.md has no version sections at all"
    newest = ".".join(heads[0])
    assert newest == ma.__version__, (
        f"CHANGELOG.md's newest section is v{newest} while metaannot.py says "
        f"__version__ = {ma.__version__!r}. One of the two has not been "
        "released: bump __version__ to the version the notes describe, or "
        "move the section.")
    nums = [tuple(int(n) for n in h) for h in heads]
    assert nums == sorted(nums, reverse=True), \
        f"CHANGELOG.md's sections are not newest-first: {nums}"
    assert len(set(nums)) == len(nums), \
        f"CHANGELOG.md sections a version twice: {nums}"


# --- documentation that had drifted from the code ---------------------
def test_the_changelog_and_the_readme_disagreed_on_the_uc_stage_count(ma):
    # symptom: CHANGELOG's v0.3.0 entry said sixteen applicable stages ran on
    # the UC dataset while README said fifteen and enumerated fifteen, so the
    # answer to "how much of this tool has been run for real" depended on
    # which file you opened.
    readme = _norm(_text(README))
    m = re.search(r"(\w+) of the twenty-one stages ran on that input: (.+?)\.",
                  readme)
    assert m, "the 'what has actually been run' stage list is gone from README"
    named = set(re.findall(r"`([a-z_]+)`", m.group(2)))
    unknown = named - set(ma.STAGE_NAMES)
    assert not unknown, f"README names stages that do not exist: {unknown}"
    # the enumeration is what settles it: the count in the prose, the count of
    # the list, and the stages left over must all be the same arithmetic.
    assert m.group(1) == "Fifteen" and len(named) == 15
    assert len(set(ma.STAGE_NAMES) - named) == 6
    assert "all fifteen applicable stages" in _norm(_text(CHANGELOG))


def test_the_readme_eggnog_pfam_paragraph_was_truncated_mid_sentence(ma):
    # symptom: the paragraph broke off at "with `run." in v0.3, taking with it
    # the one thing it exists to say — that it is the pfam stage being OFF
    # that used to hide a domain eggNOG had already assigned.
    txt = _norm(_text(README))
    m = re.search(r"\*\*eggNOG's own `PFAMs` column counts as domain "
                  r"evidence\.\*\*(.+?)Bins 2 to 4", txt)
    assert m, "the eggNOG PFAMs paragraph is gone from README.md"
    para = m.group(1)
    assert para.count("`") % 2 == 0, \
        "an unclosed backtick: the paragraph is truncated again"
    assert "`run.pfam: false`" in para, \
        "the paragraph no longer names the condition it is about"
    # the figure is the code comment's, so the two cannot drift apart
    assert "6,271 of 17,377 dark proteins (36%)" in para
    assert "6,271 of 17,377 dark proteins (36%)" in _norm(_text(METAANNOT_PY))


def test_the_tutorial_free_form_list_named_a_removed_config_block(ma):
    # symptom: TUTORIAL told the reader to proof-read `effector_predictions`
    # by hand — a block v0.3.0 removed — and omitted two blocks that really
    # are free-form, so a typo in one of those went on being invisible.
    m = re.search(r"cannot see inside are the free-form blocks \((.+?)\),",
                  _norm(_text(TUTORIAL)))
    assert m, "the free-form block list is gone from TUTORIAL.md"
    assert set(re.findall(r"`([\w.]+)`", m.group(1))) == ma.FREEFORM


def test_the_tutorial_quoted_an_esmfold_message_v0_3_cannot_emit():
    # symptom: the troubleshooting table quoted "were skipped (OOM) or never
    # folded", the one hedged sentence v0.3 split into three separate causes.
    # The row sent a reader hunting for a string the tool no longer prints,
    # and blamed OOM for proteins that were never submitted to the card.
    tut, code = _norm(_text(TUTORIAL)), _text(METAANNOT_PY)
    assert "were skipped (OOM) or never folded" not in tut
    for q in ("esmfold has not finished", "although esmfold has finished",
              "were never submitted", "attempted and failed twice",
              "absent from that list"):
        assert _norm(q) in tut, f"the esmfold rows lost '{q}'"
        assert q in code, f"TUTORIAL quotes '{q}', which the code cannot emit"


def test_the_tutorial_still_said_reporter_ions_are_not_read():
    # symptom: phase 3a told the reader "reporter-ion quantification is not
    # supported at all" and pointed at a README section v0.3.0 deleted, and
    # the troubleshooting table read the isobaric refusal as a dead end rather
    # than as "the quant_format is wrong". Both predate the fragpipe_tmt
    # reader, and between them they talk a TMT user out of a supported route.
    tut, code = _norm(_text(TUTORIAL)), _text(METAANNOT_PY)
    for gone in ("reporter-ion quantification is not supported at all",
                 'See README, "FragPipe TMT output is not supported"',
                 "Reporter-ion channels are not read"):
        assert gone not in tut, f"TUTORIAL still says: {gone}"
    assert "`quant_format: fragpipe_tmt`" in tut
    # and the way out it offers is the one the refusal itself prints
    assert "quant_format: fragpipe_tmt" in code


def test_the_example_runbook_named_a_stage_that_no_longer_exists(ma):
    # symptom: the run plan's "what is off, and why" still listed `effectors`,
    # removed in v0.3.0 and absent from all eight configs, so the prose and
    # the files it describes disagreed about what the runs turn off.
    import glob
    import yaml as _yaml
    assert "effectors" not in ma.DEFAULT_CONFIG["run"]
    m = re.search(r"## What is off, and why\n+(.+?)\n\n",
                  _text(os.path.join(EXAMPLE, "README.md")), re.S)
    assert m, "the 'what is off' section is gone from the example README"
    named = set(re.findall(r"`([a-z_]+)`", m.group(1)))
    off = None
    for c in sorted(glob.glob(os.path.join(EXAMPLE, "*", "config.yaml"))):
        run = _yaml.safe_load(open(c, encoding="utf-8"))["run"]
        this = {k for k, v in run.items() if not v}
        assert off is None or this == off, f"{c} turns off a different set"
        off = this
    assert off, "no example config was read"
    assert named == off


def test_claude_md_told_the_next_session_to_refuse_tmt(ma):
    # symptom: CLAUDE.md is the standing instruction file an agent reads
    # first, and after v0.3.0 shipped the fragpipe_tmt reader it still said
    # "FragPipe TMT output is not supported ... do not read a number out of
    # one". As written it instructed the next session to refuse exactly what
    # the release was built to do.
    assert "fragpipe_tmt" in ma.ALL_FORMATS
    txt = _norm(_text(CLAUDE))
    for gone in ("FragPipe TMT output is not supported",
                 "Reporter-ion channels are not read",
                 "do not read a number out of one",
                 "Isobaric support is a separate piece of work",
                 "Isobaric input is **refused, not read**",
                 "Do not weaken it to get a TMT run through",
                 "Only the pre-search half has been run on real data",
                 "the report and R object have only ever run on synthetic "
                 "data"):
        assert gone not in txt, f"CLAUDE.md still says: {gone}"
    assert "quant_format: fragpipe_tmt" in txt, \
        "CLAUDE.md does not name the format that reads a TMT run"
    # and what replaces rule 8 has to be the real division
    assert "Fifteen of the twenty-one stages" in txt
    for never in ("smorf", "context", "hhblits", "jackhmmer", "unipept",
                  "taxonomy"):
        assert f"`{never}`" in txt, f"CLAUDE.md no longer names {never}"


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
    assert "Report and R object: built from the real run" in txt


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


def test_the_documented_stage_workers_default_matches_the_config(ma):
    """The README said "with the default 3" six lines above its own example
    writing 4, in the paragraph explaining the scheduling behaviour a reader
    goes there to understand. test_documented_numeric_defaults_match_default_config
    pins the config VALUE and never looked at the prose."""
    import re
    # the raw text, not _norm's: paragraph breaks are the boundary here.
    txt = _text(README)
    i = txt.find("Concurrency is capped by `stage_workers`")
    assert i != -1, "the stage_workers paragraph moved; find it again"
    para = _norm(txt[i:txt.index(chr(10) * 2, i)])
    m = re.search(r"with the default (\d+)", para)
    assert m, f"no stated default in the stage_workers paragraph: {para!r}"
    assert int(m.group(1)) == ma.DEFAULT_CONFIG["stage_workers"]


# ----------------------------------------------------------------------
# what a results directory says about itself
# ----------------------------------------------------------------------
def test_the_documented_run_record_matches_the_one_a_run_writes(ma, tmp_path):
    # this JSON block is what a console author reads to build a parser, so a
    # key it does not have - or a number of a different type - is a parser
    # written against a file that does not exist. `heartbeat_s: 30.0` in the
    # real record against `30` in the example was exactly that.
    m = re.search(r'```json\n(\s*"_run":.*?)\n```', _text(README), re.S)
    assert m, "the README no longer shows a _run example"
    example = json.loads("{" + m.group(1) + "}")["_run"]
    rec = ma.RunRecord(str(tmp_path / "state.json"), {}, ["metaannot", "run"],
                       None, ma.DEFAULT_CONFIG["heartbeat_s"]).rec
    assert set(example) == set(rec), "the documented record has drifted"
    assert example["heartbeat_s"] == ma.DEFAULT_CONFIG["heartbeat_s"]
    assert type(example["heartbeat_s"]) is type(rec["heartbeat_s"]), \
        "a reader parses this example and then meets the real file"
    assert example["final_status"] == "running"


def test_the_documented_heartbeat_default_is_the_real_one(ma):
    # the same pinning tmbed_max_len gets: a documented default that nothing
    # compares to DEFAULT_CONFIG drifts silently.
    txt = _norm(_text(README))
    assert f"heartbeat_s: {ma.DEFAULT_CONFIG['heartbeat_s']} " in txt, \
        "the documented heartbeat_s default has drifted from DEFAULT_CONFIG"


def test_the_docs_never_promise_that_a_heartbeat_reclaims_a_lock(ma):
    """The heartbeat is advisory, and the docs have to say so, because the
    opposite promise is what a shared-filesystem user would act on. A stopped
    heartbeat is not a stopped process - one failed write ends the timer - so
    nothing may reclaim on it, and _holder_is_alive is where that is enforced.
    """
    src = open(os.path.join(ROOT, "metaannot.py"), encoding="utf-8").read()
    assert "_heartbeat_says_dead" not in src, \
        "the reclaim-on-a-stale-heartbeat path is back"
    for doc in (README, TUTORIAL):
        txt = _norm(_text(doc))
        assert "advisory" in txt or "does not decide for you" in txt, \
            f"{os.path.basename(doc)} does not say the heartbeat is advisory"
        assert "is treated as a corpse and removed" not in txt
        assert "proof of death" not in txt
    assert "nothing in metaannot reclaims a lock because a" in \
        _norm(_text(README))


def test_the_documented_signal_exit_codes_are_the_conventional_ones(ma):
    # 128 + the signal, which is what the shell and systemd both expect: a
    # supervisor reading 130 for a `systemctl stop` is told a user pressed
    # Ctrl-C, and SuccessExitStatus=143 never matches.
    import signal
    txt = _norm(_text(README))
    assert f"`{128 + int(signal.SIGINT)}` for Ctrl-C" in txt
    assert f"`{128 + int(signal.SIGTERM)}` for `SIGTERM`" in txt
    assert "SuccessExitStatus=143" in txt
    assert f"exit status {128 + int(signal.SIGTERM)}" in _norm(_text(TUTORIAL))


def test_the_readme_does_not_sell_requirements_as_all_of_doctor(ma):
    # measured: on one project `describe --json` reported 3 not-ok
    # requirements while `doctor` additionally reported the missing input, the
    # emapper_precomputed files, a retired config key, the DIAMOND content
    # checks and the CUDA probe. A preflight built on `requirements` alone is
    # green for a config doctor fails, and the README said "exactly what
    # doctor checks, as data".
    txt = _norm(_text(README))
    assert "exactly what `doctor` checks" not in txt
    assert "the tool-and-database half of `doctor`" in txt
    cfg = ma.load_config(None)
    ids = {r["id"] for r in ma.requirements(cfg, ma.Paths(cfg))}
    for missing in ("proteins_faa", "quant_table"):
        assert missing not in ids, \
            f"requirements now covers {missing}; the README paragraph is stale"
        assert missing in txt, f"the README does not name {missing} as a gap"
    assert "CUDA probe" in txt


def test_the_readme_says_what_the_new_files_disclose(ma):
    # both artefacts widen what a shared results directory tells a reader, and
    # .metaannot_state.json travels inside the .rds people publish as
    # supplementary data.
    txt = _norm(_text(README))
    assert "`sources.*` and `tool_args` are free-form" in txt
    assert "presigned URL" in txt
    for freeform in ("sources", "tool_args"):
        assert freeform in ma.FREEFORM or freeform in ma.DEFAULT_CONFIG, \
            f"{freeform} is not a config block any more"


# `doctor`'s own section headings, minus the two that requirements() feeds
# (`== tools ==` and `== databases ==`, printed from a variable). The README
# paragraph that tells a console author what doctor checks BEYOND
# `requirements` has to name every one of them, because it reads as a closed
# list and a preflight screen gets built from it.
DOCTOR_SECTIONS = ("config", "inputs", "precomputed emapper", "gpu", "tmt",
                   "manifest", "taxonomy", "resources", "R")


def test_the_readme_names_every_section_doctor_prints(ma):
    # symptom: the list omitted `== manifest ==`, `== resources ==` and
    # `== R ==` outright, so a preflight built from it is green for a config
    # doctor fails on the manifest-to-column mapping, on a memory split that
    # cannot give eggNOG --dbmem, or on a missing required R package.
    src = _text(os.path.join(ROOT, "metaannot.py"))
    printed = set(re.findall(r"== ([A-Za-z][A-Za-z ]*) ==", src))
    assert printed, "doctor no longer prints section headings"
    assert printed == set(DOCTOR_SECTIONS), \
        f"doctor's sections have changed: {printed ^ set(DOCTOR_SECTIONS)}"
    txt = _norm(_text(README))
    for name in sorted(printed) + ["tools", "databases"]:
        assert f"`== {name} ==`" in txt, \
            f"the README's list of what doctor checks omits == {name} =="


def test_the_readme_says_the_run_records_config_path_can_be_null(ma, tmp_path):
    # symptom: `run` with no --config writes `"config_path": null`, but the
    # README's JSON example shows a string and the prose never said otherwise,
    # so a console author writing a strict parser fails on the real file the
    # first time somebody runs metaannot without a config.
    txt = _norm(_text(README))
    assert "`config_path` is the absolute path" in txt and "`null`" in txt, \
        "the README does not say config_path is nullable"
    # ...and the file really does say null, so the prose is not the drifting
    # half. The run dies on the missing FASTA, but the record is written well
    # before that check.
    proc = subprocess.run([sys.executable, METAANNOT_PY, "run"],
                          capture_output=True, text=True, cwd=str(tmp_path),
                          timeout=300)
    assert proc.returncode == 1 and "proteins_faa not found" in proc.stderr
    raw = _text(os.path.join(str(tmp_path), "results",
                             ".metaannot_state.json"))
    assert '"config_path": null' in raw
    assert json.loads(raw)["_run"]["config_path"] is None


def test_the_docs_say_what_a_superseded_run_stops_doing(ma):
    """`--force-unlock` is used on a directory whose holder may still be
    unwinding - SIGTERM makes a killed run unwind, and unwinding writes - so
    what the old run does next is operator-facing and has to be written down.

    Pinned as BEHAVIOUR, not as call syntax: an earlier version of this test
    asserted the literal `if self.is_still_ours():`, which broke the moment the
    gate grew its third answer and said nothing about whether the two callers
    still asked it. What matters is that both writes consult the one gate and
    that the three answers stay distinguishable, because collapsing "vacant"
    into "somebody else's" cost an unsuperseded run its own final verdict."""
    src = _text(os.path.join(ROOT, "metaannot.py"))
    assert "def is_still_ours(self):" in src, "the ownership gate is gone"
    assert "self.is_still_ours() is True" in src, \
        "the lock release must remove only on a definite yes"
    assert "self.owner.is_still_ours() is False" in src, \
        "the state write must refuse only on a definite no, not on a vacancy"
    assert "this run no longer holds" in src

    # The three answers, exercised rather than grepped.
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        lk = ma.ResultsLock(os.path.join(d, "x.lock"))
        lk.__enter__()
        assert lk.is_still_ours() is True, "our own lock must read as ours"
        with open(lk.path, "w", encoding="utf-8") as fh:
            fh.write('{"pid": 1, "host": "somewhere-else", "started": "z"}')
        assert lk.is_still_ours() is False, "another run's lock must read False"
        os.remove(lk.path)
        assert lk.is_still_ours() is None, "a vacant path must be its own answer"

    txt = _norm(_text(README))
    assert "A run that is unwinding stops writing when it is superseded" in txt
    assert "this run no longer holds" in txt, \
        "the README does not quote the line the operator will actually see"
    # The honest version: the stage output is NOT protected, and saying it is
    # was wrong - the record in the file by then belongs to the replacement.
    assert "under the new run's valid signature" in txt, \
        "the README does not say the stage output can land under a live record"
    assert "remains unsupported" in txt, \
        "the README does not say --force-unlock on a live run is unsupported"
    assert "this run no longer holds" in _norm(_text(TUTORIAL))


# ----------------------------------------------------------------------
# the two signal paths, which the docs described as one
# ----------------------------------------------------------------------
def test_the_docs_describe_the_sigterm_that_actually_happens(ma):
    """A `kill`ed run does NOT take Ctrl-C's path, and both docs said it did.

    `main()` installs an unwinding SIGTERM handler and then `cmd_run`
    registers a different one over it, which releases the lock, writes one
    pre-encoded line to fd 2 and `os._exit(128 + N)`s. So for `run` and `all`
    there is no `interrupted` message, no `_run` stamp, and no waiting for the
    stage in flight - and an operator who has read "the same `interrupted`
    message, the same `_run` stamp" goes looking in the log for a line that
    was never written, or waits for an exit that already happened.

    Pinned from the handler rather than from a sentence about it: what makes
    the prose true is that this handler does not unwind, so that is what is
    read off the source here.
    """
    src = _text(METAANNOT_PY)
    i = src.index("def _release_lock_on_signal(sig, _frame):")
    handler = src[i:src.index("os._exit(128 + int(sig))", i)]
    assert "release_results_lock()" in handler, \
        "the handler no longer releases the lock; the docs below say it does"
    for unwinds in ("raise KeyboardInterrupt", "sys.exit(", "stamp_run("):
        assert unwinds not in handler, (
            f"{unwinds!r} in the signal handler: it unwinds or stamps after "
            "all, and README/TUTORIAL describe one that does neither")

    for doc in (README, TUTORIAL):
        txt = _norm(_text(doc))
        for gone in ("the same `interrupted` message, the same `_run` stamp",
                     "unwind the run exactly as Ctrl-C does",
                     "the stage already running has to finish first",
                     "a stage already running has to finish before the "
                     "process exits"):
            assert gone not in txt, \
                f"{os.path.basename(doc)} still says: {gone}"

    readme = _norm(_text(README))
    assert "does not unwind" in readme, \
        "the README does not say the SIGTERM handler skips the unwind"
    assert "no `interrupted` message, no `_run` stamp" in readme, \
        "the README does not say which two traces a `kill` does not leave"
    # and the exit-status convention has to survive the rewrite, because it is
    # the only channel a supervisor reads. See
    # test_the_documented_signal_exit_codes_are_the_conventional_ones.
    assert "`143` for `SIGTERM`" in readme


def test_a_run_record_left_at_running_is_documented_as_an_ordinary_kill(ma):
    # symptom: the README read `"final_status": "running"` with a stale
    # `last_seen` as "a run that was SIGKILLed or lost its machine". That is
    # the ORDINARY trace of `kill`: the handler os._exit()s, so nothing stamps
    # a verdict on the way out and `running` is simply the last thing the
    # record was told. Read the old way, every `systemctl stop` looks like a
    # dead machine, and the natural next move is --force-unlock on a directory
    # that never needed it.
    txt = _norm(_text(README))
    assert "is a run that was `SIGKILL`ed or lost its machine" not in txt, \
        "the README still reads a `running` record as proof of a hard death"
    assert "ordinary** trace of a `kill`" in txt, \
        "the README does not say what a `running` record usually means"
    # the code fact that makes it so: nothing between the release and the exit
    # writes a verdict, and `finished` therefore stays null.
    src = _text(METAANNOT_PY)
    i = src.index("def _release_lock_on_signal(sig, _frame):")
    assert "stamp_run(" not in src[i:src.index("os._exit(128 + int(sig))", i)]
    # ...and the TUTORIAL has to tell an operator what to do with one.
    tut = _norm(_text(TUTORIAL))
    assert "normal appearance of a `kill`" in tut, \
        "the TUTORIAL does not say a `running` record after a kill is normal"


# ----------------------------------------------------------------------
# the test counts, which are the only yardstick a reader has
# ----------------------------------------------------------------------
def _collected(*args):
    """(selected, deselected) for one pytest selection of this suite.

    A collection, not a run: the point is the size of each selection, which is
    what the README quotes, and collecting is seconds where running is
    minutes. `-p no:cacheprovider` so a nested pytest does not write over the
    outer one's cache.

    PYTEST_ADDOPTS is emptied for the child, and that is not tidiness. It is how
    CI passes `-m slow` or `-n auto` to a run, and pytest applies it to THIS
    subprocess as much as to the outer one - so under such a job the inner
    collection is a different selection from the one the README describes, and
    the count test either fails on a machine where the suite is fine or, worse,
    quietly measures a selection nobody documented. The README's numbers are
    about the DEFAULT selection, so the default selection is what is collected,
    whatever the job that invoked us wanted for itself. Everything else in the
    environment is inherited: the interpreter has to find pandas and pyyaml.
    """
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q",
         "-p", "no:cacheprovider", *args],
        capture_output=True, text=True, cwd=ROOT, timeout=600,
        env={**os.environ, "PYTEST_ADDOPTS": ""})
    m = re.search(r"(\d+)/(\d+) tests collected \((\d+) deselected\)",
                  proc.stdout)
    if not m:
        m2 = re.search(r"(\d+) tests collected", proc.stdout)
        assert m2, f"could not read a collection count from:\n{proc.stdout[-2000:]}"
        return int(m2.group(1)), 0
    return int(m.group(1)), int(m.group(3))


def _near(documented, real, what):
    """Assert a documented count is within 10% of the real one.

    A band, not equality, and the width is the whole argument. The defect this
    exists for was a README quoting 589 passed / 35 skipped / 8 xfailed / 33
    deselected against a suite that had grown by roughly 490 tests - a fifth of
    the suite missing from the one yardstick a reader has for deciding whether
    their checkout is sound. A reader who runs the suite and sees 1155 where
    1111 was promised shrugs; one who sees 1155 where 589 was promised goes
    looking for the four hundred tests they think they are missing.

    Equality would catch that too, and would also make every commit that adds a
    single test a two-file commit for ever. A pin that expensive gets deleted,
    and a deleted pin is how the number drifted by a fifth in the first place.
    Ten percent is wide enough that a release's worth of new tests does not
    trip it and narrow enough that the drift above fails by a factor of four.
    """
    assert abs(documented - real) <= max(2, real * 0.10), (
        f"the README says {documented} {what}; this suite has {real}. "
        "Rerun the numbers in the Tests section and write down what you saw.")


def test_the_readme_test_counts_are_the_counts_this_suite_really_has():
    """The Tests section's numbers, against a real collection of this suite.

    Collection rather than a run, because the point is the SIZE of each
    selection, which is what the README quotes, and collecting takes seconds
    where running takes minutes. Every test is collected whether it goes on to
    pass, skip or xfail, so the README's three outcome counts have to sum to
    the selected total and its `deselected` has to be the `slow` selection.

    The last pair - what the same run reports on a machine with no R - is
    checked as exact arithmetic on the others rather than against a band,
    because it is not a measurement at all: it is the same suite counted twice,
    and if the two readings do not subtract to each other one of them was
    guessed.
    """
    txt = _norm(_text(README))
    m = re.search(r"\*\*(\d+) passed, (\d+) skipped, (\d+) xfailed, "
                  r"(\d+) deselected\*\*", txt)
    assert m, "the README no longer states what a healthy default run reports"
    passed, skipped, xfailed, deselected = (int(g) for g in m.groups())

    selected, really_deselected = _collected()
    _near(passed + skipped + xfailed, selected, "tests in the default run")
    _near(deselected, really_deselected, "deselected (the `slow` marker)")

    r = re.search(r"selects the (\d+) R tests", txt)
    assert r, "the README no longer says how many R tests there are"
    r_selected, _ = _collected("-m", "R")
    _near(int(r.group(1)), r_selected, "R tests")

    no_r = re.search(r"reports (\d+) passed and (\d+) skipped", txt)
    assert no_r, "the README no longer gives the counts on a machine with no R"
    assert int(no_r.group(1)) == passed - int(r.group(1)), \
        "the no-R passed count is not the default one minus the R tests"
    assert int(no_r.group(2)) == skipped + int(r.group(1)), \
        "the no-R skipped count is not the default one plus the R tests"


# ----------------------------------------------------------------------
# the console, which shipped in this release
# ----------------------------------------------------------------------
CONSOLE_PY = os.path.join(ROOT, "console", "console.py")
GUI_DESIGN = os.path.join(ROOT, "docs", "gui-design.md")


def _readme_console_section():
    """The README's console section, heading to next heading."""
    txt = _text(README)
    i = txt.find("## The console")
    assert i != -1, "README.md has no console section"
    j = txt.find("\n## ", i + 1)
    return txt[i:j if j != -1 else len(txt)]


def test_the_readme_documents_the_console_this_release_ships():
    """symptom: `console/console.py` is a top-level program of several
    thousand lines shipping in v0.5.0, and README, TUTORIAL and CLAUDE.md
    mentioned it nowhere - so the only way to find out it existed was to list
    the repository. Everything asserted here is a property the console's own
    code has, not a slogan: the three refusals are what make it safe to point
    at a running job, and a reader who does not know about them will not point
    it at one.
    """
    section = _norm(_readme_console_section())
    src = _text(CONSOLE_PY)

    # what it is
    assert "console/console.py" in section
    for claim in ("read-only", "stdlib-only", "`scp`"):
        assert claim in section, f"the README's console section omits {claim}"

    # it writes nothing into a results directory, and there is nothing that acts
    assert "never writes a byte into a results directory" in section
    assert "def do_POST" not in src, \
        "the console grew a write route; the README says it has none"
    assert "no `do_POST`" in section
    assert 'add_argument("--force-unlock"' not in src
    assert "no button that acts" in section
    assert "`--force-unlock`, no launch" in section

    # it does not import metaannot; the contract is where it gets everything
    assert "never imports metaannot" in section
    assert not re.search(r"^\s*(import|from)\s+metaannot\b", src, re.M), \
        "the console imports metaannot; the README says it never does"
    assert "describe --json" in section

    # A UNIX socket, never a TCP port - read off the CODE, with every comment
    # and string literal blanked out first. `127.0.0.1` is a string this
    # console has to be able to TALK about: the reason a loopback listener was
    # rejected is the single most load-bearing sentence in docs/gui-design.md,
    # and the day somebody repeats it in a comment here is not the day the
    # console grew a TCP port. A raw scan over the whole source would have
    # failed on the explanation of the very property it is checking. What
    # would really be a TCP port is an AF_INET socket or a TCPServer /
    # HTTPServer base, so that is what is asserted against, alongside the
    # address literals themselves.
    code = _code_only(CONSOLE_PY)
    assert "UnixStreamServer" in code, \
        "the console no longer serves on a UNIX socket"
    for tcp in ("AF_INET", "TCPServer", "HTTPServer", "127.0.0.1", "0.0.0.0"):
        assert tcp not in code, \
            f"{tcp} is in the console's code, not just in its prose: the " \
            "README says it binds a mode-0700 UNIX socket and never a TCP port"
    assert "mode-0700 UNIX socket" in section and "never a TCP port" in section
    assert "ssh -N -L" in section, \
        "the README does not show how to reach the socket"

    # and it hands over evidence rather than a verdict
    assert "will not tell you a run is dead" in section
    assert "advisory" in section
    assert "ps -p" in section, \
        "the README does not say the console hands over the `ps` line"


def test_every_flag_the_console_has_is_named_in_the_readme():
    # the same rule test_every_documented_command_and_flag_exists_in_the_cli
    # applies to the engine, run the other way: a flag nobody documented is a
    # flag nobody uses, and --root is the difference between watching one
    # project and watching the eight on the machine.
    help_text = subprocess.run([sys.executable, CONSOLE_PY, "--help"],
                               capture_output=True, text=True,
                               timeout=120).stdout
    flags = set(re.findall(r"(--[\w-]+)", help_text)) - {"--help", "--version"}
    assert flags, "console.py --help printed no flags"
    section = _norm(_readme_console_section())
    for flag in sorted(flags):
        assert f"`{flag}" in section, \
            f"the console has {flag} and the README does not mention it"
    # the bare positional too: it is how most people will name one directory
    assert "bare arguments" in section


def test_the_console_version_is_documented_as_separate_from_the_tool_version(
        ma):
    # the same decoupling SIGNATURE_VERSION has, for the same reason: one
    # repository, two programs. A console release must not imply an engine
    # release, and an engine release must not silently renumber the console.
    src = _text(CONSOLE_PY)
    m = re.search(r'^CONSOLE_VERSION = "([^"]+)"', src, re.M)
    assert m, "console.py no longer declares a CONSOLE_VERSION"
    assert m.group(1) != ma.__version__, \
        "CONSOLE_VERSION now tracks __version__; the docs say it does not"
    txt = _norm(_text(README)) + " " + _norm(_text(CLAUDE))
    assert "CONSOLE_VERSION" in txt, "CONSOLE_VERSION is documented nowhere"
    assert f"`{m.group(1)}`" in _norm(_readme_console_section()), \
        "the README states a CONSOLE_VERSION the console does not have"


def test_the_tutorial_meets_the_console_where_an_operator_would():
    # phase 5b is where the tutorial tells you how to watch a run, and the two
    # things it offered were `tail -f` and a hand-rolled json.dumps. The
    # console belongs next to those, not in a document an operator reads once.
    tut = _norm(_text(TUTORIAL))
    assert "console/console.py" in tut, "the TUTORIAL never names the console"
    i = tut.index("console/console.py")
    near = tut[max(0, i - 1200):i + 1600]
    assert "read-only" in near and "ssh -N -L" in near, \
        "the TUTORIAL names the console without saying how to reach it"


def test_the_design_doc_no_longer_says_the_console_does_not_exist():
    """symptom: docs/gui-design.md opened "Status: **design, not built.**
    Nothing here exists yet" in a tree that ships the thing, and described a
    front end serving to `127.0.0.1` - the one decision the build reversed. A
    design record is allowed to be superseded; it is not allowed to say the
    product does not exist.
    """
    txt = _norm(_text(GUI_DESIGN))
    for gone in ("Status: **design, not built.**", "Nothing here exists yet"):
        assert gone not in txt, f"gui-design.md still says: {gone}"
    assert "shipped" in txt, "gui-design.md does not say what shipped"
    # A design record may keep a superseded decision, but not unanswered: the
    # `127.0.0.1` sentences are the ones a reader would otherwise act on, so
    # every section that still carries one has to carry the correction too.
    raw = _text(GUI_DESIGN)
    sections = re.split(r"\n(?=## )", raw)
    for sec in sections:
        if "127.0.0.1" not in sec:
            continue
        head = sec.splitlines()[0]
        assert "**Built:**" in sec, \
            f"'{head}' still describes a 127.0.0.1 listener and does not say " \
            "what was built instead"
    assert "UNIX socket" in txt
    assert "never a TCP port" in txt or "no TCP mode" in txt
    # M2 is not in this release and the document has to say so, or the
    # preflight screen becomes a thing people look for.
    assert "Not in v0.5.0" in txt


def test_the_changelog_says_which_console_milestone_shipped(ma):
    # symptom risk: "the console shipped" reads as "the console is finished".
    # M1 is a watcher; the preflight checklist, the directory picker and
    # config authoring are later milestones, and a changelog that does not say
    # so sends people looking for screens that do not exist.
    txt = _norm(_text(CHANGELOG))
    m = re.search(r"## v0\.5\.0 — (\d{4}-\d{2}-\d{2})(.*?)## v0\.4\.0", txt,
                  re.S)
    assert m, "there is no v0.5.0 section above v0.4.0 in CHANGELOG.md"
    entry = m.group(2)
    assert "milestone M1" in entry, \
        "the changelog does not say which milestone the console is"
    assert "not in this release" in entry, \
        "the changelog does not say what the console cannot yet do"
    for later in ("preflight", "M2"):
        assert later in entry, f"the changelog never mentions {later}"


# ----------------------------------------------------------------------
# the Windows paragraph, which nothing had ever read against the markers
# ----------------------------------------------------------------------
TESTS_DIR = os.path.join(ROOT, "tests")
NUMBER_WORDS = {"no": 0, "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
                "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
                "eleven": 11, "twelve": 12}


def _windows_paragraph():
    """The README's Windows note, from its first sentence to the next topic."""
    txt = _text(README)
    i = txt.find("Windows is not run from this machine")
    assert i != -1, "the README's Windows note is gone"
    j = txt.find("\nOffline, and needs none of the external tools", i)
    assert j != -1, "the README's Windows note no longer ends where it did"
    return _norm(txt[i:j])


def _skips_on_windows(fn):
    """The `reason` of this test's `skipif(os.name == "nt", ...)`, or None.

    The marker is matched on its argument rather than on its text, so a skip
    written `sys.platform == "win32"` would not be mistaken for one, and a
    reason reflowed by an editor still reads the same.
    """
    for dec in fn.decorator_list:
        if not isinstance(dec, ast.Call):
            continue
        if getattr(dec.func, "attr", "") != "skipif":
            continue
        if not any("'nt'" in ast.dump(a) for a in dec.args):
            continue
        for kw in dec.keywords:
            if kw.arg == "reason" and isinstance(kw.value, ast.Constant):
                return kw.value.value
        return ""
    return None


def _signal_tests():
    """({name: skip reason}, {name: file}) for tests that signal a child.

    DERIVED, every time, from the test files themselves. The Windows note is
    the one paragraph in this README that nobody can check by running anything
    - the machine that writes it is not the machine it describes - so the only
    thing standing between it and drift is a scan like this one. The paragraph
    it replaced said "the two tests that send a signal to a child" and gave one
    reason for both, in a tree that had three such functions carrying two
    different reasons, four collected items between them, and six more sending
    signals with no Windows marker at all. Every one of those numbers is read
    off the files here, so the next person to add a signal test cannot leave
    the paragraph stale without this failing and saying which name is missing.

    "Signals a child" is `send_signal` lexically inside the test function: that
    is what the note is about, and a helper that kills a process on the way out
    of a failure is not.
    """
    skipped, unmarked = {}, {}
    for base in sorted(os.listdir(TESTS_DIR)):
        if not (base.startswith("test_") and base.endswith(".py")):
            continue
        path = os.path.join(TESTS_DIR, base)
        for fn in ast.walk(ast.parse(_text(path))):
            if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if not fn.name.startswith("test_"):
                continue
            if not any(isinstance(n, ast.Call)
                       and isinstance(n.func, ast.Attribute)
                       and n.func.attr == "send_signal"
                       for n in ast.walk(fn)):
                continue
            reason = _skips_on_windows(fn)
            if reason is None:
                unmarked[fn.name] = base
            else:
                skipped[fn.name] = (base, reason)
    assert skipped, "no test in this suite is skipped on Windows for a signal"
    return skipped, unmarked


def _collect_ids(*paths):
    """Every node id pytest collects from `paths`, parametrization included.

    The README counts collected ITEMS, not functions, and one of the three is
    parametrized over two signals - so the arithmetic that turns three names
    into four skips has to come from a real collection rather than from an
    assumption about how many cases a decorator carries.
    """
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q",
         "-p", "no:cacheprovider", *paths],
        capture_output=True, text=True, cwd=ROOT, timeout=600,
        env={**os.environ, "PYTEST_ADDOPTS": ""})
    ids = [ln.strip() for ln in proc.stdout.splitlines() if "::" in ln]
    assert ids, f"collected nothing from {paths}:\n{proc.stdout[-2000:]}"
    return ids


def test_the_readme_names_every_signal_test_windows_really_skips():
    """The Windows note, against the markers it is describing.

    Nothing checked this paragraph before, and it showed: it named two tests
    where there were three, gave the SIGINT reason for a pair skipped over
    `TerminateProcess`, and did not mention that six further tests signal a
    child with no marker at all - so a Windows reader was told the signal tests
    were handled there when most of them are not.
    """
    para = _windows_paragraph()
    skipped, unmarked = _signal_tests()

    missing = [n for n in sorted(skipped) if f"`{n}`" not in para]
    assert not missing, (
        "the README's Windows note does not name %s, which %s skipped there "
        "for a signal reason. Add it, and correct the counts in that "
        "paragraph."
        % (", ".join(missing), "is" if len(missing) == 1 else "are"))

    m = re.search(r"(\w+) test functions that signal a child process are "
                  r"skipped there as well, (\w+) collected items", para)
    assert m, "the README no longer counts the Windows-skipped signal tests"
    assert NUMBER_WORDS.get(m.group(1).lower()) == len(skipped), (
        f"the README says {m.group(1)} test functions are skipped on Windows "
        f"for a signal; this suite has {len(skipped)}: {sorted(skipped)}")

    files = sorted({os.path.join("tests", b) for b, _r in skipped.values()})
    ids = _collect_ids(*files)
    items = [i for i in ids
             if i.split("::")[-1].split("[")[0] in skipped]
    assert NUMBER_WORDS.get(m.group(2).lower()) == len(items), (
        f"the README says {m.group(2)} collected items are skipped on "
        f"Windows; pytest collects {len(items)}: {items}")

    # Singular as well as plural: the count is what is being pinned, and it
    # reached one the moment five of the six unmarked tests were marked. A
    # regex that only reads "tests ... carry" would have forced the prose to
    # stay ungrammatical to keep the check alive, which is the wrong way round.
    u = re.search(r"\*\*(\w+) further tests? signals? a child and carr(?:y|ies) "
                  r"no Windows marker", para)
    assert u, "the README no longer says how many signal tests are NOT skipped"
    assert NUMBER_WORDS.get(u.group(1).lower()) == len(unmarked), (
        f"the README says {u.group(1)} signal tests carry no Windows marker; "
        f"this suite has {len(unmarked)}: {sorted(unmarked)}")
    # and it has to name one of the unmarked ones, not only count them: a
    # reader who is told "six others exist" and given no name cannot find out
    # which. Getting the two lists the wrong way round is the defect being
    # pinned here, so the name it gives has to come from the right one.
    named = set(re.findall(r"`(test_\w+)`", para))
    assert named & set(unmarked), (
        "the README's Windows note names no unmarked signal test, so a reader "
        "cannot tell which ones will simply run there. One of %s would do."
        % sorted(unmarked))

    # Two different reasons, and the paragraph has to carry both. Taken from
    # the markers rather than restated here: each reason's own CamelCase and
    # ALL_CAPS identifiers are what a reader would match on.
    for name, (_base, reason) in sorted(skipped.items()):
        ident = set(re.findall(r"\b(?:[A-Z][a-z]+[A-Z]\w*|[A-Z][A-Z_]{2,})\b",
                               reason or ""))
        if not ident:
            continue                  # "see above" defers to the one above it
        assert ident & set(re.findall(r"[A-Za-z_]+", para)), (
            f"{name} is skipped on Windows because of {sorted(ident)}, and "
            "the README's Windows note mentions none of those")


def test_the_retired_windows_path_test_is_really_separator_agnostic():
    # symptom: the note used to list "a path test asserting forward slashes"
    # among the Windows failures, and the rewrite dropped the sentence rather
    # than settling it. v0.4.0 fixed that test; the note now says so, and says
    # which test and how - so the claim is checkable, which is the point.
    para = _windows_paragraph()
    assert "normcase" in para, \
        "the README no longer says how the old path failure was retired"
    m = re.search(r"v0\.4\.0 rewrote `(test_\w+)`", para)
    assert m, "the README does not name the test it says was rewritten"
    name = m.group(1)
    found = _named_test(name)
    assert found, f"the README names {name}, which does not exist"
    for base, fn in found:
        body = ast.get_source_segment(_text(os.path.join(TESTS_DIR, base)), fn)
        assert "normcase" in body, (
            f"{name} in tests/{base} no longer normalises the separator, so "
            "the README's claim that it passes on Windows is a guess again")


def _named_test(name):
    """(basename, ast node) wherever `name` is defined under tests/."""
    out = []
    for base in sorted(os.listdir(TESTS_DIR)):
        if not (base.startswith("test_") and base.endswith(".py")):
            continue
        tree = ast.parse(_text(os.path.join(TESTS_DIR, base)))
        for fn in ast.walk(tree):
            if isinstance(fn, ast.FunctionDef) and fn.name == name:
                out.append((base, fn))
    return out


def test_the_uc_dark_rescue_numbers_reconcile_with_each_other():
    """README quoted two incompatible eggNOG-only dark figures for the same run.

    One section said 9,731 dark on eggNOG alone falling to 1,462, "8,686
    rescued" -- which does not even self-reconcile, since 9,731 - 1,462 =
    8,269. Another section said 17,377 proteins had no eggNOG KO and 6,271 of
    them carried an eggNOG PFAMs entry, which implies 11,106. Measured off the
    UC run's own outputs, 11,106 is right: 38,204 proteins, 29,084 in the
    emapper table, and 11,106 with no KO, no PFAMs and no CAZy. 1,462 is right
    too, from bin_summary.tsv and annotation_final.tsv independently.
    """
    import re
    txt = _norm(_text(README))
    m = re.search(r"took ([\d,]+) proteins dark on eggNOG alone down to "
                  r"([\d,]+) — ([\d,]+) rescued", txt)
    assert m, "the dark-rescue sentence moved; re-measure before rewording it"
    before, after, rescued = (int(g.replace(",", "")) for g in m.groups())
    assert before - after == rescued, (before, after, rescued)
    # and it has to agree with the KO/PFAMs join quoted elsewhere
    j = re.search(r"([\d,]+) of ([\d,]+) dark proteins", txt)
    assert j, "the eggNOG-PFAMs join sentence moved"
    with_pfam, no_ko = (int(g.replace(",", "")) for g in j.groups())
    assert no_ko - with_pfam == before, (no_ko, with_pfam, before)
