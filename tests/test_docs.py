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
import tempfile
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
        # Six hours, and the README argues the number rather than stating it.
        # It is here because it is the one setting whose default an operator
        # reads as a PROMISE - "a FIFO nobody writes fails by the morning" -
        # and a default that moved without the prose moving would make that
        # promise silently false.
        "fifo_wait_s": 21600,
    }
    for k, v in checks.items():
        assert d[k] == v, f"{k} changed to {d[k]}; update the docs too"
    a = d["analysis"]
    assert a["fdr"] == 0.05 and a["min_valid_per_group"] == 3
    assert a["taxon_min_proteins"] == d["taxon_min_proteins_for_factor"], \
        "the README documents these two as deliberately equal"
    assert d["thresholds"]["smorf_max_len"] == 100
    assert "hard floor of 90 nt" in _text(README)
    # ...and the FIFO wait is documented in both the seconds the config takes
    # and the hours the prose argues, because an operator plans around one and
    # writes the other.
    rd = _norm(_text(README))
    assert f"fifo_wait_s: {d['fifo_wait_s']}" in rd, \
        "the README's `fifo_wait_s` example has drifted from DEFAULT_CONFIG"
    # ...and every document that gives the wait in HOURS gives the same
    # number of them, in either spelling. None of those sentences is in the
    # count scanner's reach - a number followed by a unit is a measurement to
    # that scan, not a count - so this is the only thing holding them to the
    # config.
    n = d["fifo_wait_s"] // 3600
    word = {v: k for k, v in NUMBER_WORDS.items()}[n]
    for name, txt in (("README.md", rd), ("TUTORIAL.md", _norm(_text(TUTORIAL))),
                      ("CLAUDE.md", _norm(_text(CLAUDE)))):
        assert re.search(rf"\b(?:{word}|{n}) hours?\b", txt, re.I), \
            f"{name} does not say the FIFO wait is {n} hours, which is the " \
            "half an operator plans a working day around"


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


def test_the_readme_quotes_the_scope_statement_the_engine_authors(ma):
    # It is published as a literal string for a front end to render verbatim,
    # and the README quotes it as a block quote - which makes the README a
    # second home for the sentence, and a second home is where a claim goes
    # stale. The quote is checked against the constant rather than read.
    txt = _norm(_text(README)).replace("> ", "")
    assert _norm(ma.DOCTOR_SCOPE_STATEMENT) in txt, \
        "the README's scope quote is not the sentence metaannot.py publishes"


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


# `doctor`'s own section headings. The README paragraph that tells a console
# author what doctor checks BEYOND `requirements` has to name every one of
# them, because it reads as a closed list and a preflight screen gets built
# from it.
#
# Taken from ma.DOCTOR_SECTIONS rather than scraped out of the source with a
# regex, which is what this did until the sections became data: the scrape was
# a proxy for "the headings doctor prints", and now there is a real list that
# both the printed report and `doctor --json` are rendered from. A proxy that
# has been replaced by the thing it stood for is a test measuring the wrong
# object.
DOCTOR_SECTIONS = ("config", "inputs", "precomputed emapper", "tools",
                   "databases", "gpu", "tmt", "manifest", "taxonomy",
                   "resources", "R")


def test_the_readme_names_every_section_doctor_prints(ma):
    # symptom: the list omitted `== manifest ==`, `== resources ==` and
    # `== R ==` outright, so a preflight built from it is green for a config
    # doctor fails on the manifest-to-column mapping, on a memory split that
    # cannot give eggNOG --dbmem, or on a missing required R package.
    printed = [title for _id, title in ma.DOCTOR_SECTIONS]
    assert printed, "doctor no longer prints section headings"
    assert tuple(printed) == DOCTOR_SECTIONS, \
        f"doctor's sections have changed: {set(printed) ^ set(DOCTOR_SECTIONS)}"
    txt = _norm(_text(README))
    for name in printed:
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


def test_the_tutorial_says_what_a_fifo_really_does_to_a_run():
    """symptom, twice over, in the same paragraph.

    The first version said "each stage dies as it opens the path. A FIFO or a
    socket behaves the same way", of a state where nothing died at all: the
    open BLOCKED, and `run` was measured on one at every configured input with
    no output, no traceback and no exit status. The second version said it
    HANGS, which was true when it was written and is not true now - opener()
    waits `fifo_wait_s` and then dies naming the path.

    The THIRD version is the one this pins, and it is a fourth fact rather
    than a rewording of the other three: a verifier drove the workflow the
    second version recommended and it WORKS AT ONE of the operator-supplied
    inputs. At the rest the run opens the path again after the pipe has been
    drained, so the paragraph was selling a six-hour wait for a failure. A
    paragraph that says "a FIFO is read" without saying WHERE is the version
    that sends somebody to pipe their FASTA in.

    So this asserts the paragraph carries every fact an operator has to have
    before leaving a run overnight: that a pipe works only where the run reads
    that input once, that it is refused at once where it does not, that a
    writer with nothing to say ends the run, and that `doctor` never waits at
    all.
    """
    tut = _norm(_text(TUTORIAL))
    assert "behaves the same way" not in tut, \
        "the TUTORIAL again lumps a FIFO in with the directory beside it"
    i = tut.index("A **FIFO**")
    near = tut[i:i + 2600]
    assert "live writer" in near and "mkfifo" in near, \
        "the TUTORIAL does not say a FIFO with a writer is read at all"
    assert re.search(r"reads exactly once|read(s)? (it )?exactly once", near), \
        "the TUTORIAL does not say a pipe works only at a single-read input"
    assert "refused immediately" in near or "refused at once" in near, \
        "the TUTORIAL does not say a multi-read input refuses a pipe at once"
    assert "proteins_faa" in near and "quant_table" in near, \
        "the TUTORIAL does not name which inputs are which"
    assert "no writer" in near and "fifo_wait_s" in near, \
        "the TUTORIAL does not say what ends the wait, or what sets it"
    assert re.search(r"\bdies\b", near), \
        "the TUTORIAL does not say the run ends on an unwritten FIFO"
    assert "O_NONBLOCK" in near, \
        "the TUTORIAL does not say why `doctor` itself is immune"


def test_no_published_sentence_claims_a_fifo_still_hangs_the_run():
    """The class the sentence above belongs to, over every published surface,
    and it is the exact reverse of the scan it replaces.

    That scan read: a sentence naming a FIFO and claiming a death is wrong.
    It was right while the open blocked, and it is the wrong way round now -
    a FIFO with no writer DOES kill the stage, after `fifo_wait_s`, and the
    dangerous sentence is the one that still promises a hang, because it sends
    an operator to wait for a run that has already exited. Verbs go stale in
    whichever direction the code moves; what does not go stale is that the
    published verb has to be the one the code has.
    """
    bad = []
    for name, path in (("README.md", README), ("TUTORIAL.md", TUTORIAL),
                       ("CHANGELOG.md", CHANGELOG)):
        for sent in re.split(r"(?<=[.!?]) ", _norm(_text(path))):
            if "FIFO" not in sent:
                continue
            # QUOTED text is exempt, and has to be: this file's own entries
            # exist to reproduce the wrong sentences - `said "it dies reading
            # this"` - and a scan that cannot tell a quotation from a claim
            # would make the record of the defect impossible to write down.
            claim = re.sub(r'"[^"]*"|\u201c[^\u201d]*\u201d|`[^`]*`', " ",
                           sent).lower()
            if not re.search(r"\b(hangs?|hung|hanging|blocks? forever|"
                             r"never returns?|indefinitely)\b", claim):
                continue
            # A sentence saying it does NOT hang, or saying it USED to, is the
            # record of the fix rather than a live claim.
            if re.search(r"\b(not|never|no longer|rather than|instead of|nor|"
                         r"used to|before|until|would have)\b", claim):
                continue
            bad.append(f"{name}: {sent}")
    assert not bad, "\n".join(bad)


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
# NUMBER_WORDS used to be defined HERE as well, by hand and up to twelve, and
# the generated one two hundred lines below shadowed it - so the three
# assertions above that call NUMBER_WORDS.get() were running against a table
# that had never had the `"no": 0` this one was written for. A test whose
# table is not the table it thinks it is pins whatever the other table
# happens to say. There is one now, it is generated, and `"no"` is in it.


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


# ----------------------------------------------------------------------
# the changelog's Fixed section, which claimed a count it did not have
# ----------------------------------------------------------------------
_ONES = ("", "one", "two", "three", "four", "five", "six", "seven", "eight",
         "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
         "sixteen", "seventeen", "eighteen", "nineteen")
_TENS = {20: "twenty", 30: "thirty", 40: "forty", 50: "fifty", 60: "sixty",
         70: "seventy", 80: "eighty", 90: "ninety"}
# Generated rather than typed out, which is this section's own rule applied to
# itself: a table of number words that stops at twenty plus one hand-added
# "thirty-three" is a table that fails the day a count passes it, and the
# failure reads as a KeyError rather than as a drifted number.
NUMBER_WORDS = {w: i for i, w in enumerate(_ONES) if w}
# "no counts as zero" for the sentences that say "no unmarked items", which is
# the one word this table needs that counting cannot produce.
NUMBER_WORDS["no"] = 0
NUMBER_WORDS.update({w: n for n, w in _TENS.items()})
NUMBER_WORDS.update({f"{w}-{o}": n + i
                     for n, w in _TENS.items()
                     for i, o in enumerate(_ONES) if o})


def _unreleased_fixed_groups():
    """(heading, entry count) for each `#### ` group under Unreleased/Fixed."""
    txt = _text(CHANGELOG)
    start = txt.index("## Unreleased")
    fixed = txt.index("### Fixed", start)
    end = txt.index("\n### Changed", fixed)
    out = []
    for part in re.split(r"^#### ", txt[fixed:end], flags=re.M)[1:]:
        head = part.split("\n", 1)[0].strip()
        out.append((head, len(re.findall(r"^\*\*", part, flags=re.M))))
    return out


def test_the_changelog_fixed_section_is_countable(ma):
    """symptom: the Unreleased section said the count of false verdicts "is
    now the number of entries under Fixed". There were twenty-one bolded
    entries and fourteen false verdicts, and the two were never going to be
    the same number: a correction to what the document says about ITSELF
    changes no config's exit status and is not a verdict at all. The first
    wave also printed EIGHT entries under a paragraph saying "Seven were in
    the shipped doctor", with nothing marking the eighth as an attribution
    rather than a verdict.

    The section is grouped now, each heading leads with its own count, and
    this counts them.
    """
    groups = _unreleased_fixed_groups()
    assert groups, "the Fixed section is no longer grouped under headings"
    total = 0
    for head, n in groups:
        word = head.split()[0].lower()
        assert word in NUMBER_WORDS, \
            f"the heading {head!r} does not lead with a count"
        assert NUMBER_WORDS[word] == n, \
            f"{head!r} claims {word} entries and has {n}"
        total += n
    txt = _norm(_text(CHANGELOG))
    m = re.search(r"### Fixed (\w+(?:-\w+)?) entries, in (\w+) groups", txt)
    assert m, "the Fixed section no longer says how many entries it has"
    # _count_value() and not NUMBER_WORDS: the table of number words is
    # generated up to ninety-nine and this total has passed it, so the count
    # is written as a NUMERAL - which is what the document already does for
    # every number that big, and what the count scanner already reads.
    assert _count_value(m.group(1)) == total, \
        f"the section claims {m.group(1)} entries and has {total}"
    assert NUMBER_WORDS[m.group(2).lower()] == len(groups)


def test_the_changelog_false_verdict_count_is_the_verdict_groups(ma):
    # The OTHER number, which is the one the prose above the section is
    # actually about: the groups whose headings say "false verdicts". The two
    # counts are different on purpose and the section now says so.
    groups = _unreleased_fixed_groups()
    verdicts = sum(n for head, n in groups if "false verdict" in head.lower())
    waves = [head for head, _ in groups if "false verdict" in head.lower()]
    txt = _norm(_text(CHANGELOG))
    m = re.search(r"turned up \*\*([\w-]+)\*\* live false verdicts, in ([\w-]+) "
                  r"waves", txt)
    assert m, "the changelog no longer states a false-verdict count"
    assert NUMBER_WORDS[m.group(1).lower()] == verdicts, \
        f"the prose claims {m.group(1)} false verdicts; the groups hold " \
        f"{verdicts}"
    assert NUMBER_WORDS[m.group(2).lower()] == len(waves)
    # ...and it no longer claims that number IS the size of the section. The
    # sentence that says they are NOT the same is allowed to name it.
    assert "is now the number of entries under Fixed" not in txt
    assert "they are not the number of entries under fixed" in txt.lower()


def test_the_changelog_does_not_name_functions_that_do_not_exist(ma):
    """symptom: four sentences added by this change set named a
    `read_protein_table` that has never existed in this tool, and two named a
    `taxon_map()` where the function is `resolve_taxonomy()`. A changelog is
    read by someone going to the code next, so a name it invents costs a grep.
    """
    src = _text(METAANNOT_PY) + "".join(
        _text(os.path.join(ROOT, "tests", f))
        for f in sorted(os.listdir(os.path.join(ROOT, "tests")))
        if f.endswith(".py"))
    txt = _text(CHANGELOG)
    start = txt.index("## Unreleased")
    end = txt.index("\n## ", start + 1)
    section = txt[start:end]
    named = set(re.findall(r"`([a-z_][a-z_0-9]{3,})\(\)`", section))
    # Names this tool defines, or attribute calls on something it imports.
    for name in sorted(named):
        assert f"def {name}(" in src or f".{name}(" in src, \
            f"the changelog names {name}(), which is defined nowhere in " \
            "metaannot.py or the suite"


# ======================================================================
# V: every count in prose is DERIVED from the thing it counts, or pinned
# ======================================================================
# Three consecutive rounds of this change set shipped a wrong hand-written
# count, and twice it was inside the comment on the very constant that
# enforces the set being counted:
#
#   * "NULL IS AN EIGHTH VALUE" beside a tuple of nine
#   * the caveat on every TMT annotation row said "three checks" in PUBLISHED
#     text while the scope statement beside it said four
#   * "the 768-config sweep", of a sweep that has eleven configs
#   * TUTORIAL's "that block probes four packages", of a probe that reads
#     fourteen and was widened two releases ago
#
# ...and a FOURTH round shipped one after this scanner was written:
#
#   * "a consumer reading the eight fields it knew about is unaffected by a
#     NINTH appearing", sitting above BOTH of the two keys that had just been
#     added, over a per-stage dict that emits ten
#
# That one is the interesting one, because this scanner was already running
# and did not see it - and it did not see it for three STRUCTURAL reasons,
# each of which is fixed here rather than in the sentence:
#
#   1. DIGITS WERE OUTSIDE THE SCAN. The alternation was number WORDS only, so
#      every count written 200,000 / 200 / 9 / 12 / 14 / 768 was invisible.
#      768 had been caught only by a bespoke test written for that one string.
#   2. tests/ WAS NOT A SURFACE. Two of the three stale counts found in the
#      fourth round lived in this very file, including the one in the
#      docstring of the test that counts the CHANGELOG's own entries.
#   3. THE EMITTED DOCUMENT IS WIDER THAN THE CONSTANTS THIS READ.
#      `DMND_READS` publishes "up to 200,000 records" and "up to 200 deflines"
#      as caveat text on real rows; both were hand-written, one function away
#      from the literals that decide them. Those literals are constants now
#      and the caveats are built from them, and this scan reads them.
#
# ...and a FIFTH round shipped two more, after the scanner had been through
# four readings:
#
#   * README's "thirteen test modules", of a directory holding sixteen
#   * a CHANGELOG entry describing the state sweep as six states and thirty
#     cells, of a sweep that had already outgrown both numbers
#
# NEITHER WAS INVISIBLE BY ACCIDENT, and that is what this round changes. The
# nouns the scan knew were a WHITELIST - a hand-written alternation of the
# things this codebase has a set of - so it could not see "modules", and it
# could not see a count that points at its noun with a pronoun ("all thirty of
# THEM"). Adding `modules` to the list would have fixed one sentence and left
# the class where it was; the list cannot see "vocabularies", "tests",
# "callers", "assertions", "branches", "sentences", "helpers", "functions",
# "guards" or "arms" either, and nobody is going to think of the eleventh.
#
# So the rule is NEGATIVE SPACE now. Anything shaped like a plural noun after
# a number is a count until this table says otherwise, and what is written
# down is the exceptions: units, the words that end in s without being nouns,
# and three classes of number that are not cardinalities at all. See
# _S_NOT_A_NOUN, _count_head() and _count_phrase(), each of which carries its
# own argument.
#
# This is the class: a scanner over the six surfaces a count can live in -
# metaannot.py's comments and docstrings, the document `doctor --json` emits,
# the test suite's own comments and docstrings, the README, the TUTORIAL and
# the CHANGELOG's Unreleased section - which finds every "<number> <plural
# noun>" and requires each to be CLASSIFIED here. The number of surfaces is
# itself one of the counts it scans, and it is DERIVED from _count_surfaces().
#
#   DERIVED  the number is recomputed from the source and compared. If the
#            code grows a value, this test fails and names the sentence.
#   MEASURED not recomputable - a fact about a dataset, a machine or an
#            observation. The pin is that the SENTENCE STILL EXISTS, so an
#            edit cannot strand a measurement whose subject has gone.
#   PROSE    not a cardinality of anything this codebase has ("the first value
#            was renamed", "two columns in different namespaces"). Exempt, and
#            the reason is written down so the exemption is a decision rather
#            than an oversight.
#
# `one` and `no` are outside the scan on purpose: "one row per protein" and
# "no stage dies" are articles and negations, not sums, and including them
# buries the real counts under hundreds of them.
COUNT_WORDS = dict(NUMBER_WORDS, **{
    "first": 1, "second": 2, "third": 3, "fourth": 4, "fifth": 5, "sixth": 6,
    "seventh": 7, "eighth": 8, "ninth": 9, "tenth": 10, "eleventh": 11,
    "twelfth": 12, "twenty-one": 21, "fourteen": 14, "fifteen": 15,
    "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19,
})
# The units a number can carry between itself and a noun. "a 4 MB state file"
# and "a 40 KB file" are MEASUREMENTS of one thing, not counts of four states
# or forty files, and letting the head slot swallow a unit is how a scan that
# reads digits fills up with them.
_UNITS = frozenset("""
mb kb gb tb kib mib gib byte bytes bit bits hour hours minute minutes second
seconds ms residue residues line lines char chars cpu plex aa min px nt kda
times plddt
""".split())
# `times` is in there as a unit rather than a noun because it is a RATIO:
# "four times hungrier", "ten times better", "sixteen times" - a multiplier of
# something, never a number of things.
# Words that END IN S and are not plural nouns. THE NEGATIVE SPACE, and the
# whole point of the rewrite: `_NOUNS` used to be a WHITELIST of the things
# this codebase has a set of, and a whitelist cannot see a noun nobody thought
# of. It could not see modules, vocabularies, tests, callers, assertions,
# branches, sentences, helpers, functions, guards or arms - and two counts
# shipped wrong through the gap, "thirteen test modules" of sixteen and a cell
# count two states stale. Listing the nouns it happened to miss would have
# fixed those sentences and left the class exactly where it was; nobody is
# going to think of the next noun either.
#
# So the rule is inverted. A head is anything that LOOKS like a plural noun,
# and what has to be written down is the exceptions - verbs, possessives and
# adverbs that happen to end in s. That list is English, not this project's
# vocabulary: it does not grow when the code does, which is the property the
# whitelist did not have.
_S_NOT_A_NOUN = frozenset("""
is was has does its this thus his hers theirs ours yours says keeps means
gives leaves goes takes makes runs needs wants reads writes calls names
carries costs holds knows lives looks matters moves passes puts raises reaches
reports returns sends sets shows stays stops turns uses works less across
always perhaps yes plus versus unless whereas nevertheless else
becomes covers exits ships measures satisfies previous various serious dies
""".split())
# A count may point at its noun with a pronoun - "all thirty of them" - and
# that is still a count. It is in fact the shape of one of the two counts this
# rewrite exists for, which is why the pronouns are heads rather than noise:
# the classification such a phrase forces is either a recompute or, better, a
# rewritten sentence that names the thing.
_PRONOUN_HEADS = frozenset(("them", "these", "those"))
# Ordinals are POSITIONS - "the first column", "the fifth reading", "the
# second command above" - and there are hundreds of them. They are out of the
# negative space and INTO the table only when a sentence opts one in, which is
# what the "one past the end" idiom does: "a NINTH format", "a TENTH value",
# "the FIFTH entry" are claims about the size of a set and are recomputed as
# such below. An ordinal the table does not name is a position and is not
# scanned.
_ORDINALS = frozenset("""first second third fourth fifth sixth seventh eighth
ninth tenth eleventh twelfth""".split())
# Generated from COUNT_WORDS, longest first so "twenty-seven" wins over
# "seven". Typing the alternation out by hand is how "twenty-seven entries"
# was read as "seven entries" the first time this test ran.
#
# DIGITS as well as words, which is the fix for the fourth wrong count: a
# number spelled 200,000 was not a number as far as this was concerned. The
# lookbehind is what keeps that from drowning the scan - it refuses a digit
# that is part of a larger token, so `2.4 MB` is not "4 MB", `mode-000` is not
# "000", `cost-3 stages` is not "3 stages" and `#5` is not "5".
_NUMBER = ("|".join(sorted((w for w in COUNT_WORDS if w not in ("one", "no")),
                           key=len, reverse=True))
           + r"|\d+(?:,\d{3})*")
# The number, and the two words after it. TWO, because a count is written
# either bare ("twenty-one stages") or with one adjective ("twenty-four LIVE
# false verdicts", "thirteen TEST modules") and a wider window starts finding
# the plural at the far end of the next clause instead.
_COUNT_RE = re.compile(r"(?<![\w.,#-])(" + _NUMBER + r")((?:[ \-]+[a-z]+){0,2})",
                       re.I)


def _count_head(number, tail):
    """The noun a number counts, or None when it counts nothing nameable.

    `tail` is the two words after the number, with the separator that joined
    each. The head is the first of them that is a plural noun by SHAPE, or a
    pronoun standing in for one, or - for the first word only, and only when a
    hyphen joined it to the number - an attributive singular, because "a
    9-field row" counts fields and says so.

    Returning None is a decision and not a gap: a number with no plural within
    two words is not making a claim about the size of anything this codebase
    has, and where a sentence hides a real count that way the answer is to
    rewrite the sentence so it names what it counts. That is what happened to
    "asserts stdout parses as JSON in all thirty of them", which is now a
    number followed by the noun it is about.
    """
    # Lowercased first, so that a unit written the way units are written -
    # "4 GB resources", "94 GB to WSL2" - is seen as the unit it is. A
    # case-sensitive scan skips it as if it were punctuation and reads the
    # word after it as the head, which is how a machine's memory became a
    # count of four resources.
    for i, (sep, word) in enumerate(re.findall(r"([ \-]+)([a-z]+)",
                                               tail.lower())):
        if word in _UNITS:
            # A unit ENDS the search rather than being stepped over: "2 cpu,
            # 4 GB resources" is a machine's size twice, and a scan that walks
            # past the unit to the next plural reads it as four resources.
            return None
        # A TIGHT hyphen, and a word long enough to be one: "a 9-field row"
        # is attributive, while "median length 15 - and returned 0 hits" and
        # "on one job of seven -- a race" are a number, a dash and the next
        # clause. Spaced dashes read as hyphens is how a scan fills up with
        # "15 and" and "seven a".
        if i == 0 and sep == "-" and len(word) >= 3:
            return word
        if word in _PRONOUN_HEADS:
            return word
        if len(word) >= 3 and word.endswith("s") and word not in _S_NOT_A_NOUN:
            return word
    return None


def _count_phrase(m):
    """"<number> <noun>" for one match, or None when it is out of the scan.

    What is out of the scan is out by CLASS, and every class has its reason
    written here rather than left to be inferred from a table with nothing in
    it:

    * an ORDINAL, unless the table names the phrase - see _ORDINALS.
    * a number with a LEADING ZERO, which in this codebase is a file mode and
      nothing else - 0700 on a socket directory, 0755 on a results directory -
      and a mode is one number about one thing.
    * the digits 0 and 1, for exactly the reason "one" and "no" are out of the
      alternation: "0 disables", "1 if fails else 0", "1 row per protein" are
      values, articles and negations, and including them buries every real
      count under them.
    * TWO, and it is the only cardinal that is out. It is the number this
      codebase writes most and the one it writes least dangerously: "the two
      stage lists", "two files on purpose", "the difference of the two group
      means" - a pair is SPELLED OUT in the sentence that counts it, so the
      number and the list cannot part company the way "thirteen test modules"
      and a directory of sixteen can. Twenty-five of the exemptions this table
      carried before the rule changed said nothing but "a named pair", which
      is a rule being retyped rather than written down. A `two` that IS the
      cardinality of a set - DOCTOR_FAIL_REASONS has two values - is opted
      back in by naming it below, and the table always wins over a class.
    * a number of a thousand or more. Nothing in this codebase has a thousand
      of anything - the largest set it counts is twenty-one stages - so a
      number that big is a measurement of a dataset, a file or a machine, and
      its pin is that the sentence still exists, which is what MEASURED
      already means. The one big number that IS derived, the DIAMOND record
      cap, is named in the table and checked there.
    """
    tok = m.group(1).lower()
    head = _count_head(tok, m.group(2))
    if head is None:
        return None
    phrase = f"{tok} {head}"
    if phrase in COUNT_PROSE:
        return phrase
    if tok in _ORDINALS or tok == "two":
        return None
    if tok[0].isdigit():
        # A LEADING ZERO is a file mode in this codebase and nothing else -
        # 0700 on a socket directory, 0755 on a results directory, 0777 on a
        # parent - and a mode is one number about one thing, never a count.
        if tok.startswith("0"):
            return None
        n = int(tok.replace(",", ""))
        if n in (0, 1) or n >= 1000:
            return None
    return phrase


def _phrase_re(key):
    """Where a classified phrase is allowed to be written.

    The same shape the scan reads, run the other way: the number, an optional
    adjective, the noun. It is what keeps a SINGULAR or ORDINAL entry in the
    table alive - "a ninth format", "a 9-field row" - when the negative-space
    rule above would not have found it, and it is what makes the staleness
    check honest, because a key that matches nothing anywhere is a rule with
    nothing under it.
    """
    num, noun = key.split(" ", 1)
    return re.compile(r"(?<![\w.,#-])" + re.escape(num)
                      + r"[ \-]+(?:(?!(?:" + "|".join(sorted(_UNITS))
                      # {1,14}, not {2,14}: `_count_head` steps over a
                      # one-letter adjective and this has to reach the same
                      # phrase, or "the 35 R tests" is a count the scan finds
                      # and the table cannot locate.
                      + r")[ \-])[a-z]{1,14}[ \-]+)?"
                      + re.escape(noun).replace("\\ ", r"[ \-]+") + r"\b",
                      re.I)


def _count_value(tok):
    """The integer a matched number means, word or digits."""
    tok = tok.lower()
    return COUNT_WORDS[tok] if tok in COUNT_WORDS else int(tok.replace(",", ""))


def _py_prose(path):
    """Every comment and docstring in a Python file, and no code.

    The surface the bad counts keep turning up in - in metaannot.py, and in
    this suite, which was not scanned at all until a count in this very file's
    docstrings turned out to be two rounds stale.
    """
    src = _text(path)
    out = [l for l in src.splitlines() if l.lstrip().startswith("#")]
    for n in ast.walk(ast.parse(src)):
        if isinstance(n, (ast.Module, ast.ClassDef, ast.FunctionDef,
                          ast.AsyncFunctionDef)):
            d = ast.get_docstring(n)
            if d:
                out.append(d)
    return "\n".join(out)


def _unreleased(text):
    """The CHANGELOG section this change set is allowed to rewrite.

    A shipped release's entries are a RECORD and are not edited to match a
    later code base, so the scan stops at the first released heading.
    """
    return text[text.index("## Unreleased"):text.index("\n## v0.5.0")]


def _tests_prose():
    """Every comment and docstring in the suite, as one surface.

    One surface and not one per file, so "six surfaces" stays a number about
    KINDS of place a count can live rather than a number that moves whenever
    somebody adds a test module.
    """
    d = os.path.join(ROOT, "tests")
    return "\n".join(_py_prose(os.path.join(d, f))
                     for f in sorted(os.listdir(d)) if f.endswith(".py"))


_EMITTED = {}


def _doctor_documents(ma):
    """Real `doctor --json` documents, built in process.

    THE DOCUMENT, not a list of the constants somebody remembered to name.
    The surface used to be five hand-listed constants, and it grew by one each
    time a reading found a number in a sentence that was not on the list -
    which is the same defect as the noun whitelist, one surface over: a
    hand-kept list of where to look cannot see the place nobody thought of.
    Every string in a real document is in the scan by construction, including
    the ones built per row out of a path, a stage list or a config value.

    Three configs, because a document only carries the rows its config
    reaches. The bare one is most of them; the DIAMOND one is the only way to
    reach `db:diamond:<tag>:usable`, whose caveat is where the record and
    defline caps are published; and `fragpipe_tmt` is the only way to reach
    the `tmt` block. Memoised, because building one costs a few seconds of
    probing PATH and this file scans the surfaces more than once.
    """
    if "docs" in _EMITTED:
        return _EMITTED["docs"]
    import fixtures as F
    out = []
    for name in ("bare", "diamond", "tmt"):
        root = tempfile.mkdtemp(prefix=f"ma_doc_{name}_")
        cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
        cfg["results_dir"] = os.path.join(root, "results")
        if name == "diamond":
            db = os.path.join(root, "db")
            os.makedirs(db)
            F.write_dmnd(os.path.join(db, "vfdb.dmnd"))
            with open(os.path.join(db, "vfdb.fasta"), "w",
                      encoding="utf-8") as fh:
                fh.write(">a one\nMKV\n>b two\nMKW\n")
            cfg["db"]["diamond"] = {"vfdb": os.path.join(db, "vfdb.dmnd")}
            cfg["run"]["diamond"] = True
        if name == "tmt":
            plex = os.path.join(root, "tmt", "TMT1")
            os.makedirs(plex)
            with open(os.path.join(plex, "ion.tsv"), "w",
                      encoding="utf-8") as fh:
                fh.write("Peptide\tCharge\t126\n")
            cfg["quant_table"] = os.path.join(root, "tmt")
            cfg["quant_format"] = "fragpipe_tmt"
        # PATH EMPTIED for every one of them, because a document's text
        # depends on what is installed and a surface that changes shape with
        # the machine would make this scan pass on a laptop and fail on the
        # server. The DIAMOND row is the one that proves it: the caveat
        # carrying the record and defline caps is written only when `diamond`
        # is NOT on PATH - with the binary there the row reads the database
        # header instead and says something with no numbers in it at all, so
        # on a machine with DIAMOND installed those two counts would simply
        # stop being anywhere, and the table would call them stale.
        old_path = os.environ.get("PATH", "")
        try:
            os.environ["PATH"] = ""
            p = ma.Paths(cfg)
            reqs = ma.requirements(cfg, p)
            out.append(ma.doctor(cfg, p, reqs,
                                 ma.doctor_checks(cfg, p, reqs, None)))
        finally:
            os.environ["PATH"] = old_path
    _EMITTED["docs"] = out
    return out


def _emitted_prose(ma):
    """Every string a real `doctor --json` document publishes, as one surface.

    Every string, and not the ones that look like sentences: a `detail` is
    where the numbers are, but a `caveat`, a `finding`, a `not_counted_reason`
    and a `verdict.rule` are all published prose too, and picking which keys
    to read would be the hand-kept list again.
    """
    out = []

    def walk(v):
        if isinstance(v, str):
            out.append(v)
        elif isinstance(v, dict):
            for x in v.values():
                walk(x)
        elif isinstance(v, list):
            for x in v:
                walk(x)

    walk(_doctor_documents(ma))
    # Two things this join has to get right, both learned from the document
    # itself. A string with NO WHITESPACE is an identifier and not prose - a
    # stage name, a check id, a version, a hostname, the `generated`
    # timestamp - and none of them can carry a claim about a count. And the
    # pieces are separated by something a count cannot be read across:
    # joined with a space, the `generated` timestamp ran straight into the
    # `host` beside it, and the seconds field followed by a hostname is a
    # number followed by a plural-looking word - so the scan reported a count
    # of the machine this suite happened to run on.
    return " || ".join(x for x in out if " " in x.strip())


def _count_surfaces(ma):
    """(name, prose) for each of the six surfaces a count can live in.

    The emitted document is one of them and is generated, not typed: every
    sentence in it that carries a number is built from a constant, and this is
    what says so rather than leaving it to be assumed.
    """
    return [
        ("metaannot.py comments and docstrings", _py_prose(METAANNOT_PY)),
        ("the emitted document", _emitted_prose(ma)),
        ("tests/ comments and docstrings", _tests_prose()),
        ("README.md", _text(README)),
        ("TUTORIAL.md", _text(TUTORIAL)),
        ("CHANGELOG.md (Unreleased)", _unreleased(_text(CHANGELOG))),
    ]


def _no_manifest_formats(ma):
    """How many quant formats have no reader that opens a manifest at all."""
    n = 0
    for fmt in sorted(ma.ALL_FORMATS):
        cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
        cfg["quant_format"] = fmt
        cfg["run"].update(join=True, unipept=True, taxonomy=True)
        if not ma._manifest_readers(cfg):
            n += 1
    return n


def _example_configs():
    d = os.path.join(ROOT, "examples", "server-run-plan")
    return sum(1 for _r, _ds, fs in os.walk(d)
               for f in fs if f.endswith((".yaml", ".yml")))


def _evidence_columns():
    """The columns the README's `peptide_evidence.tsv` sentence lists.

    A count of what the sentence itself enumerates, and each name checked
    against a literal in the engine. It cannot see a TENTH column appearing in
    `rollup_features()` that nobody documented - that is what
    `test_every_bin_is_documented` and its neighbours are for - but it does
    stop the number and the list from parting company, and stops the list
    naming a column the tool does not write.
    """
    m = re.search(r"`peptide_evidence\.tsv` has one row per protein and "
                  r"(\w+) columns:\s*(.+?)\. The last two",
                  _norm(_text(README)))
    assert m, "the README no longer describes peptide_evidence.tsv's columns"
    names = re.findall(r"`([a-z_]+)`", m.group(2))
    src = _text(METAANNOT_PY)
    for name in names:
        assert f'"{name}"' in src or f"'{name}'" in src, \
            f"the README lists a peptide_evidence column the tool never " \
            f"writes: {name}"
    return COUNT_WORDS[m.group(1).lower()], len(names)


def _scope_heading_states():
    """How many states the doctor suite's scope heading enumerates.

    The heading is `exists / right kind / not empty - the three states of one
    path`, and the number in it was classified DERIVED with a recompute of
    `3`. Counting the slash-separated items is what makes the classification
    mean something: the sentence and its own number can no longer part
    company, which is the whole point of the registry.
    """
    src = _text(os.path.join(TESTS_DIR, "test_doctor_json.py"))
    m = re.search(r"^# (.+?) - the \w+ states of one path$", src, re.M)
    assert m, "the doctor suite no longer carries that scope heading"
    return len([x for x in m.group(1).split("/") if x.strip()])


def _miss_cols_slice():
    """How many unmatched quant columns the `manifest:mapping` row prints.

    The README and the TUTORIAL both quote this number, and it is a literal
    slice in `_manifest_checks` - `miss_cols[:N]` - so it is derivable rather
    than a thing to remember.
    """
    src = _text(METAANNOT_PY)
    tree = ast.parse(src)
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "_manifest_checks")
    body = ast.get_source_segment(src, fn)
    # Scoped to the doctor row on purpose: stage_join's own WARN truncates the
    # same list at a different length, and the sentence in the docs is about
    # what `doctor` prints.
    m = re.search(r"miss_cols\[:(\d+)\]", body)
    assert m, "the mapping row no longer truncates its unmatched columns"
    return int(m.group(1))


def _describe_per_stage_keys():
    """The keys describe() emits for one stage, off the dict literal itself.

    The count that was wrong for the fourth consecutive round. Read from the
    source rather than from a document, because the sentence being checked is
    a comment sitting on that literal and the two have to be the same thing.
    """
    src = _text(METAANNOT_PY)
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "describe")
    for node in ast.walk(fn):
        if isinstance(node, ast.ListComp) and isinstance(node.elt, ast.Dict):
            return [k.value for k in node.elt.keys
                    if isinstance(k, ast.Constant)]
    raise AssertionError("describe() no longer builds its stage dicts inline")


def _suite_tuple_len(module, name):
    """How many entries a module-level tuple in a test module has.

    tests/ is a surface now, so a count in a test's own comment about a list
    in a test module is derivable exactly as one about an engine constant is.
    Read by AST rather than by import: importing a test module to count one
    tuple pulls its whole fixture stack in for nothing.
    """
    src = _text(os.path.join(ROOT, "tests", module))
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name
                for t in node.targets):
            return len(node.value.elts)
    raise AssertionError(f"tests/{module} no longer defines {name}")


def _fn_list_len(fn_name, var):
    """How many elements a list literal assigned inside a function has."""
    src = _text(METAANNOT_PY)
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == fn_name)
    for node in ast.walk(fn):
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == var
                        for t in node.targets)
                and isinstance(node.value, ast.List)):
            return len(node.value.elts)
    raise AssertionError(f"{fn_name}() no longer assigns a list to {var}")


# phrase -> (bucket, note, recompute or None)
COUNT_PROSE = {
    # ---- DERIVED: recomputed from the source, right here -------------
    "four checks": ("DERIVED", "the checks that read INTO a file",
                    lambda ma: len(ma.DOCTOR_DEEP_CHECKS)),
    "fifth entry": ("DERIVED", "the next DOCTOR_DEEP_CHECKS entry",
                    lambda ma: len(ma.DOCTOR_DEEP_CHECKS) + 1),
    "tenth value": ("DERIVED", "null, one past DOCTOR_FOUND_KINDS",
                    lambda ma: len(ma.DOCTOR_FOUND_KINDS) + 1),
    "twenty-one stages": ("DERIVED", "STAGES", lambda ma: len(ma.STAGES)),
    "eight formats": ("DERIVED", "ALL_FORMATS", lambda ma: len(ma.ALL_FORMATS)),
    "eight format": ("DERIVED", "ALL_FORMATS", lambda ma: len(ma.ALL_FORMATS)),
    "ninth format": ("DERIVED", "the next quant format",
                     lambda ma: len(ma.ALL_FORMATS) + 1),
    "four formats": ("DERIVED", "formats with no manifest reader",
                     _no_manifest_formats),
    "fourteen packages": ("DERIVED", "RNEED + ROPT",
                          lambda ma: len(ma.RNEED) + len(ma.ROPT)),
    "three values": ("DERIVED", "TMBED_GPU_MODES",
                     lambda ma: len(ma.TMBED_GPU_MODES)),
    "three command": ("DERIVED", "TMBED_GPU_MODES' command lines",
                      lambda ma: len(ma.TMBED_GPU_MODES)),
    "fourth mode": ("DERIVED", "the next TMBED_GPU_MODES entry",
                    lambda ma: len(ma.TMBED_GPU_MODES) + 1),
    "fourth value": ("DERIVED", "the next TMBED_GPU_MODES entry",
                     lambda ma: len(ma.TMBED_GPU_MODES) + 1),
    "two values": ("DERIVED", "DOCTOR_FAIL_REASONS",
                   lambda ma: len(ma.DOCTOR_FAIL_REASONS)),
    "five values": ("DERIVED", "ASSIGNMENT_CLASSES",
                    lambda ma: len(ma.ASSIGNMENT_CLASSES)),
    "eight configs": ("DERIVED", "examples/server-run-plan",
                      lambda ma: _example_configs()),
    "nine columns": ("DERIVED", "the columns the README's own sentence lists",
                     lambda ma: _evidence_columns()[1]),
    # The self-correction groups are the ones whose headings are about what
    # the document CLAIMS - "Six claims the first draft made about ITSELF",
    # "Six claims the third reading corrected". The attribution group and the
    # engine-defect group are each their own noun in the same sentence and are
    # counted beside this one.
    "twenty-eight corrections": ("DERIVED", "the groups that correct the document",
                           lambda ma: sum(
                               n for h, n in _unreleased_fixed_groups()
                               if "claims" in h.lower())),

    "four columns": ("DERIVED", "the miss_cols slice in _manifest_checks",
                     lambda ma: _miss_cols_slice()),

    # ---- DERIVED, added by the FIFTH reading, when the scan grew -----
    # Digits, tests/ and the wider emitted document all arrived together, and
    # these are the counts that only became visible because of it.
    "ten keys": ("DERIVED", "the per-stage dict describe() emits",
                 lambda ma: len(_describe_per_stage_keys())),
    "200,000 records": ("DERIVED", "DMND_FASTA_RECORD_CAP",
                        lambda ma: ma.DMND_FASTA_RECORD_CAP),
    "200 deflines": ("DERIVED", "DMND_DEFLINE_CAP",
                     lambda ma: ma.DMND_DEFLINE_CAP),
    "six surfaces": ("DERIVED", "the surfaces this scan reads",
                     lambda ma: len(_count_surfaces(ma))),
    "eleven configs": ("DERIVED", "CONFIGS in tests/test_doctor_json.py",
                       lambda ma: _suite_tuple_len("test_doctor_json.py",
                                                   "CONFIGS")),
    "four rows": ("DERIVED", "UNREADABLE_FALSE_PASSES",
                  lambda ma: _suite_tuple_len("test_doctor_json.py",
                                              "UNREADABLE_FALSE_PASSES")),
    "four keys": ("DERIVED", "UNREADABLE_FALSE_PASSES",
                  lambda ma: _suite_tuple_len("test_doctor_json.py",
                                              "UNREADABLE_FALSE_PASSES")),
    "21 stages": ("DERIVED", "STAGES", lambda ma: len(ma.STAGES)),
    "9 field": ("DERIVED", "the custom --outfmt 6 columns parse_diamond reads",
                lambda ma: _fn_list_len("parse_diamond", "cols")),
    # C4: this was ("DERIVED", ..., lambda ma: 3) - a DERIVED entry whose
    # recompute was the literal it was checking, so it could only ever agree
    # with itself. That is the "rule with nothing under it" the staleness half
    # of this test exists to prevent, one level up: a reader meets DERIVED and
    # believes the number is pinned to something. It is pinned now, to the
    # thing the sentence itself enumerates - the slash-separated states in the
    # scope heading it sits in - which is the same shape as "nine columns"
    # above and fails if somebody adds a fourth state to the list and leaves
    # the word alone.
    "three states": ("DERIVED", "the states the scope heading itself lists",
                     lambda ma: _scope_heading_states()),

    # ---- MEASURED added by the fifth reading --------------------------
    # A fact about a dataset, a host, someone else's file format, or an
    # observation of what a tool or a browser really did.
    # Not recomputable; the pin is that the SENTENCE STILL EXISTS.
    "seven paths": ("MEASURED", "the paths `run` was measured hanging on, "
                                "one FIFO at a time"),
    # Was DERIVED, as "one input key through every state", and the sweep has
    # grown two states since - so the recompute moved to "nine cells" above
    # and this phrase is left pointing at what it always meant in the sentence
    # it appears in: the cells where `doctor` and `run` actually disagreed,
    # counted once, on a machine, at a time.
    "seven cells": ("MEASURED", "the FIFO cells where the differential sweep "
                                "found doctor and run disagreeing"),
    "1 row": ("MEASURED", "how often a structure column is non-empty, "
                          "observed on the real run"),
    "1,237,468 values": ("MEASURED", "the real 8-plex run's matrix"),
    "6 groups": ("MEASURED", "the real run's design"),
    "4 stages": ("MEASURED", "stage_workers in a worked sizing example"),
    "5 databases": ("MEASURED", "the DIAMOND databases of the real run"),
    "262 entries": ("MEASURED", "BAGEL4's motif seed set"),
    "12 columns": ("MEASURED", "DIAMOND's own default --outfmt 6"),
    "10 field": ("MEASURED", "a legacy foldseek row's width, observed"),
    "10 column": ("MEASURED", "a legacy row's width before a column was "
                              "appended, observed"),
    "11 column": ("MEASURED", "the same row with target_db appended"),
    "tenth column": ("PROSE", "'a TENTH column appearing in "
                              "rollup_features()' - a hypothetical about a "
                              "column nobody has documented"),
    "11 field": ("MEASURED", "a foldseek row's width with target_db appended"),
    "ten fields": ("MEASURED", "the foldseek output fields both lists ask "
                               "for, observed against the binary"),
    "18 columns": ("MEASURED", "hmmsearch --tblout's fixed columns"),
    "23 column": ("MEASURED", "hmmsearch --domtblout's fixed columns"),
    "14 column": ("MEASURED", "InterProScan's TSV with -iprlookup -goterms"),
    "9 column": ("MEASURED", "the DIAMOND custom format, in a fixture's "
                             "docstring and in this file's own example"),
    "30 columns": ("MEASURED", "the hhsearch hit table's name column"),
    "46 columns": ("MEASURED", "the annotation columns one heuristic swept "
                               "into an assay, observed"),
    "1000 rows": ("MEASURED", "readr's type-guess window"),
    "3 values": ("MEASURED", "limma's minimum valid values per group"),
    "21 rows": ("MEASURED", "the console rows one distinction governs"),
    "ten rows": ("MEASURED", "the console rows a caveat used to sit on"),
    "thirteen rows": ("MEASURED", "the console rows one reproduction showed"),
    "eight rows": ("MEASURED", "the console rows that all read 'results'"),
    "five columns": ("MEASURED", "what a browser did to a table whose widths "
                                 "were style attributes"),
    "four stages": ("MEASURED", "the auditor's reproduction of the queue"),
    "nine stages": ("MEASURED", "the stages the README's worked example "
                                "turns on"),
    "nine fields": ("MEASURED", "the documented describe fields one commit "
                                "deleted, listed beside the number"),

    # ---- PROSE added by the fifth reading -----------------------------
    "0 row": ("PROSE", "'a 0-row matrix' - an empty one, not a count"),
    "3 file": ("PROSE", "'the chunk 3 whose file is on disk' - an index"),
    "1 stage": ("PROSE", "a quotation of a rendered console line"),
    "3 stages": ("PROSE", "'`cost-3 stages` is not \"3 stages\"' - this "
                          "scanner's own worked example of what it refuses"),
    "forty files": ("PROSE", "'not counts of four states or forty files' - "
                             "this scanner's own example of a unit"),
    "768 config": ("PROSE", "the sweep that never existed, quoted as the "
                            "record of the defect and asserted nowhere; "
                            "test_the_changelog_names_the_config_sweep_this_"
                            "suite_really_runs pins that it stays a quote"),
    "eight fields": ("PROSE", "a quotation of the fourth wrong count, kept "
                              "as the record of the defect"),
    "seven fields": ("PROSE", "the per-stage fields before `cost` - a "
                              "history, and the rule it illustrates"),
    "eight entries": ("PROSE", "a history: what the first wave printed under "
                               "a paragraph that said seven"),
    "twenty-one entries": ("PROSE", "a history: the bolded entries at the "
                                    "time, quoted"),
    "fourteen verdicts": ("PROSE", "a history: the false verdicts at the "
                                   "time, quoted"),
    "thirty-three entries": ("PROSE", "a quotation of the sentence that "
                                      "accounted for thirty-two of them"),
    "14 entries": ("PROSE", "a history: the non-verdict entries the section "
                            "held when the sentence was written"),
    "sixteen stages": ("PROSE", "a quotation of what a shipped CHANGELOG "
                                "entry used to say"),
    "seven entries": ("PROSE", "'was read as \"seven entries\"' - the "
                               "misreading this alternation exists to stop"),
    "twenty-seven entries": ("PROSE", "the phrase that was misread, quoted"),
    "fifth status": ("PROSE", "'a fifth status cannot appear' - a "
                              "hypothetical about a set that has not grown"),
    "first flag": ("PROSE", "'bounded by the first flag' - a position"),
    "second key": ("PROSE", "'a SECOND key the row never mentioned' - an "
                            "ordinal in a narrative, not a cardinality"),
    "three files": ("PROSE", "a named trio, listed in the same sentence"),
    "three reasons": ("PROSE", "'for three STRUCTURAL reasons', listed "
                               "immediately and numbered"),
    "two reasons": ("PROSE", "'for two reasons', and both follow"),
    "three rows": ("PROSE", "'three rows in the engine also had' - the "
                            "sentences one round corrected, named"),
    "two stage": ("PROSE", "'the two stage lists' - `blocks` and `degrades`, "
                           "a named pair"),
    "two file": ("PROSE", "'a two-file commit' - what a pin would cost"),
    "two files": ("PROSE", "'Two files on purpose' - a named pair"),
    "two flags": ("PROSE", "a named pair of booleans"),
    "two group": ("PROSE", "'the difference of the two group means' - a "
                           "contrast between two named groups"),
    "two keys": ("PROSE", "'the LAST of two identical keys' - a named pair"),
    "two paths": ("PROSE", "'two different paths' - a named pair"),
    "two record": ("PROSE", "'the last two record how' - `record` is the "
                            "verb here"),
    "two records": ("PROSE", "'writes two records for every protein' - an "
                             "arity, not a total"),
    "two rows": ("PROSE", "'the two usability rows' - a named pair"),

    # ---- MEASURED added by the fourth reading -------------------------
    "three columns": ("MEASURED", "the eggNOG 2.0.x columns spelled "
                                  "differently from 2.1.x, observed"),
    "twelve stages": ("MEASURED", "a worked example of the scheduler queue"),

    # ---- PROSE added by the fourth reading ----------------------------
    "two verdicts": ("PROSE", "'the difference is two false verdicts' - what "
                              "one wrong helper cost, named"),
    "two column": ("PROSE", "'two different column names' - a named pair"),
    "four check": ("PROSE", "'it used to probe four and check neither' - a "
                            "history, and the four is pinned as four packages"),
    "first entries": ("PROSE", "'a table whose first three entries' - a "
                               "position in a runbook's own timing table, not "
                               "a count of anything this codebase has"),
    "seven verdicts": ("PROSE", "a group heading, counted by the group test"),
    "five false verdicts": ("PROSE", "a group heading, counted by the group "
                                     "test"),
    # ---- MEASURED: a fact about a dataset, a host or an observation ---
    "five databases": ("MEASURED", "the DIAMOND databases of the real run"),
    "four tiers": ("MEASURED", "the real search database's prefixes"),
    "three key spaces": ("MEASURED", "the real search database's namespaces"),
    "two tiers": ("MEASURED", "two prefixes of one namespace, measured"),
    "six groups": ("MEASURED", "the real run's design"),
    "three stages": ("MEASURED", "the timing table's concurrency, and "
                      "'adopted three stages later', a distance"),

    # ---- PROSE: not a cardinality of anything this codebase has -------
    "first value": ("PROSE", "'the first value was called stage_dies' - a "
                             "history of one rename, not a count"),
    "first wave": ("PROSE", "the scheduler's first dispatch, not a count"),
    "first verdict": ("PROSE", "'the first verdict wins' - precedence"),
    "first column": ("PROSE", "a column position in someone else's file"),
    "second command": ("PROSE", "'the second command above' - a reference"),
    "third state": ("PROSE", "'a directory is a third state again' - one more"),
    "third kind": ("PROSE", "'a third kind of failure' - a hypothetical"),
    "third wave": ("PROSE", "'the third wave' names a heading, which is "
                            "counted by the group test"),
    "sixth kind": ("PROSE", "'or a sixth kind to the other' - a hypothetical "
                            "about a set that has not grown"),
    "seven values": ("PROSE", "a quotation of what the README USED to say"),
    "three checks": ("PROSE", "a quotation of what the document USED to say"),
    "eighth value": ("PROSE", "a quotation of the wrong count, kept as the "
                              "record of the defect"),
    "two checks": ("PROSE", "'two checks about one file' - a pair, named"),
    "two stages": ("PROSE", "a named pair - esmfold and tmbed, or the two "
                            "that fall back"),
    "two columns": ("PROSE", "a named pair of columns being compared"),
    "two commands": ("PROSE", "'the two commands above' - a reference"),
    "two database": ("PROSE", "'two database downloads' - a named pair"),
    "two states": ("PROSE", "'two states of one file' - a named pair"),
    "two depths": ("PROSE", "'anything comparing two depths' - an arity"),
    "three depths": ("PROSE", "'the other three depths', immediately listed"),
    "two sections": ("PROSE", "'two sections later' - a distance"),
    "three verdicts": ("PROSE", "'THREE verdicts on one key', listed"),
    "four waves": ("PROSE", "counted by the false-verdict group test"),
    "fourth wave": ("PROSE", "'the fourth wave' names a group, which the "
                             "group test counts"),
    "two waves": ("PROSE", "'two waves ago' - a distance in this file"),
    "twelve corrections": ("PROSE", "a quotation of the sentence that left "
                                    "the engine defect out"),
    "twenty-four false verdicts": (
        "DERIVED", "the entries under the false-verdict groups",
        lambda ma: sum(n for h, n in _unreleased_fixed_groups()
                       if "false verdict" in h.lower())),
    "four states": ("MEASURED", "the four path states measured taking the "
                                "whole document down"),
    "four packages": ("PROSE", "a quotation of what the TUTORIAL used to say"),
    "ten columns": ("PROSE", "'asked for the same ten columns again' - a "
                             "foldseek output format, quoted"),
    "twelve tier": ("PROSE", "'the twelve-tier budget' - a cap, not a count "
                             "of things that exist"),
    "two tier": ("PROSE", "'under two tier prefixes' - a named pair"),
    "four package": ("PROSE", "'the old four-package probe' - the record of "
                              "what it used to be"),
    # ---- DERIVED, added by the SIXTH reading, when the whitelist went ----
    # The negative-space rule surfaced a pile of counts that had never been
    # visible. These are the ones that are recomputable, and several of them
    # are the sweeps' own sizes, which is what the two wrong counts this round
    # fixed were about.
    "sixteen modules": ("DERIVED", "the test modules in tests/",
                        lambda ma: len([f for f in os.listdir(
                            os.path.join(ROOT, "tests"))
                            if f.startswith("test_") and f.endswith(".py")])),
    "nine states": ("DERIVED", "PATH_STATES in tests/test_doctor_json.py",
                     lambda ma: _suite_tuple_len("test_doctor_json.py",
                                                 "PATH_STATES")),
    "fifty-four cells": ("DERIVED", "the input sweep: every key in every state",
                        lambda ma: (_suite_tuple_len("test_doctor_json.py",
                                                     "INPUT_PATH_KEYS")
                                    * _suite_tuple_len("test_doctor_json.py",
                                                       "PATH_STATES"))),
    "81 cases": ("DERIVED", "the whole sweep: the input keys and the TMT tree",
                 lambda ma: ((_suite_tuple_len("test_doctor_json.py",
                                               "INPUT_PATH_KEYS")
                              + _suite_tuple_len("test_doctor_json.py",
                                                 "TMT_READ_PATHS"))
                             * _suite_tuple_len("test_doctor_json.py",
                                                "PATH_STATES"))),
    "nine cells": ("DERIVED", "one input key through every state",
                    lambda ma: _suite_tuple_len("test_doctor_json.py",
                                                "PATH_STATES")),
    "104 entries": ("DERIVED", "the Unreleased/Fixed groups",
                            lambda ma: sum(
                                n for _h, n in _unreleased_fixed_groups())),
    "fifteen groups": ("DERIVED", "the Unreleased/Fixed groups",
                      lambda ma: len(_unreleased_fixed_groups())),
    "eighty entries": ("DERIVED", "the entries that are NOT false verdicts",
                           lambda ma: sum(
                               n for h, n in _unreleased_fixed_groups()
                               if "false verdict" not in h.lower())),
    "three configs": ("DERIVED", "the documents _emitted_prose() is built from",
                      lambda ma: len(_doctor_documents(ma))),
    "seven vocabularies": ("PROSE", "'five of the seven closed vocabularies' "
                                    "- a history: what there were before "
                                    "`section` turned out to be the eighth"),

    # ---- MEASURED, added by the sixth reading ---------------------------
    # Facts about a dataset, a run, a machine or somebody else's file format.
    # Not recomputable; the pin is that the SENTENCE STILL EXISTS, so an edit
    # cannot strand a measurement whose subject has gone.
    "30 rows": ("MEASURED", "what a taxonomy run over a CHAR-DEVICE manifest "
                            "really returned - the measurement that proved "
                            "the manifest row's 'every stage goes with it' "
                            "false and its `blocks: [join]` right"),
    "141 proteins": ("MEASURED", "the real run's `3s_structure_only` bin"),
    "604 proteins": ("MEASURED", "the first real run's effector shortlist"),
    "40 proteins": ("MEASURED", "a worked example of the min_plexes message"),
    "5 proteins": ("MEASURED", "the ~10^5 proteins a real database holds"),
    "10 proteins": ("MEASURED", "the fixture taxon, 7 of its 10 in one plex"),
    "10 members": ("MEASURED", "the same fixture taxon, in the README"),
    "thirty proteins": ("MEASURED", "the TMT fixture's protein set"),
    "twelve proteins": ("MEASURED", "the taxonomy fixture's T1 taxon"),
    "fifty proteins": ("MEASURED", "the dark-fold example in the README"),
    "fifty singletons": ("MEASURED", "the same example's comparison"),
    "30 features": ("MEASURED", "the assignment sweep's fixture"),
    "four features": ("MEASURED", "the four kinds of feature one fixture has"),
    "four peptides": ("MEASURED", "a finding's worked example"),
    "three peptides": ("MEASURED", "the TMT fixture, per protein"),
    "36 samples": ("MEASURED", "a worked example of the centring problem"),
    "75 samples": ("MEASURED", "the second real dataset"),
    "four samples": ("MEASURED", "the TMT fixture's design"),
    "8 plexes": ("MEASURED", "the real 8-plex TMT run"),
    "six plexes": ("MEASURED", "where the pool sat in the real 8-plex run"),
    "88 channels": ("MEASURED", "the same run, 11 channels x 8 plexes"),
    "16 channels": ("MEASURED", "a TMTpro plex, in the mixed-sizes note"),
    "16 channel": ("MEASURED", "the annotation a 16-channel plex writes"),
    "four channels": ("MEASURED", "the TMT fixture's plex"),
    "three channels": ("MEASURED", "the one-plex fixture"),
    "262 sequences": ("MEASURED", "BAGEL4's motif seed set"),
    "91 sequences": ("MEASURED", "the ESMFold run's long tail"),
    "7 sequences": ("MEASURED", "the sequences one OOM killed a stage over"),
    "5 taxids": ("MEASURED", "eggNOG 5's seed taxids, from a 2018 taxonomy"),
    "965 hits": ("MEASURED", "the largest VFDB category on a real gut set"),
    "22 cores": ("MEASURED", "the workstation the timing table was run on"),
    "4 threads": ("MEASURED", "DIAMOND's sublinear thread scaling, worked"),
    "4 jobs": ("MEASURED", "the same worked example"),
    "four jobs": ("MEASURED", "the README's peak-memory example"),
    "90 chunks": ("MEASURED", "a resumed run's chunk count"),
    "33 requests": ("MEASURED", "the console's concurrency reproduction"),
    "30 attempts": ("MEASURED", "the SIGTERM race, 8 of 30"),
    "30 tries": ("MEASURED", "the same race, the other sentence about it"),
    "256 parts": ("MEASURED", "the chunk planner's overshoot, worked"),
    "4 points": ("MEASURED", "what a zeroed vfdb weight still scored"),
    "10 directories": ("MEASURED", "the banner one console reproduction read"),
    "eleven directories": ("MEASURED", "the fixture tree behind that banner"),
    "eight directories": ("MEASURED", "the directories one cache bug hid"),
    "fifteen directories": ("MEASURED", "what `describe` used to leave behind"),
    "eight datasets": ("MEASURED", "the worked example's real plan"),
    "eight functions": ("MEASURED", "the Windows-skipped signal tests"),
    "eight items": ("MEASURED", "what those eight functions collect as"),
    "five tests": ("MEASURED", "the tests the open xfail markers sit on"),
    "35 tests": ("MEASURED", "the R selection, counted exactly by "
                             "test_the_readme_test_counts_are_the_counts_"
                             "this_suite_really_has"),
    "seven strings": ("MEASURED", "the class names three outputs share"),
    "seven markers": ("MEASURED", "the console's status markers"),
    "four endpoints": ("MEASURED", "the console's JSON routes"),
    "fourteen routes": ("MEASURED", "the console's routes, in a test's note"),
    "six predicates": ("MEASURED", "the predicates one fixture combines"),
    "three genus": ("MEASURED", "the taxonomy fixture's tree"),
    "21 claims": ("MEASURED", "what the console page asserted about a "
                              "directory it knew nothing about"),
    "500 characters": ("MEASURED", "the stage-error cap"),
    "500 character": ("MEASURED", "the same cap, attributively"),
    "eight failures": ("MEASURED", "what fits on the console page at that cap"),
    "five failures": ("MEASURED", "what pushed the vitals off screen at "
                                  "1440x860"),
    "three failures": ("MEASURED", "the failures one `.dead` block held"),
    "five misses": ("MEASURED", "the heartbeat reproduction, 5 x 0.5s"),
    "three days": ("MEASURED", "how long a first full run takes"),
    "3 days": ("MEASURED", "the same estimate, in the runbook's table"),
    "three day": ("MEASURED", "the run a reversed kill order would hit"),
    "six months": ("MEASURED", "when somebody opens the directory again"),
    "three tools": ("MEASURED", "the concurrency the timing table really used"),
    "three orders": ("MEASURED", "the fixture's plex effect, in log2"),
    "four orders": ("MEASURED", "the size-factor example, in log2"),
    "2 orders": ("MEASURED", "how much slower tmbed is on a CPU"),
    "four ratios": ("MEASURED", "what a median can rest on at the threshold"),
    "three spaces": ("MEASURED", "the real search database's namespaces"),
    "4 mer": ("MEASURED", "the LP.TG motif's length"),
    "3 letter": ("MEASURED", "the `Xre` substring that matched everything"),
    "7 digit": ("MEASURED", "the accession widths one database uses"),
    "six character": ("MEASURED", "the mark the human report prints"),
    "three levels": ("MEASURED", "how deep `--root` scans"),
    "eleven configs": ("DERIVED", "CONFIGS in tests/test_doctor_json.py",
                       lambda ma: _suite_tuple_len("test_doctor_json.py",
                                                   "CONFIGS")),

    # ---- PROSE, added by the sixth reading ------------------------------
    # Not a cardinality of anything this codebase has: a pair named in the
    # sentence, a position, a history quoted as the record of a defect, or
    # this scanner's own worked example of something it refuses.
    "143 matches": ("PROSE", "an exit status, not a count of anything"),
    "2 others": ("PROSE", "'minus 2 if others share the machine' - a margin"),
    "2 settings": ("PROSE", "'your phase 2 settings' - a phase, not a count"),
    "50 these": ("PROSE", "'50 for these three because' - a threshold"),
    "25 them": ("PROSE", "'25 of them is hours between log lines' - a rate"),
    "264 them": ("PROSE", "'264 of them' - the parts the planner overshot to"),
    "993 them": ("PROSE", "the models that passed the pLDDT gate, measured "
                          "in the same sentence as the 1,821 built"),
    "eight them": ("PROSE", "'All eight of them' - the refused matrices, "
                            "listed immediately"),
    "eight these": ("PROSE", "'a real run writes eight of these' - the "
                             "refused files, named in the same comment"),
    "five them": ("PROSE", "'five of them went unguarded' - a history: the "
                           "vocabularies before `section` became the eighth"),
    "five those": ("PROSE", "'Five of those seven markers are new' - a "
                            "subset of a count pinned beside it"),
    "seven them": ("PROSE", "'all SEVEN of them have now been measured' - "
                            "the seven paths, pinned as `seven paths`"),
    "three these": ("PROSE", "'All three of these were real' - a named trio"),
    "four these": ("PROSE", "'false for four of these twelve' - a subset"),
    "thirty them": ("PROSE", "'all thirty of them' - the sentence this rule "
                             "was written for, quoted in the comment that "
                             "explains why a pronoun is a head"),
    "three ones": ("PROSE", "'the three security ones' - a named trio"),
    "five ones": ("PROSE", "'like the five stock ones' - the stock DIAMOND "
                           "databases, pinned as `five databases`"),
    "six others": ("PROSE", "'six others exist' - a quotation of the banner "
                            "that counted without naming"),
    "3 aminomutase": ("PROSE", "'3-aminomutase' - part of an identifier"),
    "3 tuple": ("PROSE", "'a 3-tuple' - an arity"),
    "four step": ("PROSE", "'a four-step runbook' - the steps are listed"),
    "three way": ("PROSE", "'the three-way split' - the branches follow"),
    "768 way": ("PROSE", "'the 768-way product this sentence imagined' - the "
                         "sweep that never existed, quoted"),
    "four constants": ("PROSE", "a history: what `_emitted_prose()` read "
                                "before it was built from a document"),
    "four readings": ("PROSE", "a history: the readings before this one"),
    "four rounds": ("PROSE", "'four verification rounds' - a history"),
    "three rounds": ("PROSE", "'for three rounds' - a history"),
    "six reproductions": ("PROSE", "a heading over the reproductions it lists"),
    "six claims": ("PROSE", "a group heading, counted by the group test"),
    "ten claims": ("PROSE", "a group heading, counted by the group test"),
    "ten defects": ("PROSE", "a group heading, counted by the group test"),
    "twelve defects": ("PROSE", "a group heading, counted by the group test"),
    "seventeen defects": ("PROSE", "a group heading, counted by the group "
                                   "test"),
    "seven defects": ("PROSE", "a group heading, counted by the group test"),
    "six tests": ("PROSE", "'five of the six unmarked tests were marked' - a "
                           "history of one count reaching one"),
    "490 tests": ("PROSE", "a history: how far the README's numbers had "
                           "drifted, quoted"),
    "thirteen modules": ("PROSE", "the wrong count, kept as the record of "
                                  "the defect"),
    "six states": ("PROSE", "a quotation of what the sweep entry used to say"),
    "thirty cells": ("PROSE", "the same quotation's other half"),
    "twelve cells": ("PROSE", "the DIFFERENCE between two counts in one "
                              "sentence, both of which are pinned above"),
    "four cells": ("PROSE", "a history: the cells the differential sweep "
                            "turned up before the other three were measured"),
    "seven hangs": ("PROSE", "'all seven hangs' - the seven paths, pinned "
                             "as `seven paths`"),
    "four spurious": ("PROSE", "a history: what one duplicated row cost"),
    "three previous rounds": ("PROSE", "a history, in the CHANGELOG"),
    "three assertions": ("PROSE", "a history: what a shadowed table broke"),
    "three cases": ("PROSE", "'the three cases say three different things' - "
                             "the cases are named in the sentence"),
    "three sources": ("PROSE", "'Three sources, two answers' - a named trio"),
    "five instances": ("PROSE", "'The class, not the five instances' - the "
                                "same rule, in the helper it names"),
    "four resources": ("PROSE", "'a count of four resources' - this "
                                "scanner's own example of what a unit does"),
    "four callers": ("PROSE", "the callers of one helper, listed beside it"),
    "four copies": ("PROSE", "'one row rather than four copies' - the "
                             "copies that would have been"),
    "four hmmpress": ("PROSE", "'the four hmmpress siblings' - the files "
                               "hmmpress writes, pinned as `four siblings`"),
    "four siblings": ("MEASURED", "what hmmpress writes beside a library"),
    "four items": ("PROSE", "'four collected items between them' - the "
                            "tests one tree held, named"),
    "four places": ("PROSE", "'written by hand in four places' - a history"),
    "four sentences": ("PROSE", "'four sentences added by this change set' - "
                                "a history"),
    "four shapes": ("PROSE", "'the four falsy shapes' - listed in the test"),
    "four skips": ("PROSE", "'turns three names into four skips' - an "
                            "arithmetic worked in the sentence"),
    "four slots": ("PROSE", "'compete for four slots' - the default "
                            "stage_workers, quoted as an example"),
    "four statements": ("PROSE", "'three of the four statements of the "
                                 "verdict rule' - the statements are listed"),
    "four terms": ("PROSE", "'an OR of four terms' - the terms are listed"),
    "four tests": ("PROSE", "'all four describe tests' - a named set"),
    "four things": ("PROSE", "'an OR of four things' - listed in the doc"),
    "three answers": ("PROSE", "'Three answers, not two' - each is given"),
    "three arms": ("PROSE", "'exactly three arms' - the arms are the code "
                            "below the sentence"),
    "three behaviours": ("PROSE", "'Three behaviours on one key' - listed"),
    "three branches": ("PROSE", "'the three branches that used to write' - a "
                                "history"),
    "three buckets": ("PROSE", "'the three buckets it happened to name' - a "
                               "history of one sentence"),
    "three causes": ("PROSE", "'The three causes are not the same' - listed"),
    "three claims": ("PROSE", "'Three claims, and each of them' - listed"),
    "three classes": ("PROSE", "'the three utility classes' - a section "
                               "heading over the three it holds"),
    "three conditions": ("PROSE", "'as three separate conditions' - the "
                                  "conditions are written out"),
    "three copies": ("PROSE", "'three copies of it had already' - a history"),
    "three counts": ("PROSE", "'let three other counts drift' - a history"),
    "three decisions": ("PROSE", "'three decisions matter more' - listed"),
    "three facts": ("PROSE", "'the three urgent facts' - named in the test"),
    "three functions": ("PROSE", "'a NoneType traceback three functions "
                                 "later' - a distance"),
    "three instances": ("PROSE", "'the STRUCTURE and not the three "
                                 "instances' - the same rule this scanner is"),
    "three limits": ("PROSE", "'Three limits are worth knowing' - listed"),
    "three lists": ("PROSE", "'Three different lists come' - named"),
    "three outcomes": ("PROSE", "'Three outcomes:' - a table follows"),
    "three places": ("PROSE", "'in this list in three places' - a history"),
    "three probes": ("PROSE", "'Three probes, because no one' - listed"),
    "three refusals": ("PROSE", "'the three refusals' - the console's, named "
                                "in CLAUDE.md and in the test"),
    "three renderings": ("PROSE", "'three renderings of one rule' - named"),
    "three sentences": ("PROSE", "'Three published sentences said' - the "
                                 "history this verb's docstring records"),
    "three settings": ("PROSE", "'the three settings cmd_run validates' - "
                                "listed immediately"),
    "three tests": ("PROSE", "'the three tests beside it' - named"),
    "three things": ("PROSE", "'the same three things as ion.tsv' - named"),
    "three ways": ("PROSE", "'The verdict, three ways' - each is written"),
}


def test_every_count_in_prose_is_derived_from_the_thing_it_counts_or_pinned():
    """The class, not the five instances. See the block comment above.

    Two halves, and they do different jobs. The NEGATIVE SPACE half reads the
    six surfaces and refuses any "<number> <plural noun>" the table does not
    classify - that is the half that would have caught "thirteen test modules"
    and "all thirty of them", and it is the half that is new. The TABLE half
    goes the other way: it looks for each classified phrase where it is
    allowed to be written, recomputes the DERIVED ones, and fails on a key
    nothing in the tree says any more. That half is what keeps a singular or
    ordinal pin - "a ninth format", "a 9-field row" - alive, since the shape
    rule above would not have found either.
    """
    import metaannot as _unused                              # noqa: F401
    ma = _load_ma()
    surfaces = [(name, _norm(text).replace("**", ""))
                # Markdown emphasis stripped: `**twenty-four** live false
                # verdicts` put two asterisks between the number and the
                # thing, and the scan read straight past the one count in this
                # file that had already been wrong once.
                for name, text in _count_surfaces(ma)]

    unclassified = []
    for name, flat in surfaces:
        for m in _COUNT_RE.finditer(flat):
            phrase = _count_phrase(m)
            if phrase is None or phrase in COUNT_PROSE:
                continue
            unclassified.append(
                f"{name}: {phrase!r} is a count of something this codebase "
                f"has and nothing pins it. Context: "
                f"...{flat[max(0, m.start() - 90):m.end() + 60]}...")
    assert not unclassified, "\n".join(unclassified)

    wrong, stale = [], []
    for phrase, entry in sorted(COUNT_PROSE.items()):
        rx = _phrase_re(phrase)
        seen = [(name, m) for name, flat in surfaces for m in rx.finditer(flat)]
        if not seen:
            stale.append(phrase)
            continue
        if entry[0] != "DERIVED":
            continue
        want, got = entry[2](ma), _count_value(phrase.split(" ", 1)[0])
        if got != want:
            name, m = seen[0]
            flat = dict(surfaces)[name]
            wrong.append(f"{name}: {phrase!r} says {got}, but {entry[1]} is "
                         f"{want}. Context: "
                         f"...{flat[max(0, m.start() - 90):m.end() + 60]}")
    assert not wrong, "\n".join(wrong)
    # A registry entry for a sentence nobody writes any more is a rule with
    # nothing under it, and the next reader believes it is still enforced.
    assert not stale, \
        f"COUNT_PROSE classifies phrases that appear in none of the six " \
        f"surfaces any more: {stale}"


def _load_ma():
    import importlib.util
    spec = importlib.util.spec_from_file_location("_ma_docs", METAANNOT_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_readme_peptide_evidence_column_count_is_the_list_beside_it():
    # The narrow half of the rule above, spelled out: the number and the list
    # in one sentence cannot part company, and no name in the list is one the
    # tool never writes.
    said, listed = _evidence_columns()
    assert said == listed, \
        f"the README says {said} columns and lists {listed}"


def test_the_changelog_non_verdict_entries_add_up_to_the_section(ma):
    """symptom: the sentence accounted for thirty-two of thirty-three entries.

    "one attribution that was never a verdict and twelve corrections to what
    the document says about ITSELF" came to 13 against a section that then
    held 14 entries that are not false verdicts - the fourteenth was the
    ENGINE defect a verdict uncovered, which changes no config's `doctor`
    verdict either and was in no bucket the sentence named.

    BOTH NUMBERS HAVE MOVED SINCE, and this docstring said them in the present
    tense for two rounds after they stopped being true - which is exactly the
    class the scanner above is for, and exactly why tests/ is one of its
    surfaces now. They are written here as the history they are; what the
    section holds TODAY is not written down anywhere, because this test sums
    EVERY clause of the sentence against the groups rather than checking the
    three buckets it happened to name, so a bucket left out of the sentence
    fails here whatever the totals have grown to.
    """
    groups = _unreleased_fixed_groups()
    total = sum(n for _h, n in groups)
    verdicts = sum(n for h, n in groups if "false verdict" in h.lower())
    txt = _norm(_text(CHANGELOG))
    m = re.search(r"That section also carries (.+?) — ([\w-]+) entries that "
                  r"change no config", txt)
    assert m, "the changelog no longer says what the non-verdict entries are"
    # Every clause of the sentence, summed, rather than three named ones: the
    # defect this pins is a bucket left OUT of it, so the test may not know in
    # advance how many buckets there are.
    named = sum(NUMBER_WORDS[w.lower()]
                for w in re.findall(r"\b([\w-]+)\b", m.group(1))
                if w.lower() in NUMBER_WORDS)
    assert NUMBER_WORDS[m.group(2).lower()] == named, \
        f"the sentence lists {named} and then totals {m.group(2)}"
    assert named == total - verdicts, \
        f"the sentence accounts for {named} non-verdict entries; the section " \
        f"has {total - verdicts}"


def test_the_changelog_names_the_config_sweep_this_suite_really_runs(ma):
    # symptom: "The 768-config sweep could not have caught this". There is no
    # 768-config sweep: `CONFIGS` in tests/test_doctor_json.py is the list the
    # doctor tests parametrise over, and it is a named handful.
    txt = _norm(_text(CHANGELOG))
    assert "There is no 768-config sweep." in txt, \
        "the entry that corrects the number is the only place it may appear"
    assert txt.count("768-config sweep") == 2, \
        "768 is quoted twice, in the entry that retires it, and asserted " \
        "nowhere"
    assert "`CONFIGS` in `tests/test_doctor_json.py`" in txt, \
        "the changelog should name the sweep rather than count it"
    src = _text(os.path.join(ROOT, "tests", "test_doctor_json.py"))
    assert re.search(r"^CONFIGS = \(", src, re.M), \
        "the changelog names a list this suite no longer has"


def test_no_derived_count_recomputes_a_literal():
    """C4, as a rule rather than as the one entry that broke it.

    `COUNT_PROSE["three states"]` was classified DERIVED with a recompute of
    `lambda ma: 3` - a rule with nothing under it, which is exactly the defect
    the staleness half of the scan above exists to prevent, one level up. A
    reader meeting DERIVED believes the number is pinned to the thing it
    counts; a recompute that is the literal can only ever agree with itself,
    so the classification is a claim about the registry that the registry does
    not keep.

    The test is structural, not a list of exceptions: a recompute that reads
    NOTHING - no global, no attribute, no closure - cannot be deriving
    anything. Every honest entry names at least one thing (`ma.STAGES`,
    `_example_configs()`, `ma.DMND_FASTA_RECORD_CAP`), so this costs the real
    ones nothing and catches the next literal the day it is written.
    """
    dead = []
    for phrase, entry in sorted(COUNT_PROSE.items()):
        if entry[0] != "DERIVED":
            continue
        code = entry[2].__code__
        if not code.co_names and not code.co_freevars:
            dead.append(f"{phrase!r} ({entry[1]}) recomputes "
                        f"{code.co_consts!r} and reads nothing")
    assert not dead, (
        "a DERIVED count derives nothing - classify it MEASURED or PROSE "
        "with the reason, or give it a real recompute:\n" + "\n".join(dead))
