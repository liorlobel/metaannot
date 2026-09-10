"""Documentation against the code.

Every number and every flag in README.md and TUTORIAL.md is an instruction
someone will follow on a shared server, so a doc that has drifted from the
code is a defect, not a typo.
"""
from __future__ import annotations

import json
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
    assert f"## v{ma.__version__}" in _text(CHANGELOG)


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
