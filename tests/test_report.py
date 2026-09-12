"""The embedded R report and the R object.

Everything that needs R is skipped cleanly when Rscript or a package is
missing; the static checks below run everywhere.
"""
from __future__ import annotations

import json
import os
import random
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest

from html import unescape

import fixtures as F
from conftest import METAANNOT_PY, ROOT, build_project, needs_r, r_has, \
    run_metaannot

R_CORE = ("readr", "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "purrr",
          "limma", "knitr", "rmarkdown")
R_BIOC = ("SummarizedExperiment", "S4Vectors")


def _chunks(rmd):
    """(label, code) for every R chunk in an Rmd."""
    return [(m.group(1), m.group(2)) for m in
            re.finditer(r"^```\{r([^}]*)\}\n(.*?)^```", rmd, re.S | re.M)]


def _rscript(code, *args):
    return subprocess.run(["Rscript", "-e", code, *[str(a) for a in args]],
                          capture_output=True, text=True, timeout=900)


# ----------------------------------------------------------------------
# static checks: no R needed
# ----------------------------------------------------------------------
# --- finding 55 -------------------------------------------------------
def test_the_analysis_block_covers_every_param_the_report_reads(ma):
    # symptom: a param the document reads but the generated header omits is
    # NULL at knit time, and the comparison it reaches fails obscurely.
    used = set(ma.template_params())
    have = set(ma.DEFAULT_CONFIG["analysis"]) | {"results_dir"}
    assert used <= have, f"params with no analysis: key: {sorted(used - have)}"


def test_a_missing_analysis_param_is_refused_before_the_knit(ma, project):
    cfg = ma.load_config(project.config_path)
    cfg["analysis"].pop("fdr")
    p = ma.Paths(cfg)
    p.mkdirs()
    with pytest.raises(ma.StageError) as e:
        ma.write_report_rmd(cfg, p)
    assert "missing parameters the report declares" in str(e.value)
    assert "fdr" in str(e.value)


def test_the_report_reads_the_per_protein_dominance_flag(ma):
    # symptom: taxon_unique_dominated is written into annotated_quant.tsv and
    # into the R object's rowData, and NOTHING read it back - the only place
    # the finding existed was one line of stderr, which the reader of an HTML
    # report never sees. The feature-support chunk is where it belongs,
    # because `aq` there is already this document's own protein set.
    code = next(c for lab, c in _chunks(ma.RMD_TEMPLATE)
                if lab.strip().startswith("feature-support"))
    assert "taxon_unique_dominated" in code
    assert "col_or_na(aq, \"taxon_unique_dominated\")" in code


def test_the_report_states_its_dominance_rate_over_its_own_protein_set(ma):
    # the denominator is nrow(aq) AFTER analysis.min_features has run, not
    # anything the join stage logged: the join stage reports over the proteins
    # its assignment rule decided about, and this document reports over the
    # proteins it is about. Pinned because the whole point of the change that
    # introduced this line was that a rate must name its own population.
    code = next(c for lab, c in _chunks(ma.RMD_TEMPLATE)
                if lab.strip().startswith("feature-support"))
    dom = code[code.index("dom <- as.logical"):]
    assert "nrow(aq)" in dom
    assert "params$min_features" not in dom, \
        "the rate must be computed after the filter, not alongside it"
    # and it carries the same denominator floor as every other coverage claim
    # in this document, so the two halves of the tool call the same size
    # 'too few'.
    assert "COVERAGE_MIN_N" in dom
    assert "Too few for that to be a rate" in dom


def test_the_two_halves_of_the_tool_agree_on_what_counts_as_too_few(ma):
    # ASSESSABLE_MIN_N gates the join stage's percentage and COVERAGE_MIN_N
    # gates the report's. They are deliberately the same number, and a docs
    # pin is the only thing that keeps one from being tuned without the other.
    r = re.search(r"^COVERAGE_MIN_N <- (\d+)$", ma.RMD_TEMPLATE, re.M)
    assert r, "COVERAGE_MIN_N is no longer declared in the report"
    assert ma.ASSESSABLE_MIN_N == int(r.group(1))


def test_the_analysis_block_has_no_key_the_report_ignores(ma):
    # the other direction: analysis: is checked for typos precisely because
    # its keys ARE the Rmd's params, so a key nothing reads is dead config.
    used = set(ma.template_params())
    unused = set(ma.DEFAULT_CONFIG["analysis"]) - used
    assert unused == set(), f"analysis keys the Rmd never reads: {unused}"


# --- finding 54 -------------------------------------------------------
def test_the_generated_rmd_carries_absolute_paths(ma, tmp_path):
    # symptom: rmarkdown::render() evaluates a document with the working
    # directory set to the Rmd's own folder, so anything relative in the
    # header resolves inside results/analysis/ and is not found.
    proj = build_project(tmp_path / "p")
    proj.write_config(results_dir="results")
    cfg = ma.load_config(proj.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    out = ma.write_report_rmd(cfg, p)
    header = open(out, encoding="utf-8").read().split("---")[1]
    for key in ("results_dir", "metadata"):
        m = re.search(rf"^  {key}: !r '\"(.*)\"'", header, re.M)
        assert m, f"{key} is missing from the generated header"
        assert os.path.isabs(m.group(1)), f"{key} is relative: {m.group(1)}"


def test_out_subdir_stays_relative(ma, tmp_path):
    # it is joined onto results_dir inside the document and passed to
    # build_object.R as a bare name; absolutising it breaks both.
    proj = build_project(tmp_path / "p")
    cfg = ma.load_config(proj.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    out = ma.write_report_rmd(cfg, p)
    header = open(out, encoding="utf-8").read().split("---")[1]
    assert re.search(r"^  out_subdir: !r '\"analysis\"'", header, re.M)


def test_a_path_shaped_param_added_later_is_absolutised_too(ma):
    assert ma._is_report_path_param("some_new_dir")
    assert ma._is_report_path_param("metadata")
    assert not ma._is_report_path_param("out_subdir")


def test_report_params_are_emitted_as_r_expressions(ma):
    # symptom: plain YAML gave every value YAML's escaping on the way into a
    # document that reads them as R data — a Windows path lost its
    # backslashes and an empty numeric arrived as the character "".
    assert ma.r_literal(True) == "TRUE"
    assert ma.r_literal(3) == "3L"
    assert ma.r_literal(0.05) == "0.05"
    assert ma.r_literal(None) == "NULL"
    assert ma.r_literal(r"C:\data\x") == '"C:\\\\data\\\\x"'
    assert ma.r_literal(["a", "b"]) == 'c("a", "b")'


def test_a_mapping_parameter_is_refused(ma):
    with pytest.raises(ma.StageError) as e:
        ma.r_literal({"a": 1})
    assert "the Rmd's params are scalars" in str(e.value)


def test_the_params_block_round_trips(ma):
    block = ma._render_params_block({"fdr": 0.05, "contrasts": "a - b",
                                     "quote": 'he said "no"',
                                     "win": r"C:\x\y"})
    assert "!r" in block


# --- finding 52 -------------------------------------------------------
def _calls(code, fname):
    """Every `fname(...)` call body in `code`, with balanced parentheses."""
    out = []
    for m in re.finditer(rf"\b{fname}\s*\(", code):
        i, depth = m.end(), 1
        while i < len(code) and depth:
            if code[i] == "(":
                depth += 1
            elif code[i] == ")":
                depth -= 1
            i += 1
        out.append(code[m.end():i - 1])
    return out


def _split_args(body):
    args, depth, cur = [], 0, ""
    for ch in body:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            args.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        args.append(cur)
    return args


def _strip_r_comments(code):
    """Drop `#` comments, keeping quoted strings intact."""
    out = []
    for line in code.splitlines():
        res, in_str, quote, i = "", False, "", 0
        while i < len(line):
            c = line[i]
            if in_str:
                res += c
                if c == "\\":
                    if i + 1 < len(line):
                        res += line[i + 1]
                    i += 2
                    continue
                if c == quote:
                    in_str = False
            elif c in "\"'":
                in_str, quote = True, c
                res += c
            elif c == "#":
                break
            else:
                res += c
            i += 1
        out.append(res)
    return "\n".join(out)


def _mentions(expr, name):
    """`name` used as a bare column reference in `expr`.

    Not `name(...)` (a function of the same name, e.g. dplyr's n()), and not
    `df$name` / `pkg::name` (a column of ANOTHER frame, e.g. sf$method).
    """
    return bool(re.search(rf"(?<![$@:\w.]){re.escape(name)}\b(?!\s*\()",
                          expr))


def test_no_summarise_call_reuses_a_name_it_is_redefining(ma):
    # symptom: `summarise(kept = sum(kept), pct = 100 * mean(kept))` — the
    # second argument sees the ALREADY-REDEFINED column, so the percentage
    # came out as 400%. Only a SELF-referential assignment counts as a
    # redefinition; `total = a + b` followed by `pct = x / total` is the
    # ordinary, intended use of sequential evaluation.
    bad = []
    for source in (ma.RMD_TEMPLATE, ma.ROBJECT_SCRIPT):
        code = _strip_r_comments(source)
        for verb in ("summarise", "summarize", "mutate", "transmute"):
            for body in _calls(code, verb):
                redefined = []
                for arg in _split_args(body):
                    lhs = re.match(r"\s*([A-Za-z._][\w.]*)\s*=(?!=)", arg)
                    rhs = arg.split("=", 1)[1] if lhs else arg
                    for d in redefined:
                        if _mentions(rhs, d):
                            bad.append((verb, d, arg.strip()[:70]))
                    if lhs and _mentions(rhs, lhs.group(1)):
                        redefined.append(lhs.group(1))
    assert bad == [], f"a redefined name is read later in the same call: {bad}"


# --- finding 58 -------------------------------------------------------
def _r_vector(source, name):
    # comments first: both lists carry explanatory comments that themselves
    # quote strings ("816"), which would otherwise be read as members.
    code = _strip_r_comments(source)
    m = re.search(rf"{name} <- c\((.*?)\)\n", code, re.S)
    assert m, f"{name} is gone"
    return [x for x in re.findall(r'"([^"]+)"', m.group(1))]


def test_the_report_and_the_object_agree_on_what_is_not_a_sample(ma):
    # symptom: the fallback heuristic "which columns look numeric" could sweep
    # 46 annotation columns into the assay. Both components carry the same
    # exclusion list, and if the two disagree the object and the model are
    # built on different matrices.
    a = _r_vector(ma.ROBJECT_SCRIPT, "ANNOT_NUM")
    b = _r_vector(ma.RMD_TEMPLATE, "ANNOT_NUMERIC")
    assert a == b


def test_sample_columns_are_identified_positively_before_any_heuristic(ma):
    # the order must be: recorded list, then the design, then size factors,
    # then the feature table, and only then column types.
    for source in (ma.ROBJECT_SCRIPT, ma.RMD_TEMPLATE):
        # the Rmd seeds the variable with a default before the if-chain, so
        # compare the positions of the branches rather than the first write.
        pos = {m.group(1): m.start() for m in re.finditer(
            r'sample_s(?:ource|rc) <- "([^"]+)"', source)}
        named = {k: v for k, v in pos.items() if "types" not in k}
        assert any(k.endswith("sample_columns.txt") for k in named)
        first = min(named, key=named.get)
        assert first.endswith("sample_columns.txt"), first
        last = max(pos, key=pos.get)
        assert "types" in last, last


def test_the_join_stage_writes_the_sample_column_list(tmp_path):
    # symptom: both build_object.R and the Rmd read
    # quant/sample_columns.txt as the most authoritative statement of which
    # columns are samples, and nothing wrote it — so every run fell through to
    # design_from_input.tsv, or with no manifest to the column-type heuristic
    # the file exists to avoid.
    proj = build_project(tmp_path / "p")
    proj.run()
    path = proj.rpath("quant", "sample_columns.txt")
    assert os.path.exists(path)
    names = open(path, encoding="utf-8").read().splitlines()
    assert names == proj.samples, "one sample name per line, in column order"


def test_the_recorded_sample_columns_are_columns_of_annotated_quant(tmp_path):
    # the file is only worth reading if it names the columns of the table it
    # sits beside: the R side intersects the two.
    proj = build_project(tmp_path / "p")
    proj.run()
    names = open(proj.rpath("quant", "sample_columns.txt"),
                 encoding="utf-8").read().splitlines()
    aq = pd.read_csv(proj.rpath("quant", "annotated_quant.tsv"), sep="\t",
                     nrows=0)
    assert names, "the list must not be empty for a run that quantified"
    assert set(names) <= set(aq.columns)
    # and it must name every column that actually carries intensities
    assert set(names) == set(proj.samples)


def test_the_sample_column_list_survives_a_manifest_rename(tmp_path):
    # the quant table's columns are 'A_1 Intensity'; the manifest renames them
    # to 'A_1'. The recorded list has to be the names as WRITTEN, not as read.
    proj = build_project(tmp_path / "p")
    proj.run()
    names = open(proj.rpath("quant", "sample_columns.txt"),
                 encoding="utf-8").read().splitlines()
    assert not any(n.endswith("Intensity") for n in names)
    assert names == ["A_1", "A_2", "B_1", "B_2"]


def test_the_design_is_used_when_no_sample_column_list_exists(tmp_path):
    # the documented second choice: it must actually be produced by a run.
    proj = build_project(tmp_path / "p")
    proj.run()
    d = pd.read_csv(proj.rpath("quant", "design_from_input.tsv"), sep="\t")
    assert sorted(d["sample"]) == sorted(proj.samples)


# --- the TMT design: the plex is a batch, not a hypothesis ------------
def _tmt_project(ma, tmp_path, rows, name="results", **analysis):
    """A results directory holding only a TMT design_from_input.tsv.

    `rows` is [(sample, plex, group?)]. Returns (cfg, Paths).
    """
    import json
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["results_dir"] = str(tmp_path / name)
    cfg["quant_format"] = "fragpipe_tmt"
    cfg["analysis"].update(analysis)
    p = ma.Paths(cfg)
    p.mkdirs()
    os.makedirs(p.quant_dir, exist_ok=True)
    cols = ["sample", "plex"] + (["group"] if len(rows[0]) > 2 else [])
    with open(os.path.join(p.quant_dir, "design_from_input.tsv"), "w",
              encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return cfg, p


def _header(out):
    return open(out, encoding="utf-8").read().split("---")[1]


def _chunk(ma, label):
    """One chunk of the report by name; the label carries its options."""
    for lab, code in _chunks(ma.RMD_TEMPLATE):
        if lab.split(",")[0].strip() == label:
            return code
    raise AssertionError(f"the report has no chunk named {label}")


TMT_CROSSED = [("s1", "TMT1", "a"), ("s2", "TMT1", "b"),
               ("s3", "TMT2", "a"), ("s4", "TMT2", "b")]


def test_a_tmt_run_puts_the_plex_in_the_model_and_contrasts_the_condition(
        ma, tmp_path):
    # symptom: the input names the PLEX per sample, not the condition, so the
    # ordinary default would contrast plex against plex - a batch effect
    # presented as a hypothesis.
    cfg, p = _tmt_project(ma, tmp_path, TMT_CROSSED)
    out = ma.write_report_rmd(cfg, p)
    h = _header(out)
    assert re.search(r'design_formula: !r \'"~ 0 \+ group \+ plex"\'', h)
    assert re.search(r'factor_cols: !r \'"group,plex"\'', h)
    # the contrast is over the condition, and names no plex
    m = re.search(r'^  contrasts: !r \'"(.*)"\'', h, re.M)
    assert m and m.group(1) == "b_vs_a = groupb - groupa", h


def test_a_hand_written_tmt_formula_is_left_alone(ma, tmp_path):
    # a formula the user wrote is their statement of the model.
    cfg, p = _tmt_project(ma, tmp_path, TMT_CROSSED,
                          design_formula="~ 0 + group", factor_cols="group")
    cfg["analysis"]["design_formula"] = "~ group + plex"
    out = ma.write_report_rmd(cfg, p)
    assert re.search(r'design_formula: !r \'"~ group \+ plex"\'', _header(out))


def test_a_formula_that_does_not_model_the_plex_keeps_it_out_of_factor_cols(
        ma, tmp_path):
    # symptom: the formula and factor_cols were replaced independently, so a
    # user who wrote "~ group" kept their formula and still got factor_cols
    # 'group,plex'; the knit then stopped on the report's own
    # "factor_cols names a missing column" over a column their metadata had no
    # reason to carry, naming a config key they never set.
    cfg, p = _tmt_project(ma, tmp_path, TMT_CROSSED)
    meta = tmp_path / "meta.tsv"          # the plex is not modelled, so it is
    meta.write_text("sample\tgroup\ns1\ta\ns2\tb\ns3\ta\ns4\tb\n",  # not here
                    encoding="utf-8")
    cfg["analysis"]["metadata"] = str(meta)
    cfg["analysis"]["design_formula"] = "~ group"
    h = _header(ma.write_report_rmd(cfg, p))
    assert re.search(r'design_formula: !r \'"~ group"\'', h)
    assert re.search(r'factor_cols: !r \'"group"\'', h)
    for cc in re.search(r'factor_cols: !r \'"(.*)"\'', h).group(1).split(","):
        assert cc in list(pd.read_csv(meta, sep="\t", nrows=0).columns), \
            "the Rmd stops the knit on a factor_cols column the metadata lacks"
    # a formula that DOES model the plex still types it, interaction included
    cfg["analysis"]["design_formula"] = "~ 0 + group * plex"
    cfg["analysis"]["factor_cols"] = "group"
    cfg["analysis"]["metadata"] = ""
    assert re.search(r'factor_cols: !r \'"group,plex"\'',
                     _header(ma.write_report_rmd(cfg, p)))


def test_a_single_plex_run_does_not_get_a_plex_term(ma, tmp_path, capsys):
    # model.matrix() cannot make a contrast for a one-level factor, and one
    # plex is not a batch effect.
    cfg, p = _tmt_project(ma, tmp_path, [("s1", "TMT1", "a"),
                                         ("s2", "TMT1", "b")])
    out = ma.write_report_rmd(cfg, p)
    assert re.search(r'design_formula: !r \'"~ 0 \+ group"\'', _header(out))
    assert "one level is not a batch effect" in capsys.readouterr().err


def test_a_tmt_design_with_no_condition_demands_the_metadata_file(ma,
                                                                  tmp_path):
    # symptom: the annotation carries MF#### codes and nothing else, so there
    # is no condition anywhere. Guessing one from the plex would report a
    # batch effect as biology; this says exactly what to write instead.
    cfg, p = _tmt_project(ma, tmp_path, [("MF0030", "TMT1"),
                                         ("MF0071", "TMT2")])
    with pytest.raises(ma.StageError) as e:
        ma.write_report_rmd(cfg, p)
    msg = str(e.value)
    assert "no 'group' column" in msg
    assert "plex is a TMT batch, not a condition" in msg
    assert "cp " in msg and "design_from_input.tsv metadata.tsv" in msg
    assert "sample\tplex\tgroup" in msg          # the header to write
    assert "MF0030\tTMT1\t<condition>" in msg    # a filled-in row
    assert "metadata: metadata.tsv" in msg
    assert "tmt.condition_from_name" in msg


def test_a_tmt_metadata_file_without_the_plex_is_refused_by_name(ma,
                                                                 tmp_path):
    # the formula models the plex; a metadata file that dropped the column
    # otherwise fails inside the knit, naming the formula and not the fix.
    cfg, p = _tmt_project(ma, tmp_path, TMT_CROSSED)
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tgroup\ns1\ta\ns2\tb\ns3\ta\ns4\tb\n",
                    encoding="utf-8")
    cfg["analysis"]["metadata"] = str(meta)
    with pytest.raises(ma.StageError) as e:
        ma.write_report_rmd(cfg, p)
    msg = str(e.value)
    assert "no 'plex' column" in msg
    assert "design_from_input.tsv" in msg


def test_tmt_contrasts_come_from_the_metadata_the_user_wrote(ma, tmp_path):
    # the recovered design has the plex and often no condition at all; the
    # condition is in analysis.metadata, so that is where the contrast is.
    cfg, p = _tmt_project(ma, tmp_path, [("s1", "TMT1"), ("s2", "TMT1"),
                                         ("s3", "TMT2"), ("s4", "TMT2")])
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tplex\tgroup\ns1\tTMT1\tresp\ns2\tTMT1\tnon\n"
                    "s3\tTMT2\tresp\ns4\tTMT2\tnon\n", encoding="utf-8")
    cfg["analysis"]["metadata"] = str(meta)
    out = ma.write_report_rmd(cfg, p)
    m = re.search(r'^  contrasts: !r \'"(.*)"\'', _header(out), re.M)
    assert m and m.group(1) == "resp_vs_non = groupresp - groupnon"


def test_the_label_free_default_design_is_untouched(ma, tmp_path):
    # the TMT branch must be reachable only from quant_format fragpipe_tmt.
    cfg, p = _tmt_project(ma, tmp_path, TMT_CROSSED)
    cfg["quant_format"] = "fragpipe_peptide"
    out = ma.write_report_rmd(cfg, p)
    h = _header(out)
    assert re.search(r'design_formula: !r \'"~ 0 \+ group"\'', h)
    assert re.search(r'factor_cols: !r \'"group"\'', h)


def test_the_report_says_the_plex_is_a_batch_where_it_prints_the_design(ma):
    # the coefficient list is where a reader decides what the model means.
    code = _chunk(ma, "design")
    assert "plex" in code and "BATCH" in code
    rec = _chunk(ma, "export")
    assert "design_notes.txt" in rec, \
        "design_record.txt must carry how the design was made"


# --- finding: plex confounded with the condition ----------------------
def _plex_guard_source(ma):
    """The report's own plex check, lifted out so a test can run it."""
    code = _chunk(ma, "design-diagnostics")
    m = re.search(r"^plex_confounding <- function.*?^\}$", code, re.S | re.M)
    assert m, "the design-diagnostics chunk no longer defines plex_confounding"
    return m.group(0)


def test_the_report_names_plex_when_it_is_confounded_with_the_condition(ma):
    # static half: the message has to name the term, because limma's own
    # symptom is "coefficient plexTMT8 is not estimable".
    src = _plex_guard_source(ma)
    assert "perfectly confounded" in src
    assert "spread over several plexes" in src


@needs_r()
@pytest.mark.parametrize("rows,expect", [
    # every plex holds one condition: the batch and the biology are one vector
    ([("s1", "TMT1", "a"), ("s2", "TMT1", "a"),
      ("s3", "TMT2", "b"), ("s4", "TMT2", "b")], "perfectly confounded"),
    # each condition crosses both plexes: adjustable, and must not stop
    ([("s1", "TMT1", "a"), ("s2", "TMT1", "b"),
      ("s3", "TMT2", "a"), ("s4", "TMT2", "b")], "ok"),
])
def test_the_plex_guard_stops_a_confounded_tmt_design(ma, tmp_path, rows,
                                                      expect):
    # symptom: a plex nested inside the condition fits, or fails deep in
    # makeContrasts naming a coefficient rather than the design decision.
    md = tmp_path / "design.tsv"
    md.write_text("sample\tplex\tgroup\n"
                  + "".join(f"{s}\t{p}\t{g}\n" for s, p, g in rows),
                  encoding="utf-8")
    code = (_plex_guard_source(ma) +
            '\nmd <- read.delim(commandArgs(TRUE)[1], stringsAsFactors = TRUE)'
            '\nplex_confounding(md)\ncat("ok")\n')
    r = _rscript(code, str(md))
    out = r.stdout + r.stderr
    assert expect in out, out
    if expect != "ok":
        assert "plex" in out and "TMT1" in out


# ----------------------------------------------------------------------
# R
# ----------------------------------------------------------------------
# --- finding 50 -------------------------------------------------------
@needs_r()
def test_every_r_chunk_in_the_report_parses(ma, tmp_path):
    # symptom: a syntax error in a chunk halfway down the document is only
    # found by a knit, which needs a finished run and several minutes.
    chunks = _chunks(ma.RMD_TEMPLATE)
    assert len(chunks) > 10, "chunk extraction is wrong"
    names = []
    for i, (label, code) in enumerate(chunks):
        f = tmp_path / f"chunk_{i:02d}.R"
        f.write_text(code, encoding="utf-8")
        names.append((label.strip() or f"#{i}", str(f)))
    # one R process for the lot: starting R per chunk is most of the cost.
    r = _rscript(
        'for (f in commandArgs(TRUE)) { '
        '  e <- tryCatch({parse(f); NULL}, error = function(c) c); '
        '  if (!is.null(e)) cat(f, ":", conditionMessage(e), "\n") }',
        *[path for _, path in names])
    assert r.returncode == 0, r.stderr
    if r.stdout.strip():
        label_of = {path: lab for lab, path in names}
        lines = [f"{label_of.get(l.split(' :')[0], l)}: {l}"
                 for l in r.stdout.strip().splitlines()]
        raise AssertionError("chunks that do not parse:\n" + "\n".join(lines))


@needs_r()
def test_the_object_script_parses(ma, tmp_path):
    f = tmp_path / "build_object.R"
    f.write_text(ma.ROBJECT_SCRIPT, encoding="utf-8")
    r = _rscript('invisible(parse(commandArgs(TRUE)[1]))', str(f))
    assert r.returncode == 0, r.stderr


@needs_r()
def test_the_generated_rmd_header_parses_as_yaml_and_r(ma, tmp_path):
    proj = build_project(tmp_path / "p")
    cfg = ma.load_config(proj.config_path)
    p = ma.Paths(cfg)
    p.mkdirs()
    out = ma.write_report_rmd(cfg, p)
    r = _rscript('cat(length(rmarkdown::yaml_front_matter('
                 'commandArgs(TRUE)[1])$params))', out)
    assert r.returncode == 0, r.stderr
    assert int(r.stdout.strip() or 0) > 5


# --- finding 51 -------------------------------------------------------
@needs_r("dplyr", "tidyr")
def test_no_report_helper_shadows_a_dplyr_or_tidyr_export(ma, tmp_path):
    # symptom: a helper named pick() was intercepted by dplyr inside mutate()
    # and never called, so the report failed on its second data chunk.
    helpers = sorted(set(re.findall(r"^([A-Za-z._][\w.]*)\s*<-\s*function",
                                    ma.RMD_TEMPLATE, re.M)))
    assert helpers, "no helpers found; the extraction is wrong"
    script = tmp_path / "probe.R"
    script.write_text(
        'a <- commandArgs(TRUE)\n'
        'ex <- unique(c(getNamespaceExports("dplyr"), '
        'getNamespaceExports("tidyr")))\n'
        'cat(paste(intersect(a, ex), collapse=" "))\n', encoding="utf-8")
    r = subprocess.run(["Rscript", str(script), *helpers],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr
    clashes = r.stdout.split()
    assert clashes == [], f"report helpers shadow tidyverse exports: {clashes}"


def _dominance_probe(ma, tmp_path, aq_expr, name="dom.R"):
    """Run the report's OWN dominance lines over a constructed `aq`.

    The helpers and the floor are lifted from RMD_TEMPLATE rather than
    restated here, so the probe cannot pass against an R snippet the report no
    longer contains.
    """
    t = ma.RMD_TEMPLATE
    helpers = t[t.index("note <- function"):t.index("need <- function")]
    floor = re.search(r"^COVERAGE_MIN_N <- \d+$", t, re.M).group(0)
    code = next(c for lab, c in _chunks(t)
                if lab.strip().startswith("feature-support"))
    snippet = code[code.index("  dom <- as.logical"):code.index("  ev_path <-")]
    script = tmp_path / name
    script.write_text(f"{helpers}\n{floor}\naq <- {aq_expr}\n{snippet}\n",
                      encoding="utf-8")
    r = subprocess.run(["Rscript", str(script)], capture_output=True,
                       text=True, timeout=600)
    assert r.returncode == 0, r.stderr
    return r.stdout


def _aq(n_dom, n_clean):
    return (f"data.frame(taxon_unique_dominated = c(rep(TRUE, {n_dom}), "
            f"rep(FALSE, {n_clean})))")


@needs_r()
def test_the_reports_dominance_rate_is_denominated_on_the_rows_it_kept(
        ma, tmp_path):
    # the number a reader of the RESULTS stands on, which the join stage
    # cannot state honestly because the report filters again after it. Its
    # denominator is nrow(aq) at the point the line runs - this document's own
    # protein set, after min_features_per_protein and analysis.min_features.
    out = _dominance_probe(ma, tmp_path, _aq(7, 5))
    assert out.startswith("GATE")
    assert "7/12 (58.3%)" in out
    assert "taxon_unique_dominated" in out
    assert "protein_unique" in out, \
        "the line must name the run that shows what survives the assumption"


@needs_r()
def test_the_reports_dominance_rate_withholds_a_percentage_below_the_floor(
        ma, tmp_path):
    # same floor as every other coverage claim in this document: over a
    # handful of proteins a percentage is an anecdote wearing a decimal point.
    out = _dominance_probe(ma, tmp_path, _aq(3, 4), name="few.R")
    assert out.startswith("NOTE")
    assert "3 of the 7 protein(s)" in out
    assert "%" not in out


@needs_r()
def test_the_reports_dominance_line_is_silent_when_nothing_is_flagged(
        ma, tmp_path):
    # a line that fires on every protein_unique run is the line a reader
    # learns to skip, and the zero is already one column of the table.
    assert _dominance_probe(ma, tmp_path, _aq(0, 40), name="none.R") == ""
    # and a protein-level input, which carries no evidence columns at all,
    # must produce a rate rather than an error - col_or_na fills with NA.
    out = _dominance_probe(ma, tmp_path, "data.frame(bin = rep('x', 12))",
                           name="na.R")
    assert "rest more" not in out


# --- a full knit ------------------------------------------------------
@pytest.fixture(scope="module")
def knitted(tmp_path_factory):
    """One finished run with enough replication to fit the model, knitted
    once. Rendering is the expensive step, so every report assertion reads
    the same rendered output rather than knitting again."""
    root = tmp_path_factory.mktemp("knit")
    samples = ["A_1", "A_2", "A_3", "A_4", "B_1", "B_2", "B_3", "B_4"]
    proj = build_project(root / "p", samples=samples, n_extra=40)
    proj.run()
    proj.rendered = False
    if shutil.which("Rscript") and shutil.which("pandoc") and r_has(*R_CORE):
        run_metaannot("report", "--config", proj.config_path, cwd=proj.root,
                      timeout=1800)
        proj.rendered = True
    return proj


@needs_r(*R_CORE, *R_BIOC)
def test_the_report_knits_to_html(knitted, ma):
    # symptom: the README says the report has only ever knitted by hand, on
    # synthetic data, with no test. This is that check, in the suite.
    if not knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    html = knitted.rpath("analysis", "analyse_metaannot.html")
    assert os.path.exists(html) and os.path.getsize(html) > 10000
    txt = open(html, encoding="utf-8").read()
    # it must have reached the end, not stopped after the setup chunk
    for marker in ("1_ko_pathway", "invisible to KEGG", "GATE"):
        assert marker in txt, f"the rendered report never mentions {marker}"


@needs_r(*R_CORE, *R_BIOC)
def test_the_report_writes_the_tables_the_object_folds_back_in(knitted):
    if not knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    out = knitted.rpath("analysis")
    written = set(os.listdir(out))
    assert any(f.startswith("differential_abundance_") for f in written)
    assert "taxon_normalisation_risk.tsv" in written
    assert any(f.startswith("effector_candidates_") for f in written)


@needs_r(*R_CORE, *R_BIOC)
def test_the_risk_table_only_uses_the_documented_verdicts(knitted):
    if not knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    risk = pd.read_csv(knitted.rpath("analysis",
                                     "taxon_normalisation_risk.tsv"), sep="\t")
    assert "assumption" in risk.columns
    assert set(risk["assumption"]) <= {
        "AT RISK", "ok", "check (median in a sparse region)",
        "fragile (too few proteins)"}


# --- finding 53 -------------------------------------------------------
@needs_r("readr")
def test_a_zero_row_group_conflicts_file_does_not_crash_the_gates_chunk(
        knitted, tmp_path):
    # symptom: a zero-row file still has a header, readr types empty columns
    # as character, and `|` fails on them — so the GOOD case crashed.
    gp = knitted.rpath("quant", "group_conflicts.tsv")
    df = pd.read_csv(gp, sep="\t")
    assert len(df) == 0, "the fixture already has no conflicts"
    code = (
        'suppressPackageStartupMessages({library(readr)})\n'
        'gcf <- read_tsv(commandArgs(TRUE)[1], show_col_types = FALSE)\n'
        'lg <- function(x) if (is.null(x)) logical(0) else as.logical(x)\n'
        'cat(sum(lg(gcf$bin_conflict) | lg(gcf$ko_conflict), na.rm = TRUE),\n'
        '    sum(lg(gcf$unannotated_group), na.rm = TRUE))\n')
    r = _rscript(code, gp)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "0 0"


# --- finding 56 -------------------------------------------------------
DESIGN_GUARD = """
suppressPackageStartupMessages({library(limma)})
a <- commandArgs(TRUE)
md <- read.delim(a[1], stringsAsFactors = TRUE)
design <- model.matrix(as.formula(a[2]), data = md)
r <- qr(design)$rank
if (r < ncol(design)) stop("rank deficient")
if (nrow(design) - r < 1) stop("0 residual degrees of freedom")
fac <- names(md)[vapply(md, is.factor, logical(1))]
one <- fac[vapply(fac, function(f) nlevels(md[[f]]) < 2, logical(1))]
if (length(one)) stop("factor(s) with only one level: ",
                      paste(one, collapse = ", "))
cat("ok")
"""


@needs_r("limma")
@pytest.mark.parametrize("rows,formula,expect", [
    # one sample per group: no residual degrees of freedom
    ([("s1", "a"), ("s2", "b")], "~ 0 + group", "0 residual degrees"),
    # a factor with a single level
    # a factor with a single level: nothing can be contrasted against
    # anything, whether the Rmd's own guard catches it or model.matrix does
    ([("s1", "a"), ("s2", "a"), ("s3", "a")], "~ 0 + group", "level"),
    # a usable design
    ([("s1", "a"), ("s2", "a"), ("s3", "b"), ("s4", "b")], "~ 0 + group",
     "ok"),
])
def test_the_design_guards_stop_before_the_model_is_fitted(tmp_path, rows,
                                                           formula, expect):
    # symptom: rank deficiency fails deep inside makeContrasts with an
    # unhelpful error, or worse, fits and returns nonsense.
    md = tmp_path / "design.tsv"
    md.write_text("sample\tgroup\n"
                  + "".join(f"{s}\t{g}\n" for s, g in rows), encoding="utf-8")
    r = _rscript(DESIGN_GUARD, str(md), formula)
    out = (r.stdout + r.stderr)
    assert expect in out, out


@needs_r(*R_CORE, *R_BIOC)
def test_zero_variance_proteins_are_counted_and_listed(tmp_path):
    # symptom: limma MODERATES a constant row rather than failing, which turns
    # it into an enormous t-statistic — usually one shared peptide or an
    # imputed constant, not biology.
    if not shutil.which("pandoc"):
        pytest.skip("pandoc not on PATH")
    samples = ["A_1", "A_2", "A_3", "A_4", "B_1", "B_2", "B_3", "B_4"]
    proj = build_project(tmp_path / "zv", samples=samples, n_extra=30)
    q = proj.path("input", "combined_peptide.tsv")
    d = pd.read_csv(q, sep="\t")
    flat = d["Protein"] == "P_ko_path"
    for c in [f"{s} Intensity" for s in samples]:
        d.loc[flat, c] = 50000.0            # identical in every sample
    d.to_csv(q, sep="\t", index=False)
    # normalise: none, or the per-sample median shift makes a protein that is
    # constant in raw intensity non-constant by the time the model sees it.
    proj.write_config(analysis=dict(proj.cfg.get("analysis") or {},
                                    normalise="none"))
    proj.run()
    run_metaannot("report", "--config", proj.config_path, cwd=proj.root,
                  timeout=1800)
    zv = proj.rpath("analysis", "zero_variance.tsv")
    assert os.path.exists(zv), "a constant protein was not reported"
    listed = pd.read_csv(zv, sep="\t")
    assert listed.columns.tolist() == ["group_id"]
    assert "P_ko_path" in set(listed["group_id"])
    html = open(proj.rpath("analysis", "analyse_metaannot.html"),
                encoding="utf-8").read()
    assert "zero variance within every group" in html


# --- finding 40 -------------------------------------------------------
RISK_RULE = """
suppressPackageStartupMessages({library(dplyr)})
a <- commandArgs(TRUE)
d <- read.delim(a[1])
p <- list(min_lfc = 0.585, assumption_frac_deviating = 0.5,
          assumption_opposing = 0.6, assumption_min_deviating = 3,
          assumption_opposing_frac = 0.3)
out <- d %>% group_by(taxid, factor_log2FC) %>%
  summarise(n_tested = n(),
            n_dev = sum(deviates),
            n_opposing = sum(deviates & sign(lfc) != sign(factor_log2FC)),
            method = dplyr::first(method),
            .groups = "drop") %>%
  mutate(frac_deviating = n_dev / n_tested,
         opposing = ifelse(n_dev > 0, n_opposing / n_dev, 0),
         opposing_frac = n_opposing / n_tested,
         assumption = case_when(
           method == "sum_fallback" ~ "fragile (too few proteins)",
           abs(factor_log2FC) > p$min_lfc &
             n_dev >= p$assumption_min_deviating &
             opposing > p$assumption_opposing &
             opposing_frac >= p$assumption_opposing_frac ~ "AT RISK",
           frac_deviating > p$assumption_frac_deviating ~
             "check (median in a sparse region)",
           TRUE ~ "ok"))
write.csv(out[, c("taxid", "assumption")], row.names = FALSE)
"""


def test_the_risk_probe_matches_the_rule_in_the_report(ma):
    # the probe above reproduces the Rmd's case_when; if the Rmd's rule
    # changes, the probe must change with it or it protects nothing.
    rmd = _strip_r_comments(ma.RMD_TEMPLATE)
    for frag in ('method == "sum_fallback" ~ "fragile (too few proteins)"',
                 'n_dev >= params$assumption_min_deviating',
                 'opposing > params$assumption_opposing',
                 'opposing_frac >= params$assumption_opposing_frac ~ "AT RISK"',
                 'frac_deviating > params$assumption_frac_deviating'):
        assert frag in re.sub(r"\s+", " ", rmd), frag
    probe = re.sub(r"\s+", " ", RISK_RULE).replace("p$", "params$")
    for frag in ('method == "sum_fallback"',
                 'n_dev >= params$assumption_min_deviating',
                 'opposing > params$assumption_opposing'):
        assert frag in probe, frag


@needs_r("dplyr")
def test_normalisation_risk_flags_a_taxon_whose_majority_moved_together(
        tmp_path):
    # symptom: the discriminator is deviations one-sided AGAINST a factor that
    # moved — not simply "many proteins deviate", which flags every taxon that
    # merely changed abundance.
    rows = ["taxid\tfactor_log2FC\tdeviates\tlfc\tmethod"]
    # AT RISK: the size factor moved by +1, and most members deviate the OTHER
    # way, which is what a whole-taxon shift absorbed into the reference looks
    # like.
    for i in range(10):
        rows.append(f"risky\t1.0\t{'TRUE' if i < 8 else 'FALSE'}\t-1.2"
                    "\tmedian_of_ratios")
    # ok: the taxon changed abundance and its members agree with it, so the
    # few that deviate do so in the SAME direction as the factor.
    for i in range(10):
        rows.append(f"clean\t1.0\t{'TRUE' if i < 2 else 'FALSE'}\t1.4"
                    "\tmedian_of_ratios")
    # ok: nothing moved.
    for i in range(10):
        rows.append("still\t0.05\tFALSE\t0.02\tmedian_of_ratios")
    f = tmp_path / "dev.tsv"
    f.write_text("\n".join(rows) + "\n", encoding="utf-8")
    r = _rscript(RISK_RULE, str(f))
    assert r.returncode == 0, r.stderr
    got = dict(l.strip('"').split('","') for l in
               r.stdout.strip().splitlines()[1:])
    assert got["clean"] == "ok"
    assert got["still"] == "ok"
    assert got["risky"] == "AT RISK"


@needs_r("dplyr")
def test_a_single_opposing_protein_does_not_flag_a_taxon(tmp_path):
    # symptom: opposing = n_opposing / n_dev gives 1 for one-out-of-one, which
    # used to flag the whole taxon. assumption_min_deviating is 3 for that
    # reason.
    rows = ["taxid\tfactor_log2FC\tdeviates\tlfc\tmethod",
            "t\t1.0\tTRUE\t-1.2\tmedian_of_ratios"]
    rows += ["t\t1.0\tFALSE\t1.0\tmedian_of_ratios"] * 9
    f = tmp_path / "dev.tsv"
    f.write_text("\n".join(rows) + "\n", encoding="utf-8")
    r = _rscript(RISK_RULE, str(f))
    assert r.returncode == 0, r.stderr
    assert "AT RISK" not in r.stdout


# --- finding 57 -------------------------------------------------------
@needs_r("SummarizedExperiment", "S4Vectors", "QFeatures")
def test_the_object_builds_as_qfeatures_with_a_real_assay_link(tmp_path):
    # symptom: the QFeatures assay-link branch — the headline deliverable of
    # the object script — sits inside two nested try(..., silent = TRUE) and
    # had never been exercised.
    proj = build_project(tmp_path / "p")
    proj.run()
    run_metaannot("object", "--config", proj.config_path, cwd=proj.root,
                  timeout=900)
    rds = proj.rpath("metaannot.rds")
    assert os.path.exists(rds)
    r = _rscript(
        'suppressPackageStartupMessages(library(QFeatures)); '
        'q <- readRDS(commandArgs(TRUE)[1]); '
        'al <- assayLink(q, "proteins"); '
        'cat(class(q)[1], paste(names(q), collapse=","), al@from, al@fcol, '
        'length(al@hits))', rds)
    assert r.returncode == 0, r.stderr
    cls, assays, frm, fcol, hits = r.stdout.split()
    assert cls == "QFeatures", f"got {cls}: the link branch was swallowed"
    assert assays == "peptides,proteins"
    assert frm == "peptides"
    # the link records the assignment metaannot MADE, not one re-derived in R
    assert fcol == "assigned_protein"
    assert int(hits) > 0, "the assay link matched no feature"


@needs_r("SummarizedExperiment", "S4Vectors", "QFeatures")
def test_a_failing_assay_link_reports_the_real_condition(tmp_path):
    # symptom: it printed a hard-coded "QFeatures version differs" whatever
    # failed, which sent people after the wrong problem. addAssayLink fails
    # for several reasons; "No match found" is one of them.
    proj = build_project(tmp_path / "p")
    proj.run()
    fq = proj.rpath("quant", "feature_quant.tsv")
    d = pd.read_csv(fq, sep="\t")
    d["assigned_protein"] = "NOT_A_GROUP_ID"
    d.to_csv(fq, sep="\t", index=False)
    run_metaannot("object", "--config", proj.config_path, "--no-run",
                  cwd=proj.root)
    script = proj.rpath("analysis", "build_object.R")
    r = subprocess.run(["Rscript", script, proj.results,
                        proj.rpath("metaannot.rds"), "analysis"],
                       capture_output=True, text=True, timeout=900)
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    assert "QFeatures version differs" not in out
    assert "assay link not added" in out, out[-800:]
    # and the message must carry the REAL condition, not a guess
    assert re.search(r"assay link not added \((?!\)).+\)", out)


@needs_r("SummarizedExperiment", "S4Vectors")
def test_the_object_falls_back_rather_than_trapping_the_data(tmp_path):
    # the data is never trapped behind a missing dependency: without
    # QFeatures the object is a SummarizedExperiment with the peptide assay
    # in metadata()$peptides.
    proj = build_project(tmp_path / "p")
    proj.run()
    script = proj.rpath("analysis", "build_object.R")
    run_metaannot("object", "--config", proj.config_path, "--no-run",
                  cwd=proj.root)
    assert os.path.exists(script)
    # simulate QFeatures being absent by shadowing requireNamespace
    patched = tmp_path / "no_qf.R"
    body = open(script, encoding="utf-8").read().replace(
        'have <- function(pkg) requireNamespace(pkg, quietly = TRUE)',
        'have <- function(pkg) pkg != "QFeatures" && '
        'requireNamespace(pkg, quietly = TRUE)')
    patched.write_text(body, encoding="utf-8")
    out = tmp_path / "se.rds"
    r = subprocess.run(["Rscript", str(patched), proj.results, str(out),
                        "analysis"], capture_output=True, text=True,
                       timeout=900)
    assert r.returncode == 0, r.stderr
    assert out.exists()
    probe = _rscript('x <- readRDS(commandArgs(TRUE)[1]); '
                     'cat(class(x)[1], !is.null(S4Vectors::metadata(x)$peptides))',
                     str(out))
    cls, has_pep = probe.stdout.split()
    assert cls == "SummarizedExperiment"
    assert has_pep == "TRUE"


@needs_r("SummarizedExperiment", "S4Vectors")
def test_an_rscript_that_writes_nothing_is_reported_as_a_failure(ma, tmp_path):
    # symptom: an R script that stops inside a tryCatch, or saves to a path it
    # could not create, can still exit 0 — and the caller was told nothing.
    with pytest.raises(ma.StageError) as e:
        ma._run_rscript(["Rscript", "-e", "stop('boom')"], "a test", "hint")
    assert "Rscript exited" in str(e.value)
    assert "boom" in str(e.value)


def test_rscript_is_launched_by_the_path_that_was_resolved_for_it(
        ma, tmp_path, monkeypatch):
    # symptom: every tool is meant to be launched by the absolute path PATH
    # resolves to, because CreateProcess ignores PATHEXT and an Rscript.bat
    # earlier on PATH loses to an Rscript.exe later on it — but _run_rscript
    # handed the OS the bare name, so the R that ran need not be the one
    # have() reported and the logged command line named.
    d = tmp_path / "rbin"
    d.mkdir()
    (d / "Rscript").write_text("#!/usr/bin/env python3\n", encoding="utf-8")
    os.chmod(d / "Rscript", 0o755)
    if os.name == "nt":
        # PATHEXT decides what is executable there, so an extension-less stub
        # is never found — the same reason conftest's stub_bin writes a shim.
        (d / "Rscript.cmd").write_text(
            f'@echo off\r\n"{sys.executable}" "%~dp0Rscript" %*\r\n',
            encoding="utf-8")
    monkeypatch.setenv("PATH", str(d) + os.pathsep + os.environ["PATH"])
    # The spy is on Popen because that is what _run_rscript calls now, and it
    # would see the launch either way: subprocess.run calls Popen itself, so
    # nothing here distinguishes the two. THIS TEST IS ABOUT THE PATH AND ONLY
    # THE PATH, which is what its name says. What covers the spawn site is
    # test_rscript_goes_through_the_spawn_helper_so_a_kill_reaches_it below.
    launched = []
    real = ma.subprocess.Popen

    def spy(argv, *a, **k):
        launched.append([str(c) for c in argv])
        return real(argv, *a, **k)
    monkeypatch.setattr(ma.subprocess, "Popen", spy)
    ma._run_rscript(["Rscript", "-e", "invisible(NULL)"], "a test")
    assert launched, "no Rscript process was launched"
    assert os.path.isabs(launched[0][0]), \
        f"Rscript was launched as {launched[0][0]!r}, not the resolved path"
    assert os.path.dirname(launched[0][0]) == str(d)


def test_rscript_goes_through_the_spawn_helper_so_a_kill_reaches_it(
        ma, monkeypatch):
    """The Mac's whole exposure to the tool-group kill, actually pinned.

    `all` calls cmd_run, then cmd_report, then cmd_object in ONE process and
    signal handlers are never unregistered, so _release_lock_on_signal is
    still installed during the R phase: a `kill` there used to leave Rscript
    holding a 455k-row data frame in RAM with nothing able to stop it, because
    subprocess.run registers nothing and the handler exits through os._exit().

    THE TEST ABOVE DOES NOT COVER THAT, and it was cited as though it did.
    Measured: with metaannot.py reverted to HEAD -- _run_rscript still on
    subprocess.run -- and that test file unchanged, it PASSES. Its spy is on
    subprocess.Popen and subprocess.run calls Popen itself, so moving the spy
    from one to the other changed what the test watches and not what it
    proves.

    What actually changed at that call site is the SESSION and the
    REGISTRATION, so those are what this asserts. Only _tool_process claims a
    slot, so a claim for the launched pid is proof the launch went through it;
    start_new_session is what makes one killpg reach whatever R forks. A stub
    that exits immediately is enough: neither property is about R.
    """
    kw, claimed = {}, []
    real_popen = ma.subprocess.Popen
    real_claim = ma._claim_tool_slot

    def popen_spy(argv, *a, **k):
        kw.update(k)
        proc = real_popen(argv, *a, **k)
        kw["pid"] = proc.pid
        return proc

    def claim_spy(pgid):
        claimed.append(pgid)
        return real_claim(pgid)

    monkeypatch.setattr(ma.subprocess, "Popen", popen_spy)
    monkeypatch.setattr(ma, "_claim_tool_slot", claim_spy)
    ma._run_rscript([sys.executable, "-c", "pass"], "a test")
    assert kw, "no process was launched"
    if os.name == "nt":
        # No POSIX sessions there, so the helper passes none and claims no
        # slot - which is the gap the README and the CHANGELOG both state.
        assert "start_new_session" not in kw
        assert claimed == []
    else:
        assert kw.get("start_new_session") is True, (
            "Rscript was not launched in a session of its own, so one killpg "
            "cannot reach what R forks and the run's own group kill would "
            "have to signal itself to try")
        assert claimed == [kw["pid"]], (
            "no tool slot was claimed for the Rscript pid, so a `kill` during "
            f"the R phase of `all` cannot see it: claimed {claimed}, pid "
            f"{kw['pid']}")
    # stdin on /dev/null, for the reason the helper gives: a child in a
    # session of its own that reads the inherited terminal gets SIGTTIN and
    # STOPS, and a three-day run does not come back from that.
    assert kw.get("stdin") == subprocess.DEVNULL


@needs_r("SummarizedExperiment", "S4Vectors")
def test_the_object_identifies_samples_from_the_recorded_list(tmp_path):
    # the whole point of writing the file: the object must say it used it,
    # rather than falling back to the design or to column types.
    proj = build_project(tmp_path / "p")
    proj.run()
    run_metaannot("object", "--config", proj.config_path, "--no-run",
                  cwd=proj.root)
    r = subprocess.run(["Rscript", proj.rpath("analysis", "build_object.R"),
                        proj.results, proj.rpath("metaannot.rds"), "analysis"],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr
    assert "samples (quant/sample_columns.txt):" in r.stdout, r.stdout[-600:]
    assert "inferred from column types" not in r.stdout


@needs_r("SummarizedExperiment", "S4Vectors")
def test_the_object_needs_no_heuristic_when_there_is_no_manifest(tmp_path):
    # symptom this fixes: with no manifest there is no design_from_input.tsv,
    # so the sample columns used to be guessed. The recorded list covers it.
    proj = build_project(tmp_path / "p")
    proj.write_config(manifest="")
    proj.run()
    assert not os.path.exists(proj.rpath("quant", "design_from_input.tsv"))
    run_metaannot("object", "--config", proj.config_path, "--no-run",
                  cwd=proj.root)
    r = subprocess.run(["Rscript", proj.rpath("analysis", "build_object.R"),
                        proj.results, proj.rpath("metaannot.rds"), "analysis"],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr
    assert "samples (quant/sample_columns.txt):" in r.stdout
    assert "WARNING: no sample-column list" not in r.stdout


@needs_r(*R_CORE, *R_BIOC)
def test_the_report_identifies_samples_from_the_recorded_list(knitted):
    if not knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    html = open(knitted.rpath("analysis", "analyse_metaannot.html"),
                encoding="utf-8").read()
    # The rendered NOTE, not the chunk source: code_folding echoes the whole
    # document into the page, so any phrase that also appears in a comment
    # there would match whatever the report actually did. This line carries
    # the resolved path and count, which the source cannot.
    assert re.search(r"sample columns identified from "
                     r"quant/sample_columns\.txt: 8", html), \
        "the report did not use the recorded sample-column list"


# --- phase 3: the two filters, the printed design, ratio compression --
def _prep_filter_source(ma):
    """The report's own two-filter block, lifted out so a test can run it."""
    code = _chunk(ma, "prep")
    m = re.search(r"^keep_valid <- apply.*?^cat\(sprintf\(\"%d/%d groups "
                  r"retained by both filters.*?\n", code, re.S | re.M)
    assert m, "the prep chunk no longer contains the two-filter block"
    return m.group(0)


FILTER_PREAMBLE = """
gate <- function(fmt, ...) cat("GATE:", sprintf(fmt, ...), "\\n")
note <- function(fmt, ...) cat("NOTE:", sprintf(fmt, ...), "\\n")
a <- commandArgs(TRUE)
params <- list(min_valid_per_group = as.integer(a[1]),
               min_plexes = as.integer(a[2]),
               group_col_for_filtering = "group")
design_path0 <- "design_from_input.tsv"
IS_ISOBARIC <- as.logical(a[3])
# four samples per plex, two plexes, two conditions crossed over both
X <- matrix(NA_real_, nrow = 4, ncol = 8,
            dimnames = list(paste0("P", 1:4), paste0("s", 1:8)))
PLEX_OF <- setNames(rep(c("TMT1", "TMT2"), each = 4), colnames(X))
fgrp <- factor(rep(c("a", "a", "b", "b"), 2))
X["P1", ] <- 1:8                      # every sample, both plexes
X["P2", 1:4] <- 1:4                   # plex TMT1 only, but 2 per group
X["P3", c(1, 2, 5, 6)] <- 1           # both plexes, only group 'a'
X["P4", ] <- 1:8
if (!IS_ISOBARIC) PLEX_OF <- NULL
"""


def _rscript_file(tmp_path, code, *args):
    """Run R code from a FILE, not from `Rscript -e`.

    `-e` takes the program as a command-line argument, and on Windows a
    program this long crashes Rscript with an access violation before it runs
    a line. A file has no such limit and behaves the same everywhere.
    """
    path = tmp_path / "snippet.R"
    path.write_text(code, encoding="utf-8")
    return subprocess.run(["Rscript", str(path), *[str(a) for a in args]],
                          capture_output=True, text=True, timeout=900)


@needs_r()
@pytest.mark.parametrize("mv,mp,iso,expect", [
    # min_valid 2, no plex filter: P3 fails the group count, P2 survives on a
    # count one batch supplied on its own - which is the thing being reported
    (2, 1, True, ["min_valid_per_group >= 2 in every level of group: removes 1 of 4",
                  "min_plexes >= 1: removes 0 of 4",
                  "GATE: 1 protein(s) pass min_valid_per_group but are "
                  "quantified in a single plex",
                  "3/4 groups retained by both filters"]),
    # min_plexes 2 as well: P2 goes too, and each filter is counted against
    # the same starting set rather than in sequence
    (2, 2, True, ["min_valid_per_group >= 2 in every level of group: removes 1 of 4",
                  "min_plexes >= 2: removes 1 of 4, 1 of which "
                  "min_valid_per_group would have kept",
                  "2/4 groups retained by both filters"]),
    # label-free: there are no plexes and the filter says so instead of
    # silently passing everything
    (2, 1, False, ["min_plexes: not applied (no per-sample plex; this is not "
                   "an isobaric run)",
                   "3/4 groups retained by both filters"]),
])
def test_the_report_applies_both_filters_and_reports_them_separately(
        ma, tmp_path, mv, mp, iso, expect):
    # symptom: min_valid_per_group counts SAMPLES, and an isobaric run's
    # missingness is shaped by the plex, so "3 valid values in every group"
    # can be satisfied entirely inside one batch and the difference the model
    # then reports is that batch.
    r = _rscript_file(tmp_path, FILTER_PREAMBLE + _prep_filter_source(ma),
                      mv, mp, "TRUE" if iso else "FALSE")
    out = r.stdout + r.stderr
    for want in expect:
        assert want in out, out


@needs_r()
def test_min_plexes_without_a_plex_stops_instead_of_passing_everything(
        ma, tmp_path):
    # a filter the user asked for that cannot be applied must not report
    # itself as satisfied.
    r = _rscript_file(tmp_path, FILTER_PREAMBLE + _prep_filter_source(ma),
                      2, 2, "FALSE")
    out = r.stdout + r.stderr
    assert "min_plexes is 2" in out and "not isobaric" in out


def test_the_report_prints_the_design_it_used_and_the_reference_treatment(ma):
    # symptom: the log scrolls away, and a report that does not say which
    # reference treatment produced its numbers cannot be checked afterwards.
    code = _chunk(ma, "isobaric-design")
    assert "DESIGN_NOTES" in code and "readLines(DESIGN_NOTES" in code
    assert "reference treatment" in code
    assert "COVARIATE" in code and "never a condition" in code
    assert "samples per plex" in code and "condition x plex" in code


def test_the_isobaric_flag_comes_from_the_file_only_the_tmt_reader_writes(ma):
    # design_notes.txt is written by the isobaric reader and by nothing else,
    # so a label-free run cannot accidentally take the isobaric branches.
    code = _chunk(ma, "read-quant")
    assert 'DESIGN_NOTES <- file.path(RD, "quant", "design_notes.txt")' in code
    assert "IS_ISOBARIC  <- file.exists(DESIGN_NOTES)" in code


def test_ratio_compression_is_stated_as_a_limitation_and_not_corrected(ma):
    code = _chunk(ma, "caveat-ratio-compression")
    assert "if (IS_ISOBARIC)" in code            # label-free must not see it
    assert "Ratio compression" in code
    assert "LOWER BOUNDS" in code
    assert "does **not** correct" in code
    # and nothing anywhere claims to undo it
    assert "compression_correction" not in ma.RMD_TEMPLATE
    assert "correct_compression" not in ma.RMD_TEMPLATE


def test_min_plexes_is_recorded_beside_the_numbers_it_produced(ma):
    code = _chunk(ma, "export")
    assert 'paste("min_plexes:    ", params$min_plexes' in code
    assert "not an isobaric run; not applied" in code


ISOBARIC_PREAMBLE = """
suppressPackageStartupMessages(library(readr))
note <- function(fmt, ...) cat("NOTE:", sprintf(fmt, ...), "\\n")
read_tsv_full <- function(p, ...) readr::read_tsv(p, show_col_types = FALSE)
a <- commandArgs(TRUE)
RD <- a[1]
meta_path <- file.path(RD, "meta.tsv")
design_path0 <- file.path(RD, "quant", "design_from_input.tsv")
DESIGN_NOTES <- file.path(RD, "quant", "design_notes.txt")
IS_ISOBARIC <- file.exists(DESIGN_NOTES)
meta <- as.data.frame(read_tsv_full(meta_path))
meta$.sample <- as.character(meta$sample)
meta$.col <- meta$.sample
"""


@needs_r("readr")
@pytest.mark.parametrize("with_plex", [True, False])
def test_the_isobaric_design_chunk_runs_and_prints_what_was_done(
        ma, tmp_path, with_plex):
    # symptom: a chunk that only parses can still die at knit time on a
    # missing column, and this one is the record of every choice the numbers
    # were made under.
    q = tmp_path / "quant"
    q.mkdir()
    (q / "design_notes.txt").write_text(
        "input:            FragPipe TMT, 2 plex(es), 4 sample column(s)\n"
        "condition source: analysis.metadata\n"
        "within-plex norm: median centring per channel\n"
        "reference:        covariate (dropped from the design; plex stays "
        "in the model)\n"
        "reference channel: TMT1=131C/Pool01, TMT2=131N/Pool02\n",
        encoding="utf-8")
    (q / "design_from_input.tsv").write_text(
        "sample\tplex\tchannel\ns1\tTMT1\t126\ns2\tTMT1\t127N\n"
        "s3\tTMT2\t126\ns4\tTMT2\t127N\n", encoding="utf-8")
    cols = "sample\tgroup\tplex\n" if with_plex else "sample\tgroup\n"
    rows = [("s1", "a", "TMT1"), ("s2", "b", "TMT1"),
            ("s3", "a", "TMT2"), ("s4", "b", "TMT2")]
    (tmp_path / "meta.tsv").write_text(
        cols + "".join(("\t".join(r if with_plex else r[:2])) + "\n"
                       for r in rows), encoding="utf-8")
    r = _rscript_file(tmp_path, ISOBARIC_PREAMBLE + _chunk(ma, "isobaric-design"),
                      str(tmp_path))
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    # the record, verbatim, including which reference treatment was applied
    assert "reference:        covariate" in out
    assert "NOTE: reference treatment: covariate" in out
    assert "within-plex norm: median centring per channel" in out
    # the design itself, and the plex named as what it is
    assert "samples per plex" in out and "TMT1" in out and "TMT2" in out
    assert "condition x plex" in out
    assert "COVARIATE here: a TMT batch, and never a condition" in out


@needs_r("readr")
def test_the_isobaric_chunk_is_inert_without_design_notes(ma, tmp_path):
    # design_notes.txt is written by the isobaric reader and by nothing else,
    # so a label-free run must take the other branch rather than fail on a
    # missing file.
    (tmp_path / "quant").mkdir()
    (tmp_path / "meta.tsv").write_text("sample\tgroup\ns1\ta\ns2\tb\n",
                                       encoding="utf-8")
    r = _rscript_file(tmp_path, ISOBARIC_PREAMBLE + _chunk(ma, "isobaric-design"),
                      str(tmp_path))
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    assert "Not an isobaric run" in out


# --- a whole isobaric run, knitted ------------------------------------
@pytest.fixture(scope="module")
def tmt_knitted(tmp_path_factory):
    """One finished TMT run, knitted once.

    Two plexes with the condition crossed over both, a reference channel in
    each, and a fifth of the proteins seen in TMT1 only - which is the
    plex-shaped missingness min_plexes exists for, and the thing
    min_valid_per_group cannot see.
    """
    import random
    root = tmp_path_factory.mktemp("tmtknit") / "p"
    os.makedirs(os.path.join(str(root), "input"), exist_ok=True)
    proteins = F.protein_set(n_extra=40)
    faa = F.write_fasta(os.path.join(str(root), "input", "proteins.faa"),
                        proteins)
    emp = F.write_emapper(
        os.path.join(str(root), "input", "cat.emapper.annotations"), proteins)
    plexes = {
        "TMT1": [("126", "s1"), ("127N", "s2"), ("128N", "s3"),
                 ("129N", "s4"), ("131C", "Pool01")],
        "TMT2": [("126", "s5"), ("127N", "s6"), ("128N", "s7"),
                 ("129N", "s8"), ("131N", "Pool02")],
    }
    group = {"s1": "a", "s2": "a", "s3": "b", "s4": "b",
             "s5": "a", "s6": "a", "s7": "b", "s8": "b"}
    run_dir = os.path.join(str(root), "input", "run")
    rng = random.Random(3)
    for pi, (plex, chans) in enumerate(plexes.items()):
        rows = []
        for k, p in enumerate(proteins):
            if plex != "TMT1" and k % 5 == 0:
                continue                     # confined to the first plex
            for j in range(2):
                base = 20000 * (1 + rng.random())
                rows.append({
                    "peptide": f"{p.pid}PEP{j}K".upper().replace("_", ""),
                    "razor": p.pid,
                    # (1 + pi) is a plain plex effect for the batch term to
                    # absorb; the lift is the only real difference
                    "values": {s: round(base * (1 + pi)
                                        * (2.0 if (group.get(s) == "b"
                                                   and k % 7 == 0) else 1.0)
                                        * (0.8 + 0.4 * rng.random()))
                               for _c, s in chans}})
        F.write_tmt_plex(run_dir, plex, chans, rows, seed=7 + pi)
    meta = os.path.join(str(root), "input", "metadata.tsv")
    with open(meta, "w", encoding="utf-8") as fh:
        fh.write("sample\tplex\tgroup\n")
        for plex, chans in plexes.items():
            for _c, s in chans:
                if not s.startswith("Pool"):
                    fh.write(f"{s}\t{plex}\t{group[s]}\n")
    proj = build_project(root, proteins=proteins)
    proj.write_config(quant_table=run_dir, quant_format="fragpipe_tmt",
                      proteins_faa=faa, emapper_precomputed=[emp],
                      manifest="",
                      tmt=dict(reference_name="Pool*",
                               condition_from_name=""),
                      analysis=dict(proj.cfg.get("analysis") or {},
                                    metadata=meta, min_valid_per_group=2,
                                    min_plexes=2, fdr=0.2, min_lfc=0.2))
    proj.run()
    proj.rendered = False
    if shutil.which("Rscript") and shutil.which("pandoc") and r_has(*R_CORE):
        run_metaannot("report", "--config", proj.config_path, cwd=proj.root,
                      timeout=1800)
        proj.rendered = True
    return proj


@needs_r(*R_CORE)
def test_an_isobaric_run_knits_and_reports_both_filters(tmt_knitted):
    # symptom: everything about the TMT path had only ever been checked one
    # piece at a time; a chunk that parses can still die at knit time.
    if not tmt_knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    txt = open(tmt_knitted.rpath("analysis", "analyse_metaannot.html"),
               encoding="utf-8").read()
    # the design it used, and what was done to get there
    assert "The isobaric design" in txt
    assert "reference treatment: covariate" in txt
    assert "samples per plex" in txt and "condition x plex" in txt
    assert "COVARIATE here: a TMT batch, and never a condition" in txt
    # both filters, separately, and the 10 plex-confined proteins are exactly
    # the ones the sample count could not see
    assert "min_valid_per_group &gt;= 2 in every level of group: removes 0 of 50" in txt
    assert ("min_plexes &gt;= 2: removes 10 of 50, 10 of which "
            "min_valid_per_group would have kept") in txt
    assert "proteins by number of plexes" in txt
    assert "40/50 groups retained by both filters" in txt
    # and the limitation that is stated rather than corrected
    assert "Ratio compression" in txt and "LOWER BOUNDS" in txt


@needs_r(*R_CORE)
def test_the_isobaric_run_records_its_choices_beside_the_numbers(tmt_knitted):
    if not tmt_knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    rec = open(tmt_knitted.rpath("analysis", "design_record.txt"),
               encoding="utf-8").read()
    assert "design_formula: ~ 0 + group + plex" in rec
    assert "plex:           in the model as a batch term" in rec
    assert "within-plex norm: median centring per channel" in rec
    assert "reference:        covariate" in rec
    assert "reference channel: TMT1=131C/Pool01, TMT2=131N/Pool02" in rec
    assert "min_plexes:     2 (protein level, counted across plexes)" in rec


# --- phase 4: the planted effect, and the reference in the object ------
def _planted_project(tmp_path, name="p", n_proteins=8, **over):
    """A whole project over a two-plex run with a planted effect."""
    proteins = F.protein_set()
    root, truth = F.tmt_planted_run(str(tmp_path / "run"), proteins=proteins,
                                    n_proteins=n_proteins,
                                    layout=(("a", "a", "b", "b"),
                                            ("a", "a", "b", "b")))
    cfg = dict(quant_table=root, quant_format="fragpipe_tmt", manifest="",
               tmt={"reference_name": "Pool*"})
    cfg.update(over)
    proj = build_project(tmp_path / name, proteins=proteins, **cfg)
    proj.truth = truth
    return proj


LIMMA_RECOVERY = """
suppressPackageStartupMessages(library(limma))
a <- commandArgs(TRUE)
X <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
md <- read.delim(a[2], stringsAsFactors = TRUE)
md <- md[match(colnames(X), as.character(md$sample)), ]
group <- factor(md$group); plex <- factor(md$plex)
fit_one <- function(design, coef) {
  f <- lmFit(X, design)
  cm <- makeContrasts(contrasts = coef, levels = design)
  topTable(eBayes(contrasts.fit(f, cm)), number = Inf, sort.by = "none")
}
with_plex <- fit_one(model.matrix(~ 0 + group + plex), "groupb - groupa")
without   <- fit_one(model.matrix(~ 0 + group), "groupb - groupa")
out <- data.frame(protein = rownames(X),
                  with_plex = with_plex$logFC, without = without$logFC)
write.table(out, a[3], sep = "\\t", quote = FALSE, row.names = FALSE)
"""


@needs_r("limma")
def test_limma_recovers_the_planted_effect_only_with_the_plex_in_the_model(
        ma, tmp_path):
    # symptom: the report's default formula for an isobaric run gained
    # '+ plex', and the case for it cannot be made by reading the formula.
    # This runs the real model, in limma, on a matrix with a KNOWN 2x effect
    # and a KNOWN 3x plex loading that are deliberately not orthogonal, and
    # asks what each formula reports.
    root, truth = F.tmt_planted_run(str(tmp_path / "run"))
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg.update(quant_table=root, quant_format="fragpipe_tmt")
    cfg["tmt"]["reference_name"] = "Pool*"
    feats, int_cols, design = ma.read_feature_table(root, "fragpipe_tmt", cfg)
    prot = ma.rollup_features(feats, int_cols, {}, "razor", 0)[0]

    import numpy as np
    mat = tmp_path / "prot.tsv"
    out = prot[["group_id"] + int_cols].copy()
    out[int_cols] = np.log2(out[int_cols].astype(float))
    out.to_csv(mat, sep="\t", index=False)
    md = tmp_path / "design.tsv"
    design.to_csv(md, sep="\t", index=False)
    res = tmp_path / "fits.tsv"
    r = _rscript_file(tmp_path, LIMMA_RECOVERY, mat, md, res)
    assert r.returncode == 0, r.stdout + r.stderr
    fit = pd.read_csv(res, sep="\t").set_index("protein")

    up = [p for p in truth.regulated if truth.effect_of[p] > 0]
    down = [p for p in truth.regulated if truth.effect_of[p] < 0]
    # the model the report writes recovers +1 and -1 log2, as planted
    assert fit.loc[up, "with_plex"].mean() == pytest.approx(1.0, abs=0.15)
    assert fit.loc[down, "with_plex"].mean() == pytest.approx(-1.0, abs=0.15)
    assert abs(fit.loc[sorted(truth.null), "with_plex"]).max() < 0.15
    # and without the plex term every protein carries half the batch, so the
    # proteins that truly went DOWN are reported as barely moving and the
    # ones that did not move at all are reported as up
    assert fit.loc[sorted(truth.null), "without"].mean() == pytest.approx(
        0.5 * truth.plex_log2["TMT2"], abs=0.15)
    assert fit.loc[down, "without"].mean() > -0.4


@pytest.fixture(scope="module")
def tmt_object(tmp_path_factory):
    """One finished TMT run with the R object built from it."""
    root = tmp_path_factory.mktemp("tmtobj")
    proj = _planted_project(root, "p")
    proj.run()
    proj.built = False
    if shutil.which("Rscript") and r_has(*R_BIOC):
        run_metaannot("object", "--config", proj.config_path, cwd=proj.root,
                      timeout=900)
        proj.built = os.path.exists(proj.rpath("metaannot.rds"))
    return proj


@needs_r("SummarizedExperiment", "S4Vectors")
def test_the_tmt_reference_channel_is_in_neither_the_design_nor_the_coldata(
        tmt_object):
    # symptom: colData is what every downstream model reads its samples from,
    # so a pooled bridge that reached it would acquire a condition, join a
    # group's mean and be modelled as biology. Under both reference
    # treatments the pool stops being a sample in the reader; this is the
    # check that nothing put it back.
    if not tmt_object.built:
        pytest.skip("the object was not built")
    r = _rscript('x <- readRDS(commandArgs(TRUE)[1]); '
                 'cd <- SummarizedExperiment::colData(x); '
                 'cat(paste(rownames(cd), collapse=","), "|", '
                 'paste(colnames(cd), collapse=","))',
                 tmt_object.rpath("metaannot.rds"))
    assert r.returncode == 0, r.stderr
    rows, cols = [s.strip() for s in r.stdout.split("|")]
    samples = rows.split(",")
    assert sorted(samples) == sorted(tmt_object.truth.samples)
    assert not any(s.startswith("Pool") for s in samples)
    # the plex IS there, because the model needs it as a batch term
    assert "plex" in cols.split(",")


@needs_r("SummarizedExperiment", "S4Vectors", "QFeatures")
def test_the_tmt_object_links_its_peptide_assay_to_the_proteins(tmt_object):
    # the peptide assay is the layer that the tmt-report matrices would have
    # deleted, and the reason this reads the per-plex tables at all. It has
    # to be present and linked on isobaric input exactly as on label-free.
    if not tmt_object.built:
        pytest.skip("the object was not built")
    if not r_has("QFeatures"):
        pytest.skip("QFeatures is not installed")
    r = _rscript('suppressPackageStartupMessages(library(QFeatures)); '
                 'q <- readRDS(commandArgs(TRUE)[1]); '
                 'al <- assayLink(q, "proteins"); '
                 'cat(class(q)[1], paste(names(q), collapse=","), al@fcol, '
                 'length(al@hits), nrow(q[["peptides"]]), '
                 'ncol(q[["peptides"]]))',
                 tmt_object.rpath("metaannot.rds"))
    assert r.returncode == 0, r.stderr
    cls, assays, fcol, hits, npep, ncol = r.stdout.split()
    assert cls == "QFeatures"
    assert assays == "peptides,proteins"
    assert fcol == "assigned_protein"          # the assignment metaannot made
    assert int(hits) > 0
    # two peptides per protein in the fixture, and the pools are not columns
    assert int(npep) == 2 * len(tmt_object.truth.null
                                | tmt_object.truth.regulated)
    assert int(ncol) == len(tmt_object.truth.samples)


# ----------------------------------------------------------------------
# issue #5: a bin that contributes nothing to the model has to say so
# ----------------------------------------------------------------------
def _coverage_source(ma):
    """The report's coverage block, with the constants it is written against.

    Returned as (constants, block): the preamble between them builds its
    fixture out of BIN_LEVELS, so the two cannot simply be concatenated here.
    Lifted out of the document rather than retyped, so a threshold or a
    sentence that changes there changes here too. The constants come from the
    setup chunk because the rule is meaningless without them: BIN_KOLESS
    decides which bins the headline is about, and the two COVERAGE_ values are
    the whole of when it is raised.

    Only BIN_LEVELS has to be there, and it is older than this feature. The
    coverage constants and the block come back EMPTY when they are absent,
    deliberately: a case whose whole assertion is that a line was not printed
    has to fail on the report going silent, and with an assertion here it
    would instead fail on a regex in the test file - which is a statement
    about this harness and not about the report. Empty, the script still runs,
    prints nothing, and every case fails on the retention table it did not
    find.
    """
    setup = _chunk(ma, "setup")
    levels = re.search(r"^BIN_LEVELS <- c\(.*?\)$", setup, re.S | re.M)
    assert levels, "the setup chunk no longer defines BIN_LEVELS"
    # `\d+(?:\.\d+)?`, not `\d+`: a threshold edited to 0.5 made this regex
    # miss, which returned EMPTY constants, which failed 17 of the cases below
    # with `object 'BIN_KOLESS' not found` and told the reader the escalation
    # was not in the document. A harness that blames the report for a change to
    # a number is worse than no harness.
    consts = re.search(r"^BIN_KOLESS <- .*?^COVERAGE_MIN_PCT <- \d+(?:\.\d+)?$",
                       setup, re.S | re.M)
    block = re.search(r"^TESTED <- aq\$group_id .*", _chunk(ma, "prep"),
                      re.S | re.M)
    return (levels.group(0) + "\n" + (consts.group(0) if consts else ""),
            block.group(0) if block else "")


# aq, keep and X as the document has them by the time the block runs: one row
# per quantified group, `keep` the two retention filters, and X already cut
# down to the rows the model is fitted on. The first two arguments are the
# per-bin quantified and tested counts in BIN_LEVELS order, which is exactly
# the shape of the `retention by bin` table the issue quotes. The third is the
# same pair for the quantified groups with NO annotation row - no bin, so
# dropped from every per-bin number and still fitted - and the fourth is
# min_plexes and the number of plexes, off unless a case asks for them.
COVERAGE_PREAMBLE = """
suppressPackageStartupMessages({library(dplyr); library(tibble)})
gate <- function(fmt, ...) cat("GATE:", sprintf(fmt, ...), "\\n")
note <- function(fmt, ...) cat("NOTE:", sprintf(fmt, ...), "\\n")
a <- commandArgs(TRUE)
nq <- as.integer(strsplit(a[1], ",")[[1]])
nt <- as.integer(strsplit(a[2], ",")[[1]])
nna <- if (length(a) >= 3) as.integer(strsplit(a[3], ",")[[1]]) else c(0L, 0L)
px <- if (length(a) >= 4) as.integer(strsplit(a[4], ",")[[1]]) else c(1L, 0L)
stopifnot(length(nq) == length(BIN_LEVELS), all(nt <= nq), nna[2] <= nna[1])
params <- list(min_valid_per_group = 3L, min_features = 0L,
               min_plexes = px[1], group_col_for_filtering = "group")
fgrp <- factor(rep(c("a", "b"), each = 18))
pbatch <- if (px[2] > 0) rep(paste0("P", seq_len(px[2])), length.out = 36) else
          NULL
ntot <- sum(nq) + nna[1]
aq <- tibble(group_id = if (ntot) paste0("g", seq_len(ntot)) else character(0),
             bin = factor(c(rep(BIN_LEVELS, nq), rep(NA_character_, nna[1])),
                          levels = BIN_LEVELS))
keep <- c(unlist(Map(function(q, t) c(rep(TRUE, t), rep(FALSE, q - t)), nq, nt),
                 use.names = FALSE),
          rep(TRUE, nna[2]), rep(FALSE, nna[1] - nna[2]))
X <- matrix(0, nrow = sum(keep), ncol = 36,
            dimnames = list(aq$group_id[keep], paste0("s", seq_len(36))))
"""

# The issue's own run, bin by bin: 16,070 quantified protein groups, 212 in
# the model, and not one of them from 3d_duf_only or 4_dark.
ISSUE_5_QUANTIFIED = "7023,5153,2914,170,0,0,810"
ISSUE_5_TESTED = "26,175,11,0,0,0,0"


def _loose(n):
    """`n` as the printed table may break it up, thin space, comma or not."""
    return r"[\s,]*".join(str(n))


@needs_r("dplyr", "tibble")
@pytest.mark.parametrize("quantified,tested,unbinned,plex,expect,forbid", [
    # The run in the issue. Two bins quantified in the thousands contribute
    # nothing, and the KO-less fraction as a whole is in the model only just.
    (ISSUE_5_QUANTIFIED, ISSUE_5_TESTED, "0,0", "1,0",
     ["GATE: 0/170 3d_duf_only and 0/810 4_dark quantified group(s) reach "
      "the model",
      "not the effector shortlist - is about them",
      "GATE: the statistics below cover 0.28% of the KO-less groups this run "
      "quantified: 11 of 3894 are tested, and they are 11 of the 212 "
      "group(s) in the model",
      "Two knobs decide it",
      "analysis.min_valid_per_group is 3",
      "min_features_per_protein"],
     # the two bins that held nothing are not named anywhere: a gate that
     # fires on an empty bin every run is the noise this rule exists to avoid
     ["3s_structure_only", "3p_profile_only",
      # and the knob is named where the config really keeps it: `join` is a
      # stage name, so a reader who writes the block this implies is told
      # "unrecognised key 'join'" and their setting is silently ignored
      "join.min_features_per_protein"]),
    # The same run with the quantified groups that have no annotation row at
    # all put back, every one of them in the model, because an unbinned group
    # is a well-covered one and passes the filter the sparse dark ones fail.
    # They are not a bin and are rightly absent from every per-bin number -
    # but the model they are fitted in is 48 and not 34, and the fraction is
    # read against the model.
    (ISSUE_5_QUANTIFIED, "20,13,1,0,0,0,0", "14,14", "1,0",
     ["they are 1 of the 48 group(s) in the model"],
     ["1 of the 34 group(s)"]),
    # A healthy run: every bin keeps most of what it brought, and the block
    # says nothing at all.
    ("31,12,14,11,10,0,15", "30,12,13,10,9,0,14", "0,0", "1,0",
     [], ["GATE:", "NOTE:"]),
    # 3p_profile_only as it is on every run ever made - fed only by hhblits
    # and jackhmmer, so quantified-nothing, tested-nothing. Silence.
    ("31,12,14,11,10,0,15", "30,12,13,10,9,0,14", "0,0", "1,0",
     [], ["3p_profile_only"]),
    # Below the denominator floor the same zero is a NOTE, not a GATE: "none
    # of 7" is an anecdote about a handful of proteins.
    ("400,300,200,0,0,0,7", "380,290,190,0,0,0,0", "0,0", "1,0",
     ["NOTE: 0/7 4_dark reach the model", "fewer than 10 quantified group(s)"],
     ["GATE:"]),
    # Every KO-less bin at zero, with enough behind it to be a population:
    # the headline sentence the issue asked for, and it names only the bins
    # that actually held proteins.
    ("400,300,0,0,0,0,60", "380,290,0,0,0,0,0", "0,0", "1,0",
     ["GATE: 0/60 4_dark quantified group(s) reach the model",
      "is about it.",
      "GATE: no KO-less protein is tested in any contrast: 0 of 60 "
      "quantified KO-less group(s), in 4_dark"],
     ["3_annotated_no_ko"]),
    # The same run isobaric, with min_plexes in force. `keep` is keep_valid &
    # keep_plex, so a third filter is behind the counts the NOTE explains, and
    # min_valid_per_group is not on its own what they measure.
    ("400,300,0,0,0,0,60", "380,290,0,0,0,0,0", "0,0", "2,3",
     ["Three knobs decide it",
      "analysis.min_plexes is 2, over 3 plex(es); the counts above measure "
      "both"],
     ["Two knobs", "and it is what the counts above measure"]),
    # A bin that goes to zero is gated whether or not it is KO-less: nothing
    # is tested at all here, and that is not a KO-less-only problem.
    ("120,300,200,0,0,0,190", "0,290,190,0,0,0,180", "0,0", "1,0",
     ["GATE: 0/120 1_ko_pathway quantified group(s) reach the model"], []),
])
def test_a_bin_filtered_out_is_escalated_but_an_empty_one_is_not(
        ma, tmp_path, quantified, tested, unbinned, plex, expect, forbid):
    # symptom: the report printed `retention by bin` with a 0 in it and
    # nothing escalated it, while a GATE was raised for findings orders of
    # magnitude smaller. The zeros of a bin that HELD nothing and of a bin
    # that lost everything looked identical.
    consts, block = _coverage_source(ma)
    # the constants first, because the preamble builds its bins out of
    # BIN_LEVELS, and the block last, exactly as the document orders them
    r = _rscript_file(tmp_path, consts + COVERAGE_PREAMBLE + block,
                      quantified, tested, unbinned, plex)
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    # First, that the block RAN. Every case below asserts some line was not
    # printed, and a block that printed nothing at all satisfies all of them
    # - so the table it always prints is checked before anything else, with
    # the counts each row is supposed to carry.
    for col in ("bin", "n_kept", "n_tested", "pct_tested"):
        assert col in out, f"the retention table has no {col} column:\n{out}"
    for b, q, t in zip(_BIN_LEVELS, quantified.split(","), tested.split(",")):
        if int(q) == 0:      # dropped by group_by(); its absence is the rule
            continue
        row = rf"^\s*\d+\s+{b}\s+{_loose(q)}\s+{_loose(t)}\s+{_loose(t)}\s"
        assert re.search(row, out, re.M), f"no {b} row of {q}/{t}:\n{out}"
    nq_na, nt_na = unbinned.split(",")
    if int(nq_na):
        row = (rf"^\s*\d+\s+(?:NA|<NA>)\s+{_loose(nq_na)}\s+{_loose(nt_na)}"
               rf"\s+{_loose(nt_na)}\s")
        assert re.search(row, out, re.M), f"no unbinned row:\n{out}"
    for want in expect:
        assert want in out, out
    for no in forbid:
        assert no not in out, out


# The bins in the order the document lists them, for the row check above.
_BIN_LEVELS = ("1_ko_pathway", "2_ko_orphan", "3_annotated_no_ko",
               "3d_duf_only", "3s_structure_only", "3p_profile_only", "4_dark")


def test_the_bin_levels_this_file_checks_rows_against_are_the_documents(ma):
    # a literal list here would go stale the first time a bin is added, and
    # the row check above would then silently stop checking that bin.
    consts, _ = _coverage_source(ma)
    assert re.findall(r'"([0-9a-z_]+)"', consts.split("BIN_KOLESS")[0]) \
        == list(_BIN_LEVELS)


def test_the_report_names_config_keys_where_the_config_keeps_them(ma):
    # symptom: the coverage NOTE told the reader to set
    # `join.min_features_per_protein`. The key is TOP-LEVEL; `join` is a stage
    # name and `run.join` is a boolean, so the config that line describes is
    # refused with "unrecognised key 'join' - did you mean 'join'?" and the
    # setting is silently not in effect. Every dotted config path the report
    # prints is checked against DEFAULT_CONFIG, in both directions.
    tpl = ma.RMD_TEMPLATE
    top = {k for k, v in ma.DEFAULT_CONFIG.items() if not isinstance(v, dict)}
    sections = {k: v for k, v in ma.DEFAULT_CONFIG.items()
                if isinstance(v, dict)}
    for m in re.finditer(r"(?<![\w.$])([a-z_]{2,})\.([a-z_]{4,})\b", tpl):
        prefix, key = m.groups()
        if key in top:
            raise AssertionError(
                f"the report writes {prefix}.{key}, but {key} is a top-level "
                f"config key: a reader who writes that block is told the "
                f"'{prefix}' key is unrecognised and their setting is "
                f"silently ignored")
        if prefix in sections:
            assert key in sections[prefix], \
                f"the report names {prefix}.{key}, which {prefix} has no key for"
    # and the one this was found on, named in full rather than left to the
    # scan: it is in the document, and it is in the config, unprefixed
    assert "min_features_per_protein decided which proteins" in tpl
    assert "min_features_per_protein" in top


def _ratio_koless_source(ma):
    """The ratio model's KO-less escalation, lifted whole.

    Empty when it is absent, for the reason _coverage_source is: the cases
    that assert silence must fail on a silent report rather than on a regex.
    """
    block = re.search(r"^    if \(sum\(usable\) > 0 .*?^    \}$",
                      _chunk(ma, "taxon-adjusted"), re.S | re.M)
    return block.group(0) if block else ""


RATIO_PREAMBLE = """
suppressPackageStartupMessages({library(dplyr); library(tibble)})
gate <- function(fmt, ...) cat("GATE:", sprintf(fmt, ...), "\\n")
note <- function(fmt, ...) cat("NOTE:", sprintf(fmt, ...), "\\n")
a <- as.integer(commandArgs(TRUE))
# a[1] KO-less groups in the model, a[2] of them with a usable taxon, a[3]
# mapped groups with one. aqk is the retained set, `usable` its taxon filter.
aqk <- tibble(bin = factor(c(rep("4_dark", a[1]), rep("1_ko_pathway", a[3])),
                           levels = BIN_LEVELS))
usable <- c(rep(TRUE, a[2]), rep(FALSE, a[1] - a[2]), rep(TRUE, a[3]))
KOLESS_TESTED <- a[1]
"""

# The input the escalation is loudest on, used to prove the block is there
# when a case's whole assertion is that it said nothing: a report that CANNOT
# speak satisfies every such case, and nothing else in a silent run would tell
# the two apart.
RATIO_LOUD = (30, 0, 30)
RATIO_LOUD_SAYS = "GATE: and none of them is KO-less: all 30 KO-less group(s)"


@needs_r("dplyr", "tibble")
@pytest.mark.parametrize("koless,koless_usable,mapped,expect,forbid", [
    # The healthy shape, and the one this was found on: KO-less proteins ARE
    # in the model, their taxa are too small for a size factor, and there are
    # seven of them. A GATE on that is a GATE on a document with nothing
    # wrong with it.
    (7, 0, 30, ["NOTE: and none of them is KO-less: the 7 KO-less group(s) in "
                "the model are without a usable taxon",
                "fewer than 10 group(s)"], ["GATE:"]),
    # The same zero over a population: now it is the model's whole argument
    # that has nothing to stand on, and it is loud.
    (30, 0, 30, ["GATE: and none of them is KO-less: all 30 KO-less group(s) "
                 "in the model are without a usable taxon",
                 "the population it exists for"], ["NOTE:"]),
    # A different cause with the same symptom: nothing KO-less reached the
    # model at all. That is the coverage failure the retention block reports,
    # it is not a taxonomy limit, and sending the reader to
    # taxon_min_proteins would send them to the wrong knob.
    (0, 0, 30, ["NOTE: and none of them is KO-less, because no KO-less group "
                "reached the model at all",
                "not a limit of the taxonomy"],
     ["GATE:", "without a usable taxon"]),
    # KO-less proteins in the usable set: nothing to say.
    (30, 3, 30, [], ["GATE:", "NOTE:"]),
    # No usable taxon for anything. The count above it already says so, and
    # the block does not speak.
    (0, 0, 0, [], ["GATE:", "NOTE:"]),
])
def test_the_ratio_model_escalates_only_over_a_koless_population(
        ma, tmp_path, koless, koless_usable, mapped, expect, forbid):
    # symptom: the GATE fired whenever no KO-less protein had a usable taxon,
    # with no floor under the population and no way to tell "none reached the
    # model" from "their taxa are too small" - so it fired on a knit with
    # every bin fully retained, which is the noise the rule above it exists
    # to keep out of the report.
    consts, _ = _coverage_source(ma)
    code = consts + RATIO_PREAMBLE + _ratio_koless_source(ma)
    r = _rscript_file(tmp_path, code, koless, koless_usable, mapped)
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    if not expect:
        # this case says only that nothing was printed, which a document with
        # no escalation in it also satisfies. So the same block is run on the
        # input it is loudest on, and the silence counts only once it has
        # been shown to be a choice.
        loud = _rscript_file(tmp_path, code, *RATIO_LOUD)
        assert RATIO_LOUD_SAYS in loud.stdout + loud.stderr, \
            "the escalation is not in the document, so silence proves nothing"
    for want in expect:
        assert want in out, out
    for no in forbid:
        assert no not in out, out


def _shortlist_source(ma):
    """The empty-shortlist branch, lifted whole. Empty when it is absent."""
    block = re.search(r"^if \(nrow\(short\) == 0 && KOLESS_N > 0\) \{.*?^\}$",
                      _chunk(ma, "shortlist"), re.S | re.M)
    return block.group(0) if block else ""


# KOLESS_N, KOLESS_TESTED and a shortlist of the given length: base R only,
# because the branch is arithmetic and three sentences.
SHORTLIST_PREAMBLE = """
gate <- function(fmt, ...) cat("GATE:", sprintf(fmt, ...), "\\n")
note <- function(fmt, ...) cat("NOTE:", sprintf(fmt, ...), "\\n")
a <- as.integer(commandArgs(TRUE))
KOLESS_N <- a[1]
KOLESS_TESTED <- a[2]
short <- data.frame(group_id = seq_len(a[3]))
"""

SHORTLIST_LOUD = (19, 0, 0)
SHORTLIST_LOUD_SAYS = "GATE: and none was possible: no KO-less group reached"


@needs_r()
@pytest.mark.parametrize("koless_n,koless_tested,candidates,expect,forbid", [
    # Nothing KO-less reached the model, over a population: the coverage
    # failure, and the empty list is not a result at all.
    (19, 0, 0, [SHORTLIST_LOUD_SAYS, "(0 of 19 quantified)"], ["NOTE:"]),
    # The same zero under the floor. The tier drops to a NOTE - "none of 7"
    # is an anecdote - but the CLAIM does not change with it: there was still
    # nothing to rank, and the list is still not a negative result. Saying it
    # was "a statement about those 0" is the reading this block exists to
    # prevent, printed by the block itself.
    (7, 0, 0, ["NOTE: and nothing was rankable: no KO-less group reached the "
               "model, of 7 quantified",
               "not a negative result either"],
     ["GATE:", "statement about those 0", "drawn from the 0"]),
    # KO-less proteins WERE tested and none was significant. That is a
    # result, and it says what it was drawn from.
    (7, 7, 0, ["NOTE: the list was drawn from the 7 KO-less group(s) that "
               "reached the model, of 7 quantified",
               "a statement about those 7"], ["GATE:"]),
    # A list with candidates on it explains nothing: there is nothing to
    # explain.
    (7, 7, 3, [], ["GATE:", "NOTE:"]),
    # No KO-less group quantified at all - 3p_profile_only on every run there
    # has been. Silent, like the bin rule above it.
    (0, 0, 0, [], ["GATE:", "NOTE:"]),
])
def test_an_empty_shortlist_never_reads_as_a_result_it_is_not(
        ma, tmp_path, koless_n, koless_tested, candidates, expect, forbid):
    # symptom: under the denominator floor the branch printed "the list was
    # drawn from the 0 KO-less group(s) that reached the model ... an empty
    # list is a statement about those 0" - no population, presented as a weak
    # negative result, which is the exact confusion the block was added to
    # remove. Reachable on any small run.
    consts, _ = _coverage_source(ma)
    code = consts + SHORTLIST_PREAMBLE + _shortlist_source(ma)
    r = _rscript_file(tmp_path, code, koless_n, koless_tested, candidates)
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    if not expect:
        loud = _rscript_file(tmp_path, code, *SHORTLIST_LOUD)
        assert SHORTLIST_LOUD_SAYS in loud.stdout + loud.stderr, \
            "the branch is not in the document, so silence proves nothing"
    for want in expect:
        assert want in out, out
    for no in forbid:
        assert no not in out, out


def _rendered_output(path):
    """Only what the report PRINTED, never its own source.

    `code_folding: hide` echoes every chunk into the page, so a sentence that
    lives in a format string or a comment is in the HTML whether or not the
    line it belongs to ever fired - which makes a raw search on the file
    useless for asserting that something was NOT said. knitr prefixes each
    line of a chunk's output with `## `, and that prefix is the only thing
    that tells the two apart.
    """
    txt = _text_of(path)
    out = []
    for m in re.finditer(r"<pre><code>(.*?)</code></pre>", txt, re.S):
        for line in unescape(m.group(1)).splitlines():
            if line.startswith("## "):
                out.append(line[3:])
    return "\n".join(out)


def _text_of(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


@pytest.fixture(scope="module")
def sparse_koless(tmp_path_factory):
    """A finished run shaped like the run issue #5 was found on, knitted once.

    Rendering is the expensive step, so the assertions about what it says are
    spread over several tests reading one rendered document.

    Twelve extra dark proteins so `4_dark` clears the denominator floor, and
    every KO-less protein quantified in too few samples per group to pass
    `min_valid_per_group`. They are QUANTIFIED - a real intensity in half of
    the samples - so they are in annotated_quant.tsv and in the retention
    table, and then the filter takes all of them. That is the shape the issue
    is about: not absent, filtered out.
    """
    rng = random.Random(11)
    proteins = F.protein_set(n_extra=30)
    proteins += [F.Protein(f"P_darkx{i:02d}",
                           "".join(rng.choice(F.AA) for _ in range(110)),
                           in_emapper=False, expect_bin="4_dark")
                 for i in range(12)]
    samples = ["A_1", "A_2", "A_3", "A_4", "B_1", "B_2", "B_3", "B_4"]
    proj = build_project(tmp_path_factory.mktemp("sparse") / "p",
                         proteins=proteins, samples=samples)
    koless = [p.pid for p in proteins
              if p.expect_bin not in ("1_ko_pathway", "2_ko_orphan")]
    q = proj.path("input", "combined_peptide.tsv")
    d = pd.read_csv(q, sep="\t")
    sparse = d["Protein"].isin(koless)
    # 0 is what FragPipe writes for "not quantified", and metaannot reads it
    # as missing: two valid values in each group of four, against a
    # min_valid_per_group of three.
    for c in [f"{s} Intensity" for s in ("A_3", "A_4", "B_3", "B_4")]:
        d.loc[sparse, c] = 0
    d.to_csv(q, sep="\t", index=False)
    proj.run()
    proj.rendered = False
    if shutil.which("Rscript") and shutil.which("pandoc") and r_has(*R_CORE):
        run_metaannot("report", "--config", proj.config_path, cwd=proj.root,
                      timeout=1800)
        proj.rendered = True
    proj.koless = koless
    return proj


@needs_r(*R_CORE, *R_BIOC)
def test_a_knitted_report_escalates_the_bins_that_reach_the_model_with_nothing(
        sparse_koless):
    # symptom: issue #5. Every number was already in the report - the
    # retention table printed 0 for two whole bins, and the join stage had
    # warned about the proteins it dropped - but nothing escalated any of it,
    # so "your statistics cover none of the fraction this tool exists for"
    # was left for the reader to work out from a cell in a tibble.
    if not sparse_koless.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    html = _rendered_output(sparse_koless.rpath("analysis",
                                                "analyse_metaannot.html"))
    # the bin that is big enough for its zero to be a claim about a population
    assert "GATE  0/15 4_dark quantified group(s) reach the model" in html
    # the headline: the KO-less fraction as a whole is in no statistic below
    assert ("GATE  no KO-less protein is tested in any contrast: 0 of 19 "
            "quantified KO-less group(s)") in html
    # and the knobs that decided it, named where the consequence is - and
    # named as the config really keeps them: min_features_per_protein is
    # top-level, and this run is not isobaric, so there are two of them
    assert "Two knobs decide it" in html
    assert "analysis.min_valid_per_group is 3" in html
    assert "min_features_per_protein decided which proteins" in html
    assert "join.min_features_per_protein" not in html
    # the small bins are reported without being escalated
    assert ("NOTE  0/2 3_annotated_no_ko and 0/2 3d_duf_only reach the model"
            in html)


@needs_r(*R_CORE, *R_BIOC)
def test_the_ratio_model_and_the_shortlist_say_it_too_rather_than_ranking_nothing(
        sparse_koless):
    # symptom: with no KO-less protein in the model the ratio model still runs
    # - on the enzymes it was never the argument for - and the effector
    # shortlist prints "0 candidates", which reads as "nothing was
    # significant" when the truth is that it had nobody to rank.
    if not sparse_koless.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    html = _rendered_output(sparse_koless.rpath("analysis",
                                                "analyse_metaannot.html"))
    # the ratio model says WHICH of the two things happened: nothing KO-less
    # reached the model at all, which is the coverage failure gated above and
    # not a limit of the taxonomy, so it is not escalated a second time here
    assert ("NOTE  and none of them is KO-less, because no KO-less group "
            "reached the model at all") in html
    assert "GATE  and none of them is KO-less" not in html
    assert "0 candidates (significant, no KO, secreted or surface-exposed)" in html
    assert ("GATE  and none was possible: no KO-less group reached the model "
            "at all (0 of 19 quantified)") in html


@needs_r(*R_CORE, *R_BIOC)
def test_a_bin_with_nothing_to_lose_raises_nothing_on_a_healthy_knit(knitted):
    # the other half of the rule, and the half that keeps the gate readable:
    # 3p_profile_only is fed only by hhblits and jackhmmer, so no run has ever
    # put a protein in it. A gate that fired on its 0 every time would teach
    # the reader to skip the line they most need to read.
    #
    # The absences below are only worth something against a document that
    # really has the feature, which is what the first half of this test is
    # for: asserting that nothing was said is also satisfied by a report that
    # cannot say it. So the table's new columns are checked, and the one line
    # this run DOES earn - the ratio model's sub-floor NOTE, its KO-less
    # population being too small to gate over - is required to be there. That
    # NOTE was a GATE once, on this very knit, which is the bug the absence
    # list below could not see.
    if not knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    html = _rendered_output(knitted.rpath("analysis",
                                          "analyse_metaannot.html"))
    assert "retention by bin" in html, "the retention table is gone"
    for col in ("n_kept", "n_tested", "pct_tested"):
        assert col in html, f"the retention table has no {col} column"
    assert ("NOTE  and none of them is KO-less: the 7 KO-less group(s) in the "
            "model are without a usable taxon") in html
    for absent in ("quantified group(s) reach the model",
                   "no KO-less protein is tested in any contrast",
                   "the statistics below cover",
                   # the ratio model's version of the same rule: its floor is
                   # the KO-less population in the model, and seven is under
                   # it. Without one, this GATE fired on a healthy document.
                   "GATE  and none of them is KO-less",
                   # every knob-naming NOTE belongs to a line that fired
                   "knobs decide it"):
        assert absent not in html, \
            f"a coverage line fired on a healthy run: {absent}"
    # and no GATE anywhere in the document is one of this rule's, however it
    # was worded. Scoped to the rule's own vocabulary rather than to the word
    # KO-less, which a legitimate gate two sections up uses for the
    # identification-rate odds ratio.
    for line in html.splitlines():
        if not line.startswith("GATE"):
            continue
        for phrase in ("KO-less group", "KO-less protein", "reach the model",
                       "statistics below cover"):
            assert phrase not in line, \
                f"a coverage GATE on a healthy run: {line}"


@needs_r(*R_CORE, *R_BIOC)
def test_an_empty_shortlist_says_what_population_it_was_drawn_from(knitted):
    # symptom: "0 candidates" is a negative RESULT when the KO-less proteins
    # were tested and none was significant, and a COVERAGE failure when none
    # was tested at all. The printed line was the same either way. This is the
    # first of those two, and it must not read as the second - the tool has
    # been wrong in that direction before (see the v0.2.0 CHANGELOG entry).
    if not knitted.rendered:
        pytest.skip("the report was not rendered (pandoc or a package is absent)")
    html = _rendered_output(knitted.rpath("analysis",
                                          "analyse_metaannot.html"))
    assert "0 candidates (significant, no KO, secreted or surface-exposed)" in html
    assert "NOTE  the list was drawn from the 7 KO-less group(s) that reached" \
        in html
    assert "none was possible" not in html
