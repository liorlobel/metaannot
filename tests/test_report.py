"""The embedded R report and the R object.

Everything that needs R is skipped cleanly when Rscript or a package is
missing; the static checks below run everywhere.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest

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


@pytest.mark.xfail(reason="LIVE DEFECT: both build_object.R and the Rmd treat "
                          "quant/sample_columns.txt as the most authoritative "
                          "statement of which columns are samples, and NOTHING "
                          "in metaannot.py ever writes it. Every run therefore "
                          "falls through to design_from_input.tsv, or — with "
                          "no manifest — to the column-type heuristic the file "
                          "exists to avoid.",
                   strict=True)
def test_the_join_stage_writes_the_sample_column_list(tmp_path):
    proj = build_project(tmp_path / "p")
    proj.run()
    assert os.path.exists(proj.rpath("quant", "sample_columns.txt"))


def test_the_design_is_used_when_no_sample_column_list_exists(tmp_path):
    # the documented second choice: it must actually be produced by a run.
    proj = build_project(tmp_path / "p")
    proj.run()
    d = pd.read_csv(proj.rpath("quant", "design_from_input.tsv"), sep="\t")
    assert sorted(d["sample"]) == sorted(proj.samples)


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
