"""`doctor --json`: the preflight contract, and the claims it makes about code.

Two halves, and the second is the one that matters. The first pins the SHAPE -
the keys, the closed vocabularies, the version rule - the way
tests/test_config.py pins `describe --json`. The second pins the CLAIMS: every
`blocks` entry says "this stage's function dies on this state", and the only
honest way to test that is to drive the stage function into that state and
watch it die. An earlier attempt at "doctor agrees with run" failed precisely
because the agreement lived in prose, so a test that merely asserts a string is
present would reproduce it.
"""
from __future__ import annotations

import ast
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time

import pandas as pd
import pytest

import fixtures as F
from conftest import METAANNOT_PY, build_project, run_metaannot


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def no_r_env():
    """PATH without Rscript.

    Two reasons, and the second is the important one. It keeps these tests off
    a 180-second Rscript query per subprocess, and - because the R block's
    verdict depends on which packages happen to be installed on the machine
    running the suite - it keeps the EXIT STATUS of every config below a
    property of the config rather than of the laptop. The R rows get their own
    tests, in process, with the answer supplied.

    IT IS ALSO A BLIND SPOT, and worth naming here because it hid a real
    defect. No R failure is ever in a swept document, and the R rows are where
    `blocks_commands` lives: a whole sweep of configs could pass the invariants
    without one row that fails while killing no stage, which is exactly the
    class `verdict.rule` and the TUTORIAL had both dropped from the rule. Two
    things cover it now - CONFIGS carries `no_faa`, whose `input:proteins_faa`
    row is a `blocks_commands` failure with no R involved, and
    test_the_r_rows_obey_the_same_invariants_the_swept_configs_do runs the
    sweep's own invariants over the R block in process.
    """
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep)
             if p and not os.path.exists(os.path.join(p, "Rscript"))]
    return {"PATH": os.pathsep.join(parts)}


def doctor_json(config_path=None, *extra, expect=None, cwd=None, env=None,
                timeout=300):
    """Run `doctor --json` in a subprocess and return (returncode, document).

    A subprocess rather than an in-process call, because half of what this file
    pins is about the PROCESS: that stdout parses as JSON and carries nothing
    else, and that the exit status and the document agree.

    `timeout` because one of the states below does not make the command fail -
    it makes it NEVER RETURN. A FIFO at a path doctor reads was a blocking
    open, so the process printed nothing and exited never; a suite that called
    this without a timeout would have hung on it rather than gone red, and a
    test that hangs reports nothing at all.
    """
    args = ["doctor", "--json"]
    if config_path:
        args += ["--config", config_path]
    proc = run_metaannot(*args, *extra, expect=expect, cwd=cwd, timeout=timeout,
                         env=env if env is not None else no_r_env())
    return proc, json.loads(proc.stdout)


def no_orf_finder_env():
    """PATH with no ORF finder on it, for the same reason no_r_env() exists.

    A verdict that depends on whether smorf, smorfinder or macrel happens to
    be installed is a verdict this suite cannot assert without asserting
    something about the laptop it runs on - which is exactly how a row claiming
    `blocks: ["smorf"]` unconditionally got green-lit by a test on a machine
    where the claim was false. Directories holding an ORF finder are dropped
    and everything else is kept, so the OTHER `have()` probes - and therefore
    the exit status - are unchanged.
    """
    names = ("smorf", "smorfinder", "macrel")
    parts = [d for d in no_r_env()["PATH"].split(os.pathsep)
             if d and not any(os.path.exists(os.path.join(d, n))
                              for n in names)]
    return {"PATH": os.pathsep.join(parts)}


def by_id(doc):
    return {c["id"]: c for c in doc["checks"]}


def project_with(tmp_path, name, **over):
    return build_project(tmp_path / name, **over)


def run_flags(p, **on):
    """Rewrite the project's config with some run: flags flipped."""
    p.write_config(run=dict(p.cfg["run"], **on))
    return p


# ----------------------------------------------------------------------
# the stream, and the verdict
# ----------------------------------------------------------------------
def test_the_document_is_the_only_thing_on_stdout_even_when_checks_fail(
        project):
    # The console's reader does json.loads on the WHOLE buffer, the way its
    # run_describe already does. One stray print - the install-plan line, a
    # progress note, a warning that forgot to go to stderr - and a consumer
    # gets nothing at all, which is a worse failure than a red row.
    os.remove(project.cfg["proteins_faa"])
    proc, doc = doctor_json(project.config_path, expect=1)
    assert proc.stdout.lstrip().startswith("{")
    assert json.loads(proc.stdout) == doc, "stdout is not exactly one document"
    assert doc["exit_status"] == 1


def test_the_document_stays_json_when_the_config_itself_has_unknown_keys(
        project):
    # load_config reports unknown keys through log(), which writes to stderr.
    # If it ever did not, this is the test that says so: an operator's typo
    # must not be able to corrupt a machine-readable stream.
    project.write_config(quant_tabel="oops", threadz=4)
    proc, doc = doctor_json(project.config_path, expect=1)
    assert proc.stdout.lstrip().startswith("{")
    ids = by_id(doc)
    assert ids["config:unknown:quant_tabel"]["status"] == "fail"
    assert "did you mean 'quant_table'" in ids["config:unknown:quant_tabel"]["detail"]


def test_the_document_and_the_exit_status_are_the_same_verdict(project):
    # A verdict delivered twice that can disagree is worse than one delivered
    # once. `ok`, `exit_status` and `verdict.fails` all come off one
    # expression, and $? has to be that expression too.
    for broken in (False, True):
        if broken:
            os.remove(project.cfg["quant_table"])
            run_flags(project, unipept=True)
        proc, doc = doctor_json(project.config_path, expect=None)
        assert doc["exit_status"] == proc.returncode
        assert doc["ok"] == (proc.returncode == 0)
        assert bool(doc["verdict"]["fails"]) == (proc.returncode == 1)
        assert doc["verdict"]["fails"] == [c["id"] for c in doc["checks"]
                                           if c["status"] == "fail"]


def test_the_verdict_rule_names_exactly_one_field_a_consumer_can_read(project):
    # The rule is a literal sentence rather than a hint, so a console that
    # understands nothing else in a NEWER document can still compute the exit
    # status from it.
    _proc, doc = doctor_json(project.config_path, expect=None)
    rule = doc["verdict"]["rule"]
    assert "status" in rule and "'fail'" in rule and "exit_status" in rule


def test_the_scope_sentence_is_the_engines_to_author_not_the_consoles(
        ma, project):
    # "and doctor says so" is half the scope rule. A console that has to write
    # this sentence itself becomes its second home, and the two drift.
    _proc, doc = doctor_json(project.config_path, expect=None)
    assert doc["scope"]["statement"] == ma.DOCTOR_SCOPE_STATEMENT
    assert doc["scope"]["checks"] == ["exists", "kind", "non-empty"]
    assert doc["scope"]["does_not_check"] == ["contents"]
    for word in ("exists", "right kind", "not empty", "does not parse"):
        assert word in doc["scope"]["statement"]


# ----------------------------------------------------------------------
# the verdict is a property of the CONFIG, not of the file
# ----------------------------------------------------------------------
def test_a_missing_input_of_an_enabled_stage_fails_and_names_the_stage(
        tmp_path):
    # run.context on, gff set and absent: stage_context calls die().
    p = project_with(tmp_path, "ctx", gff=str(tmp_path / "absent.gff"))
    run_flags(p, context=True)
    proc, doc = doctor_json(p.config_path, expect=1)
    gff = by_id(doc)["input:gff"]
    assert gff["status"] == "fail"
    assert gff["blocks"] == ["context"]
    assert gff["fails_reason"] == "stage_or_command_dies"
    assert proc.returncode == 1


def test_a_missing_input_of_a_disabled_stage_does_not_fail_doctor(tmp_path):
    # The same absent gff with run.context off. Nothing reads it, so there is
    # nothing to fail: the check is emitted as `skip` rather than omitted, so
    # the row set is stable across two documents an operator is diffing.
    p = project_with(tmp_path, "noctx", gff=str(tmp_path / "absent.gff"))
    proc, doc = doctor_json(p.config_path, expect=0)
    gff = by_id(doc)["input:gff"]
    assert gff["status"] == "skip" and gff["finding"] == "not_enabled"
    assert gff["blocks"] == [] and gff["degrades"] == []
    assert proc.returncode == 0


def test_an_empty_gff_is_a_warning_because_stage_context_writes_an_empty_table(
        tmp_path):
    # Three behaviours on one key, and the middle one is the reason `gff` needs
    # a per-state check: run.context on with `gff: ""` is not an error at all.
    p = project_with(tmp_path, "emptygff", gff="")
    run_flags(p, context=True)
    _proc, doc = doctor_json(p.config_path, expect=0)
    gff = by_id(doc)["input:gff"]
    assert gff["status"] == "warn" and gff["degrades"] == ["context"]


def test_stage_context_really_dies_on_a_missing_gff_and_really_does_not_on_an_empty_one(
        ma, tmp_path, paths_for):
    # The claim `input:gff` makes about the code, driven rather than asserted.
    cfg, p = paths_for()
    cfg["gff"] = str(tmp_path / "nowhere.gff")
    with pytest.raises(ma.StageError) as e:
        ma.stage_context(cfg, p)
    assert "gff not found" in str(e.value)
    cfg["gff"] = ""
    ma.stage_context(cfg, p)                 # warns, writes, returns
    assert os.path.exists(p.context)


# ----------------------------------------------------------------------
# the stage_join case the whole contract exists for
# ----------------------------------------------------------------------
def test_stage_join_really_warns_and_returns_on_a_missing_quant_table(
        ma, tmp_path, paths_for, capsys):
    # Verified before anything is claimed about it: stage_join's first act is
    # `if not os.path.exists(qpath): log(...); return`. run.join is true in
    # DEFAULT_CONFIG, so this is the ordinary annotate-without-MS-quant case
    # and doctor must not fail it.
    cfg, p = paths_for()
    cfg["quant_table"] = str(tmp_path / "nowhere.tsv")
    ma.stage_join(cfg, p)                    # no StageError
    assert "quant table not found, skipping join" in capsys.readouterr().err
    assert ma.DEFAULT_CONFIG["run"]["join"] is True, \
        "run.join stopped being on by default; this test's premise moved"


def test_an_annotate_without_ms_quant_config_passes_doctor(tmp_path):
    p = project_with(tmp_path, "nq")
    os.remove(p.cfg["quant_table"])
    os.remove(p.cfg["manifest"])
    proc, doc = doctor_json(p.config_path, expect=0)
    q = by_id(doc)["input:quant_table"]
    assert q["status"] == "warn"
    assert q["blocks"] == [] and q["degrades"] == ["join"]
    assert "stage_join logs a WARN and returns" in q["detail"]
    assert proc.returncode == 0, "doctor failed a config `run` exits 0 on"


def test_the_same_missing_quant_table_is_fatal_once_unipept_is_on(tmp_path):
    # Same file, same id, opposite verdict. Severity is a property of the
    # config: stage_unipept calls peptide_features(), which reads the table and
    # has no skip branch.
    p = project_with(tmp_path, "nq2")
    os.remove(p.cfg["quant_table"])
    run_flags(p, unipept=True)
    proc, doc = doctor_json(p.config_path, expect=1)
    q = by_id(doc)["input:quant_table"]
    assert q["status"] == "fail"
    assert q["blocks"] == ["unipept"]
    assert q["degrades"] == ["join"], "join still only skips"
    assert proc.returncode == 1


def test_unipept_with_a_result_file_does_not_claim_the_quant_table(tmp_path):
    # stage_unipept returns after ingesting unipept.result and never reaches
    # peptide_features, so the missing table stops blocking it. A `blocks` that
    # ignored this would be a fatal doctor cannot justify.
    p = project_with(tmp_path, "nq3")
    res = tmp_path / "pept2lca.tsv"
    res.write_text("peptide\tgenus_id\nPEPTIDE\t1\n", encoding="utf-8")
    os.remove(p.cfg["quant_table"])
    p.write_config(run=dict(p.cfg["run"], unipept=True),
                   unipept={"result": str(res)})
    _proc, doc = doctor_json(p.config_path, expect=0)
    assert by_id(doc)["input:quant_table"]["blocks"] == []


# ----------------------------------------------------------------------
# exists / right kind / not empty - the three states of one path
# ----------------------------------------------------------------------
def test_a_zero_byte_quant_table_fails_where_a_missing_one_only_warns(tmp_path):
    # The scope's first sentence, and the asymmetry that makes it necessary:
    # os.path.exists() is true for an empty file, so stage_join never takes its
    # skip branch and dies in the reader instead. Missing is survivable; empty
    # is not, and only a per-state check can say that.
    p = project_with(tmp_path, "zero")
    open(p.cfg["quant_table"], "w", encoding="utf-8").close()
    proc, doc = doctor_json(p.config_path, expect=1)
    q = by_id(doc)["input:quant_table"]
    assert q["status"] == "fail" and q["finding"] == "empty"
    assert q["found"]["kind"] == "empty_file" and q["found"]["bytes"] == 0
    assert q["blocks"] == ["join"]
    assert proc.returncode == 1


def test_a_zero_byte_proteins_faa_blocks_integrate_rather_than_the_run_command(
        tmp_path):
    # cmd_run's own test is os.path.exists, so an empty FASTA gets past it,
    # every search stage runs happily against nothing, and build_annotation
    # dies on "no sequences in ..." hours later. Same id, different finding,
    # different victim.
    p = project_with(tmp_path, "zfaa")
    open(p.cfg["proteins_faa"], "w", encoding="utf-8").close()
    _proc, doc = doctor_json(p.config_path, expect=1)
    faa = by_id(doc)["input:proteins_faa"]
    assert faa["status"] == "fail" and faa["finding"] == "empty"
    assert "integrate" in faa["blocks"]
    assert faa["blocks_commands"] == [], "`run` itself does not refuse this one"


def test_a_missing_proteins_faa_blocks_the_run_command_not_a_list_of_stages(
        tmp_path):
    # cmd_run refuses before it schedules anything, so naming nine stages would
    # overstate what happens. `blocks` stays a pure list of stage names that
    # joins to describe --json's stage_names; the command goes in its own list.
    p = project_with(tmp_path, "nofaa")
    os.remove(p.cfg["proteins_faa"])
    _proc, doc = doctor_json(p.config_path, expect=1)
    faa = by_id(doc)["input:proteins_faa"]
    assert faa["blocks"] == [] and faa["blocks_commands"] == ["run"]
    # ...and the reason names BOTH halves, because this row is the half the
    # old name left out: nothing in `blocks`, a command in `blocks_commands`,
    # and a value called `stage_dies` describing it.
    assert faa["fails_reason"] == "stage_or_command_dies"


def test_a_directory_where_a_file_is_expected_is_reported_as_a_wrong_kind(
        tmp_path):
    # os.path.exists() is true for a directory, which is how a path typo used
    # to pass every check doctor had.
    p = project_with(tmp_path, "dir")
    os.remove(p.cfg["quant_table"])
    os.makedirs(p.cfg["quant_table"])
    open(os.path.join(p.cfg["quant_table"], "x.tsv"), "w",
         encoding="utf-8").write("x\n")
    _proc, doc = doctor_json(p.config_path, expect=1)
    q = by_id(doc)["input:quant_table"]
    assert q["status"] == "fail" and q["finding"] == "wrong_kind"
    assert q["expect"]["kind"] == "file" and q["found"]["kind"] == "dir"
    assert q["expect"]["derived_from"] == ["quant_format"], \
        "the expected kind must say which config key decided it"


def test_a_file_where_fragpipe_tmt_expects_a_directory_is_a_wrong_kind(
        tmp_path):
    # "Right kind" is config-dependent and must never be hard-coded to "file":
    # of the quant formats only fragpipe_tmt makes quant_table a directory.
    p = project_with(tmp_path, "tmtfile", quant_format="fragpipe_tmt")
    _proc, doc = doctor_json(p.config_path, expect=1)
    q = by_id(doc)["input:quant_table"]
    assert q["status"] == "fail" and q["finding"] == "wrong_kind"
    assert q["expect"]["kind"] == "dir" and q["found"]["kind"] == "file"
    assert "per-plex" in q["expect"]["because"]
    assert q["blocks"] == ["join"]


def test_a_dangling_symlink_is_named_rather_than_reported_as_merely_absent(
        tmp_path):
    # os.path.exists() follows the link and answers False, so every check has
    # always rejected one "for the right reason by accident" and could not say
    # why - and an unmounted volume is the common way this happens here.
    p = project_with(tmp_path, "sym")
    faa = p.cfg["proteins_faa"]
    os.remove(faa)
    os.symlink(str(tmp_path / "nowhere.faa"), faa)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["input:proteins_faa"]
    assert c["finding"] == "dangling_symlink"
    assert c["found"]["kind"] == "symlink_broken" and c["found"]["symlink"]
    assert "symlink with no target" in c["detail"]


# ----------------------------------------------------------------------
# config keys
# ----------------------------------------------------------------------
def test_a_retired_config_key_warns_and_does_not_fail_the_document(tmp_path):
    # A retired key is not a failure: the config predates a removal, the
    # setting is inert, and exiting non-zero over it would block a run that is
    # otherwise correct.
    p = project_with(tmp_path, "retired")
    proc, doc = doctor_json(p.config_path, expect=0)
    ids = by_id(doc)
    retired = [c for c in doc["checks"] if c["id"].startswith("config:retired:")]
    assert retired, "the fixture no longer carries a retired key"
    for c in retired:
        assert c["status"] == "warn" and c["fails_reason"] is None
        assert c["blocks"] == [] and c["blocks_commands"] == []
    assert proc.returncode == 0


def test_an_unrecognised_key_fails_and_names_itself_as_the_one_exception(
        tmp_path):
    # report_unknown_keys only logs and load_config carries on, so no stage
    # dies - and yet the setting the operator believes is in effect is not.
    # This is the single declared exception to "fails exactly when a stage
    # dies", and it says so in a field rather than in prose nobody can switch
    # on.
    p = project_with(tmp_path, "typo")
    p.write_config(emapper_min_coverag=0.5)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["config:unknown:emapper_min_coverag"]
    assert c["status"] == "fail"
    assert c["fails_reason"] == "setting_ignored"
    assert c["blocks"] == [] and c["blocks_commands"] == []
    others = [x for x in doc["checks"]
              if x["status"] == "fail"
              and x["fails_reason"] != "stage_or_command_dies"]
    assert [x["id"] for x in others] == [c["id"]], \
        "a second failure that kills nothing appeared; the exception is " \
        "meant to be the only one"


# ----------------------------------------------------------------------
# requirements, and what they really block
# ----------------------------------------------------------------------
def test_a_missing_database_of_an_enabled_stage_fails_and_names_that_stage(
        tmp_path):
    p = project_with(tmp_path, "pfam")
    p.write_config(run=dict(p.cfg["run"], pfam=True),
                   db={"pfam_hmm": str(tmp_path / "absent" / "Pfam-A.hmm")})
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["req:pfam_hmm"]
    assert c["status"] == "fail" and c["blocks"] == ["pfam"]
    assert c["requirement_id"] == "pfam_hmm"
    reqs = {r["id"]: r for r in doc["requirements"]}
    assert reqs["pfam_hmm"]["cmds"], \
        "the install commands live in requirements[], reached by requirement_id"
    assert "cmds" not in c and "size_gb" not in c, \
        "a check must not carry a second copy of what requirements[] has"


def test_the_same_database_missing_with_that_stage_off_is_not_even_listed(
        tmp_path):
    # requirements() never produces a row for a disabled stage, and doctor does
    # not synthesise one: a table showing every stage's databases as "-" would
    # be mostly noise and the count would mislead.
    p = project_with(tmp_path, "nopfam",
                     db={"pfam_hmm": str(tmp_path / "absent" / "Pfam-A.hmm")})
    proc, doc = doctor_json(p.config_path, expect=0)
    assert "req:pfam_hmm" not in by_id(doc)
    assert proc.returncode == 0


def test_a_missing_diamond_database_warns_because_stage_diamond_skips_it(
        tmp_path):
    # The counter-intuitive one, and the opposite of what doctor has always
    # said: an ABSENT diamond database is a warning, because stage_diamond logs
    # "diamond database missing, skipping" and searches the rest.
    p = project_with(tmp_path, "dia")
    p.write_config(run=dict(p.cfg["run"], diamond=True),
                   db={"diamond": {"vfdb": str(tmp_path / "absent.dmnd")}})
    proc, doc = doctor_json(p.config_path, expect=None)
    c = by_id(doc)["req:diamond:vfdb"]
    assert c["status"] == "warn"
    assert c["blocks"] == [] and c["degrades"] == ["diamond"]
    assert "req:diamond:vfdb" not in doc["verdict"]["fails"]


def test_stage_diamond_really_skips_an_absent_database_and_really_dies_on_a_broken_one(
        ma, tmp_path, paths_for, stub_bin, capsys):
    # Both halves of the asymmetry, driven. This is what `blocks` is claiming.
    cfg, p = paths_for()
    cfg["proteins_faa"] = F.write_fasta(str(tmp_path / "p.faa"),
                                        F.protein_set()[:2])
    cfg["db"]["diamond"] = {"vfdb": str(tmp_path / "absent.dmnd")}
    ma.stage_diamond(cfg, p)                 # no StageError
    assert "diamond database missing, skipping" in capsys.readouterr().err

    broken = tmp_path / "broken.dmnd"
    broken.write_bytes(b"")                  # a failed `diamond makedb`
    cfg["db"]["diamond"] = {"vfdb": str(broken)}
    with pytest.raises(ma.StageError) as e:
        ma.stage_diamond(cfg, p)
    assert "cannot answer" in str(e.value)


def test_a_present_but_unanswerable_diamond_database_is_a_second_check(
        tmp_path):
    # Existence and usability are two questions about one file, with different
    # depths, different marks and different verdicts - which is the clearest
    # argument for keying the document on the CHECK rather than on the subject.
    p = project_with(tmp_path, "dia2")
    broken = tmp_path / "broken.dmnd"
    broken.write_bytes(b"")
    p.write_config(run=dict(p.cfg["run"], diamond=True),
                   db={"diamond": {"vfdb": str(broken)}})
    _proc, doc = doctor_json(p.config_path, expect=1)
    ids = by_id(doc)
    assert ids["req:diamond:vfdb"]["status"] == "warn", "existence is separate"
    usable = ids["db:diamond:vfdb:usable"]
    assert usable["status"] == "fail" and usable["blocks"] == ["diamond"]
    # ...and the depth is the one this refusal actually reached. A zero-byte
    # .dmnd is refused on os.path.getsize alone - diamond_db_check returns
    # before it profiles anything - so the row that used to hard-code
    # `header` and "read from the DIAMOND header" was describing a read that
    # never happened.
    assert usable["depth"] == "kind"
    assert "size" in usable["caveat"]
    assert usable["depends_on"] == ["req:diamond:vfdb"]


def test_a_missing_smorf_tool_only_degrades_while_its_input_is_fatal(tmp_path):
    # The inversion stage_smorf makes deliberately: the INPUT is fatal, but a
    # missing smorf/macrel binary warns and writes an empty file, because
    # "every other stage is independent of this one".
    p = project_with(tmp_path, "smorf")
    run_flags(p, smorf=True)
    _proc, doc = doctor_json(p.config_path, expect=1)
    ids = by_id(doc)
    assert ids["req:smorf"]["status"] == "warn"
    assert ids["req:smorf"]["degrades"] == ["smorf"]
    assert ids["input:contigs_fna"]["status"] == "fail"
    assert ids["input:contigs_fna"]["blocks"] == ["smorf"]


def test_stage_smorf_really_dies_without_contigs_and_really_does_not_without_its_tools(
        ma, tmp_path, paths_for, capsys):
    cfg, p = paths_for()
    cfg["contigs_fna"] = ""
    with pytest.raises(ma.StageError) as e:
        ma.stage_smorf(cfg, p)
    assert "needs contigs_fna" in str(e.value)
    cfg["contigs_fna"] = F.write_fasta(str(tmp_path / "c.fna"),
                                       F.protein_set()[:1])
    ma.stage_smorf(cfg, p)                   # neither smorf nor macrel is here
    assert "no small ORFs" in capsys.readouterr().err
    assert os.path.exists(p.smorf_faa)


def test_emapper_precomputed_is_not_a_failure_when_nothing_reads_it(tmp_path):
    # Folded into `ok` with no run.eggnog guard until now, so a stale path
    # failed doctor even when the only stage that reads it was off.
    p = project_with(tmp_path, "pre")
    p.write_config(run=dict(p.cfg["run"], eggnog=False),
                   emapper_precomputed=[str(tmp_path / "gone.tsv")])
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["precomputed_emapper:0"]
    assert c["status"] == "skip" and c["finding"] == "not_enabled"
    assert proc.returncode == 0

    p.write_config(run=dict(p.cfg["run"], eggnog=True))
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["precomputed_emapper:0"]
    assert c["status"] == "fail" and c["blocks"] == ["emapper"]


# ----------------------------------------------------------------------
# the tmt block, and the false-fail it carried
# ----------------------------------------------------------------------
def test_the_tmt_block_is_skipped_when_no_enabled_stage_reads_the_quant_tree(
        tmp_path):
    # Until now the whole `== tmt ==` block ran on quant_format alone and every
    # MISS inside it set ok=False, so a TMT config with run.join, run.unipept
    # and run.taxonomy all off exited 1 - on a config where nothing opens the
    # quant tree at all.
    p = project_with(tmp_path, "tmtoff", quant_format="fragpipe_tmt")
    p.write_config(run=dict(p.cfg["run"], join=False, unipept=False,
                            taxonomy=False))
    proc, doc = doctor_json(p.config_path, expect=0)
    ids = by_id(doc)
    assert ids["tmt:applies"]["status"] == "skip"
    assert ids["input:quant_table"]["status"] == "skip"
    assert [c["id"] for c in doc["checks"] if c["section"] == "tmt"] \
        == ["tmt:applies"]
    assert proc.returncode == 0


def test_a_tmt_config_that_does_read_the_quant_tree_still_checks_it(tmp_path):
    p = project_with(tmp_path, "tmton", quant_format="fragpipe_tmt")
    _proc, doc = doctor_json(p.config_path, expect=1)
    assert "tmt:applies" not in by_id(doc)
    assert by_id(doc)["input:quant_table"]["blocks"] == ["join"]


# ----------------------------------------------------------------------
# the invariants the contract rests on
# ----------------------------------------------------------------------
CONFIGS = ("healthy", "no_quant", "unipept", "empty_faa", "no_faa",
           "pfam_on", "tmt_file", "smorf_on", "typo", "man_dir", "fifo_faa")


def _configured(tmp_path, which):
    p = build_project(tmp_path / which)
    if which == "no_quant":
        os.remove(p.cfg["quant_table"])
    elif which == "unipept":
        os.remove(p.cfg["quant_table"])
        p.write_config(run=dict(p.cfg["run"], unipept=True))
    elif which == "empty_faa":
        open(p.cfg["proteins_faa"], "w", encoding="utf-8").close()
    elif which == "no_faa":
        # The one configuration in this sweep whose failure kills a COMMAND
        # and no stage. Without it the invariants below never see a
        # blocks_commands row, because no_r_env() keeps the R block - the
        # other place they come from - out of every swept document.
        os.remove(p.cfg["proteins_faa"])
    elif which == "fifo_faa":
        os.remove(p.cfg["proteins_faa"])
        os.mkfifo(p.cfg["proteins_faa"])
    elif which == "man_dir":
        # A DIRECTORY at the manifest, which used to end the command with an
        # IsADirectoryError before it printed a byte.
        os.remove(p.cfg["manifest"])
        os.makedirs(p.cfg["manifest"])
    elif which == "pfam_on":
        p.write_config(run=dict(p.cfg["run"], pfam=True, dbcan=True,
                                ncbifam=True, cluster=True, interpro=True))
    elif which == "tmt_file":
        p.write_config(quant_format="fragpipe_tmt")
    elif which == "smorf_on":
        p.write_config(run=dict(p.cfg["run"], smorf=True, context=True))
    elif which == "typo":
        p.write_config(thresholdz={"diamond_evalue": 1})
    return p


@pytest.mark.parametrize("which", CONFIGS)
def test_every_row_obeys_the_invariants_the_exit_status_is_derived_from(
        ma, tmp_path, which):
    # The scope's second sentence, mechanised, over a spread of real configs:
    #   fail            <=> fails_reason is set
    #   fails_reason == "stage_or_command_dies"
    #                   <=> blocks or blocks_commands is non-empty
    #   warn / skip     =>  neither list is non-empty
    #   skip            =>  degrades is empty too
    # Without these, `blocks` is prose in a JSON field.
    p = _configured(tmp_path, which)
    _proc, doc = doctor_json(p.config_path, expect=None)
    on = set(ma.enabled_stages(ma.load_config(p.config_path)))
    for c in doc["checks"]:
        where = f"{which}/{c['id']}"
        named = bool(c["blocks"] or c["blocks_commands"])
        assert (c["status"] == "fail") == (c["fails_reason"] is not None), where
        assert (c["fails_reason"]
                == "stage_or_command_dies") == named, where
        if c["status"] in ("warn", "skip", "ok"):
            assert not named, where
        if c["status"] == "skip":
            assert not c["degrades"], where
        for name in c["blocks"] + c["degrades"]:
            assert name in ma.STAGE_NAMES, f"{where}: {name} is not a stage"
            assert name in on, f"{where}: {name} is not enabled in this config"
        for name in c["blocks_commands"]:
            assert name in ma.DOCTOR_COMMANDS, where
        assert c["status"] in ma.DOCTOR_STATUSES, where
        assert c["remedy"] in ma.DOCTOR_REMEDIES, where
        assert c["depth"] in ma.DOCTOR_DEPTHS, where
        assert c["fails_reason"] in (None,) + ma.DOCTOR_FAIL_REASONS, where


@pytest.mark.parametrize("which", CONFIGS)
def test_the_ids_are_unique_and_the_sections_are_the_ones_declared(ma, tmp_path,
                                                                   which):
    # A console diffs two documents on `id`, so a repeated id makes the diff
    # meaningless, and a section nobody declared has no printed heading.
    p = _configured(tmp_path, which)
    _proc, doc = doctor_json(p.config_path, expect=None)
    ids = [c["id"] for c in doc["checks"]]
    assert len(ids) == len(set(ids)), f"{which}: duplicate check id"
    declared = [i for i, _t in ma.DOCTOR_SECTIONS]
    assert [s["id"] for s in doc["sections"]] == \
        [i for i in declared if i in {c["section"] for c in doc["checks"]}]
    # ...and the rows arrive grouped, in that order, so the printed report and
    # a table rendered from the document read the same way.
    order = [declared.index(c["section"]) for c in doc["checks"]]
    assert order == sorted(order), f"{which}: sections are interleaved"


# ----------------------------------------------------------------------
# the attribution that makes `blocks` possible at all
# ----------------------------------------------------------------------
def test_every_stage_declares_what_it_dies_without(ma):
    # requirements() gates on run: flags with inline boolean logic, so it knows
    # that SOMETHING wants hmmsearch and not which stage dies without it. This
    # is the tripwire: a stage added with no attribution gets `blocks` silently
    # wrong for every requirement it needs.
    for st in ma.STAGES:
        assert "requires" in st and "degraded_by" in st, \
            f"{st['name']} declares neither requires nor degraded_by"
        assert isinstance(st["requires"], tuple), st["name"]
        assert isinstance(st["degraded_by"], tuple), st["name"]


def test_every_requirement_this_tool_can_produce_is_claimed_by_some_stage(
        ma, tmp_path):
    # The other half: a requirement nobody claims reaches the document with an
    # empty `blocks`, which reads as "nothing dies without this" - the exact
    # false pass the second scope sentence was written against.
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["run"] = {k: True for k in cfg["run"]}
    cfg["emapper_precomputed"] = ""
    cfg["results_dir"] = str(tmp_path / "r")
    reqs = ma.requirements(cfg, ma.Paths(cfg))
    assert len(reqs) > 15, "the maximal config produces suspiciously little"
    unclaimed = []
    for r in reqs:
        blocks, degrades = ma.requirement_effect(cfg, r["id"])
        if not blocks and not degrades:
            unclaimed.append(r["id"])
    assert unclaimed == [], \
        f"no stage claims {unclaimed}; add them to a stage's requires/degraded_by"


def test_the_taxdump_blocks_join_only_when_taxon_rank_asks_for_a_rank(
        ma, tmp_path):
    # collapse_taxon_rank dies on a missing nodes.dmp under taxon_rank='genus'
    # and returns the raw seed taxids untouched when taxon_rank is empty. A
    # static attribution could not express that, so it is applied where the
    # condition is, with the derivation written beside it.
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["results_dir"] = str(tmp_path / "r")
    cfg["run"].update(taxonomy=True, unipept=True, join=True)
    # A quant table that is THERE, because a join that skips for want of one
    # reaches nothing it requires - see the test below, which is the other
    # half of this condition.
    quant = tmp_path / "q.tsv"
    quant.write_text("Peptide Sequence\n", encoding="utf-8")
    cfg["quant_table"] = str(quant)
    cfg["taxon_rank"] = ""
    blocks, degrades = ma.requirement_effect(cfg, "ncbi_taxonomy")
    assert "join" not in blocks and degrades == ["taxonomy"]
    cfg["taxon_rank"] = "genus"
    blocks, _ = ma.requirement_effect(cfg, "ncbi_taxonomy")
    assert blocks == ["join"]

    # ...and the condition is only ever consulted if requirements() PRODUCED
    # the row. With run.taxonomy off it did not, so the whole derivation above
    # ran over an entry that was not in the document and doctor exited 0 on a
    # config stage_join dies on. Both directions, at the place the row is made:
    rows = lambda c: [r["id"] for r in ma.requirements(c, ma.Paths(c))]
    cfg["run"].update(taxonomy=False)
    assert "ncbi_taxonomy" in rows(cfg), "a rank needs the taxdump for join"
    blocks, degrades = ma.requirement_effect(cfg, "ncbi_taxonomy")
    assert blocks == ["join"] and degrades == []
    cfg["taxon_rank"] = ""
    assert "ncbi_taxonomy" not in rows(cfg), \
        "with no rank and no comparison stage, nothing opens a taxdump"


def test_the_taxdump_row_promises_a_restored_quant_table_only_where_it_helps(
        ma, tmp_path):
    """symptom: `depends_on: ["input:quant_table"]` was a promise that
    restoring the quant table would make this row fatal, and with
    `taxon_rank: ""` it would not.

    TWO conditions take `join` out of `blocks` - an absent quant table, and an
    empty `taxon_rank` - and the contingency sentence is about only one of
    them. Asking it as "join is on, the quant table is gone, and join's
    requires tuple names this id" gets the same answer either way, so the row
    told an operator to restore a file that would change nothing.

    It is one question now, asked of `requirement_effect()` itself: would
    restoring the quant table put join in `blocks`.
    """
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["results_dir"] = str(tmp_path / "r")
    cfg["run"].update(taxonomy=True, unipept=True, join=True)
    cfg["quant_table"] = str(tmp_path / "gone.tsv")      # absent on purpose
    reqs = lambda c: ma.requirements(c, ma.Paths(c))
    row = lambda c: {x["requirement_id"]: x
                     for x in ma._requirement_checks(c, reqs(c))}

    cfg["taxon_rank"] = "genus"
    c = row(cfg)["ncbi_taxonomy"]
    assert c["blocks"] == [], "the quant table is gone, so join dies on nothing"
    assert c["depends_on"] == ["input:quant_table"], \
        "restoring it WOULD make this fatal for join, so the row says so"

    cfg["taxon_rank"] = ""
    c = row(cfg)["ncbi_taxonomy"]
    assert c["blocks"] == []
    assert c["depends_on"] == [], \
        "with no rank, collapse_taxon_rank never needs a taxdump - restoring " \
        "the quant table changes nothing and the row may not say it does"
    # ...and the derivation is the function's own answer, not a second copy
    # of half of it.
    assert ma.requirement_effect(cfg, "ncbi_taxonomy",
                                 join_reaches=True)[0] == []
    cfg["taxon_rank"] = "genus"
    assert ma.requirement_effect(cfg, "ncbi_taxonomy",
                                 join_reaches=True)[0] == ["join"]


# ----------------------------------------------------------------------
# the version constant, and the rule it carries
# ----------------------------------------------------------------------
# What `doctor_version` promises not to break silently. Changing either set
# means REMOVING or renaming something a reader outside this file depends on,
# which is exactly what DOCTOR_VERSION is for: bump it in the same commit, and
# update this test. Adding a key is not a bump and not a failure here, which is
# why both are compared exactly rather than as supersets - the pin has to be
# able to see a deletion.
DOCTOR_TOP_LEVEL = {
    "doctor_version", "describe_version", "metaannot_version",
    "signature_version", "generated", "host", "config_path",
    "ok", "exit_status", "counts", "verdict", "scope", "sections",
    "requirements", "checks", "totals", "install_plan",
}
DOCTOR_PER_CHECK = {
    "id", "section", "label", "status", "finding", "target", "expect", "found",
    "depth", "caveat", "config_keys", "requirement_id", "blocks",
    "blocks_commands", "degrades", "depends_on", "fails_reason", "remedy",
    "remedy_reason", "detail",
}


def test_doctor_version_pins_the_shape_it_versions(ma, project):
    # The sibling of test_describe_version_pins_the_shape_it_versions, and it
    # exists for the same reason: a version constant that pins nothing let nine
    # documented fields be deleted from `describe --json` while the whole suite
    # stayed green.
    _proc, doc = doctor_json(project.config_path, expect=None)
    assert set(doc) == DOCTOR_TOP_LEVEL, \
        "the top-level contract changed; bump DOCTOR_VERSION if a key went"
    for c in doc["checks"]:
        assert set(c) == DOCTOR_PER_CHECK, \
            f"{c['id']}: the per-check contract changed"
    assert doc["doctor_version"] == ma.DOCTOR_VERSION
    assert doc["describe_version"] == ma.DESCRIBE_VERSION, \
        "one call has to tell a consumer the level of BOTH contracts"


def test_adding_a_key_to_the_document_does_not_move_doctor_version(ma, project):
    # The RULE, exercised rather than merely written down: a consumer reading
    # the keys it already knew is unaffected by a new one appearing, so an
    # addition must not bump the constant. A renamed or deleted key is the
    # opposite, and the exact-set assertions above are what catch that.
    before = ma.DOCTOR_VERSION
    cfg = ma.load_config(project.config_path)
    p = ma.Paths(cfg)
    reqs = ma.requirements(cfg, p)
    checks = ma.doctor_checks(cfg, p, reqs, _Args(project.config_path))
    doc = ma.doctor(cfg, p, reqs, checks, project.config_path)
    doc["a_field_a_later_release_added"] = 42
    doc["checks"][0]["a_per_check_field_a_later_release_added"] = "x"
    assert doc["doctor_version"] == before, \
        "adding a key is not a version bump; see the rule at the constant"
    # ...and a consumer that only knows the old keys still reaches the verdict.
    assert (doc["exit_status"] == 1) == bool(
        [c for c in doc["checks"] if c["status"] == "fail"])


def test_the_closed_vocabularies_are_the_ones_the_constant_protects(ma):
    # `finding`, the stage names and every `detail` are OPEN and may grow
    # without a bump. These five are CLOSED: a consumer switches on them, so a
    # new value is a meaning change. Naming them here is what makes that rule
    # checkable rather than a paragraph in a comment.
    assert ma.DOCTOR_STATUSES == ("ok", "warn", "fail", "skip")
    assert ma.DOCTOR_REMEDIES == ("auto", "manual", "config", "input", "none")
    # Two values, not three: every "config contradiction" this vocabulary can
    # express turned out to kill something outright, so it is reported as
    # stage_or_command_dies with that stage OR that command named. A value
    # nothing emits is a branch nobody exercises.
    #
    # The name carries both halves because the old one - `stage_dies` - named
    # only the first, and was a plain lie on every row that kills a command
    # and no stage: `r:package:limma` and an absent `proteins_faa` are both
    # `blocks: []`.
    assert ma.DOCTOR_FAIL_REASONS == ("stage_or_command_dies",
                                      "setting_ignored")
    assert ma.DOCTOR_DEPTHS == ("config", "existence", "kind", "header",
                                "parsed", "probe")
    assert ma.DOCTOR_COMMANDS == ("run", "report", "object")
    src = open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index("DOCTOR_VERSION = ")
    rule = src[max(0, i - 2400):i]
    assert "CLOSED ENUM IS A MEANING CHANGE" in rule.upper(), \
        "the rule the constant carries is not written beside the constant"
    assert "adding a key is not a bump" in rule
    # Seven sets, not five: `found.kind` and `expect.kind` were outside the
    # versioned set while the README told a consumer to switch on one of them,
    # so a value could have been added to either with nothing to notice. The
    # rule has to name every set it governs, or naming any of them is decoration.
    for name in ("DOCTOR_STATUSES", "DOCTOR_REMEDIES", "DOCTOR_FAIL_REASONS",
                 "DOCTOR_DEPTHS", "DOCTOR_COMMANDS", "DOCTOR_FOUND_KINDS",
                 "DOCTOR_EXPECT_KINDS"):
        assert name in rule, f"the rule does not name {name}"


class _Args:
    """The argparse namespace cmd_doctor reads, for in-process calls."""

    def __init__(self, config=None, install_plan=None):
        self.config = config
        self.install_plan = install_plan
        self.json = True
        self.fix = False
        self.yes = False


# ----------------------------------------------------------------------
# --json alongside the other flags
# ----------------------------------------------------------------------
def test_fix_and_json_are_refused_together_and_the_status_says_which(project):
    # Refused by argparse, which exits 2, rather than by die(), which exits 1 -
    # and 1 is also "problems found", so a caller could not tell a refused
    # invocation from a failed check.
    proc = subprocess.run(
        [sys.executable, METAANNOT_PY, "doctor", "--config",
         project.config_path, "--json", "--fix"],
        capture_output=True, text=True)
    assert proc.returncode == 2, proc.stderr
    assert proc.stdout == "", "a refusal must not put half a document on stdout"
    assert "not allowed with" in proc.stderr


def test_install_plan_still_writes_its_file_under_json_and_is_recorded(
        tmp_path):
    # --install-plan is a side effect the operator asked for, not output, so
    # --json suppresses the sentence about it and not the file. The one line
    # that used to print would have landed in stdout beside the document.
    p = project_with(tmp_path, "plan")
    p.write_config(run=dict(p.cfg["run"], pfam=True),
                   db={"pfam_hmm": str(tmp_path / "absent" / "Pfam-A.hmm")})
    plan = tmp_path / "install.sh"
    proc, doc = doctor_json(p.config_path, "--install-plan", str(plan),
                            expect=1)
    assert proc.stdout.lstrip().startswith("{")
    assert "install plan ->" not in proc.stdout
    assert plan.exists() and "hmmpress" in plan.read_text(encoding="utf-8")
    rec = doc["install_plan"]
    assert rec["written"] and rec["path"] == str(plan)
    assert "pfam_hmm" in rec["covers"]


def test_install_plan_is_null_when_it_was_not_asked_for(project):
    _proc, doc = doctor_json(project.config_path, expect=None)
    assert doc["install_plan"] is None


# ----------------------------------------------------------------------
# sizes
# ----------------------------------------------------------------------
def test_the_totals_name_what_they_leave_out_rather_than_hiding_it(tmp_path):
    # `missing` is `not ok and not manual`, so both totals silently omit every
    # licence-gated tool and un-URLed database - InterProScan among them. A
    # total that quietly leaves out the biggest item is worse than no total.
    p = project_with(tmp_path, "sizes")
    p.write_config(run=dict(p.cfg["run"], interpro=True, pfam=True),
                   db={"pfam_hmm": str(tmp_path / "absent" / "Pfam-A.hmm")})
    _proc, doc = doctor_json(p.config_path, expect=1)
    t = doc["totals"]
    assert "interproscan" in t["not_counted"]
    assert "interproscan" not in t["counted"]
    assert "pfam_hmm" in t["counted"]
    assert t["download_gb"] > 0 and t["accuracy"] == "indicative"
    assert t["basis"] == "hand-maintained"
    assert "AFDB50" in t["detail"] or "123 GB" in t["detail"]


def test_an_unestimated_size_is_named_instead_of_being_read_as_free(tmp_path):
    # requirements() defaults both sizes to 0.0 and uses that for "negligible"
    # AND for "not estimated", and the printed report hides anything under
    # 0.05 GB - so a consumer reading 0.0 as "free" reads a large Java
    # distribution as free. The document says which zeros mean nothing.
    p = project_with(tmp_path, "unsized")
    p.write_config(run=dict(p.cfg["run"], interpro=True))
    _proc, doc = doctor_json(p.config_path, expect=1)
    t = doc["totals"]
    assert "interproscan" in t["unsized"]
    assert "not 'free'" in t["unsized_means"] or "never 'free'" in t["unsized_means"]
    reqs = {r["id"]: r for r in doc["requirements"]}
    assert reqs["interproscan"]["size_gb"] == 0.0, \
        "requirements[] is emitted verbatim; the caveat lives in totals"


# ----------------------------------------------------------------------
# one vocabulary, two renderings
# ----------------------------------------------------------------------
def test_the_printed_report_says_exactly_what_the_document_says(tmp_path):
    # `detail` IS the printed sentence. Two wordings of one verdict drift, and
    # the day they do, an operator reading the terminal and a console reading
    # the JSON disagree about the same config.
    p = project_with(tmp_path, "both")
    os.remove(p.cfg["quant_table"])
    p.write_config(run=dict(p.cfg["run"], pfam=True))
    human = run_metaannot("doctor", "--config", p.config_path, expect=1,
                          env=no_r_env())
    _proc, doc = doctor_json(p.config_path, expect=1)
    for c in doc["checks"]:
        assert c["detail"] in human.stdout, \
            f"{c['id']}'s sentence is not in the printed report"
    titles = {s["title"] for s in doc["sections"]}
    for title in titles:
        assert f"== {title} ==" in human.stdout


def test_manual_is_a_remedy_and_never_a_severity(tmp_path):
    # The printed mark fuses the two for the sake of one line, and that is the
    # ONLY place they are fused: a MANUAL item that an enabled stage dies
    # without is an ordinary failure everywhere else in the document. Treating
    # MANUAL as satisfied is what once let doctor say "all checks passed" with
    # no InterProScan and the run die at that stage hours later.
    p = project_with(tmp_path, "manual")
    p.write_config(run=dict(p.cfg["run"], interpro=True))
    human = run_metaannot("doctor", "--config", p.config_path, expect=1,
                          env=no_r_env())
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["req:interproscan"]
    assert c["status"] == "fail" and c["blocks"] == ["interpro"]
    assert c["remedy"] == "manual" and c["remedy_reason"]
    assert "MANUAL InterProScan" in human.stdout.replace("  ", " ")


def test_an_unfilled_database_path_asks_for_a_config_edit_not_a_licence(
        tmp_path):
    # requirements()[].manual is overloaded: dbentry() sets it to "no path
    # configured" for an unfilled key, while signalp6 and interproscan set it
    # to a licence or a version-specific distribution. The first needs an edit
    # box, the second a link-out and structurally no button.
    p = project_with(tmp_path, "unset_db")
    p.write_config(run=dict(p.cfg["run"], kofam=True),
                   db={"kofam_profiles": "", "kofam_ko_list": ""})
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["req:kofam_db"]
    assert c["remedy"] == "config", "an unfilled path is not a manual download"
    assert c["remedy_reason"] == "no path configured"
    assert c["blocks"] == ["kofam"]
    reqs = {r["id"]: r for r in doc["requirements"]}
    assert reqs["kofam_db"]["manual"] == "no path configured", \
        "requirements[].manual is left exactly as it was; the split is in remedy"


# ----------------------------------------------------------------------
# the R block, in process, with the answer supplied
# ----------------------------------------------------------------------
class _RAnswer:
    """A canned `Rscript -e requireNamespace(...)` reply."""

    def __init__(self, missing=()):
        self.returncode = 0
        self.stderr = ""
        self.stdout = ""
        self.missing = set(missing)

    def for_packages(self, names):
        self.stdout = " ".join(
            f"{n}={'FALSE' if n in self.missing else 'TRUE'}" for n in names)
        return self


def _r_checks_with(ma, monkeypatch, missing=()):
    want = list(ma.RNEED) + list(ma.ROPT)
    answer = _RAnswer(missing).for_packages(want)
    monkeypatch.setattr(ma, "have", lambda tool: True)
    monkeypatch.setattr(ma, "resolve_tool", lambda tool: tool)
    monkeypatch.setattr(ma.subprocess, "run",
                        lambda *a, **kw: answer)
    return {c["id"]: c for c in ma._r_checks(ma.DEFAULT_CONFIG)}


def test_each_r_package_is_its_own_row_with_its_own_install_line(ma,
                                                                 monkeypatch):
    # One row per package rather than one row for the set: the remedy differs
    # per package - BiocManager::install against install.packages, decided by
    # Bioconductor membership - and a preflight checklist offers its affordance
    # per row. Collapsing fourteen printed lines into two checks would make one
    # missing package a member of a check rather than a row of its own.
    ids = _r_checks_with(ma, monkeypatch, missing=["limma", "readr"])
    assert ids["r:package:limma"]["remedy_reason"] == \
        'BiocManager::install("limma")'
    assert ids["r:package:readr"]["remedy_reason"] == 'install.packages("readr")'
    assert len([i for i in ids if i.startswith("r:package:")]) == \
        len(ma.RNEED) + len(ma.ROPT)


def test_a_missing_required_r_package_blocks_the_commands_not_a_stage(
        ma, monkeypatch):
    # `report` and `object` are metaannot SUBCOMMANDS, not stages, so
    # describe --json's stage_names will never contain them - which is exactly
    # why they go in blocks_commands and never in blocks, where every name is
    # meant to join to that list.
    ids = _r_checks_with(ma, monkeypatch, missing=["limma"])
    c = ids["r:package:limma"]
    assert c["status"] == "fail"
    assert c["blocks"] == [] and c["blocks_commands"] == ["report", "object"]
    for name in c["blocks_commands"]:
        assert name not in ma.STAGE_NAMES


def test_a_missing_optional_r_package_only_warns(ma, monkeypatch):
    # The difference between RNEED and ROPT, made machine-readable: the report
    # is written either way, it just carries less.
    ids = _r_checks_with(ma, monkeypatch, missing=["patchwork"])
    c = ids["r:package:patchwork"]
    assert c["status"] == "warn"
    assert c["blocks_commands"] == [] and c["fails_reason"] is None


def test_the_r_rows_say_what_they_actually_asked(ma, monkeypatch):
    # requireNamespace answers whether the package can be FOUND, not whether it
    # loads. The caveat is per row rather than a footnote, because this is one
    # of the places the check goes past "the file is there".
    ids = _r_checks_with(ma, monkeypatch)
    c = ids["r:package:limma"]
    assert "requireNamespace" in c["caveat"]
    assert c["depth"] == "probe"


def test_no_rscript_is_a_warning_and_not_a_failure(ma, monkeypatch):
    # run.report/run.object still write the scripts, so they can be knitted on
    # the MacBook that has the R stack. A machine without R is not a broken
    # preflight.
    monkeypatch.setattr(ma, "have", lambda tool: False)
    checks = ma._r_checks(ma.DEFAULT_CONFIG)
    assert [c["id"] for c in checks] == ["r:rscript"]
    assert checks[0]["status"] == "warn"


# ----------------------------------------------------------------------
# the contingent verdict, and the boundary the scope sentence draws
# ----------------------------------------------------------------------
def test_a_missing_manifest_is_only_fatal_once_the_quant_table_is_there(
        tmp_path):
    # read_manifest is reached from read_feature_table and directly from
    # stage_join - there is no read_protein_table anywhere in this tool; the
    # protein-level path is inline in stage_join - and both are downstream of
    # stage_join's `if not os.path.exists(quant_table): return`. With the quant
    # table gone, nothing ever opens the manifest, so an unconditional failure
    # here is the same family of false-fail as the quant_table one.
    p = project_with(tmp_path, "man")
    os.remove(p.cfg["manifest"])
    os.remove(p.cfg["quant_table"])
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "warn" and c["degrades"] == ["join"]
    assert c["depends_on"] == ["input:quant_table"]
    assert "becomes fatal for join" in c["detail"]
    assert proc.returncode == 0

    # Restore the quant table and the same row turns fatal, exactly as the
    # detail said it would.
    F.write_peptide_table(p.cfg["quant_table"], p.proteins, p.samples)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "fail" and c["blocks"] == ["join"]


def test_every_check_that_looked_past_existence_says_what_it_read(tmp_path):
    # `depth` is a RATCHET rather than a menu, and it is what keeps "doctor
    # does not parse" auditable in a document that already contains a DIAMOND
    # header read, a plex annotation read, a read_manifest and a
    # pd.read_csv(nrows=0) on the quant header. A row deeper than a file's
    # existence has to name what it actually opened, per row, rather than
    # leaving the scope sentence to cover for it. How MANY of them there are is
    # DOCTOR_DEEP_CHECKS's business and nobody else's: this comment used to
    # count them by hand and say three.
    p = project_with(tmp_path, "depth")
    broken = tmp_path / "broken.dmnd"
    broken.write_bytes(b"")
    p.write_config(run=dict(p.cfg["run"], diamond=True),
                   db={"diamond": {"vfdb": str(broken)}})
    _proc, doc = doctor_json(p.config_path, expect=1)
    deep = [c for c in doc["checks"] if c["depth"] in ("header", "parsed")]
    assert deep, "the fixture no longer exercises a check past existence"
    for c in deep:
        assert c["caveat"], f"{c['id']} went past existence and said nothing"
    # ...and nothing in the document has anywhere to put a parse RESULT. No
    # rows, no columns, no record counts, no protein ids: the shape itself is
    # what stops this command growing into a different program.
    forbidden = {"rows", "columns", "records", "ids", "sequences", "n_rows"}
    for c in doc["checks"]:
        assert not (set(c) & forbidden), c["id"]
        assert not (set(c["found"] or {}) & forbidden), c["id"]


def test_doctor_json_works_with_no_config_at_all(tmp_path):
    # `describe --json` with no --config is a documented invocation and the
    # console's startup call; doctor's has to parse the same way, because a
    # front end asking "what would this build need by default?" gets there
    # first.
    proc, doc = doctor_json(None, expect=None, cwd=str(tmp_path))
    assert doc["config_path"] is None
    assert doc["checks"] and doc["doctor_version"] == 1
    assert proc.returncode == doc["exit_status"]


def test_a_real_plex_tree_is_read_and_the_rows_say_they_read_it(tmp_path):
    # The `== tmt ==` block against the layout this format is most particular
    # about: a run directory of per-plex folders, each with its own annotation
    # file. These rows are `depth: "parsed"` - they open every annotation - and
    # that is one of the checks in DOCTOR_DEEP_CHECKS, so each one names what
    # it opened.
    F.tmt_planted_run(str(tmp_path / "run"))
    p = project_with(tmp_path, "tmtreal", quant_format="fragpipe_tmt",
                     quant_table=str(tmp_path / "run"), manifest="")
    proc, doc = doctor_json(p.config_path, expect=0)
    ids = by_id(doc)
    assert ids["input:quant_table"]["status"] == "ok"
    assert ids["tmt:plex_dirs"]["status"] == "ok"
    counts = ids["tmt:channel_counts"]
    assert counts["status"] == "ok" and counts["depth"] == "parsed"
    assert counts["caveat"], "a row that read every annotation must say so"
    # No reference channel named is a warning about what the join will be able
    # to do, not a failure: the reader handles it.
    assert ids["tmt:reference"]["status"] == "warn"
    assert ids["tmt:reference"]["degrades"] == ["join"]
    assert proc.returncode == 0


def test_a_plex_missing_its_level_file_blocks_the_stage_that_reads_the_tree(
        tmp_path):
    F.tmt_planted_run(str(tmp_path / "run"))
    victim = tmp_path / "run" / "TMT2" / "ion.tsv"
    assert victim.exists(), "the fixture's layout moved"
    victim.unlink()
    p = project_with(tmp_path, "tmtbroken", quant_format="fragpipe_tmt",
                     quant_table=str(tmp_path / "run"), manifest="")
    proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["tmt:TMT2:level_file"]
    assert c["status"] == "fail" and c["blocks"] == ["join"]
    assert proc.returncode == 1


# ----------------------------------------------------------------------
# group A: rows that disagreed with the stage they describe
#
# Every test here drives the STAGE (or the command) into the state the row
# claims, and then reads the row. A test that only asserted the row would have
# passed against every one of these defects, because each of them was a
# consistent, well-formed, confidently wrong document.
# ----------------------------------------------------------------------
def test_a_taxon_rank_needs_the_taxdump_even_with_run_taxonomy_off(tmp_path):
    # The row was gated on run.taxonomy alone, so on this config it did not
    # EXIST - and a requirement with no row cannot be wrong about anything,
    # which is why this was a silent exit 0 rather than a visible bad verdict.
    # requirement_effect()'s taxon_rank condition was reading an entry nothing
    # had produced.
    p = project_with(tmp_path, "rank", taxon_rank="genus",
                     db={"ncbi_taxonomy": ""})
    proc, doc = doctor_json(p.config_path, expect=1)
    row = by_id(doc)["req:ncbi_taxonomy"]
    assert row["status"] == "fail" and row["blocks"] == ["join"]
    assert row["remedy"] == "config", "an unfilled path, not a licence"
    assert proc.returncode == 1


def test_the_taxdump_is_not_required_when_no_rank_is_asked_for(tmp_path):
    # The other direction, and it is what keeps the fix from being "always
    # emit it": with no rank, collapse_taxon_rank returns the raw seed taxids
    # and nothing opens a taxdump, so there is nothing to report and nothing
    # to fail.
    p = project_with(tmp_path, "norank", taxon_rank="")
    proc, doc = doctor_json(p.config_path, expect=0)
    assert "req:ncbi_taxonomy" not in by_id(doc)
    assert proc.returncode == 0


def test_resolve_taxonomy_really_dies_on_a_rank_with_no_taxdump(
        ma, tmp_path, paths_for):
    # The claim, driven: this is the call stage_join makes on every run -
    # resolve_taxonomy -> collapse_taxon_rank - so a config that sets a rank
    # and no taxdump dies here whether or not the comparison stage is on.
    cfg, p = paths_for()
    cfg["db"]["ncbi_taxonomy"] = ""
    ann = pd.DataFrame({"protein_id": ["P1"], "seed_taxid": ["1280"]})
    cfg["taxon_rank"] = "genus"
    with pytest.raises(ma.StageError) as e:
        ma.resolve_taxonomy(cfg, p, ann)
    assert "needs db.ncbi_taxonomy" in str(e.value)
    cfg["taxon_rank"] = ""
    assert ma.resolve_taxonomy(cfg, p, ann) == {"P1": "1280"}


def test_a_zero_byte_unipept_result_is_not_a_result(tmp_path):
    # The row asserts `depth: "kind"` - "also asked what kind of thing it is
    # and whether it is empty" - and then branched on `found.present`, which a
    # zero-byte file and a directory both satisfy. It computed the answer and
    # threw it away.
    p = project_with(tmp_path, "ures")
    res = tmp_path / "pept2lca.tsv"
    res.write_bytes(b"")
    p.write_config(run=dict(p.cfg["run"], unipept=True),
                   unipept={"result": str(res)})
    proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["taxonomy:unipept_result"]
    assert c["status"] == "fail" and c["finding"] == "empty"
    assert c["found"]["kind"] == "empty_file" and c["depth"] == "kind"
    assert c["blocks"] == ["unipept"]
    assert proc.returncode == 1


def test_a_directory_at_the_unipept_cache_is_not_a_cache(tmp_path):
    p = project_with(tmp_path, "ucachedir")
    p.write_config(run=dict(p.cfg["run"], unipept=True),
                   unipept={"result": "", "allow_http": False})
    cache = os.path.join(p.results, "unipept", "pept2lca_cache.tsv")
    os.makedirs(cache)
    proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["taxonomy:unipept_cache"]
    assert c["status"] == "fail" and c["finding"] == "wrong_kind"
    assert c["found"]["kind"] == "empty_dir" and c["depth"] == "kind"
    assert c["blocks"] == ["unipept"]
    assert proc.returncode == 1


def test_the_unipept_result_rows_quote_the_refusal_the_reader_really_raises(
        ma, tmp_path, paths_for):
    """The reader is driven, and the doctor row is held against what it said.

    os.path.exists() is stage_unipept's whole test, so a zero-byte export is
    ingested rather than queried and read_unipept_result() is what meets the
    file next - which is the fact both halves of this assert.

    THE QUOTE IS THE POINT. The row used to promise pandas' "Could not
    determine delimiter", which was true only while read_unipept_result()
    reached the file through `pd.read_csv(sep=None, engine="python")`. It goes
    through opener() now, and pandas says "No columns to parse from file"
    instead - so the row was quoting an error nothing in this program could
    produce, and nothing noticed, because the row's sentence and the reader's
    sentence were two independent strings. They are one constant now, and this
    drives the reader to prove the constant is really what comes out.
    """
    cfg, p = paths_for()
    empty = tmp_path / "pept2lca.tsv"
    empty.write_bytes(b"")
    cfg["unipept"] = dict(cfg["unipept"], result=str(empty))
    with pytest.raises(ma.StageError) as e:
        ma.stage_unipept(cfg, p)
    assert ma.EMPTY_PEPT2LCA in str(e.value), \
        "the reader, not the existence check"
    assert str(empty) in str(e.value), "the refusal has to name the file"
    assert not os.path.exists(p.unipept_lca)

    proj = project_with(tmp_path, "emptyres")
    proj.write_config(run=dict(proj.cfg["run"], unipept=True),
                      unipept={"result": str(empty)})
    _proc, doc = doctor_json(proj.config_path, expect=1)
    row = by_id(doc)["taxonomy:unipept_result"]
    assert ma.EMPTY_PEPT2LCA in row["detail"], \
        "the row quotes a refusal the reader does not raise"


def test_stage_unipept_really_dies_on_a_zero_byte_cache(tmp_path):
    # One step further in: a zero-byte cache is not adopted (getsize > 0), so
    # every peptide is still to do, and with allow_http false that is fatal -
    # the message this row exists to pre-empt.
    p = project_with(tmp_path, "ucache")
    p.write_config(run=dict(p.cfg["run"], unipept=True),
                   unipept={"result": "", "allow_http": False})
    cache = os.path.join(p.results, "unipept", "pept2lca_cache.tsv")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    open(cache, "wb").close()
    proc = p.run(expect=1)
    assert "the cache is incomplete" in proc.stderr


def test_a_manifest_is_not_claimed_by_a_unipept_that_reads_a_result_file(
        tmp_path):
    # The same species as the stage_join false verdict this build was written
    # to fix, in the block that had no equivalent of it: stage_unipept returns
    # after ingesting unipept.result and never opens a quant table, so it
    # never reaches read_manifest either.
    p = project_with(tmp_path, "manres")
    res = tmp_path / "pept2lca.tsv"
    res.write_text("peptide\tgenus_id\nPEPTIDE\t1\n", encoding="utf-8")
    os.remove(p.cfg["manifest"])
    p.write_config(run=dict(p.cfg["run"], unipept=True, join=False,
                            taxonomy=False),
                   unipept={"result": str(res)})
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "skip" and c["blocks"] == []
    assert proc.returncode == 0, "doctor failed a config `run` exits 0 on"


def test_stage_unipept_really_returns_before_it_could_open_a_manifest(
        ma, tmp_path, paths_for):
    # The claim above, driven: the manifest and the quant table are both
    # nowhere, and the stage still finishes.
    cfg, p = paths_for()
    res = tmp_path / "pept2lca.tsv"
    res.write_text("peptide\ttaxon_id\tgenus_id\nPEPTIDE\t1280\t1279\n",
                   encoding="utf-8")
    cfg["manifest"] = str(tmp_path / "nowhere.fp-manifest")
    cfg["quant_table"] = str(tmp_path / "nowhere.tsv")
    cfg["unipept"] = dict(cfg["unipept"], result=str(res))
    ma.stage_unipept(cfg, p)                    # no StageError, no read
    assert os.path.exists(p.unipept_lca)


def test_a_manifest_is_not_claimed_on_a_format_whose_reader_never_opens_one(
        tmp_path):
    # read_fragpipe_tmt logs "`manifest` is ignored for quant_format
    # 'fragpipe_tmt'" and takes its sample names from each plex's annotation
    # file. Nothing opens the manifest, so a missing one kills nothing.
    F.tmt_planted_run(str(tmp_path / "run"))
    p = project_with(tmp_path, "tmtman", quant_format="fragpipe_tmt",
                     quant_table=str(tmp_path / "run"))
    os.remove(p.cfg["manifest"])
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "skip" and c["finding"] == "not_used"
    assert c["blocks"] == [] and "not used" in c["detail"]
    assert proc.returncode == 0


def test_a_tmt_reference_conflict_blocks_join_and_not_the_peptide_readers(
        ma, tmp_path):
    # read_fragpipe_tmt dies on the pair; read_fragpipe_tmt_peptides never
    # looks at a reference, and peptide_features() falls back to it. So the
    # row may name join and may not name the two stages that survive.
    root, _truth = F.tmt_planted_run(str(tmp_path / "run"))
    p = project_with(tmp_path, "tmtref", quant_format="fragpipe_tmt",
                     quant_table=root, manifest="",
                     tmt={"reference_name": "Pool*", "reference_channel": "126"})
    # allow_http so the unipept cache row is a warning: the verdict under test
    # is the reference conflict's, not that stage's.
    p.write_config(run=dict(p.cfg["run"], unipept=True, taxonomy=True,
                            join=True),
                   unipept=dict(p.cfg.get("unipept") or {}, allow_http=True))
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["tmt:reference_conflict"]
    assert c["status"] == "fail" and c["blocks"] == ["join"]
    assert "peptides" in c["detail"], "the row says who survives it"

    # ...and both halves of that, driven.
    cfg = ma.load_config(p.config_path)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(root, "fragpipe_tmt", cfg)
    assert "reference_name" in str(e.value)
    feats = ma.peptide_features(cfg, "unipept")
    assert len(feats) and "peptide" in feats.columns


def test_a_tmt_refusal_is_a_warning_when_join_is_the_stage_that_is_off(
        tmp_path):
    # The same tree and the same conflict with run.join off: nothing dies, so
    # the row must not say `fail` - and a `fail` with an empty `blocks` would
    # trip _check()'s own invariant and take the command down.
    root, _truth = F.tmt_planted_run(str(tmp_path / "run"))
    p = project_with(tmp_path, "tmtref2", quant_format="fragpipe_tmt",
                     quant_table=root, manifest="",
                     tmt={"reference_name": "Pool*", "reference_channel": "126"})
    p.write_config(run=dict(p.cfg["run"], unipept=True, taxonomy=True,
                            join=False),
                   unipept=dict(p.cfg.get("unipept") or {}, allow_http=True))
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["tmt:reference_conflict"]
    assert c["status"] == "warn" and c["blocks"] == []
    assert proc.returncode == 0


def test_esmfold_is_satisfied_by_the_backend_the_stage_falls_back_to(
        ma, tmp_path):
    # STAGES gives esmfold `requires=("esmfold",)`, and stage_esmfold tries
    # fair-esm and THEN transformers, dying only when both fail - so probing
    # "esm" alone failed a host where ESMFold runs. On an sm_120 card that is
    # the ordinary host: fair-esm's esmfold extra needs an openfold whose CUDA
    # kernels do not compile there at all.
    #
    # The stage half of this is already driven, in tests/test_esmfold.py: its
    # `folding` fixture puts `esm` in sys.modules as None and a fake
    # transformers beside it - "fair-esm must not be importable, so the
    # transformers branch is taken" - and every test in that file then folds
    # through it. This is the half that was missing: the requirement agreeing
    # with the stage that already works that way.
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["run"] = dict(cfg["run"], structure=True)
    cfg["results_dir"] = str(tmp_path / "r")
    got = {}
    for have_esm, have_hf in ((True, False), (False, True), (False, False)):
        seen = {"esm": have_esm, "transformers": have_hf}
        orig = ma._pyhas
        try:
            ma._pyhas = lambda mod: seen.get(mod, orig(mod))
            reqs = ma.requirements(cfg, ma.Paths(cfg))
        finally:
            ma._pyhas = orig
        got[(have_esm, have_hf)] = next(r["ok"] for r in reqs
                                        if r["id"] == "esmfold")
    assert got[(True, False)] and got[(False, True)]
    assert not got[(False, False)], "both absent is what the stage dies on"


def test_a_directory_at_proteins_faa_does_not_claim_that_run_refuses(tmp_path):
    # cmd_run's guard is os.path.exists(), which a directory passes: `run`
    # schedules everything and each stage dies as it opens the path. The row
    # folded absent, dangling and directory into one branch whose wording was
    # written for "absent", so the machine-readable claim - the one a
    # preflight screen gates the `run` button on - was false.
    p = project_with(tmp_path, "faadir")
    faa = p.cfg["proteins_faa"]
    os.remove(faa)
    os.makedirs(faa)
    proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["input:proteins_faa"]
    assert c["blocks_commands"] == [], "`run` does not refuse a directory"
    assert c["blocks"] == ["emapper", "integrate"]
    assert c["finding"] == "wrong_kind" and c["depth"] == "kind"
    assert proc.returncode == 1


def test_run_really_refuses_a_missing_proteins_faa_and_not_a_directory(
        tmp_path):
    # The COMMAND, driven both ways. `run --dry-run` is past the guard at the
    # top of cmd_run and before anything executes, so it answers exactly the
    # question the row makes a claim about.
    p = project_with(tmp_path, "faacmd")
    faa = p.cfg["proteins_faa"]
    os.remove(faa)
    os.makedirs(faa)
    ok = run_metaannot("run", "--dry-run", "--config", p.config_path,
                       expect=0, cwd=p.root)
    assert "stage" in ok.stdout, "the plan printed: `run` did not refuse"
    os.rmdir(faa)
    gone = run_metaannot("run", "--dry-run", "--config", p.config_path,
                         expect=1, cwd=p.root)
    assert "proteins_faa not found" in gone.stderr


def test_an_empty_directory_is_a_wrong_kind_rather_than_an_empty_file(
        tmp_path):
    # `empty` and `wrong_kind` carry different remedies, and _finding() used
    # to answer `empty` for an empty DIRECTORY - indistinguishable in the
    # document from a zero-byte file, which on this key is a warning and not a
    # failure at all.
    p = project_with(tmp_path, "gffdir", gff=str(tmp_path / "g"))
    os.makedirs(tmp_path / "g")
    run_flags(p, context=True)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["input:gff"]
    assert c["found"]["kind"] == "empty_dir"
    assert c["finding"] == "wrong_kind" and c["status"] == "fail"
    assert c["depth"] == "kind" and c["blocks"] == ["context"]


def test_stage_context_really_dies_on_a_directory_it_was_given_as_a_gff(
        ma, tmp_path, paths_for):
    cfg, p = paths_for()
    d = tmp_path / "gffdir"
    os.makedirs(d)
    cfg["gff"] = str(d)
    with pytest.raises((ma.StageError, OSError)):
        ma.stage_context(cfg, p)             # past os.path.exists, into parse


def test_a_protein_level_format_blocks_taxonomy_whatever_unipept_result_says(
        tmp_path):
    # The guard was `and not u.get("result")`, which is right for unipept and
    # wrong for taxonomy: stage_unipept returns before its own format check,
    # stage_taxonomy's is guarded on nothing at all.
    p = project_with(tmp_path, "protfmt", quant_format="fragpipe")
    res = tmp_path / "pept2lca.tsv"
    res.write_text("peptide\tgenus_id\nPEPTIDE\t1\n", encoding="utf-8")
    p.write_config(run=dict(p.cfg["run"], unipept=True, taxonomy=True),
                   unipept={"result": str(res)})
    proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["taxonomy:format"]
    assert c["status"] == "fail" and c["blocks"] == ["taxonomy"]
    assert "unipept.result does not answer this one" in c["detail"]
    assert proc.returncode == 1


def test_stage_taxonomy_really_dies_on_a_protein_level_format(
        ma, tmp_path, paths_for):
    # Driven past the check before it - a unipept_lca is on disk, so the
    # format check is what this dies on.
    cfg, p = paths_for()
    cfg["quant_format"] = "fragpipe"
    with open(p.unipept_lca, "w", encoding="utf-8") as fh:
        fh.write("peptide\ttaxon_id\tgenus_id\nPEPTIDE\t1280\t1279\n")
    with pytest.raises(ma.StageError) as e:
        ma.stage_taxonomy(cfg, p)
    assert "peptide-level input" in str(e.value)


# ----------------------------------------------------------------------
# group B: rows that described themselves wrongly
# ----------------------------------------------------------------------
def test_the_diamond_row_says_which_file_it_actually_read(
        ma, tmp_path, monkeypatch):
    # With diamond off PATH - the ordinary state the first time anyone runs
    # doctor - the check falls through to the source FASTA and reads up to
    # 200,000 records for a median, and motif_seed_evidence reads 200 deflines
    # after it. The hard-coded caveat said "the sequences themselves were not
    # parsed", which is the opposite of what just happened.
    p = project_with(tmp_path, "diadepth")
    db = F.write_dmnd(tmp_path / "bagel.dmnd")
    # BAGEL's motif seed set, which is the database that prompted this check:
    # 15-residue entries, and the only way to see them is to read the FASTA.
    with open(tmp_path / "bagel.fas", "w", encoding="utf-8") as fh:
        for i in range(5):
            fh.write(">ggmotif-%d\n%s\n" % (i, "M" * (14 + i)))
    p.write_config(run=dict(p.cfg["run"], diamond=True),
                   db={"diamond": {"bagel": str(db)}})
    cfg = ma.load_config(p.config_path)
    monkeypatch.setattr(ma, "have", lambda t: t != "diamond")
    rows = ma._diamond_checks(cfg, set(ma.enabled_stages(cfg)))
    c = next(r for r in rows if r["id"] == "db:diamond:bagel:usable")
    assert c["status"] == "warn" and "motif or seed set" in c["detail"]
    assert c["depth"] == "parsed", "it read the FASTA, records and all"
    assert "200,000 records" in c["caveat"] and "deflines" in c["caveat"]
    assert "DIAMOND header" not in c["caveat"], \
        "the .dmnd was never opened: diamond is not on PATH"


def test_the_interproscan_row_names_the_key_it_is_really_tested_with(
        tmp_path):
    # requirements() tests _exists(db.interproscan_sh) and never looks at
    # PATH, so "found on PATH by have()" sent an operator to copy a shell
    # script onto PATH. And an unfilled key is a CONFIG remedy, not a licence:
    # _remedy_for exists to make that split, and this row went round it by
    # setting manual= directly.
    p = project_with(tmp_path, "ipr", db={"interproscan_sh": ""})
    run_flags(p, interpro=True)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["req:interproscan"]
    assert c["expect"]["kind"] == "file"
    assert "db.interproscan_sh" in c["expect"]["because"]
    assert "never looks at PATH" in c["expect"]["because"]
    assert c["config_keys"] == ["db.interproscan_sh"]
    assert c["remedy"] == "config" and c["remedy_reason"] == "no path configured"
    assert c["depth"] == "config" and c["found"] is None


def test_a_set_interproscan_path_that_is_not_there_is_the_licence_case(
        tmp_path):
    p = project_with(tmp_path, "ipr2",
                     db={"interproscan_sh": str(tmp_path / "nowhere.sh")})
    run_flags(p, interpro=True)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["req:interproscan"]
    assert c["remedy"] == "manual"
    assert c["depth"] == "kind" and c["found"]["present"] is False


def test_the_esmfold_row_says_it_probed_an_import_and_not_a_path(tmp_path):
    p = project_with(tmp_path, "esm")
    run_flags(p, structure=True)
    _proc, doc = doctor_json(p.config_path, expect=None)
    c = by_id(doc)["req:esmfold"]
    assert c["expect"]["kind"] == "probe" and c["depth"] == "probe"
    assert "_pyhas" in c["expect"]["because"]
    assert "transformers" in c["expect"]["because"]


def test_an_unconfigured_database_path_claims_no_depth_past_the_config(
        tmp_path):
    # The ratchet, downward. Nothing outside the config was consulted for a
    # key nobody filled in, and `depth: "kind"` claimed a path had been
    # stat'd - which is the same kind of untrue as a caveat that names a read
    # that did not happen.
    p = project_with(tmp_path, "nopath", db={"pfam_hmm": ""})
    run_flags(p, pfam=True)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["req:pfam_hmm"]
    assert c["depth"] == "config" and c["found"] is None
    assert c["remedy"] == "config"


def test_a_row_cannot_claim_a_depth_it_did_not_reach(ma):
    # Structural, like the caveat rule beside it: the guard is in _check(), so
    # a row that claims to have stat'd something without carrying what that
    # answered is a loud failure here rather than a quiet lie in the document.
    with pytest.raises(AssertionError) as e:
        ma._check("x:y", "config", "x", "warn", "ok", "detail", depth="kind")
    assert "no `found` to show for it" in str(e.value)
    # ...and an UNSET `found` is not evidence of a stat either. This is the
    # arm that let a `skip` row claim depth "kind" over a config key nobody
    # filled in: _found("") consults nothing at all.
    with pytest.raises(AssertionError):
        ma._check("x:y", "config", "x", "skip", "not_enabled", "d",
                  depth="kind", found=ma._found(""))
    # ...and the upward half still holds.
    with pytest.raises(AssertionError):
        ma._check("x:y", "config", "x", "warn", "ok", "d", depth="parsed",
                  found=ma._found(None))


def test_a_row_that_statted_a_path_cannot_report_depth_config(ma, tmp_path):
    # The guard the ratchet did NOT have, and the half the `skip` rows were
    # getting through: "config" is a positive claim that nothing outside the
    # config was consulted, and a row carrying a real `found` has consulted
    # something. Every skipped input row - a gff with run.context off, a
    # quant_table nothing reads, an emapper_precomputed entry - was stat'ing
    # its path to fill `found` and then reporting "config".
    real = tmp_path / "f"
    real.write_text("x", encoding="utf-8")
    with pytest.raises(AssertionError) as e:
        ma._check("x:y", "inputs", "x", "skip", "not_enabled", "d",
                  depth="config", found=ma._found(str(real)))
    assert "something outside the config WAS consulted" in str(e.value)
    # ...and an unset path is still "config", because _found() really did
    # consult nothing for it.
    row = ma._check("x:y", "inputs", "x", "skip", "not_enabled", "d",
                    depth="config", found=ma._found(""))
    assert row["depth"] == "config"


def test_the_depth_of_every_skip_row_is_the_one_it_really_reached(tmp_path):
    # The guard above, driven through the rows it was written for rather than
    # through _check() directly: a `skip` row over a path that EXISTS says
    # "kind", and the same row over a config key nobody filled in says
    # "config". Both used to say "config".
    p = project_with(tmp_path, "depths")
    (tmp_path / "real.gff").write_text("##gff-version 3\n", encoding="utf-8")
    (tmp_path / "real.fna").write_text(">c\nACGT\n", encoding="utf-8")
    p.write_config(gff=str(tmp_path / "real.gff"),
                   contigs_fna=str(tmp_path / "real.fna"))
    _proc, doc = doctor_json(p.config_path, expect=None)
    ids = by_id(doc)
    for cid in ("input:gff", "input:contigs_fna"):
        c = ids[cid]
        assert c["status"] == "skip", f"{cid} is no longer the skipped case"
        assert c["depth"] == "kind", f"{cid} stat'd its path and said config"
        assert c["found"]["kind"] == "file"
    p.write_config(gff="", contigs_fna="")
    _proc, doc = doctor_json(p.config_path, expect=None)
    for cid in ("input:gff", "input:contigs_fna"):
        c = by_id(doc)[cid]
        assert c["depth"] == "config" and c["found"]["kind"] == "unset"


def test_the_two_vocabularies_that_were_outside_the_versioned_set(ma,
                                                                  tmp_path):
    # `found.kind` is what _found()'s docstring and the README both tell a
    # consumer to switch on instead of computing `bytes > 0`, and
    # `expect.kind` was enumerated in no constant, no document and no test -
    # so either could have grown a value with nothing to notice, which is the
    # difference between a contract and a suggestion.
    assert ma.DOCTOR_FOUND_KINDS == ("file", "empty_file", "dir", "empty_dir",
                                     "symlink_broken", "absent", "unset",
                                     "other", "unreadable")
    assert ma.DOCTOR_EXPECT_KINDS == ("file", "dir", "on_path", "probe",
                                      "setting", "r_package")
    # every state _found() can name is in the set, checked against the
    # function rather than against a copy of the list
    d = tmp_path / "d"
    os.makedirs(d / "full")
    open(d / "full" / "x", "w", encoding="utf-8").write("x")
    open(d / "file", "w", encoding="utf-8").write("x")
    open(d / "empty_file", "w", encoding="utf-8").close()
    os.makedirs(d / "empty_dir")
    os.symlink(str(d / "nowhere"), str(d / "dangling"))
    os.mkfifo(str(d / "pipe"))
    os.makedirs(d / "locked")
    os.chmod(d / "locked", 0o000)
    names = ("full", "file", "empty_file", "empty_dir", "dangling", "absent",
             "pipe", "locked")
    try:
        seen = {ma._found(str(d / n))["kind"] for n in names} \
            | {ma._found("")["kind"]}
    finally:
        os.chmod(d / "locked", 0o755)
    assert seen == set(ma.DOCTOR_FOUND_KINDS)


def test_a_fifo_is_present_because_os_path_exists_says_it_is(ma, tmp_path):
    # The whole reason `other` exists. Every stage in metaannot guards its
    # inputs with os.path.exists(), which is TRUE for a FIFO - so a row that
    # reported one as `absent` said "run refuses before any stage starts" and
    # `blocks_commands: ["run"]` about a path `run` is perfectly happy to
    # schedule against. `present` has to be os.path.exists()'s own answer, or
    # every sentence built on it is about a different file.
    fifo = tmp_path / "pipe"
    os.mkfifo(str(fifo))
    f = ma._found(str(fifo))
    assert f["kind"] == "other"
    assert f["present"] is True and os.path.exists(str(fifo))
    assert ma._finding(f) == "wrong_kind", "a FIFO is not a missing file"


def test_a_fifo_at_proteins_faa_does_not_claim_that_run_refuses(tmp_path):
    # The same false claim A5 fixed for directories, one path state over, and
    # driven the same way: cmd_run's only test is os.path.exists(), which a
    # FIFO passes, so `run` schedules everything rather than refusing.
    p = project_with(tmp_path, "faafifo")
    os.remove(p.cfg["proteins_faa"])
    os.mkfifo(p.cfg["proteins_faa"])
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["input:proteins_faa"]
    assert c["found"]["kind"] == "other" and c["found"]["present"] is True
    assert c["blocks_commands"] == [], "`run` does not refuse a FIFO"
    assert c["blocks"], "the stages that open it are what dies"
    assert "refuses before any stage starts" not in c["detail"]


def test_a_directory_that_will_not_list_itself_is_unreadable_not_absent(
        ma, tmp_path):
    # The other state that used to come out as `absent`: os.listdir raises,
    # the old `except OSError` answered "kind: absent, present: False", and an
    # operator was sent looking for a directory sitting exactly where they put
    # it. The remedy is a permission, and `unreadable` is how the row says so.
    locked = tmp_path / "locked"
    os.makedirs(locked)
    os.chmod(locked, 0o000)
    try:
        f = ma._found(str(locked))
    finally:
        os.chmod(locked, 0o755)
    assert f["kind"] == "unreadable" and f["present"] is True
    assert ma._finding(f) == "unreadable"


@pytest.mark.parametrize("which", CONFIGS)
def test_every_kind_a_document_emits_is_in_its_closed_set(ma, tmp_path, which):
    # `None` is in the accepted set on purpose and is part of the contract: a
    # row may carry a non-null `found` whose `kind` is null, meaning "this
    # build did not compute it". This assertion used to be the only place that
    # fact was written down, where it read as a loophole rather than a value -
    # see test_found_kind_may_be_null_and_the_document_says_which_rows_do_that.
    p = _configured(tmp_path, which)
    _proc, doc = doctor_json(p.config_path, expect=None)
    for c in doc["checks"]:
        assert (c["found"] or {}).get("kind") in (None,) + ma.DOCTOR_FOUND_KINDS
        if c["expect"]:
            assert c["expect"]["kind"] in ma.DOCTOR_EXPECT_KINDS
        # `probe` means the same thing in both fields now: a row that asks the
        # HOST. A row that consulted only the config says `setting`, which is
        # what the collision was.
        if c["expect"] and c["expect"]["kind"] == "probe":
            assert c["depth"] == "probe", f"{c['id']} probes nothing"


def test_expect_members_names_the_siblings_the_engine_tests_for(ma, tmp_path):
    # A permanently-empty field in a versioned contract is worse than no
    # field: every case it was written for was being described in prose by
    # _why_requirement and left empty here.
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["run"] = {k: True for k in cfg["run"]}
    cfg["emapper_precomputed"] = ""
    cfg["db"]["ncbi_taxonomy"] = str(tmp_path / "tax")
    cfg["results_dir"] = str(tmp_path / "r")
    reqs = ma.requirements(cfg, ma.Paths(cfg))
    rows = {c["requirement_id"]: c for c in ma._requirement_checks(cfg, reqs)}
    assert rows["pfam_hmm"]["expect"]["members"] == [".h3f", ".h3i", ".h3m",
                                                     ".h3p"]
    assert rows["foldseek_target"]["expect"]["members"] == [".dbtype", ".index"]
    assert rows["ncbi_taxonomy"]["expect"]["members"] == ["nodes.dmp"]
    assert rows["ncbi_taxonomy"]["expect"]["kind"] == "dir"
    # ...and empty where the engine's own test names no sibling, with the row
    # saying so rather than leaving it to be read as "nobody filled this in".
    assert rows["hhblits_db"]["expect"]["members"] == []
    assert "members` is empty" in rows["hhblits_db"]["expect"]["because"]


def test_the_hmmpress_members_are_the_ones_pressed_really_requires(ma,
                                                                   tmp_path):
    # The tripwire that keeps the list above from being decoration: it is
    # _pressed()'s own requirement, driven.
    hmm = tmp_path / "Pfam-A.hmm"
    hmm.write_text("HMMER3/f\n", encoding="utf-8")
    members = [".h3f", ".h3i", ".h3m", ".h3p"]
    for s in members:
        (tmp_path / ("Pfam-A.hmm" + s)).write_text("x", encoding="utf-8")
    assert ma._pressed(str(hmm))
    for s in members:
        os.remove(str(hmm) + s)
        assert not ma._pressed(str(hmm)), f"{s} is not really required"
        (tmp_path / ("Pfam-A.hmm" + s)).write_text("x", encoding="utf-8")


def _fna_dir_project(tmp_path, name):
    """A project whose contigs_fna is a directory, with run.smorf on."""
    p = project_with(tmp_path, name, contigs_fna=str(tmp_path / "asm"))
    os.makedirs(tmp_path / "asm", exist_ok=True)
    open(tmp_path / "asm" / "c.fna", "w", encoding="utf-8").write(">c\nACGT\n")
    return run_flags(p, smorf=True)


def _orf_finder_on_path(tmp_path, name="macrel"):
    """PATH with a stub ORF finder that refuses whatever it is handed."""
    binp = tmp_path / "orfbin"
    os.makedirs(binp, exist_ok=True)
    exe = binp / name
    exe.write_text("#!/bin/sh\nexit 3\n", encoding="utf-8")
    os.chmod(exe, 0o755)
    env = no_r_env()
    return {"PATH": str(binp) + os.pathsep + env["PATH"]}


def test_a_directory_at_contigs_fna_blocks_smorf_only_when_an_orf_finder_is_there(
        tmp_path):
    # The verdict here is not a property of the PATH state alone, and the row
    # used to claim it was: `blocks: ["smorf"]` unconditionally, pinned by a
    # test that ran on a machine with no smorf, no smorfinder and no macrel -
    # so the suite green-lit a verdict that is false in the environment the
    # suite runs in. Measured on that machine: doctor exited 1 and `run`
    # exited 0.
    #
    # With no ORF finder installed stage_smorf never opens the assembly at
    # all: it logs, writes an empty smorf_proteins.faa and returns. The
    # assertion is therefore made against the PROBE, driven both ways, and not
    # against whatever this laptop happens to have.
    #
    # Both halves now supply their own PATH. This one used to pass `no_r_env()`
    # - the machine's real PATH minus Rscript - so the "no ORF finder" half was
    # only true while nobody running the suite had macrel installed, which is
    # the same defect one level up: a test that asserts a verdict by relying on
    # a property of the laptop.
    p = _fna_dir_project(tmp_path, "fnadir")
    _proc, doc = doctor_json(p.config_path, expect=0, env=no_orf_finder_env())
    c = by_id(doc)["input:contigs_fna"]
    assert c["finding"] == "wrong_kind" and c["depth"] == "kind"
    assert c["status"] == "warn" and c["blocks"] == []
    assert c["depends_on"] == ["req:smorf"]
    assert "no ORF finder is installed" in c["detail"]

    # ...and the same directory with one on PATH is fatal, which is what the
    # row always said.
    _proc, doc = doctor_json(p.config_path, expect=1,
                             env=_orf_finder_on_path(tmp_path))
    c = by_id(doc)["input:contigs_fna"]
    assert c["status"] == "fail" and c["blocks"] == ["smorf"]
    assert c["blocks_commands"] == []


def test_the_contigs_fna_row_does_not_promise_what_an_orf_finder_does(tmp_path):
    """The row hands the path over; it does not claim to know how the tool fails.

    This sentence was corrected because it had asserted that smorf or macrel
    "cannot read an assembly out of it" - a claim about a third-party program
    that nothing here establishes and that is false on a FIFO, where the tool
    sits blocked in its own open() while this run holds the results directory.
    The correction was then shipped with NOTHING pinning it: reverting the
    sentence left all 1521 tests green, which is the same gap that let the
    manifest promise regress, so the row is pinned here rather than trusted.

    Asserted as a rule and not as a quotation: what may not come back is a
    claim about the TOOL's behaviour, and what must come back is that the path
    is handed over by name. Rewording is free; re-promising is not.
    """
    p = _fna_dir_project(tmp_path, "promise")
    _proc, doc = doctor_json(p.config_path, expect=1,
                             env=_orf_finder_on_path(tmp_path))
    detail = by_id(doc)["input:contigs_fna"]["detail"]

    # It says whose answer it is, and that this file never opens the path.
    assert "FILENAME" in detail and "never opens it here" in detail, detail
    assert "that tool's answer" in detail, detail
    # ...and it does not describe an outcome it cannot have measured. "refuses"
    # is the retired verb; the FIFO clause may still say the tool can BLOCK,
    # which is the one thing that was actually driven.
    for promise in ("cannot read an assembly out of it",
                    "is refused by smorf", "is refused by macrel",
                    "the tool refuses"):
        assert promise not in detail, (promise, detail)


def test_the_proteins_faa_row_does_not_say_a_search_tool_dies_on_a_pipe(
        tmp_path):
    """The half of that sentence that was false, and why it is a separate row.

    `input:proteins_faa` said every stage "dies on this, in read_fasta() or in
    the search tool it was handed to". The read_fasta() half is true. The tool
    half was measured false: a FIFO there leaves hmmsearch blocked in its own
    open() and `run` waiting for a stage it has already said it cannot
    interrupt, with the results directory still locked - so the row promised a
    death where the real outcome is a run that never ends. The hang itself is
    older than this contract and is not doctor's to fix; the SENTENCE is.
    """
    p = project_with(tmp_path, "faapipe")
    os.remove(p.cfg["proteins_faa"])
    os.mkfifo(p.cfg["proteins_faa"])
    _proc, doc = doctor_json(p.config_path, expect=1)
    detail = by_id(doc)["input:proteins_faa"]["detail"]

    assert "read_fasta()" in detail, detail
    assert "FILENAME" in detail and "the tool's answer" in detail, detail
    assert "blocked in its own open()" in detail, detail
    # The retired form: a single verb covering both readers at once.
    assert "in read_fasta() or in the search tool" not in detail, detail


def test_stage_smorf_really_survives_a_directory_with_no_orf_finder_installed(
        tmp_path):
    # The stage itself, driven into the state the row above describes, both
    # ways round. This is the assertion the old test could not make, because
    # it asserted the verdict rather than the behaviour behind it.
    p = _fna_dir_project(tmp_path, "fnarun")
    p.run(env=no_orf_finder_env(), expect=0)
    faa = p.rpath("smorf", "smorf_proteins.faa")
    assert os.path.exists(faa) and os.path.getsize(faa) == 0, \
        "the stage was meant to write an EMPTY candidate list and return"
    # ...and with an ORF finder that refuses the directory, the stage dies.
    shutil.rmtree(p.results)
    proc = p.run(env=_orf_finder_on_path(tmp_path), expect=1)
    assert "smorf" in proc.stderr


def test_an_msstats_table_carries_its_own_design_and_needs_no_manifest(
        tmp_path):
    # The other half of _manifest_readers: the MSstats readers never call
    # read_manifest at all, so a manifest that is set and gone kills nothing -
    # and the row says which format decided that rather than leaving an
    # operator to infer it.
    p = project_with(tmp_path, "msstats", quant_format="msstats_csv")
    os.remove(p.cfg["manifest"])
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "skip" and c["blocks"] == []
    assert "msstats_csv" in c["detail"] and "own design" in c["detail"]
    assert "manifest:mapping" not in by_id(doc), \
        "nothing maps a manifest onto columns on this format"
    assert proc.returncode == 0


def test_peptide_only_reader_is_what_decides_who_a_refusal_kills(tmp_path):
    # The switch behind `_full_reader_refusal`, both ways round on one input.
    # A manifest whose runs match no quant column is a refusal by the FULL
    # reader: peptide_features() catches it and re-reads with the peptide-only
    # reader, so unipept lives - unless peptide_only_reader says never, which
    # is the mode that means refuse rather than fall back.
    p = project_with(tmp_path, "pormode")
    with open(p.cfg["manifest"], "w", encoding="utf-8") as fh:
        fh.write("/nowhere/ZZZ.mzML\tZZZ\t1\tDDA\n")
    common = dict(run=dict(p.cfg["run"], unipept=True, join=False),
                  unipept={"allow_http": True})
    p.write_config(peptide_only_reader="auto", **common)
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["manifest:mapping"]
    assert c["status"] == "warn" and c["blocks"] == []
    assert "peptides only" in c["detail"]
    assert proc.returncode == 0

    p.write_config(peptide_only_reader="never", **common)
    proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["manifest:mapping"]
    assert c["status"] == "fail" and c["blocks"] == ["unipept"]
    assert proc.returncode == 1


def test_the_full_reader_really_refuses_and_the_peptide_reader_really_does_not(
        ma, tmp_path):
    # That claim, driven: same table, same manifest, two readers.
    p = project_with(tmp_path, "pordrive")
    with open(p.cfg["manifest"], "w", encoding="utf-8") as fh:
        fh.write("/nowhere/ZZZ.mzML\tZZZ\t1\tDDA\n")
    cfg = ma.load_config(p.config_path)
    with pytest.raises(ma.StageError) as e:
        ma.read_feature_table(cfg["quant_table"], cfg["quant_format"], cfg)
    assert "match no column" in str(e.value)
    feats = ma.peptide_features(cfg, "unipept")      # falls back, no raise
    assert len(feats) and "peptide" in feats.columns
    cfg["peptide_only_reader"] = "never"
    with pytest.raises(ma.StageError):
        ma.peptide_features(cfg, "unipept")


# ----------------------------------------------------------------------
# group N: the verdicts a second reading against the engine still caught
#
# Same discipline as group A above: every test here drives the stage, the
# command or the reader into the state the row claims, and the assertion is
# about what THAT does, not about what the row says it does.
# ----------------------------------------------------------------------
INPUT_PATH_KEYS = ("proteins_faa", "quant_table", "manifest", "gff",
                   "contigs_fna", "emapper_precomputed")
# `emapper_precomputed` was outside this list while the entry describing the
# sweep said EVERY path `doctor` can be pointed at - which is the same shape
# of wrong claim as a stale count, one surface over. It is a LIST in the
# config rather than a path, which is the only reason it was ever awkward to
# include, and _input_target() is that awkwardness written down once.
# `unreadable_file` is the state this list did not have, and its absence is
# why four rows passed a chmod-000 input for three rounds. `unreadable_dir`
# was here from the start, so `found.kind: "unreadable"` looked covered -
# while _found() could only ever REACH it for a directory: os.path.getsize()
# is a stat, and a stat answers happily for a mode-000 regular file, so the
# commoner half of the state had no state here to drive it with.
# `socket` and `device_node` are the two this list did not have, and their
# absence is why a published sentence could be false for two thirds of the
# state it was written about. `found.kind` is `other` for a FIFO, a socket AND
# a device node, and the derived verb keyed on `kind` alone - so every row
# printed the FIFO's paragraph over all three, and the only one of the three
# the sweep could create was the one it was true of.
PATH_STATES = ("absent", "empty_file", "directory", "dangling_symlink",
               "fifo", "socket", "device_node",
               "unreadable_dir", "unreadable_file")


def _input_target(p, key):
    """The real path behind one configured input key.

    `emapper_precomputed` is a LIST of paths and every other key is one path;
    a sweep over "every configured input" has to know that in one place or it
    knows it nowhere, which is how that key stayed outside the sweep while the
    prose said the sweep covered everything.
    """
    v = p.cfg[key]
    return v[0] if isinstance(v, list) else v


def _input_row_id(key):
    """The row a configured input key writes.

    Every key but one writes `input:<key>`; `emapper_precomputed` is a LIST,
    so its rows are `precomputed_emapper:<i>` - one per entry. The sweep has
    to know which row to look for, and knowing it here is what let the key
    into the sweep at all.
    """
    return ("precomputed_emapper:0" if key == "emapper_precomputed"
            else f"input:{key}")


def _put_in_state(path, state):
    """Leave `path` in one of the states a configured input can be in."""
    _clear(path)
    if state == "absent":
        return
    if state == "empty_file":
        open(path, "w", encoding="utf-8").close()
    elif state == "directory":
        os.makedirs(path)
    elif state == "dangling_symlink":
        os.symlink(path + ".nowhere", path)
    elif state == "fifo":
        os.mkfifo(path)
    elif state == "socket":
        # Bound from INSIDE the directory: an AF_UNIX path is capped at about
        # a hundred bytes and a pytest tmp_path is most of that already. The
        # inode stays on disk when the socket object is closed, which is what
        # the sweep needs - the state is the FILE, not a listener.
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        cwd = os.getcwd()
        try:
            os.chdir(os.path.dirname(path))
            sock.bind(os.path.basename(path))
        finally:
            os.chdir(cwd)
            sock.close()
    elif state == "device_node":
        # A SYMLINK to a device node, because mknod needs root and this has to
        # run on any machine. os.stat() follows it, so _found() sees a
        # character device exactly as it would see a real one - and the
        # symlink flag on the row is the only difference, which is itself
        # worth having in the sweep.
        os.symlink("/dev/zero", path)
    elif state == "unreadable_dir":
        os.makedirs(path)
        os.chmod(path, 0o000)
    elif state == "unreadable_file":
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("whatever this file is for\n")
        os.chmod(path, 0o000)
    else:                                        # pragma: no cover - typo
        raise AssertionError(f"no such state {state}")


def _clear(path):
    """Remove whatever is at `path`, including a mode-000 directory."""
    if not (os.path.islink(path) or os.path.exists(path)):
        return
    if os.path.isdir(path) and not os.path.islink(path):
        os.chmod(path, 0o755)
        shutil.rmtree(path)
    else:
        os.remove(path)


def _restore(path, state):
    """Give a chmod-000 path its mode back, so tmp_path can be cleaned up.

    In a `finally`, always: a test that fails while a directory is mode 000
    leaves pytest unable to remove its own tmp_path, and the next failure is
    then about the harness rather than about the code.
    """
    if state in ("unreadable_dir", "unreadable_file"):
        try:
            os.chmod(path, 0o755 if state == "unreadable_dir" else 0o644)
        except OSError:                          # pragma: no cover - gone
            pass


@pytest.mark.parametrize("key", INPUT_PATH_KEYS)
@pytest.mark.parametrize("state", PATH_STATES)
def test_no_state_of_any_input_path_can_stop_the_document_being_printed(
        tmp_path, key, state):
    """symptom: a DIRECTORY at `manifest` produced NO DOCUMENT AT ALL.

    `_manifest_checks` branched on `found.present`, which a directory
    satisfies, and fell through to `read_manifest(path)` guarded only by
    `except StageError`. read_manifest opens through opener(), which is open(),
    so a directory raised IsADirectoryError: rc=1, a traceback on stderr, and
    nothing on stdout. That is the most expensive failure this command can
    have, because a consumer gets no verdict about the tools, the databases,
    the resources or R either - every other row goes down with the one that
    could not be read, and there is nothing in the document to be wrong, so
    nothing can notice.

    Parametrised over every configured input path and every state one can be
    in rather than over the one that was found, because the next such state
    will be at a different key and the point is that none of them can do this.
    """
    p = project_with(tmp_path, f"{key[:4]}_{state}")
    p.write_config(gff=p.path("input", "x.gff"),
                   contigs_fna=p.path("input", "x.fna"),
                   run=dict(p.cfg["run"], context=True, smorf=True))
    target = _input_target(p, key)
    _put_in_state(target, state)
    try:
        proc, doc = doctor_json(p.config_path, expect=None, timeout=90)
    finally:
        _restore(target, state)
    assert proc.stdout.lstrip().startswith("{"), \
        f"{key} as a {state} left stdout unparseable:\n{proc.stderr[-800:]}"
    assert json.loads(proc.stdout) == doc
    assert doc["exit_status"] == proc.returncode
    row = by_id(doc).get(_input_row_id(key))
    assert row is not None, f"{key} lost its row entirely"


def test_a_directory_at_the_manifest_fails_and_names_what_dies(tmp_path):
    # The split itself, rather than only the survival above: a directory is
    # its own branch BEFORE read_manifest is reached, with its own finding and
    # its own list of what it costs. read_manifest raises an OSError, which is
    # not the full reader REFUSING - so peptide_features() does not catch it
    # and nothing falls back to the peptide-only reader.
    p = project_with(tmp_path, "mandir")
    os.remove(p.cfg["manifest"])
    os.makedirs(p.cfg["manifest"])
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "fail" and c["finding"] == "wrong_kind"
    assert c["depth"] == "kind" and c["blocks"] == ["join"]
    assert "OSError" in c["detail"]
    assert by_id(doc)["manifest:mapping"]["status"] == "skip"


def test_stage_join_really_dies_on_a_manifest_that_is_a_directory(tmp_path):
    # The claim above, driven. read_manifest goes straight to opener(), which
    # is open(), so there is no refusal to catch and no fallback to take.
    p = project_with(tmp_path, "mandirrun")
    os.remove(p.cfg["manifest"])
    os.makedirs(p.cfg["manifest"])
    proc = p.run(env=no_r_env(), expect=1)
    assert "join" in proc.stderr and "Is a directory" in proc.stderr


def test_a_directory_at_the_manifest_costs_nothing_while_the_quant_table_is_gone(
        tmp_path):
    # The contingency is the absent branch's, for the absent branch's reason:
    # nothing opens the manifest until something has opened the quant table.
    # Driven both ways, and the run is what says which.
    p = project_with(tmp_path, "mandircontingent")
    os.remove(p.cfg["manifest"])
    os.makedirs(p.cfg["manifest"])
    os.remove(p.cfg["quant_table"])
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "warn" and c["blocks"] == []
    assert c["depends_on"] == ["input:quant_table"]
    assert proc.returncode == 0
    p.run(env=no_r_env(), expect=0)


def test_an_empty_manifest_costs_join_and_not_the_peptide_readers(tmp_path):
    # The other half of the split, and a different list: an EMPTY manifest
    # reaches read_manifest and dies inside it with a StageError ("no rows"),
    # which peptide_features() catches and re-reads around - so it costs join
    # and spares the taxonomy stages, exactly as an unreadable one does.
    p = project_with(tmp_path, "manempty")
    open(p.cfg["manifest"], "w", encoding="utf-8").close()
    run_flags(p, unipept=True)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "fail" and c["finding"] == "empty"
    assert c["blocks"] == ["join"], "unipept re-reads with the peptide reader"
    assert "unipept" in c["detail"]


def test_read_manifest_really_refuses_an_empty_one_and_raises_on_a_directory(
        ma, tmp_path):
    # The two exceptions the split is built on, from the function itself: one
    # is StageError, which the fallback catches, and the other is an OSError,
    # which nothing does.
    empty = tmp_path / "empty.fp-manifest"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(ma.StageError) as e:
        ma.read_manifest(str(empty))
    assert "no rows" in str(e.value)
    d = tmp_path / "dir.fp-manifest"
    os.makedirs(d)
    with pytest.raises(OSError) as oe:
        ma.read_manifest(str(d))
    assert not isinstance(oe.value, ma.StageError), \
        "the whole point is that this one is NOT a StageError"


def test_the_taxdump_does_not_block_a_join_that_skips_for_want_of_a_quant_table(
        tmp_path):
    """symptom: doctor exited 1 with `req:ncbi_taxonomy` on a config `run`
    completes, and run.taxonomy is off by DEFAULT so this is the more common
    half of the configs.

    stage_join's FIRST act is `if not os.path.exists(quant_table): log("quant
    table not found, skipping join"); return`, before resolve_taxonomy and
    before anything else it requires. doctor already held this contingency -
    `input:quant_table` says `degrades: ["join"]` on the same document, and
    `_manifest_checks` gates its whole verdict on it - and applied it nowhere
    near the requirement rows.
    """
    p = project_with(tmp_path, "taxskip", taxon_rank="genus")
    os.remove(p.cfg["quant_table"])
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["req:ncbi_taxonomy"]
    assert c["status"] == "warn" and c["blocks"] == []
    assert c["depends_on"] == ["input:quant_table"], \
        "the row that decided this verdict has to be named"
    assert by_id(doc)["input:quant_table"]["degrades"] == ["join"]
    # ...and `run` agrees, which is the whole claim.
    assert proc.returncode == 0
    run = p.run(env=no_r_env(), expect=0)
    assert "skipping join" in run.stderr

    # Restore the quant table and the same row turns fatal, because now join
    # gets as far as resolve_taxonomy.
    F.write_peptide_table(p.cfg["quant_table"], p.proteins, p.samples)
    _proc, doc = doctor_json(p.config_path, expect=1)
    c = by_id(doc)["req:ncbi_taxonomy"]
    assert c["status"] == "fail" and c["blocks"] == ["join"]
    assert c["depends_on"] == []


def test_stage_join_really_returns_before_it_consults_a_taxon_rank(tmp_path):
    # The engine half of the row above: with the quant table gone, join never
    # reaches collapse_taxon_rank at all, so the taxdump it would have died
    # without costs nothing.
    p = project_with(tmp_path, "joinskiprank", taxon_rank="genus")
    os.remove(p.cfg["quant_table"])
    proc = p.run(env=no_r_env(), expect=0)
    assert "skipping join" in proc.stderr
    assert "needs db.ncbi_taxonomy" not in proc.stderr


def test_a_zero_byte_gff_really_leaves_a_run_that_finishes(tmp_path):
    """symptom: doctor exited 0 on a zero-byte gff and `run` exited 1 with
    "FATAL stage 'integrate' failed: No columns to parse from file".

    The row said an empty GFF writes "a table with no rows, exactly as the
    unset case does". It did not: the unset branch wrote
    pd.DataFrame(columns=["protein_id"]), and a gff that parsed to nothing fell
    through to pd.DataFrame([]).to_csv(), which writes one newline and NO
    COLUMNS - so context.tsv was one byte, nonempty() was true, and
    parse_context's pd.read_csv raised EmptyDataError. `integrate` always runs.
    The engine is the half that was wrong, and this is the end-to-end pin.
    """
    gff = tmp_path / "empty.gff"
    gff.write_text("", encoding="utf-8")
    p = project_with(tmp_path, "gffempty", gff=str(gff))
    run_flags(p, context=True)
    proc, doc = doctor_json(p.config_path, expect=0)
    c = by_id(doc)["input:gff"]
    assert c["status"] == "warn" and c["blocks"] == []
    assert proc.returncode == 0
    p.run(env=no_r_env(), expect=0)


@pytest.mark.parametrize("body", ["", "##gff-version 3\n##sequence-region c1 1 9\n"])
def test_stage_context_writes_a_table_with_columns_when_the_gff_has_no_features(
        ma, tmp_path, body, monkeypatch):
    """The engine defect itself, at both of the two doors into it.

    A zero-byte gff is one; a HEADER-ONLY gff with no CDS records is the other,
    and that one is outside doctor's scope entirely - the file is neither
    missing nor empty, so doctor reports `ok` for it and must. Which is why the
    fix had to be the engine's: no verdict could have covered it.
    """
    gff = tmp_path / "g.gff"
    gff.write_text(body, encoding="utf-8")
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["results_dir"] = str(tmp_path / "r")
    cfg["gff"] = str(gff)
    pth = ma.Paths(cfg)
    pth.mkdirs()
    ma.stage_context(cfg, pth)
    with open(pth.context, encoding="utf-8") as fh:
        first = fh.readline().rstrip("\n")
    assert first == "protein_id", \
        f"an empty context table has no columns: {first!r}"
    # ...and the reader the run actually uses gets through it.
    assert ma.parse_context(pth.context) == {}


def test_an_unset_gff_and_one_that_parses_to_nothing_write_the_same_file(
        ma, tmp_path):
    # "Exactly as the unset case does" is the row's own sentence, so it is
    # asserted as an equality rather than as two separate shapes.
    def context_of(gff_value):
        cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
        cfg["results_dir"] = str(tmp_path / f"r{abs(hash(gff_value))}")
        cfg["gff"] = gff_value
        pth = ma.Paths(cfg)
        pth.mkdirs()
        ma.stage_context(cfg, pth)
        return open(pth.context, encoding="utf-8").read()

    empty = tmp_path / "e.gff"
    empty.write_text("", encoding="utf-8")
    assert context_of("") == context_of(str(empty))


# ----------------------------------------------------------------------
# group P: what the published contract says about itself
# ----------------------------------------------------------------------
def test_the_verdict_rule_in_the_document_names_the_command_class_too(ma):
    """symptom: `verdict.rule` - the sentence a CONSUMER parses, because it is
    in the document - said a check fails exactly when an enabled STAGE dies on
    it, or the declared `setting_ignored` exception. Neither covers
    `r:package:limma`, which is status "fail", blocks [], blocks_commands
    ["report", "object"] and no stage anywhere. doctor()'s own docstring and
    the README's closing sentence had it right, so three renderings of one rule
    disagreed and the one a consumer parses was the loose one.
    """
    rule = ma.DOCTOR_VERDICT_RULE
    assert "blocks_commands" in rule, "the rule does not name the field"
    for name in ma.DOCTOR_COMMANDS:
        assert name in rule, f"the rule does not name the {name} command"
    assert "stage_or_command_dies" in rule
    assert "setting_ignored" in rule
    # ...and the constant's own comment names both halves as well, because it
    # is where the next person writes the next copy of this sentence.
    src = open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index("DOCTOR_FAIL_REASONS = ")
    comment = src[max(0, i - 2000):i]
    assert "blocks_commands" in comment


def test_the_r_rows_obey_the_same_invariants_the_swept_configs_do(
        ma, monkeypatch):
    """The sweep's blind spot, closed where it was created.

    no_r_env() strips Rscript from PATH for every subprocess doctor call in
    this file - deliberately, so an exit status is a property of the config
    and not of the laptop - and the consequence is that NO R FAILURE IS EVER
    IN A SWEPT DOCUMENT. The R block is where `blocks_commands` lives, so a
    whole sweep could pass without one row that fails while killing no stage:
    precisely the class three of the four statements of the verdict rule had
    dropped. The invariants are therefore run over the R rows here, in
    process, with the answer supplied.
    """
    ids = _r_checks_with(ma, monkeypatch, missing=["limma", "patchwork"])
    saw_command_failure = False
    for c in ids.values():
        where = c["id"]
        named = bool(c["blocks"] or c["blocks_commands"])
        assert (c["status"] == "fail") == (c["fails_reason"] is not None), where
        assert (c["fails_reason"]
                == "stage_or_command_dies") == named, where
        assert c["status"] in ma.DOCTOR_STATUSES, where
        assert c["remedy"] in ma.DOCTOR_REMEDIES, where
        assert c["depth"] in ma.DOCTOR_DEPTHS, where
        assert c["fails_reason"] in (None,) + ma.DOCTOR_FAIL_REASONS, where
        for name in c["blocks_commands"]:
            assert name in ma.DOCTOR_COMMANDS, where
            assert name not in ma.STAGE_NAMES, where
        if c["status"] == "fail":
            assert c["blocks"] == [], f"{where}: no stage dies for R"
            saw_command_failure = True
    assert saw_command_failure, \
        "the fixture no longer produces a failure that kills only a command"


def test_a_tool_row_says_what_was_really_probed_and_what_it_really_costs(
        ma, tmp_path, monkeypatch):
    """symptom: `expect.because` said "found on PATH by have(); every stage
    that uses it dies first" for EVERY tool row - and on `req:smorf` it said
    that beside `blocks: []` and `degrades: ["smorf"]`, so one object
    contradicted itself. stage_smorf logs "smorf/smorfinder not found;
    skipping that half" and carries on.

    Both halves of the sentence are derived now: the probe from the table that
    a second test drives against requirements(), and the cost from the same
    requirement_effect() call that fills `blocks` and `degrades`.
    """
    monkeypatch.setattr(ma, "have", lambda tool: False)
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["run"] = {k: True for k in cfg["run"]}
    cfg["emapper_precomputed"] = ""
    cfg["results_dir"] = str(tmp_path / "r")
    reqs = ma.requirements(cfg, ma.Paths(cfg))
    rows = {c["requirement_id"]: c for c in ma._requirement_checks(cfg, reqs)}

    smorf = rows["smorf"]
    assert smorf["blocks"] == [] and smorf["degrades"] == ["smorf"]
    because = smorf["expect"]["because"]
    assert "dies first" not in because, "the row still contradicts itself"
    assert "nothing dies without it" in because
    for name in ("smorf", "smorfinder", "macrel"):
        assert f"`{name}`" in because, f"{name} is a name the stage accepts"

    # ...and a tool that really is fatal still says so, with the stage named.
    hmmer = rows["hmmer"]
    assert hmmer["blocks"], "the fixture no longer has a fatal tool"
    assert "dies first" in hmmer["expect"]["because"]
    for st in hmmer["blocks"]:
        assert st in hmmer["expect"]["because"]


def test_every_tool_probe_the_document_describes_is_the_one_requirements_makes(
        ma, tmp_path, monkeypatch):
    """The table against the engine, one executable name at a time.

    A hand-written sentence about a probe is exactly the thing that drifts
    from the probe, which is how `req:kofamscan` came to claim `blocks:
    ["kofam"]` on a host where stage_kofam runs: requirements() probed
    `exec_annotation` alone while the stage accepts `exec_annotation` OR
    `kofamscan`. So each name is asserted to be SUFFICIENT on its own, and the
    key set is asserted to be every tool requirements() can produce.
    """
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["run"] = {k: True for k in cfg["run"]}
    cfg["emapper_precomputed"] = ""
    cfg["results_dir"] = str(tmp_path / "r")

    def ids_ok(present):
        monkeypatch.setattr(ma, "have", lambda tool: tool in present)
        return {r["id"]: r["ok"]
                for r in ma.requirements(cfg, ma.Paths(cfg))
                if r["kind"] == "tool"}

    all_tools = set(ids_ok(set()))
    off_path = set(ma.REQUIREMENT_OFF_PATH_PROBES)
    assert set(ma.REQUIREMENT_TOOL_PROBES) == all_tools - off_path, \
        "the probe table and the tool requirements have come apart"
    for rid, names in ma.REQUIREMENT_TOOL_PROBES.items():
        for name in names:
            assert ids_ok({name})[rid], \
                f"{rid}: `{name}` is documented as enough and is not"


def test_stage_kofam_accepts_the_second_name_requirements_now_probes_for(ma,
                                                                        monkeypatch):
    # The engine half: the stage's own fallback, which is what makes the
    # second name in the table a fact rather than a hope.
    seen = {}
    monkeypatch.setattr(ma, "have", lambda tool: tool == "kofamscan")
    monkeypatch.setattr(ma, "run_cmd", lambda cmd, **kw: seen.setdefault(
        "exe", cmd[0]))
    monkeypatch.setattr(ma, "atomic_out", __import__("contextlib").contextmanager(
        lambda path: iter([str(path) + ".tmp"])))
    monkeypatch.setattr(ma.os, "makedirs", lambda *a, **kw: None)
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["results_dir"] = "/nonexistent"
    ma.stage_kofam(cfg, ma.Paths(cfg))
    assert seen["exe"] == "kofamscan", \
        "stage_kofam no longer falls back to the second name"


def test_the_scope_statement_and_every_caveat_agree_on_how_many_read_deep(ma):
    """symptom: the document said "three checks" in PUBLISHED text - the
    caveat on every `tmt:<plex>:annotation` row - while its own
    `scope.statement` said four, and "go past a header" was not the same
    predicate either, since the quant-table check reads the header line and
    nothing else. Four places were counting by hand.
    """
    n = len(ma.DOCTOR_DEEP_CHECKS)
    assert ma.DOCTOR_DEEP_COUNT == ma._count_word(n)
    assert f"{ma.DOCTOR_DEEP_COUNT.capitalize()} checks do read into a file" \
        in ma.DOCTOR_SCOPE_STATEMENT
    for name in ma.DOCTOR_DEEP_CHECKS:
        assert name in ma.DOCTOR_SCOPE_STATEMENT
    assert ma.DOCTOR_DEEP_COUNT in ma.DOCTOR_DEEP_CLAUSE
    # ...and nothing in the file writes the count as a word of its own beside
    # the word "checks", which is how the two came apart.
    src = open(METAANNOT_PY, encoding="utf-8").read()
    for word in ("three", "four", "five"):
        for phrase in (f"the {word} checks", f"{word.capitalize()} checks"):
            assert phrase not in src, \
                f"{phrase!r} is a hand-written count; derive it"


def test_the_tmt_annotation_caveat_is_the_derived_one(tmp_path):
    # The published half of the test above: the row itself, out of a real
    # document, rather than the constant it is built from.
    F.tmt_planted_run(str(tmp_path / "run"))
    bad = tmp_path / "run" / "TMT1" / "TMT1_annotation.txt"
    assert bad.exists(), "the fixture's layout moved"
    bad.write_text("", encoding="utf-8")     # read_tmt_annotation refuses it
    p = project_with(tmp_path, "tmtcaveat", quant_format="fragpipe_tmt",
                     quant_table=str(tmp_path / "run"), manifest="")
    _proc, doc = doctor_json(p.config_path, expect=None)
    rows = [c for c in doc["checks"] if c["id"].endswith(":annotation")]
    assert rows, "the fixture no longer produces an annotation row"
    for c in rows:
        assert "three checks" not in c["caveat"]
        assert doc["scope"]["statement"].count("checks do read into a file")
        assert c["caveat"].endswith(
            "one of the four checks in this command that read INTO a file")


def test_found_kind_may_be_null_and_the_document_says_which_rows_do_that(
        ma, tmp_path):
    """symptom: the README enumerated seven values for `found.kind` "so a
    consumer never has to do arithmetic on bytes", and _check() explicitly
    admits an eighth - its guard is `not in (None,) + DOCTOR_FOUND_KINDS`.
    In the DEFAULT document the null is the MAJORITY of the rows, because
    every requirement row, both DIAMOND rows, the CUDA probe and the whole R
    block answer with one boolean and cannot say what kind of thing is there.
    The null is deliberate; the defect was that nothing said so.
    """
    p = project_with(tmp_path, "nullkind")
    p.write_config(run=dict(p.cfg["run"], pfam=True, structure=True))
    _proc, doc = doctor_json(p.config_path, expect=None)
    nulls = [c for c in doc["checks"]
             if c["found"] is not None and c["found"]["kind"] is None]
    assert nulls, "the fixture no longer produces a non-null found with a " \
                  "null kind, which is the case the documentation is about"
    for c in nulls:
        assert isinstance(c["found"]["present"], bool), \
            f"{c['id']}: a null kind leaves `present` as the only answer"
    # The constant's comment is where a consumer's author looks after the
    # README, so it has to say it too.
    src = open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index("DOCTOR_FOUND_KINDS = ")
    comment = " ".join(src[max(0, i - 2500):i].replace("#", " ").split())
    assert "NULL IS" in comment.upper()
    assert "this build did not compute it" in comment


def test_an_empty_manifest_with_join_off_says_nothing_dies_and_means_it(
        tmp_path):
    # The corner the contingent sentence has to be guarded on. With join off
    # and unipept on, the full reader's refusal costs nobody - the peptide
    # readers re-read around it - so there is nothing for a restored quant
    # table to make fatal, and the row must not promise one.
    p = project_with(tmp_path, "manemptynojoin")
    open(p.cfg["manifest"], "w", encoding="utf-8").close()
    run_flags(p, join=False, unipept=True)
    _proc, doc = doctor_json(p.config_path, expect=None)
    c = by_id(doc)["input:manifest"]
    assert c["status"] == "warn" and c["blocks"] == []
    assert "becomes fatal for ." not in c["detail"]
    assert "becomes fatal" not in c["detail"]
    assert "peptide-only reader" in c["detail"]
    # (run.unipept with no cache and no allow_http fails on its own row, which
    # is not this one; the claim here is only that the manifest is not in it.)
    assert "input:manifest" not in doc["verdict"]["fails"]


# ----------------------------------------------------------------------
# group S: the path-state class, ended structurally rather than one at a time
#
# Three rounds of this change set each found another place where `doctor`
# opened something it had not established was a regular readable file, and
# each fix was local to the site that had just been caught: a directory at
# `manifest`, then a widened `except` beside it, then the identical class
# again at the TMT annotations - where a chmod-000 file printed a traceback
# and a FIFO made the command HANG INDEFINITELY.
#
# So these tests pin the STRUCTURE and not the three instances: one gate that
# every check reading into a file goes through, every path `doctor` can be
# pointed at walked through every state a path can be in, and the two
# properties that matter more than any individual verdict - a document on
# stdout, and a process that EXITS.
# ----------------------------------------------------------------------
def _tmt_project(tmp_path, name, plexes=("TMT1",)):
    """A project whose quant_table is a FragPipe TMT run directory.

    The label-free fixture cannot reach the TMT block at all, which is why
    none of its paths were ever swept: `CONFIGS`' `tmt_file` entry leaves
    `quant_table` a FILE, so `tmt:root` skips and the annotation reader is
    never called in any swept document. This builds the tree the reader wants.
    """
    p = project_with(tmp_path, name)
    root = p.path("input", "tmtrun")
    for plex in plexes:
        d = os.path.join(root, plex)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "ion.tsv"), "w", encoding="utf-8") as fh:
            fh.write("Peptide Sequence\tCharge\tProtein\t126\t127N\n")
            fh.write("PEPTIDE\t2\tP1\t100\t200\n")
        with open(os.path.join(d, f"{plex}_annotation.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write("126 S1\n127N S2\n")
    p.write_config(quant_table=root, quant_format="fragpipe_tmt")
    return p


# Every path inside the TMT tree that `doctor` is pointed at, by the name this
# file calls it. `annotation` is the one that hung; `plex_dir` is the one whose
# os.listdir() in an f-string took the document down; `level_file` had
# os.path.exists() as its only test, like every other site in this class.
TMT_READ_PATHS = ("plex_dir", "level_file", "annotation")


def _tmt_target(p, which, plex="TMT1"):
    root = p.cfg["quant_table"]
    pdir = os.path.join(root, plex)
    return {"plex_dir": pdir,
            "level_file": os.path.join(pdir, "ion.tsv"),
            "annotation": os.path.join(pdir, f"{plex}_annotation.txt")}[which]


@pytest.mark.parametrize("which", TMT_READ_PATHS)
@pytest.mark.parametrize("state", PATH_STATES)
def test_no_state_of_the_tmt_tree_stops_the_document_or_the_process(
        tmp_path, which, state):
    """symptom: a FIFO at a plex's annotation file made `doctor` HANG.

    Not fail - hang. `tmt_annotation_path()` guards only os.path.exists(),
    `read_tmt_annotation()` goes straight to `opener()`, which is `open()`,
    and a read-only open of a FIFO with no writer blocks forever. No document
    AND no exit, which is worse than any traceback: a caller waiting on the
    process has nothing to time out against either. A chmod-000 annotation
    file and a directory at the same path each printed a traceback and zero
    bytes of stdout, and a chmod-000 plex DIRECTORY did it from an
    `os.listdir()` inside an f-string building a `detail`.

    The timeout is the assertion about the hang: without it this test would
    have wedged the suite instead of failing it.
    """
    p = _tmt_project(tmp_path, f"tmt_{which[:4]}_{state[:6]}")
    target = _tmt_target(p, which)
    _put_in_state(target, state)
    try:
        proc, doc = doctor_json(p.config_path, expect=None, timeout=90)
    except subprocess.TimeoutExpired:
        pytest.fail(f"doctor never returned with a {state} at the TMT "
                    f"{which}: it is blocked on an open, not failing")
    finally:
        _restore(target, state)
    assert proc.stdout.lstrip().startswith("{"), \
        f"{which} as a {state} left stdout unparseable:\n{proc.stderr[-800:]}"
    assert json.loads(proc.stdout) == doc
    assert doc["exit_status"] == proc.returncode
    # ...and the tmt block still said SOMETHING about the tree, rather than
    # silently dropping the rows it could not compute.
    assert [c for c in doc["checks"] if c["section"] == "tmt"], \
        f"the tmt block vanished entirely for a {state} at {which}"


def test_a_fifo_at_a_tmt_annotation_is_named_rather_than_waited_on(tmp_path):
    # The hang itself, and the row that replaced it. `other` is the kind,
    # because os.path.exists() is TRUE for a FIFO and every guard in the
    # engine is built on that - so the honest answer is "present, and not a
    # thing a reader can open", not "absent".
    p = _tmt_project(tmp_path, "fifoann")
    ann = _tmt_target(p, "annotation")
    os.remove(ann)
    os.mkfifo(ann)
    t0 = time.time()
    _proc, doc = doctor_json(p.config_path, expect=1, timeout=60)
    assert time.time() - t0 < 30, "doctor is still waiting on the FIFO"
    c = by_id(doc)["tmt:TMT1:annotation"]
    assert c["status"] == "fail" and c["found"]["kind"] == "other"
    assert c["finding"] == "wrong_kind"
    assert "FIFO" in c["detail"]
    assert c["caveat"] is None, \
        "nothing was read, so the row may not claim it read the file in full"
    assert c["depth"] != "parsed"


def test_regular_readable_refuses_what_a_stat_cheerfully_accepts(ma, tmp_path):
    """The primitive, and the two lies a stat tells about a path.

    os.path.isfile() and os.path.getsize() are STATS, and a stat needs only
    search permission on the parent directory - so a mode-000 regular file
    answers its own size perfectly happily. os.path.exists() is TRUE for a
    FIFO. Both were the basis of a "this is a file I can read" decision.
    """
    good = tmp_path / "good.txt"
    good.write_text("x\n", encoding="utf-8")
    assert ma.regular_readable(str(good))

    mode000 = tmp_path / "locked.txt"
    mode000.write_text("x\n", encoding="utf-8")
    os.chmod(mode000, 0o000)
    try:
        assert os.path.isfile(str(mode000)) and os.path.getsize(str(mode000))
        assert not ma.regular_readable(str(mode000)), \
            "a stat says yes to this and an open does not"
    finally:
        os.chmod(mode000, 0o644)

    assert not ma.regular_readable(str(tmp_path))          # a directory
    assert not ma.regular_readable(str(tmp_path / "nope"))  # absent

    # The FIFO, in a thread, because the failure being tested for is that the
    # call NEVER RETURNS. Asserting the answer without asserting the return
    # would hang the suite on a regression instead of reporting one.
    fifo = tmp_path / "pipe"
    os.mkfifo(str(fifo))
    box = {}
    th = threading.Thread(target=lambda: box.update(r=ma.regular_readable(
        str(fifo))), daemon=True)
    th.start()
    th.join(10)
    assert not th.is_alive(), \
        "regular_readable BLOCKED on a FIFO - O_NONBLOCK is what stops that"
    assert box["r"] is False


def test_found_reaches_unreadable_for_a_regular_file_and_not_only_a_dir(
        ma, tmp_path):
    """symptom: `unreadable` was reachable for a DIRECTORY and never a FILE.

    `_found()` got there only when os.listdir() or os.path.getsize() raised,
    and getsize() SUCCEEDS on a mode-000 regular file - so the state had a
    name in DOCTOR_FOUND_KINDS, a sentence in `_present_kind_phrase()` and a
    `finding`, and the commoner half of it came back `kind: "file"`. Four rows
    branch on `kind == "file"` and said ok about exactly that.
    """
    d = tmp_path / "dir"
    d.mkdir()
    (d / "x").write_text("x", encoding="utf-8")
    f = tmp_path / "file.txt"
    f.write_text("hello\n", encoding="utf-8")
    os.chmod(f, 0o000)
    os.chmod(d, 0o000)
    try:
        assert ma._found(str(f))["kind"] == "unreadable"
        assert ma._found(str(f))["present"] is True
        assert ma._found(str(f))["bytes"] is None, \
            "a size nothing may read is not an answer a consumer should get"
        assert ma._found(str(d))["kind"] == "unreadable"
        assert ma._finding(ma._found(str(f))) == "unreadable"
    finally:
        os.chmod(f, 0o644)
        os.chmod(d, 0o755)
    os.chmod(f, 0o644)
    assert ma._found(str(f))["kind"] == "file"


# The four keys whose rows branch on `kind == "file"`, with the row each one
# writes and the stage `run` really dies in. Every entry was MEASURED as a
# false pass: doctor exited 0 and `run` exited 1 on "Permission denied".
UNREADABLE_FALSE_PASSES = (
    ("proteins_faa", "input:proteins_faa"),
    ("quant_table", "input:quant_table"),
    ("gff", "input:gff"),
    ("emapper_precomputed", "precomputed_emapper:0"),
)


def _unreadable_project(tmp_path, key):
    """A project with `key` present, configured, read by an enabled stage."""
    p = project_with(tmp_path, f"unread_{key[:6]}")
    gff = p.path("input", "x.gff")
    with open(gff, "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
    p.write_config(gff=gff, run=dict(p.cfg["run"], context=(key == "gff")))
    target = (p.cfg["emapper_precomputed"][0]
              if key == "emapper_precomputed" else p.cfg[key])
    return p, target


@pytest.mark.parametrize("key,row", UNREADABLE_FALSE_PASSES)
def test_a_chmod_000_input_fails_doctor_exactly_where_it_fails_run(
        tmp_path, key, row):
    """symptom: doctor exited 0 on a mode-000 input and `run` exited 1 on it.

    Measured, all four, in a differential sweep of doctor against the real
    `run`. Shared FragPipe or eggNOG output on a cluster is the ordinary way
    an input ends up in this state, and the remedy is a permission - which is
    exactly what `found.kind: "unreadable"` exists to say and could not,
    because os.path.getsize() answers for a file nobody may open.

    The claim tested here is the strong one: not merely that doctor now fails,
    but that the stage it NAMES in `blocks` is the stage `run` dies in.
    """
    p, target = _unreadable_project(tmp_path, key)
    os.chmod(target, 0o000)
    try:
        _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
        c = by_id(doc)[row]
        assert c["status"] == "fail", f"{row} still passes a mode-000 input"
        assert c["found"]["kind"] == "unreadable"
        assert c["finding"] == "unreadable"
        assert c["blocks"], f"{row} fails and names nothing that dies"
        assert c["remedy"] == "input", \
            "the remedy for a permission is not a config edit"
        run = p.run(env=no_r_env(), expect=1)
        assert "Permission denied" in run.stderr
        m = re.search(r"stage '([a-z_]+)' failed", run.stderr)
        assert m, f"run did not name a failing stage:\n{run.stderr[-800:]}"
        assert m.group(1) in c["blocks"], \
            f"{row} blocks {c['blocks']} and run died in {m.group(1)}"
    finally:
        os.chmod(target, 0o644)


# What a check that reads INTO a file calls. A function that calls any of
# these AND emits doctor rows has to go through the one gate and carry the
# backstop; nothing else in this file needs to.
DEEP_READERS = ("open", "read_manifest", "read_tmt_annotation",
                "diamond_db_check", "read_csv")


def _called_names(node):
    """Every function name called anywhere inside `node`, bare or attribute."""
    out = set()
    for n in ast.walk(node):
        if isinstance(n, ast.Call):
            fn = n.func
            if isinstance(fn, ast.Name):
                out.add(fn.id)
            elif isinstance(fn, ast.Attribute):
                out.add(fn.attr)
    return out


def test_every_doctor_row_builder_that_opens_a_path_goes_through_the_gate():
    """The structural half, which is the point of the whole exercise.

    Each of the three rounds before this one fixed the site that had just been
    caught. This asserts the rule instead: any function that emits rows (it
    calls `_check()`) and reads into a file (it calls one of DEEP_READERS)
    must go through `_deep_readable()`, which will not hand back a path it has
    not opened, and must carry a wide `except` as the backstop for the state
    the gate did not anticipate. A fifth such site cannot be added without
    either obeying this or failing here.
    """
    tree = ast.parse(open(METAANNOT_PY, encoding="utf-8").read())
    checked, offenders = [], []
    for fn in ast.walk(tree):
        if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        called = _called_names(fn)
        if "_check" not in called or not (called & set(DEEP_READERS)):
            continue
        checked.append(fn.name)
        if "_deep_readable" not in called:
            offenders.append(f"{fn.name} reads into a file without going "
                             "through _deep_readable()")
        wide = set()
        for h in ast.walk(fn):
            if isinstance(h, ast.ExceptHandler) and h.type is not None:
                for n in ast.walk(h.type):
                    if isinstance(n, ast.Name):
                        wide.add(n.id)
        if not ({"OSError", "Exception"} & wide):
            offenders.append(f"{fn.name} has no `except OSError` backstop "
                             "around the read it makes anyway")
    assert checked, "no doctor row builder reads into a file any more - if " \
                    "that is deliberate, DEEP_READERS needs updating"
    assert not offenders, "\n".join(offenders)


def test_the_gate_lets_an_empty_file_through_because_the_reader_says_it_better(
        ma, tmp_path):
    # Emptiness is not an opening failure, and three of the readers behind the
    # gate say something better about it than the gate could: read_manifest
    # refuses "no rows; expected a FragPipe .fp-manifest". A gate that
    # refused an empty file would replace that sentence with its own.
    empty = tmp_path / "empty.txt"
    empty.write_text("", encoding="utf-8")
    f, refusal = ma._deep_readable(str(empty))
    assert refusal is None and f["kind"] == "empty_file"
    fifo = tmp_path / "p"
    os.mkfifo(str(fifo))
    f, refusal = ma._deep_readable(str(fifo))
    assert refusal and "FIFO" in refusal and f["kind"] == "other"


# ----------------------------------------------------------------------
# group T: verdicts that were asserting a tool-dependent outcome
#          unconditionally
# ----------------------------------------------------------------------
def _no_gpu(ma, monkeypatch):
    monkeypatch.setattr(ma, "cuda_probe", lambda: (False, "no CUDA device"))


def _topology_cfg(ma, **over):
    cfg = json.loads(json.dumps(ma.DEFAULT_CONFIG))
    cfg["run"] = dict(cfg["run"], topology=True, structure=False)
    cfg.update(over)
    return cfg


def _gpu_row(ma, monkeypatch, **over):
    _no_gpu(ma, monkeypatch)
    cfg = _topology_cfg(ma, **over)
    rows = {c["id"]: c for c in ma._gpu_checks(cfg, set(ma.enabled_stages(cfg)))}
    return rows["gpu:topology"]


def test_the_gpu_topology_row_is_not_the_same_for_all_three_gpu_settings(
        ma, monkeypatch):
    """symptom: `gpu:topology` named `tmbed_use_gpu` and then gave one answer
    for every value of it.

    Measured warn / blocks [] / degrades ["tmbed"] under auto, true AND false,
    on a host with no CUDA. `stage_tmbed` maps the three values onto three
    different command lines, and the row's own detail said `true` makes "a
    missing GPU fatal instead of slow" - while the row it was on said nothing
    dies. It also told an operator who had already set `false` to set `false`.
    """
    seen = {w: _gpu_row(ma, monkeypatch, tmbed_use_gpu=w)
            for w in ("auto", "true", "false")}
    shapes = {w: (c["status"], tuple(c["blocks"]), tuple(c["degrades"]),
                  c["finding"]) for w, c in seen.items()}
    assert len(set(shapes.values())) == 3, \
        f"the three settings still produce {len(set(shapes.values()))} " \
        f"distinct verdicts: {shapes}"

    assert seen["true"]["status"] == "fail"
    assert seen["true"]["blocks"] == ["tmbed"]
    assert "tmbed_allow_partial" in seen["true"]["config_keys"], \
        "which of the two fatal outcomes it is depends on a second key"

    assert seen["auto"]["status"] == "warn"
    assert seen["auto"]["degrades"] == ["tmbed"]
    assert seen["auto"]["finding"] == "cpu_fallback"

    assert seen["false"]["status"] == "warn"
    assert seen["false"]["degrades"] == ["tmbed"]
    assert "Set tmbed_use_gpu: false" not in seen["false"]["detail"], \
        "this operator has already set it to false"


def test_tmbed_allow_partial_turns_the_fatal_gpu_setting_into_a_degradation(
        ma, monkeypatch):
    # The second key, which the old row never mentioned. With it on, every
    # chunk still fails but stage_tmbed logs the shortfall and returns with an
    # EMPTY prediction file instead of dying - so nothing is blocked and the
    # topology evidence is gone rather than the run.
    c = _gpu_row(ma, monkeypatch, tmbed_use_gpu="true", tmbed_allow_partial=True)
    assert c["status"] == "warn" and c["blocks"] == []
    assert c["degrades"] == ["tmbed"]
    assert "EMPTY" in c["detail"] or "empty" in c["detail"]


def test_an_unrecognised_tmbed_use_gpu_is_fatal_before_any_device_is_probed(
        ma, monkeypatch):
    # stage_tmbed die()s on the value itself, so this is not a GPU verdict and
    # the row does not pretend to be one. The choices it offers are read off
    # TMBED_GPU_MODES, which is the map the stage uses.
    c = _gpu_row(ma, monkeypatch, tmbed_use_gpu="maybe")
    assert c["status"] == "fail" and c["blocks"] == ["tmbed"]
    assert c["finding"] == "config_missing"
    for name in ma.TMBED_GPU_MODES:
        assert name in c["detail"]


def test_the_row_and_the_stage_read_the_same_tmbed_gpu_map(ma):
    # The tripwire that keeps the row honest: the flags the row prints are the
    # flags the stage passes, because there is one map. A fourth mode cannot
    # be added to one and described in the other.
    src = open(METAANNOT_PY, encoding="utf-8").read()
    assert src.count("TMBED_GPU_MODES") >= 3, \
        "the map is no longer read by both the stage and the check"
    assert '"--use-gpu", "--no-cpu-fallback"' not in \
        src[src.index("def stage_tmbed"):src.index("def stage_tmbed") + 3000], \
        "stage_tmbed has its own copy of the flags again"
    assert ma.TMBED_GPU_MODES["true"][1] is True, \
        "tmbed_use_gpu: true is the mode that makes a missing GPU fatal"
    assert ma.TMBED_GPU_MODES["auto"][1] is False
    assert ma.TMBED_GPU_MODES["false"][1] is False


def test_an_empty_contigs_fna_reads_the_same_probe_the_directory_branch_does(
        tmp_path):
    """symptom: the EMPTY branch asserted a tool-dependent outcome flatly.

    "it hands the empty assembly to whichever ORF finder is installed and
    writes an empty candidate list" is two claims welded together, and which
    one happens is the HOST's business: with a finder on PATH the path really
    is handed over, with none installed `stage_smorf` never opens the assembly
    at all. The directory branch one state over was corrected for exactly
    this; the empty branch was left asserting it.
    """
    p = project_with(tmp_path, "fnaempty",
                     contigs_fna=str(tmp_path / "empty.fna"))
    open(tmp_path / "empty.fna", "w", encoding="utf-8").close()
    run_flags(p, smorf=True)

    _proc, doc = doctor_json(p.config_path, expect=0, env=no_orf_finder_env())
    c = by_id(doc)["input:contigs_fna"]
    assert c["status"] == "warn" and c["finding"] == "empty"
    assert c["degrades"] == ["smorf"] and c["blocks"] == []
    assert c["depends_on"] == ["req:smorf"], \
        "the verdict is read off a probe, so it names the row that made it"
    assert "no ORF finder is installed" in c["detail"]
    assert "smorf_proteins.faa" in c["detail"]

    _proc, doc = doctor_json(p.config_path, expect=None,
                             env=_orf_finder_on_path(tmp_path))
    c = by_id(doc)["input:contigs_fna"]
    assert c["depends_on"] == ["req:smorf"]
    assert "hands the empty assembly to macrel" in c["detail"]
    assert "no ORF finder is installed" not in c["detail"]
    # ...and it is still only a degradation: no ORF can come out of nothing,
    # but what the tool does with an empty FASTA is the tool's business and
    # this row does not claim to know it.
    assert c["status"] == "warn" and c["blocks"] == []


# ----------------------------------------------------------------------
# group U: the tail, and the counts
# ----------------------------------------------------------------------
def test_the_rows_that_carry_a_null_found_kind_are_the_ones_named(ma, tmp_path):
    """symptom: "the DIAMOND usability rows" was in that list in three places
    for three rounds, and was never true.

    `db:diamond:<tag>:usable` is built with `found=_found(path)` and carries a
    REAL kind - it is the document's own counterexample to the sentence it was
    named in. The list is data now, and this drives a document that contains
    every family in it and holds the prose against the rows that really carry
    a null.
    """
    import fnmatch
    dmnd = tmp_path / "vfdb.dmnd"
    dmnd.write_bytes(b"x" * 4096)
    p = project_with(tmp_path, "nullrows",
                     db={"diamond": {"vfdb": str(dmnd)}},
                     diamond_weights={"vfdb": 1})
    p.write_config(run=dict(p.cfg["run"], pfam=True, diamond=True,
                            structure=True))
    _proc, doc = doctor_json(p.config_path, expect=None, env=no_r_env())
    globs = [g for _prose, g in ma.DOCTOR_NULL_KIND_ROWS]
    nulls = sorted(c["id"] for c in doc["checks"]
                   if c["found"] is not None and c["found"]["kind"] is None)
    assert nulls, "the fixture no longer produces a null kind at all"
    stray = [i for i in nulls
             if not any(fnmatch.fnmatchcase(i, g) for g in globs)]
    assert stray == [], \
        f"rows carry a null found.kind that DOCTOR_NULL_KIND_ROWS does not " \
        f"name: {stray}"
    # THE OTHER DIRECTION, which is the one the defect was in and which this
    # test did not have: it asserted only that every null row is named, so
    # putting "the DIAMOND usability rows" BACK into the constant - the exact
    # false entry that stood for three rounds - passed it unchanged. A glob
    # may only name rows that really carry a null, and it has to name some.
    lying = sorted(c["id"] for c in doc["checks"]
                   if any(fnmatch.fnmatchcase(c["id"], g) for g in globs)
                   and not (c["found"] is not None
                            and c["found"]["kind"] is None))
    assert lying == [], \
        f"DOCTOR_NULL_KIND_ROWS names rows that carry a REAL found.kind, " \
        f"which is how 'the DIAMOND usability rows' survived three rounds " \
        f"in it: {lying}"
    unmatched = [g for g in globs
                 if not any(fnmatch.fnmatchcase(c["id"], g)
                            for c in doc["checks"])]
    assert unmatched == [], \
        f"DOCTOR_NULL_KIND_ROWS names a family this document has no row " \
        f"for, so the prose describes rows nobody can check: {unmatched}"
    # ...and the DIAMOND usability row is in the document, carrying a real
    # kind, so the correction is not merely a deletion from a list nobody
    # could check.
    usable = by_id(doc)["db:diamond:vfdb:usable"]
    assert usable["found"]["kind"] in ma.DOCTOR_FOUND_KINDS
    assert usable["found"]["kind"] is not None
    assert not any(fnmatch.fnmatchcase(usable["id"], g) for g in globs)


def test_the_null_found_kind_is_counted_from_the_set_it_is_one_past(ma):
    # symptom: the comment on DOCTOR_FOUND_KINDS said "NULL IS AN EIGHTH
    # VALUE" beside a tuple of NINE. That was the third consecutive round to
    # ship a wrong hand-written count, and the second of them inside the
    # comment on the constant that enforces the very set being counted.
    # Nothing writes this number now.
    assert ma._ordinal_word(len(ma.DOCTOR_FOUND_KINDS) + 1) \
        in ma.DOCTOR_NULL_KIND_NOTE
    assert ma._ordinal_word(10) == "tenth"
    src = open(METAANNOT_PY, encoding="utf-8").read()
    i = src.index("DOCTOR_FOUND_KINDS = ")
    comment = " ".join(src[max(0, i - 3000):i].replace("#", " ").split())
    for wrong in ("AN EIGHTH VALUE", "a seventh value", "AN NINTH VALUE"):
        assert wrong not in comment
    assert "DIAMOND" not in ma.DOCTOR_NULL_KIND_NOTE, \
        "the DIAMOND usability rows carry a real kind and are not in this list"


def test_no_orf_finder_env_drops_only_the_directories_that_hold_one(tmp_path,
                                                                    monkeypatch):
    # The helper U8 turns on, driven: it has to remove a directory holding an
    # ORF finder and keep every other one, or the exit status it is used to
    # assert becomes a property of what else is on PATH.
    withorf = tmp_path / "withorf"
    plain = tmp_path / "plain"
    withorf.mkdir()
    plain.mkdir()
    (withorf / "macrel").write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    (plain / "somethingelse").write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setenv("PATH", os.pathsep.join(
        [str(withorf), str(plain)]))
    parts = no_orf_finder_env()["PATH"].split(os.pathsep)
    assert str(withorf) not in parts
    assert str(plain) in parts


def test_a_row_may_not_claim_a_depth_deeper_than_its_found_without_a_caveat(
        ma):
    """symptom: one branch hard-coded `depth="kind"` over a `found` that only
    justified `existence`, and it stayed invisible for a release.

    `_stat_depth()` derives the depth a `found` justifies, and the previous
    round's claim was that "no row writes its own". One did, and the reason
    nobody noticed is that the branch could only be reached by a DIRECTORY,
    where the derived answer and the written one agree. They part company the
    moment `unreadable` is reachable for a file. The ratchet is enforced at
    the point the row is made now, in both directions: `config` is refused
    when something WAS consulted, and a depth past what `found` justifies is
    refused unless the row carries a `caveat` saying what it read.
    """
    unreadable = {"present": True, "kind": "unreadable", "bytes": None,
                  "entries": None, "symlink": False}
    assert ma._stat_depth(unreadable) == "existence"
    row = dict(target="/p", found=unreadable, blocks=["emapper"],
               remedy="input")
    with pytest.raises(AssertionError) as e:
        ma._check("x:y", "inputs", "x", "fail", "unreadable", "d",
                  depth="kind", **row)
    assert "justifies no more than 'existence'" in str(e.value)
    # ...and the same row is fine at the depth its `found` justifies, or
    # deeper with a caveat that says what was opened.
    ma._check("x:y", "inputs", "x", "fail", "unreadable", "d",
              depth="existence", **row)
    ma._check("x:y", "inputs", "x", "fail", "unreadable", "d",
              depth="parsed", caveat="the file was read in full", **row)


def test_every_row_in_a_swept_document_derives_its_depth_from_its_found(
        ma, tmp_path):
    # The rule above, over real documents rather than a synthetic row: every
    # configured input in every state it can be in, with the invariant applied
    # to each row that comes back.
    for state in PATH_STATES:
        p = project_with(tmp_path / "depths", f"d_{state}")
        target = p.cfg["quant_table"]
        _put_in_state(target, state)
        try:
            _proc, doc = doctor_json(p.config_path, expect=None, timeout=90)
        finally:
            _restore(target, state)
        for c in doc["checks"]:
            f = c["found"]
            if f is None or f["kind"] is None:
                continue
            floor = ma._stat_depth(f)
            if c["depth"] not in ma.DOCTOR_DEPTH_ORDER:
                continue
            if ma.DOCTOR_DEPTH_ORDER.index(c["depth"]) > \
                    ma.DOCTOR_DEPTH_ORDER.index(floor):
                assert c["caveat"], \
                    f"{state}/{c['id']}: depth {c['depth']} over a " \
                    f"found.kind of {f['kind']} with nothing said about what " \
                    "it read"


# ----------------------------------------------------------------------
# the OUTCOME VERB, which the sweep does not look at
# ----------------------------------------------------------------------
# The V sweep compares counts, and the differential sweep compares VERDICTS.
# Neither has anything to say about what a row's sentence claims WILL HAPPEN,
# which is why three rows went four verification rounds saying a stage "dies"
# on a path where `run` hangs. These tests are about the verb.


def _fifo_thread(fn, seconds=4):
    """Run `fn` in a daemon thread and say whether it ever came back.

    The only honest way to test "this call does not return yet" - asserting a
    return value instead would hang the SUITE on a regression rather than
    reporting one, which is the same failure mode the rows are about. It is
    still needed with a bounded wait in place, because "waits" and "refuses at
    once" are exactly what the two halves of the FIFO contract are.
    """
    box = {}

    def run():
        # The exception is CAUGHT rather than allowed out of the thread. A
        # call left waiting by this helper goes on waiting after the test
        # returns and then raises into a daemon thread, and an uncaught one
        # there prints a traceback onto the suite's stderr from nowhere in
        # particular - which reads as a broken run rather than as the test
        # doing exactly what it meant to.
        try:
            box["r"] = fn()
        except Exception as e:                              # noqa: BLE001
            box["raised"] = e

    th = threading.Thread(target=run, daemon=True)
    th.start()
    th.join(seconds)
    return th.is_alive(), box


def test_header_columns_on_a_fifo_waits_and_then_dies(ma, tmp_path,
                                                     monkeypatch):
    """The claim the `tmt:<plex>:level_file` row makes, driven.

    Both halves, because the row now claims both: the reader WAITS on a FIFO -
    it does not refuse one, which is what would break `mkfifo p; zcat ... > p
    &` - and the wait ENDS, in a StageError naming the path, which is what
    stops it being the hang `run` was measured in.

    BOTH HALVES GO THROUGH _fifo_thread(), and the second one did not. It
    called header_columns() on the MAIN thread inside pytest.raises, which is
    the assertion that cannot fail: against an opener whose wait is not
    bounded that call does not raise, it BLOCKS, and the test wedges the whole
    suite until somebody kills it - the very failure mode this row is about,
    reproduced in the test written to rule it out. Found by reverting the
    bound (dropping O_NONBLOCK from _OPEN_FLAGS) and running the suite under a
    watchdog, which is the only way to find this class: a wedging test cannot
    be found by reading, because it looks exactly like a passing one.
    """
    fifo = tmp_path / "ion.tsv"
    os.mkfifo(str(fifo))
    monkeypatch.setattr(ma, "_FIFO_WAIT", 12.0)
    alive, box = _fifo_thread(lambda: ma.header_columns(str(fifo)))
    assert alive and not box, \
        "header_columns() came back on an unwritten FIFO before its wait was " \
        "over - a reader that refuses a FIFO on sight breaks the live-writer " \
        "workflow opener() exists to keep"
    # ...and the wait is bounded, which is the half that was missing.
    monkeypatch.setattr(ma, "_FIFO_WAIT", 1.0)
    t0 = time.time()
    alive, box = _fifo_thread(lambda: ma.header_columns(str(fifo)), 20)
    assert not alive, \
        "header_columns() is still on the open after its wait - the bound is " \
        "gone and `doctor` can hang again"
    assert time.time() - t0 < 20, "the wait is not bounded by _FIFO_WAIT"
    e = box.get("raised")
    assert isinstance(e, ma.StageError), \
        f"the wait ended in something other than a refusal: {box}"
    assert str(fifo) in str(e) and "FIFO" in str(e) \
        and "fifo_wait_s" in str(e), \
        f"the refusal names neither the path, the kind nor the setting: {e}"
    # ...and the same call on a real file comes straight back, so the test
    # above is about the FIFO and not about the helper being broken.
    real = tmp_path / "real.tsv"
    real.write_text("a\tb\n1\t2\n", encoding="utf-8")
    assert ma.header_columns(str(real)) == ["a", "b"]


def test_prepare_emapper_on_a_fifo_waits_and_then_dies(ma, tmp_path,
                                                       monkeypatch):
    """The claim the `precomputed_emapper:<i>` row makes, driven.

    prepare_emapper() is what reads that path, and a `.annotations` table is
    the one input most likely to arrive as a `.gz` - so this is also the row
    whose reader goes through opener()'s GZIP branch, which had to be made to
    read the descriptor rather than re-open the path.
    """
    ps = [F.Protein(f"P{i}", "MKV" * 20) for i in range(3)]
    faa = F.write_fasta(str(tmp_path / "p.faa"), ps)
    fifo = tmp_path / "cat.annotations"
    os.mkfifo(str(fifo))
    run = (lambda: ma.prepare_emapper(
        [str(fifo)], faa, str(tmp_path / "o.annotations"), "", "exact",
        0.5, 0.9, 100))
    monkeypatch.setattr(ma, "_FIFO_WAIT", 12.0)
    alive, box = _fifo_thread(run)
    assert alive and not box, \
        "prepare_emapper() came back on an unwritten FIFO before its wait " \
        "was over"
    monkeypatch.setattr(ma, "_FIFO_WAIT", 1.0)
    # THROUGH _fifo_thread() like the half above it, and for the same reason:
    # a bare run() here does not fail when the bound goes away, it blocks, and
    # a blocking test wedges the suite instead of reporting the regression.
    alive, box = _fifo_thread(run, 20)
    assert not alive, \
        "prepare_emapper() is still on the open after its wait - the bound " \
        "is gone and the run can hang again"
    e = box.get("raised")
    assert isinstance(e, ma.StageError), \
        f"the wait ended in something other than a refusal: {box}"
    # The message says WHICH silence this was, and here nothing was ever
    # attached: "no writer appeared" asserted that for both states, and a
    # writer that is attached but slow is the state where it is false.
    assert str(fifo) in str(e) \
        and "nothing is holding the write end" in str(e), \
        f"the refusal does not say what the silence was: {e}"


# Every path `doctor` is pointed at that a FIFO can be put at, with the row it
# writes. This is the list _fifo_clause() is measured against, and it is
# spelled out because "every row a FIFO can reach" was a claim a report made
# and nothing checked.
# The third column is how many times the project below READS that input, which
# is what decides whether the row may promise a wait at all. It is not written
# here: it is recomputed from the config by ma.set_read_plan() in the test, so
# a row and the engine cannot part company about it - which is the whole
# defect this parametrisation used to carry, one surface over.
FIFO_ROWS = (("proteins_faa", "input:proteins_faa"),
             ("quant_table", "input:quant_table"),
             ("manifest", "input:manifest"),
             ("gff", "input:gff"),
             ("emapper_precomputed", "precomputed_emapper:0"))


def _fifo_project(tmp_path, key):
    """A project with an enabled reader for `key`, and `key`'s real path."""
    p = project_with(tmp_path, f"fifo_{key[:7]}")
    gff = p.path("input", "x.gff")
    with open(gff, "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
    p.write_config(gff=gff, run=dict(p.cfg["run"], context=(key == "gff")))
    target = (p.cfg["emapper_precomputed"][0]
              if key == "emapper_precomputed" else p.cfg[key])
    return p, target


@pytest.mark.parametrize("key,row", FIFO_ROWS)
def test_a_fifo_row_says_what_that_input_really_costs(tmp_path, ma, key, row):
    """The verb, on every row a FIFO can reach, held against what `run` does.

    This assertion has been reversed twice and the second reversal is the one
    that matters. `run` was first measured HANGING on each of these, so a row
    saying the stage "dies" sent an operator to look for a failure that was
    never coming. opener() was then given a bounded wait, and the rows were
    rewritten to promise it - which was right for one of these paths and wrong
    for the rest, because a verifier drove the workflow the sentence
    recommended and found that a pipe cannot work at an input the run OPENS
    MORE THAN ONCE, whatever anything waits.

    So there is one test and two expectations, and which one applies is
    recomputed from the config rather than listed: a single-read input must
    promise the wait, because the wait is what keeps `mkfifo p; zcat
    big.faa.gz > p &` working there, and a multi-read input must refuse at
    once and say how many reads there are. A row that hard-coded either would
    pass on one path and lie on the others.
    """
    p, target = _fifo_project(tmp_path, key)
    _clear(target)
    os.mkfifo(target)
    reads = len(ma.set_read_plan(p.cfg).get(os.path.realpath(target), []))
    t0 = time.time()
    _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
    assert time.time() - t0 < 40, "doctor is still waiting on the FIFO"
    c = by_id(doc)[row]
    assert c["found"]["kind"] == "other", \
        f"{row} does not see the FIFO as `other`"
    # The bucket, split. `other` is a FIFO, a socket AND a device node, and
    # every sentence below used to be printed for all three.
    assert c["found"]["other_kind"] == "fifo", \
        f"{row} does not say which kind of `other` this is"
    assert c["status"] == "fail"
    assert "FIFO" in c["detail"], f"{row} does not name the FIFO"
    assert "dies" in c["detail"], \
        f"{row} does not say the reader dies: {c['detail']}"
    assert "never returns" not in c["detail"], \
        f"{row} still claims a hang that opener() no longer has: {c['detail']}"
    if reads > 1:
        assert "refused AT ONCE" in c["detail"], \
            f"{row} promises a wait at an input this run reads {reads} " \
            f"times, where waiting cannot help: {c['detail']}"
        assert f"{reads} times" in c["detail"], \
            f"{row} does not say how many reads there are: {c['detail']}"
        assert "not at once" not in c["detail"], \
            f"{row} says the death is delayed where it is immediate"
        assert "live writer" not in c["detail"], \
            f"{row} recommends a workflow that cannot work at this path"
    else:
        # ...that the death is not immediate, which is the thing an operator
        # watching a run has to plan around...
        assert "not at once" in c["detail"], \
            f"{row} promises a death at the open, and the open waits first"
        # ...and that a FIFO with a writer is still read, which is the
        # capability the wait exists to protect.
        assert "fifo_wait_s" in c["detail"] and "live writer" in c["detail"], \
            f"{row} carries no _fifo_clause(): {c['detail']}"


def test_the_manifest_row_refuses_a_pipe_on_a_config_that_reads_it_twice(
        tmp_path, ma):
    """THE CORRECTION THIS ROUND WAS CALLED FOR, guarded where it is published.

    `doctor --json` published "A FIFO here is NOT refused on sight, because
    this run reads this input exactly once" for a manifest that a taxonomy run
    opens once per peptide stage AND again in stage_join - read_feature_table
    calls read_manifest, and peptide_features() calls read_feature_table - so
    an operator was told a `mkfifo` workflow was supported at a path where the
    run then died at the second open, hours in. The fix was to compute the
    manifest's sites from _manifest_readers() instead of the hand-written row
    that named stage_join alone.

    Nothing in this file held that. Measured: with _manifest_sites() put back
    the way it was, every other test here stayed green while the document
    published the false sentence again - the only thing that fell was a plan
    comparison in tests/test_stages.py, which is a different file and caught
    it as a side effect. A promise is published HERE and has to be pinned
    here.

    THE COUNT IS NOT TAKEN FROM THE PLAN, which is the trap this test exists
    to avoid: recomputing it with set_read_plan() - which is what the
    parametrised FIFO test above does, for a different purpose - would move
    the expectation and the row together, and the reverted version would pass.
    It comes from _manifest_readers(), the predicate that answers who opens
    this file at all, and the sweep in tests/test_stages.py is what holds THAT
    against a driven run.
    """
    p = project_with(tmp_path, "manytimes")
    p.write_config(
        run=dict(p.cfg["run"], unipept=True, taxonomy=True, join=True),
        db=dict(p.cfg.get("db") or {},
                ncbi_taxonomy=F.toy_taxdump(str(tmp_path / "taxdump"))),
        # allow_http so the unipept cache row is a warning: the verdict under
        # test is the manifest's.
        unipept=dict(p.cfg.get("unipept") or {}, allow_http=True))
    cfg = ma.load_config(p.config_path)
    readers = ma._manifest_readers(cfg)
    assert len(readers) > 1, (
        "the premise is gone: this config is meant to open the manifest once "
        f"per stage in {readers}, and a single reader cannot show the defect")
    target = p.cfg["manifest"]
    os.remove(target)
    os.mkfifo(target)
    _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
    c = by_id(doc)["input:manifest"]
    assert c["found"]["other_kind"] == "fifo"
    # The published sentence, word for word, and not the substring "exactly
    # once" - the refusal itself says a pipe is "drained EXACTLY ONCE", which
    # is the opposite claim and shares four of its words.
    assert "reads this input exactly once" not in c["detail"], (
        "the document still promises a single-read manifest on a config that "
        f"opens it in each of {readers}:\n{c['detail']}")
    assert "refused AT ONCE" in c["detail"], \
        f"the row does not refuse the pipe it cannot support:\n{c['detail']}"
    assert f"{len(readers)} times" in c["detail"], (
        f"the row does not carry the real count ({len(readers)}), which is "
        f"what sends an operator to a real file:\n{c['detail']}")
    # ...and it names the reads, so the count can be checked rather than
    # believed. The peptide stages reach read_manifest THROUGH
    # peptide_features, which is the chain the old row had no entry for.
    for st in readers:
        assert f"stage_{st}" in c["detail"], \
            f"the row does not name the read {st} makes:\n{c['detail']}"
    assert "peptide_features" in c["detail"], \
        f"the row stops at the stage instead of naming the reader that " \
        f"opens the file:\n{c['detail']}"


def test_a_socket_row_does_not_publish_the_fifo_paragraph(tmp_path):
    """C2, driven: the derived verb was false for most of what `other` covers.

    `_raises_promptly()` was `kind != "other"`, and `other` is one bucket for a
    FIFO, a UNIX socket and a device node - so a socket at proteins_faa
    published "dies, but not at once ... A FIFO in particular is not refused
    on sight ... waits fifo_wait_s", while a socket cannot be opened as a
    file at all and os.open() itself fails on it, which is what the
    stage-side test asserts. The two
    halves of one change set contradicted each other in the emitted document,
    and PATH_STATES had no socket state to catch it with.
    """
    p, target = _fifo_project(tmp_path, "proteins_faa")
    _clear(target)
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    try:
        cwd = os.getcwd()
        os.chdir(os.path.dirname(target))
        try:
            sock.bind(os.path.basename(target))
        finally:
            os.chdir(cwd)
        _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
    finally:
        sock.close()
    c = by_id(doc)["input:proteins_faa"]
    assert c["found"]["kind"] == "other"
    assert c["found"]["other_kind"] == "socket", \
        "a socket is not distinguished from a FIFO in `found`"
    assert "socket" in c["detail"], \
        f"the row does not say a socket is what is there: {c['detail']}"
    assert "not at once" not in c["detail"], \
        f"the row claims a wait on a socket, which fails at os.open(): " \
        f"{c['detail']}"
    assert "fifo_wait_s" not in c["detail"] and "mkfifo" not in c["detail"], \
        f"the row offers a FIFO's paragraph for a socket: {c['detail']}"


def test_a_device_node_row_does_not_publish_the_fifo_paragraph(tmp_path):
    """The third of `other`'s kinds, for the same reason as the second.

    /dev/zero answers immediately and forever with bytes that are not a table,
    so nothing waits and nothing is coming that is not already here. A row
    that told an operator to start a writer, or to raise `fifo_wait_s`, would
    be describing a different object entirely.
    """
    p = project_with(tmp_path, "devnode")
    p.write_config(proteins_faa="/dev/zero")
    _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
    c = by_id(doc)["input:proteins_faa"]
    assert c["found"]["kind"] == "other"
    assert c["found"]["other_kind"] in ("char_device", "block_device"), \
        f"a device node is not distinguished: {c['found']}"
    assert "device" in c["detail"], \
        f"the row does not say a device node is what is there: {c['detail']}"
    assert "not at once" not in c["detail"] \
        and "fifo_wait_s" not in c["detail"], \
        f"the row offers a FIFO's paragraph for a device node: {c['detail']}"


def test_a_device_node_at_the_manifest_does_not_claim_an_oserror(tmp_path):
    """The row that names an EXCEPTION, over the state it named it wrong for.

    The sibling above drives `proteins_faa`, whose row names no exception type
    at all - so the new PATH_STATES entry for a device node never reached the
    row where the claim lives. `input:manifest` is that row: it is the one
    that says which exception the open ends in and therefore which stages go
    with it, and on a character device it published "it dies with an OSError
    ... every stage that opens a quant table on this format goes with it"
    while its own `blocks` field said `['join']`.

    `blocks` was the right half, and this drives the engine to prove it: a
    char-device manifest ends in a StageError, peptide_features() catches a
    StageError and re-reads with the peptide-only reader, so the taxonomy
    stage survives what the sentence said would take it.
    """
    p = project_with(tmp_path, "devman")
    dev = p.path("input", "dev-manifest")
    _clear(dev)
    os.symlink("/dev/zero", dev)
    p.write_config(manifest=dev,
                   run=dict(p.cfg["run"], unipept=True, taxonomy=True),
                   unipept={"result": "", "allow_http": False})
    _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
    c = by_id(doc)["input:manifest"]
    assert c["found"]["other_kind"] in ("char_device", "block_device"), \
        f"the state did not come out as a device node: {c['found']}"
    assert "OSError" not in c["detail"], (
        "the row claims an OSError for a state that raises a StageError, "
        f"which peptide_features() catches: {c['detail']}")
    assert "StageError" in c["detail"], \
        f"the row does not say which exception this is: {c['detail']}"
    assert "once the wait is over" not in c["detail"], \
        f"a device node does not wait for anything: {c['detail']}"
    assert c["blocks"] == ["join"], (
        "the sentence and `blocks` are two readings of one fact and have to "
        f"agree: blocks={c['blocks']}, detail={c['detail']}")
    # ...and the machine-readable half is what the engine really does.
    assert "taxonomy" not in c["blocks"], \
        "taxonomy survives a manifest the full reader refuses"


def test_the_exception_predicate_matches_what_the_open_really_raises(
        ma, tmp_path):
    """_dies_with_a_stage_error(), held against _open_for_read() itself.

    THE PREDICATE IS A CLAIM ABOUT PATH STATES, and a claim about path states
    that is written by hand is a claim that is true of the state its author
    had in mind - which is how `empty or kind == "other"` came to treat a
    socket, a block device and a character device as one thing when the first
    two end in an OSError and the third does not. So every state the sweep can
    create is driven through the real opener and the predicate is held against
    the exception that comes out.

    In a thread with a join, because one of these states is a FIFO and a FIFO
    with no writer is what the whole file is about: a test that called the
    opener directly here would not fail when the bound went away, it would
    wedge.
    """
    ma.set_fifo_wait(1.0)
    made = []
    try:
        wrong = []
        for state in PATH_STATES:
            path = str(tmp_path / f"s_{state}")
            _put_in_state(path, state)
            made.append(path)
            f = ma._found(path)
            if not f["present"]:
                continue          # absent and dangling_symlink raise at the open
            box = {}

            def go(path=path):
                try:
                    fd, _pb = ma._open_for_read(path)
                    os.close(fd)
                    box["opened"] = True
                except BaseException as e:                  # noqa: BLE001
                    box["error"] = e

            th = threading.Thread(target=go, daemon=True)
            th.start()
            th.join(20)
            assert not th.is_alive(), \
                f"_open_for_read blocked on a {state} - the wait is not bounded"
            err = box.get("error")
            # An empty regular file OPENS; what refuses it is the reader above,
            # which is a StageError, and the predicate says so for that reason.
            real = (True if state == "empty_file"
                    else isinstance(err, ma.StageError))
            said = ma._dies_with_a_stage_error(f, path)
            if said != real:
                wrong.append(
                    f"{state}: the predicate says "
                    f"{'StageError' if said else 'OSError'}, the open raised "
                    f"{type(err).__name__ if err else 'nothing'}")
        assert not wrong, (
            "_dies_with_a_stage_error() disagrees with the opener, so every "
            "row that names an exception and every `blocks` computed from it "
            "is wrong for these states:\n  " + "\n  ".join(wrong))
    finally:
        ma.set_fifo_wait(ma.DEFAULT_CONFIG["fifo_wait_s"])
        # Cleared by hand: two of these states are mode-000, and pytest's own
        # tmp_path cleanup cannot remove what it may not enter - it warns
        # instead, which reads as a broken suite rather than as a test that
        # made a directory nobody may list.
        for path in made:
            _clear(path)


@pytest.mark.parametrize("key,row", FIFO_ROWS)
def test_the_same_row_says_the_death_is_immediate_where_it_really_is(
        tmp_path, key, row):
    """The other half, and what stops the fix being 'always say it waits'.

    A DIRECTORY at the same path raises at once - measured, one reader at a
    time - so the timing has to come back different. A row that hard-coded
    the wait would pass the test above and fail this one.
    """
    p, target = _fifo_project(tmp_path, key)
    _clear(target)
    os.makedirs(target)
    _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
    c = by_id(doc)[row]
    assert c["found"]["kind"] in ("dir", "empty_dir")
    assert "not at once" not in c["detail"], \
        f"{row} claims a wait on a directory, which raises at once"
    assert "fifo_wait_s" not in c["detail"], \
        f"{row} offers a FIFO's remedy for a directory"


def test_the_two_tmt_rows_for_one_plex_agree_about_the_verb(tmp_path):
    """symptom: the correction reached one of two rows in the same function.

    `tmt:<plex>:annotation` said "or, on a FIFO with no writer, never returns
    at all" and `tmt:<plex>:level_file`, twenty lines above it in the same
    loop, said the reader "cannot read a table out of it".
    """
    for which, row in (("level_file", "tmt:TMT1:level_file"),
                       ("annotation", "tmt:TMT1:annotation")):
        p = _tmt_project(tmp_path, f"verb_{which[:5]}")
        target = _tmt_target(p, which)
        os.remove(target)
        os.mkfifo(target)
        t0 = time.time()
        _proc, doc = doctor_json(p.config_path, expect=1, timeout=90)
        assert time.time() - t0 < 40, "doctor is still waiting on the FIFO"
        c = by_id(doc)[row]
        assert c["found"]["kind"] == "other"
        assert "not at once" in c["detail"], \
            f"{row} does not say the reader waits before it dies"
        assert "never returns" not in c["detail"], \
            f"{row} claims a hang opener() no longer has"


# ----------------------------------------------------------------------
# why the fix is one open and not a gate, pinned as a measurement
# ----------------------------------------------------------------------
def test_probing_a_live_fifo_breaks_the_writer_on_the_other_side(ma, tmp_path):
    """The measurement that decided the SHAPE of the fix.

    regular_readable() is not a read-only observation of a FIFO: it opens the
    read end and closes it again, and the producer gets EPIPE. That is why
    opener() opens ONCE and keeps the descriptor rather than probing the path
    and letting the reader open it again - a gate in front of a stage's open
    destroys the stream it was meant to protect, and probe-close-reopen is a
    gate however it is spelled.
    """
    fifo = tmp_path / "probed.faa"
    os.mkfifo(str(fifo))
    err = {}

    def write():
        try:
            with open(str(fifo), "w", encoding="utf-8") as fh:
                fh.write(">P1\n" + "MKV" * 20000 + "\n")
                fh.flush()
        except BrokenPipeError as e:                        # noqa: PERF203
            err["e"] = e

    th = threading.Thread(target=write, daemon=True)
    th.start()
    time.sleep(0.5)
    assert ma.regular_readable(str(fifo)) is False, \
        "regular_readable accepted a FIFO, which would make the gate above " \
        "look harmless"
    th.join(10)
    assert not th.is_alive() and "e" in err, \
        "the probe left the writer alive - if that is now true, the second " \
        "of the two measured facts in regular_readable()'s docstring no " \
        "longer holds and opener() could have been written as a gate"


# ----------------------------------------------------------------------
# the closed vocabularies, ALL of them, at the place the rows are made
# ----------------------------------------------------------------------
# Two of these were guarded in _check() and five were not, and the five were
# enforced nowhere except by the eleven CONFIGS the tests walk - so a typo in
# a branch those configs do not reach shipped in silence. `status` is the one
# that costs something immediately: exit_status is `1 if fails else 0` counted
# over status == "fail", so a row that MEANS to fail and writes "faill" is a
# silent zero.
#
# A table, and not seven hand-written tests, because the defect was a
# vocabulary that nobody thought to guard: a new closed set added to
# DOCTOR_* without a row here is the same omission again, and
# test_every_closed_vocabulary_is_in_this_table is what notices.
NEAR_MISSES = (
    ("status", dict(status="faill"), "status 'faill'"),
    ("remedy", dict(kw={"remedy": "confgi"}), "remedy 'confgi'"),
    ("depth", dict(kw={"depth": "kindd"}), "depth 'kindd'"),
    ("fails_reason",
     dict(kw={"fails_reason": "setting_ignorred", "blocks": ["join"]}),
     "fails_reason 'setting_ignorred'"),
    ("blocks", dict(kw={"blocks": ["joinn"]}), "not a stage in STAGE_NAMES"),
    ("blocks_commands", dict(kw={"blocks_commands": ["repport"]}),
     "is not one of ['run', 'report', 'object']"),
    ("degrades", dict(kw={"degrades": ["contextt"]}),
     "not a stage in STAGE_NAMES"),
    ("found.kind", dict(kw={"found": {"kind": "flie", "present": True}}),
     "not in DOCTOR_FOUND_KINDS"),
    ("expect.kind", dict(kw={"expect": {"kind": "fille"}}),
     "not in DOCTOR_EXPECT_KINDS"),
    # The EIGHTH, and the one that had no guard of any kind. `section` decides
    # the heading a row prints under; print_doctor() looks it up in
    # dict(DOCTOR_SECTIONS) with no `.get`, so a typo was a bare KeyError and
    # NO DOCUMENT from the text command, while --json dropped the row from
    # `sections` in silence and left it in `checks`. It escaped the tripwire
    # below for a reason worth keeping in mind: DOCTOR_SECTIONS was a LIST,
    # and the tripwire collects tuples.
    ("section", dict(section="inputz"), "section 'inputz' is not one of"),
    # The two lists against EACH OTHER, which is the mistake an author really
    # makes: not inventing a name, but putting a real one in the wrong list.
    ("blocks takes a command", dict(kw={"blocks": ["run"]}),
     "which is a metaannot COMMAND and not a stage"),
    ("blocks_commands takes a stage", dict(kw={"blocks_commands": ["join"]}),
     "which is a STAGE and not a command"),
    ("degrades takes a command", dict(kw={"degrades": ["report"]}),
     "which is a metaannot COMMAND and not a stage"),
    # The NINTH, and it arrived because `other` was three things wearing one
    # name: a FIFO, a socket and a device node all answered `kind: "other"`,
    # so every sentence derived from that kind was written for the FIFO and
    # printed over the other two. `found.other_kind` is the field that tells
    # them apart, and it is guarded BOTH ways - a value outside the
    # vocabulary, and a value on a kind that may not carry one, because
    # putting "fifo" on a directory would be the same bucketing defect one
    # level down in the field added to remove it.
    ("found.other_kind",
     dict(kw={"found": {"kind": "other", "other_kind": "pipe",
                        "present": True}}),
     "not in DOCTOR_OTHER_KINDS"),
    ("other_kind on the wrong kind",
     dict(kw={"found": {"kind": "dir", "other_kind": "fifo",
                        "present": True}}),
     "only `other` carries one"),
)


@pytest.mark.parametrize("what,spec,message",
                         NEAR_MISSES, ids=[n for n, _s, _m in NEAR_MISSES])
def test_every_closed_vocabulary_refuses_a_near_miss(ma, what, spec, message):
    """symptom: five of the seven closed vocabularies were guarded nowhere.

    _check()'s own comment gave the reason the two guarded ones were guarded -
    "a typo here should be a loud failure and not a value a consumer's switch
    falls through on" - and that reason never applied to only two of them.
    """
    with pytest.raises(AssertionError) as e:
        ma._check("x:y", spec.get("section", "inputs"), "label",
                  spec.get("status", "warn"),
                  "finding", "detail", **spec.get("kw", {}))
    assert message in str(e.value), \
        f"{what} is refused, but not with a sentence naming what is wrong: " \
        f"{e.value}"


def test_a_row_that_misspells_its_status_cannot_reach_the_exit_status(ma):
    """The consequence, spelled out, because it is the expensive one.

    exit_status is `1 if fails else 0` over `status == "fail"`. A row that
    means to fail and writes "faill" is counted as a pass, the command exits
    0, and nothing anywhere objects - which is the whole reason `status` is
    the vocabulary that had to be guarded first.
    """
    ok = ma._check("x:y", "inputs", "label", "fail", "missing", "detail",
                   blocks=["join"], remedy="input")
    assert ok["status"] == "fail" and ok["fails_reason"] \
        == "stage_or_command_dies"
    with pytest.raises(AssertionError):
        ma._check("x:y", "inputs", "label", "faill", "missing", "detail",
                  blocks=["join"], remedy="input")


def test_a_misspelled_section_would_have_taken_the_text_command_down(ma):
    """The consequence, spelled out, because it is the expensive one.

    print_doctor() groups rows by `section` and prints
    `dict(DOCTOR_SECTIONS)[section]` - no `.get`, no default. A row with a
    typo there raised a bare KeyError before any heading was printed, so the
    TEXT command produced no document at all: not a wrong row, no rows, no
    tools block, no databases block, nothing. Under --json the same row stayed
    in `checks` and vanished from `sections`, which is the quieter half of the
    same defect: a consumer rendering section by section silently drops it.

    Driven through print_doctor() rather than asserted about it, because the
    KeyError was in the RENDERER and the guard is in the builder, and a test
    that only called _check() would not show what the guard buys.
    """
    for sec in [i for i, _t in ma.DOCTOR_SECTIONS]:
        row = ma._check(f"x:{sec}", sec, "label", "ok", "fine", "detail")
        assert row["section"] == sec
    with pytest.raises(AssertionError) as e:
        ma._check("x:y", "inputes", "label", "ok", "fine", "detail")
    assert "is not one of" in str(e.value) and "inputes" in str(e.value)
    # ...and the renderer is what it would have cost: a hand-made row with a
    # section nothing names still takes print_doctor() down, which is why the
    # guard is at the place rows are MADE.
    with pytest.raises(KeyError):
        ma.print_doctor([dict(ma._check("x:y", "inputs", "l", "ok", "f", "d"),
                              section="inputes")], [])


def test_a_row_whose_found_or_expect_has_no_kind_is_refused_by_name(ma):
    """symptom: `found` without a `kind` passed, and `expect` without one
    raised a bare KeyError.

    The two halves of one omission. `found` was read through `.get("kind")`,
    so a dict built by hand without one answered None - which is a LEGAL value
    of found.kind, meaning "this build did not compute it" - and the row
    shipped a `found` that every consumer switches on and that says nothing.
    `expect` was indexed directly, so the same mistake was a KeyError raised
    from inside the function whose entire job is to make a malformed row
    impossible, and it cost the whole document rather than one row.
    """
    for field, maker in (("found", "_found()"), ("expect", "_expect()")):
        with pytest.raises(AssertionError) as e:
            ma._check("x:y", "inputs", "label", "warn", "finding", "detail",
                      **{field: {"present": True, "because": "no kind here"}})
        assert f"{field} is" in str(e.value) and "no 'kind'" in str(e.value) \
            and maker in str(e.value), \
            f"{field} without a kind is refused, but not by name: {e.value}"
    # The null that IS legal stays legal: DOCTOR_NULL_KIND_ROWS is a documented
    # family, and a guard that refused it would be a different defect.
    ok = ma._check("x:y", "tools", "label", "ok", "found", "detail",
                   depth="existence",
                   found={"kind": None, "present": True, "bytes": None})
    assert ok["found"]["kind"] is None


def test_every_closed_vocabulary_is_in_the_near_miss_table(ma):
    """The tripwire for the defect itself, which was an OMISSION.

    Every DOCTOR_* tuple that _check() writes into a row is a closed set a
    consumer switches on, and the way five of them went unguarded is that
    nobody enumerated them. This enumerates them: a new one has to arrive with
    a near-miss row above, or this fails and names it.
    """
    guarded = {n.split()[0] for n, _s, _m in NEAR_MISSES}
    # The constants, and the row field each one is the vocabulary OF.
    vocabularies = {
        "DOCTOR_STATUSES": "status",
        "DOCTOR_REMEDIES": "remedy",
        "DOCTOR_DEPTHS": "depth",
        "DOCTOR_FAIL_REASONS": "fails_reason",
        "DOCTOR_FOUND_KINDS": "found.kind",
        "DOCTOR_EXPECT_KINDS": "expect.kind",
        "DOCTOR_COMMANDS": "blocks_commands",
        # The eighth. It is (id, title) pairs rather than bare ids because the
        # printed heading is named beside the id, and the vocabulary `section`
        # is checked against is the ids of exactly these pairs - so this is the
        # constant, and a ninth section arriving here with no guard is what
        # this test now notices. It took a one-character change to be seen at
        # all: a list is not a tuple, and the sweep below collects tuples.
        "DOCTOR_SECTIONS": "section",
        # The ninth: what `other` really is. See the near-miss rows above.
        "DOCTOR_OTHER_KINDS": "found.other_kind",
    }
    live = {n for n in dir(ma)
            if n.startswith("DOCTOR_") and isinstance(getattr(ma, n), tuple)
            and n not in ("DOCTOR_DEEP_CHECKS", "DOCTOR_DEPTH_ORDER",
                          "DOCTOR_NULL_KIND_ROWS")}
    assert live == set(vocabularies), \
        f"a closed DOCTOR_* vocabulary arrived or left without this table " \
        f"moving: {sorted(live ^ set(vocabularies))}"
    missing = sorted(f for f in vocabularies.values() if f not in guarded)
    assert not missing, \
        f"these vocabularies have no near-miss row and so no guard: {missing}"
    # ...and the two stage lists, which are not a DOCTOR_* tuple at all -
    # they join to STAGE_NAMES, which is why they were the pair that got
    # missed the longest.
    for field in ("blocks", "degrades"):
        assert field in guarded, f"{field} names stages and has no near-miss"
