"""stage_esmfold's failure handling, with torch and the model faked.

The stage is untestable on this machine's real terms - it wants a CUDA card
and 3 GB of weights - but the part that has actually cost time is not the
folding, it is what happens when a fold fails partway through a long run. That
part is pure control flow, so it is faked here rather than skipped.
"""
import io
import os
import sys
import types

import pytest


# ----------------------------------------------------------------------
# a torch and a transformers that do what the stage asks and nothing else
# ----------------------------------------------------------------------
class _OOM(RuntimeError):
    pass


def _fake_torch(free_gb=32.0, total_gb=32.0):
    t = types.ModuleType("torch")

    class _NoGrad:
        def __enter__(self):
            return None

        def __exit__(self, *a):
            return False

    t.no_grad = _NoGrad
    G = 1024 ** 3
    cuda = types.SimpleNamespace(
        is_available=lambda: True,
        get_device_name=lambda i: "Fake GPU",
        empty_cache=lambda: None,
        mem_get_info=lambda: (int(free_gb * G), int(total_gb * G)),
        OutOfMemoryError=_OOM,
    )
    t.cuda = cuda
    return t


class _Model:
    """Folds by script: `plan` maps sequence id -> list of outcomes.

    An outcome is either None (fold succeeds) or an exception to raise. The
    list is consumed one entry per attempt, so ['boom', None] is a sequence
    that fails once and succeeds on the retry at a smaller chunk.
    """

    def __init__(self, plan, seqs):
        self.plan = plan
        # The stage folds a CLEANED copy of anything carrying a stop codon
        # or a non-standard residue, so the string that arrives here is not
        # always the string on disk. Rather than re-implement the cleaning
        # rule this file is meant to be testing, an unrecognised sequence is
        # matched to the raw one it agrees with most.
        self.raw = list(seqs)
        self.by_seq = {s: pid for pid, s in seqs}
        self.folded = []
        self.chunks = []
        self.trunk = types.SimpleNamespace(set_chunk_size=self.chunks.append)
        self.attempts = {}

    def _pid(self, seq):
        if seq not in self.by_seq:
            agree = lambda r: sum(a == b for a, b in zip(seq, r))
            self.by_seq[seq] = max(self.raw, key=lambda kv: agree(kv[1]))[0]
        if seq not in self.folded:
            self.folded.append(seq)
        return self.by_seq[seq]

    @property
    def by_seq_inv(self):
        """pid -> the sequence actually handed to the model.

        Built from what was folded, not from by_seq, which also holds the raw
        fasta strings and would answer with those for a sequence the stage
        cleaned before folding it.
        """
        return {self.by_seq[s]: s for s in self.folded}

    def eval(self):
        return self

    def cuda(self):
        return self

    def infer_pdb(self, seq):
        pid = self._pid(seq)
        n = self.attempts.get(pid, 0)
        self.attempts[pid] = n + 1
        outcomes = self.plan.get(pid, [])
        exc = outcomes[n] if n < len(outcomes) else None
        if exc is not None:
            raise exc
        # One CA atom is enough: mean_plddt reads the B-factor column.
        return ("ATOM      1  CA  ALA A   1      "
                "0.000   0.000   0.000  1.00 88.50           C\n")


@pytest.fixture
def folding(ma, paths_for, monkeypatch):
    """Run stage_esmfold over `seqs` with a scripted `plan`."""
    def go(seqs, plan=None, reuse=None, free_gb=32.0, total_gb=32.0,
           no_mem_api=False, **cfg_over):
        # A FRESH results directory per call unless reuse is asked for: two
        # calls in one test must not see each other's PDBs, or an assertion
        # about what this call folded silently passes on the last one's output.
        go.n += 1
        cfg, p = reuse if reuse else paths_for(f"results{go.n}")
        cfg.update(cfg_over)
        with open(p.dark, "w", encoding="utf-8") as fh:
            for pid, s in seqs:
                fh.write(f">{pid}\n{s}\n")
        model = _Model(plan or {}, seqs)
        tr = types.ModuleType("transformers")
        tr.EsmForProteinFolding = types.SimpleNamespace(
            from_pretrained=lambda *a, **k: model)
        t = _fake_torch(free_gb, total_gb)
        if no_mem_api:
            # Old torch, or a build without mem_get_info.
            del t.cuda.mem_get_info
        monkeypatch.setitem(sys.modules, "torch", t)
        monkeypatch.setitem(sys.modules, "transformers", tr)
        # fair-esm must not be importable, so the transformers branch is taken.
        monkeypatch.setitem(sys.modules, "esm", None)
        go.last = (cfg, p, model)
        ma.stage_esmfold(cfg, p)
        return cfg, p, model
    go.last = None
    go.n = 0
    return go


_AA = "ACDEFGHIKLMNPQRSTVWY"


def _seqs(n, length=30):
    """Distinct sequences of EQUAL length.

    Distinct so the fake model can tell them apart; equal length so the
    stage's sort-by-length keeps them in file order and the ids in these
    assertions mean what they look like.
    """
    out = []
    for i in range(n):
        body = "".join(_AA[(i + j) % len(_AA)] for j in range(length))
        out.append((f"P{i:03d}", body))
    return out


def _folded(p):
    return sorted(f[:-4] for f in os.listdir(p.structures) if f.endswith(".pdb"))


# ----------------------------------------------------------------------
def test_a_clean_run_folds_every_sequence(folding):
    _, p, _ = folding(_seqs(3))
    assert _folded(p) == ["P000", "P001", "P002"]
    assert os.path.exists(p.struct_done)


def test_a_cuda_fault_is_retried_at_a_smaller_chunk_not_reraised(folding, ma):
    """The fault that cost this project a night was not an OOM.

    A long sequence at a large chunk runs one kernel long enough for the
    display driver to reset the card, which arrives as a plain RuntimeError
    saying 'device not ready'. The old code re-raised anything that was not an
    OOM, so one such sequence destroyed a stage that was 94% done.
    """
    err = RuntimeError("CUDA driver error: device not ready")
    _, p, model = folding(_seqs(3), plan={"P001": [err, None]})
    assert _folded(p) == ["P000", "P001", "P002"]
    assert model.attempts["P001"] == 2
    # The retry has to CHANGE something, or it is just the same fold again.
    assert min(model.chunks) < ma.DEFAULT_CONFIG["esmfold_chunk_size"]


def test_an_oom_is_still_retried_at_a_smaller_chunk(folding, ma):
    _, p, model = folding(_seqs(2), plan={"P000": [_OOM("out of memory"), None]})
    assert _folded(p) == ["P000", "P001"]
    assert min(model.chunks) < ma.DEFAULT_CONFIG["esmfold_chunk_size"]


def test_a_sequence_that_fails_twice_is_skipped_and_the_rest_still_fold(folding):
    err = RuntimeError("CUDA driver error: device not ready")
    _, p, _ = folding(_seqs(4), plan={"P001": [err, err]})
    assert _folded(p) == ["P000", "P002", "P003"]
    assert os.path.exists(p.struct_done)


def test_a_skipped_sequence_is_named_in_a_file_not_only_in_the_log(folding):
    err = RuntimeError("CUDA driver error: device not ready")
    _, p, _ = folding(_seqs(3), plan={"P002": [err, err]})
    tbl = os.path.join(p.structures, "esmfold_failed.tsv")
    body = open(tbl, encoding="utf-8").read()
    assert "P002" in body and "device not ready" in body


def test_every_sequence_failing_in_a_row_stops_rather_than_grinding(folding, ma):
    """A wedged card fails every sequence. Walking the rest produces nothing
    and takes hours, so the stage stops and says how far it got."""
    err = RuntimeError("CUDA driver error: device not ready")
    plan = {pid: [err, err] for pid, _ in _seqs(20)}
    with pytest.raises(ma.StageError) as e:
        folding(_seqs(20), plan=plan, esmfold_max_consecutive_failures=3)
    assert "consecutive failures" in str(e.value)
    assert "esmfold_allow_partial" in str(e.value)


def test_stopping_early_keeps_the_structures_already_written(folding, ma):
    """Stopping is not the same as losing the work.

    Every PDB is written as it is folded, so a card that dies at protein 1800
    of 1900 costs the tail and nothing else. If this ever regressed, a long
    run would have to start over.
    """
    err = RuntimeError("CUDA driver error: device not ready")
    plan = {f"P{i:03d}": [err, err] for i in range(3, 20)}
    with pytest.raises(ma.StageError):
        folding(_seqs(20), plan=plan, esmfold_max_consecutive_failures=3)
    _, p, _ = folding.last
    assert _folded(p) == ["P000", "P001", "P002"]
    # No .done marker: the stage did not finish, so foldseek must not run.
    assert not os.path.exists(p.struct_done)


def test_allow_partial_goes_on_with_what_folded(folding):
    err = RuntimeError("CUDA driver error: device not ready")
    plan = {f"P{i:03d}": [err, err] for i in range(3, 20)}
    _, p, _ = folding(_seqs(20), plan=plan,
                      esmfold_max_consecutive_failures=3,
                      esmfold_allow_partial=True)
    assert _folded(p) == ["P000", "P001", "P002"]
    assert os.path.exists(p.struct_done),         "with allow_partial on, foldseek must be allowed to search what folded"


def test_a_stopped_run_resumes_instead_of_refolding(folding, ma):
    """The second run folds only what the first one did not."""
    err = RuntimeError("CUDA driver error: device not ready")
    seqs = _seqs(8)
    plan = {f"P{i:03d}": [err, err] for i in range(3, 8)}
    with pytest.raises(ma.StageError):
        folding(seqs, plan=plan, esmfold_max_consecutive_failures=3)
    cfg, p, first = folding.last
    assert _folded(p) == ["P000", "P001", "P002"]

    # Same results directory, a card that now works.
    _, p2, second = folding(seqs, plan={}, reuse=(cfg, p))
    assert _folded(p2) == [pid for pid, _ in seqs]
    assert set(second.attempts) == {f"P{i:03d}" for i in range(3, 8)},         "a resumed run must not refold what is already on disk"


def test_the_plddt_table_gains_the_resumed_rows_and_keeps_the_old_ones(folding,
                                                                      ma):
    err = RuntimeError("CUDA driver error: device not ready")
    seqs = _seqs(8)
    plan = {f"P{i:03d}": [err, err] for i in range(3, 8)}
    with pytest.raises(ma.StageError):
        folding(seqs, plan=plan, esmfold_max_consecutive_failures=3)
    cfg, p, _ = folding.last
    folding(seqs, plan={}, reuse=(cfg, p))
    rows = ma.read_plddt(p.structures)
    assert sorted(rows) == [pid for pid, _ in seqs]


def test_a_stop_codon_is_folded_as_a_cleaned_copy(folding, ma):
    """ESMFold refuses '*' outright, and Prokka emits it.

    The whole stage used to die on protein 350 of 1912 because of this.
    """
    seqs = _seqs(3)
    seqs[1] = (seqs[1][0], seqs[1][1] + "*")
    _, p, model = folding(seqs)
    assert _folded(p) == ["P000", "P001", "P002"]
    assert "*" not in model.by_seq_inv["P001"]


def test_a_non_standard_residue_is_replaced_rather_than_refused(folding):
    seqs = _seqs(2)
    seqs[0] = (seqs[0][0], "U" + seqs[0][1][1:])
    _, p, model = folding(seqs)
    assert _folded(p) == ["P000", "P001"]
    assert model.by_seq_inv["P000"].startswith("X")


def test_a_sequence_over_max_len_structure_is_not_attempted(folding):
    seqs = [("short", "ACDEFGHIKL"), ("long", "A" * 900)]
    _, p, model = folding(seqs, max_len_structure=100)
    assert _folded(p) == ["short"]
    assert "long" not in model.attempts


# ----------------------------------------------------------------------
# the VRAM cap
# ----------------------------------------------------------------------
def test_a_sequence_that_will_not_fit_in_vram_is_not_attempted(folding, capsys):
    """ESMFold raises no OOM when it stops fitting.

    The driver pages device memory to host RAM instead, so the fold does not
    fail - it becomes one to two orders of magnitude slower. Measured on a
    16 GB card: 22 s at 478 aa, 140 s at 481 aa, 2,053 s at 491 aa. A run that
    crosses that line does not stop, it stops being finishable, so the cap has
    to come from the memory actually free rather than from a fixed number.
    """
    seqs = [("fits", "A" * 100), ("huge", "A" * 900)]
    # 1 GB free, minus the 0.5 GB reserve, at 21000 bytes/pair -> ~159 aa
    _, p, model = folding(seqs, free_gb=1.0, total_gb=16.0)
    assert _folded(p) == ["fits"]
    assert "huge" not in model.attempts
    err = capsys.readouterr().err
    assert "sequences up to" in err
    assert "will NOT be folded" in err


def test_the_sequences_the_cap_skipped_are_kept_where_they_can_be_folded(
        folding, ma):
    """A count in a log is not a work-list.

    The advice has always been "fold them elsewhere and drop the models in",
    and it always left the reader to reconstruct WHICH ones from a number: on
    the first full real run the cap excluded 1,462 of 2,000 and the list
    existed nowhere at all.
    """
    seqs = [("fits", "A" * 60), ("huge", "C" * 900), ("bigger", "D" * 950)]
    _, p, _ = folding(seqs, free_gb=1.0, total_gb=16.0)
    rows = [l.split("\t") for l in
            open(f"{p.structures}/not_folded.tsv", encoding="utf-8")
            .read().splitlines()]
    assert rows[0] == ["protein_id", "length", "limit_aa", "limit_from"]
    assert {r[0] for r in rows[1:]} == {"huge", "bigger"}
    assert [r[1] for r in rows[1:]] == ["900", "950"]
    assert {r[3] for r in rows[1:]} == {"the card's free VRAM"}
    # and the sequences themselves, byte for byte, so the file can go
    # straight to another card rather than being rebuilt from dark.faa
    faa = dict(ma.read_fasta(f"{p.structures}/not_folded.faa"))
    assert faa == {"huge": "C" * 900, "bigger": "D" * 950}


def test_nothing_skipped_still_writes_the_list_so_absence_means_one_thing(
        folding):
    # An absent file that means either "nothing was skipped" or "this stage
    # never ran" is the ambiguity esmfold_failed.tsv used to have.
    seqs = [("a", "A" * 60), ("b", "C" * 70)]
    _, p, _ = folding(seqs, free_gb=40.0, total_gb=48.0)
    assert open(f"{p.structures}/not_folded.tsv", encoding="utf-8").read() \
        == "protein_id\tlength\tlimit_aa\tlimit_from\n"
    assert os.path.getsize(f"{p.structures}/not_folded.faa") == 0


def test_the_skipped_warning_prices_the_vram_the_longest_one_needed(
        folding, ma, capsys):
    """A count of skipped proteins is a complaint; a number of GB is a step.

    And it is the datum that settles what this cap actually is: the
    coefficient reproduces both measured points exactly, so a low cap means
    the card had little free VRAM, not that the constant is a guess.
    """
    seqs = [("fits", "A" * 60), ("huge", "C" * 900)]
    cfg, _, _ = folding(seqs, free_gb=1.0, total_gb=16.0)
    need = (900 ** 2 * cfg["esmfold_bytes_per_residue_pair"]
            + cfg["esmfold_vram_reserve_gb"] * 1024 ** 3) / 1024 ** 3
    err = capsys.readouterr().err
    assert f"about {need:.1f} GB free" in err, err
    assert "not the coefficient that decides this cap" in err


def test_max_len_structure_skipping_does_not_price_vram_it_did_not_decide(
        folding, capsys):
    seqs = [("fits", "A" * 60), ("huge", "C" * 900)]
    folding(seqs, free_gb=40.0, total_gb=48.0, max_len_structure=100)
    err = capsys.readouterr().err
    assert "will NOT be folded" in err and "max_len_structure" in err
    assert "GB free with the weights resident" not in err, \
        "the VRAM had nothing to do with this cap"


def test_the_failure_record_is_written_even_when_nothing_failed(folding):
    # Its absence used to mean either "nothing failed" or "this results
    # directory predates the record", and the shortfall message had to hedge
    # across both of those at once. See the have_tbl test in test_binning.py.
    seqs = [("a", "A" * 60), ("b", "C" * 70)]
    _, p, _ = folding(seqs)
    assert open(f"{p.structures}/esmfold_failed.tsv", encoding="utf-8").read() \
        == "protein_id\tlength\terror\n"


def test_the_warning_says_it_is_memory_and_says_what_to_do(folding, capsys):
    seqs = [("fits", "A" * 60), ("huge", "A" * 900)]
    folding(seqs, free_gb=1.0, total_gb=16.0)
    err = capsys.readouterr().err
    # It must not read as a failure: these were never attempted.
    assert "never attempted, not as" in err
    assert "more memory" in err
    # And it must name the real hazard, which is silent slowness, not an OOM.
    assert "no OOM raised" in err


def test_a_roomy_card_is_not_capped_below_max_len_structure(folding, capsys):
    seqs = [("a", "A" * 600), ("b", "A" * 650)]
    _, p, _ = folding(seqs, free_gb=40.0, total_gb=48.0, max_len_structure=700)
    assert _folded(p) == ["a", "b"],         "a card with headroom must not be throttled by the guard"
    assert "will NOT be folded" not in capsys.readouterr().err


def test_max_len_structure_still_wins_when_it_is_the_tighter_limit(folding,
                                                                  capsys):
    seqs = [("short", "A" * 80), ("long", "A" * 300)]
    _, p, _ = folding(seqs, free_gb=40.0, total_gb=48.0, max_len_structure=100)
    assert _folded(p) == ["short"]
    assert "max_len_structure" in capsys.readouterr().err


def test_the_cap_can_be_switched_off(folding, capsys):
    seqs = [("huge", "A" * 900)]
    _, p, _ = folding(seqs, free_gb=1.0, total_gb=16.0,
                      max_len_structure=1000, esmfold_vram_cap=False)
    assert _folded(p) == ["huge"]
    assert "sequences up to" not in capsys.readouterr().err


def test_the_coefficient_is_configurable(folding):
    seqs = [("mid", "A" * 300)]
    # 8.5 GB usable at 21000 b/pair -> ~636 aa, so 300 aa folds
    _, p, _ = folding(seqs, free_gb=9.0, total_gb=16.0)
    assert _folded(p) == ["mid"]
    # the same card declared four times hungrier -> ~318 aa, still folds
    _, p2, _ = folding(seqs, free_gb=9.0, total_gb=16.0,
                       esmfold_bytes_per_residue_pair=84000)
    assert _folded(p2) == ["mid"]
    # sixteen times hungrier -> ~159 aa, now excluded
    _, p3, _ = folding(seqs, free_gb=9.0, total_gb=16.0,
                       esmfold_bytes_per_residue_pair=336000)
    assert _folded(p3) == []


def test_a_torch_without_mem_get_info_folds_anyway_and_says_the_guard_is_off(
        folding, capsys):
    """The guard must degrade to the old behaviour, not fail the stage."""
    _, p, _ = folding([("a", "A" * 200)], no_mem_api=True, free_gb=1.0,
                      total_gb=16.0, max_len_structure=700)
    assert _folded(p) == ["a"],         "a torch that cannot report free VRAM must not stop the fold"
    err = capsys.readouterr().err
    assert "cannot read free VRAM" in err
    assert "max_len_structure alone" in err


def test_the_real_measurement_that_prompted_this_is_reproduced(folding):
    """The exact case that took the host down: a 16 GB card holding an 11.2 GB
    model, asked to fold 505-684 aa. Every one of those must be refused."""
    seqs = ([(f"ok{i}", "A" * 470) for i in range(2)]
            + [(f"over{i}", "A" * L) for i, L in enumerate((505, 600, 684))])
    # 16 GB card, 11.2 GB of weights resident -> 4.8 GB free
    _, p, model = folding(seqs, free_gb=4.8, total_gb=16.0,
                          max_len_structure=700)
    assert _folded(p) == ["ok0", "ok1"]
    assert not any(k.startswith("over") for k in model.attempts),         "the sequences that bugchecked the host must never be submitted"


def test_a_finished_stage_does_not_touch_the_gpu_at_all(folding, capsys):
    """Loading the weights is itself a risk on a virtualised GPU.

    A resumed run whose work-list is already complete used to upload ~11 GB to
    the card before discovering it had nothing to do. That upload is the
    operation a host bugcheck landed in, so a finished stage must not perform
    it.
    """
    seqs = _seqs(3)
    _, p, first = folding(seqs)
    assert _folded(p) == ["P000", "P001", "P002"]
    cfg = folding.last[0]
    capsys.readouterr()
    _, p2, second = folding(seqs, reuse=(cfg, p))
    assert second.attempts == {}, "nothing should have been folded"
    err = capsys.readouterr().err
    assert "the GPU is not touched at all" in err
    assert "backend on" not in err, "the model must not have been loaded"


def test_a_work_list_entirely_over_max_len_structure_skips_the_gpu(folding,
                                                                  capsys):
    _, p, model = folding([("huge", "A" * 900)], max_len_structure=100)
    assert model.attempts == {}
    err = capsys.readouterr().err
    assert "the GPU is not touched at all" in err
    assert "backend on" not in err


# --- interrupted writers ----------------------------------------------
# symptom: the resume rule is "the .pdb is there, so it is folded", and the
# file was written in place. A run killed between open() and the flush — a
# wedged card, a bugcheck, a Ctrl-C, all of which this stage has seen — left a
# 0-byte structure that every later run counted as done. The protein was then
# permanently absent from Foldseek with nothing anywhere to say why.
def test_a_structure_is_renamed_into_place_never_written_in_place(
        folding, ma, monkeypatch):
    seen = []
    real = ma.atomic_out
    import contextlib

    @contextlib.contextmanager
    def spy(path):
        with real(path) as tmp:
            seen.append((path, tmp))
            yield tmp
    monkeypatch.setattr(ma, "atomic_out", spy)
    cfg, p, _ = folding(_seqs(3))
    pdbs = [(dst, tmp) for dst, tmp in seen if dst.endswith(".pdb")]
    assert len(pdbs) == 3, seen
    for dst, tmp in pdbs:
        assert tmp != dst, "the structure was written straight to its name"
        assert os.path.basename(tmp).startswith("."), \
            f"{tmp} is not hidden from glob('*.pdb')"
        assert tmp.endswith(".pdb")
    assert _folded(p) == ["P000", "P001", "P002"]
    assert not [f for f in os.listdir(p.structures)
                if ma.ATOMIC_SUFFIX in f], "a temp was left behind"


def test_an_empty_structure_left_by_a_killed_writer_is_refolded(
        folding, ma, capsys):
    cfg, p, model = folding(_seqs(2))
    open(f"{p.structures}/P001.pdb", "w", encoding="utf-8").close()
    capsys.readouterr()
    folding(_seqs(2), reuse=(cfg, p))
    err = capsys.readouterr().err
    assert "holds no atom" in err
    assert os.path.getsize(f"{p.structures}/P001.pdb") > 0
    assert "ATOM" in io.open(f"{p.structures}/P001.pdb",
                             encoding="utf-8").read()


def test_a_structure_with_no_atom_records_is_refolded(folding, ma, capsys):
    cfg, p, model = folding(_seqs(2))
    io.open(f"{p.structures}/P000.pdb", "w", encoding="utf-8").write(
        "HEADER    truncated before the coordinates\n")
    capsys.readouterr()
    folding(_seqs(2), reuse=(cfg, p))
    assert "ATOM" in io.open(f"{p.structures}/P000.pdb",
                             encoding="utf-8").read()


def test_a_complete_structure_is_still_skipped(folding, ma, capsys):
    # the guard must not turn a resume into a refold: that is the whole point
    # of the stage keeping its output.
    cfg, p, model = folding(_seqs(3))
    before = dict(model.attempts)
    capsys.readouterr()
    _, _, model2 = folding(_seqs(3), reuse=(cfg, p))
    assert model2.attempts == {}, "a finished structure was folded again"
    assert before and "holds no atom" not in capsys.readouterr().err


def test_a_torn_structure_does_not_make_the_stage_decide_it_is_finished(
        folding, ma, capsys):
    # the early exit runs BEFORE the loop and before the GPU is touched, so a
    # guard in the loop alone is unreachable: the stage would announce that
    # everything was already folded and return, leaving the damaged file for
    # Foldseek.
    cfg, p, _ = folding(_seqs(2))
    open(f"{p.structures}/P000.pdb", "w", encoding="utf-8").close()
    capsys.readouterr()
    _, _, model = folding(_seqs(2), reuse=(cfg, p))
    err = capsys.readouterr().err
    assert "the GPU is not touched at all" not in err
    assert list(model.attempts) == ["P000"], model.attempts
