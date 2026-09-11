# CLAUDE.md — metaannot on the lab server

Project context and standing rules. Read `TUTORIAL.md` for the runbook.

## What this is

`metaannot.py` is a single-file metaproteomics annotation pipeline. It takes
FragPipe, DIA-NN or MSstats **label-free** quantification — or FragPipe **TMT**
since v0.3.0 — plus a protein FASTA, annotates every protein with whatever
evidence exists for it, and bins proteins by that evidence so the fraction
KEGG-based analysis discards stops being invisible. It then writes an R
Markdown report and a QFeatures object.

One file. No package to install. `python metaannot.py --help`.

**FragPipe TMT is read, from the per-plex folders and nowhere else.** v0.3.0
added `quant_format: fragpipe_tmt`: `quant_table` becomes the FragPipe run
directory holding the `TMTn/` folders, and the reader takes the reporter
channels from each plex's `ion.tsv` or `peptide.tsv`, joins the plexes at
feature level, and carries the plex into the design as a covariate. Those
intensities are linear, so the log2 path, the roll-up and the
median-of-ratios size factor apply unchanged, and the `tmt:` block configures
the rest. Report the plex-exposure warning the join stage prints.

Every **other** TMT file is still refused by name, and that is a result rather
than a gap: the `tmt-report/` matrices are already log2, median-centred and
protein level; a per-plex `protein.tsv` would report one plex as the whole
experiment; and the TMT `msstats.csv` holds the channels in `Channel <mass>`
columns that no format here maps to samples. Every one of those refusals names
`quant_format: fragpipe_tmt` and the directory to point at instead. Apply what
it names; do not force the file through another format.

## Standing rules

1. **Run `doctor` before any long run, and read its output.** It checks tools,
   databases, the manifest-to-column mapping, and the CPU/RAM split. A missing
   database found by `doctor` costs a second; found at hour six it costs six
   hours.

2. **Install via `doctor --install-plan`, not `--fix`, unless Lior says
   otherwise.** The plan is a script to read before anything touches the
   filesystem. `--fix` downloads immediately; on a shared server with a
   Foldseek AFDB50 database (about 123 GB to download, ~200 GB on disk) that is
   not a decision to make on someone's behalf. `doctor` prints two totals — to
   download and on disk — before it does anything; report both either way, and
   say so if one looks off, since the per-item sizes it prints are
   hand-maintained.

3. **Never pass a bare `--force` without asking.** On its own it discards every
   cached stage, which on a real dataset is days of compute. Scoped, it is
   safe and supported: `--force` now discards the state records of exactly the
   stages selected by `--only`/`--from` and leaves every other stage's record
   intact, so `--force --only <stage>` redoes one stage and the rest still get
   a real signature check on the next run. Ask before the bare form; the scoped
   form is the normal way to redo something.

4. **One run per results directory.** A second run refuses with
   `another metaannot is already running here`. Do not pass `--force-unlock`
   to get past it — find out what the other process is first.

5. **Never delete anything under a results directory.** Stages are cached by
   signature, so deleting an output to "clean up" triggers a silent
   recomputation of that stage with no record of why. **Redo it instead:**
   `python metaannot.py run --config config.yaml --force --only <stage>`.
   `--only` alone reports the stage as cached and skips it, and `--force`
   scoped this way touches only the named stages, so every other stage keeps
   its signature and is still checked properly afterwards. The one thing that
   may be deleted by hand is scratch the run never records —
   `results/foldseek/tmp*`, `results/foldseek/tmpc`, `results/cluster/tmp` —
   which is never cleaned up and can reach hundreds of GB.

6. **Long runs go in tmux.** An SSH drop mid-run leaves partial state. It
   resumes correctly, but only if the process was allowed to record what it
   finished.

7. **Check ID overlap before committing to a run.** The single most common
   failure is protein identifiers in the FASTA not matching those in the
   eggNOG table or the quant file. metaannot aborts below 50% overlap and
   names the transform that fixes it; where the mismatch is a prefix the
   search database added (`uhgpSM_`, `HUMANHOST_`) it names the prefix to put
   in `emapper_strip_id_prefix` instead. Apply what it names. Do not raise
   `emapper_min_coverage` to get past it.

8. **Fifteen of the twenty-one stages have been run on real data. Six have
   not.** One label-free dataset went end to end — FragPipe
   `combined_peptide.tsv` + manifest + precomputed eggNOG, 38,204 proteins —
   through `emapper`, `pfam`, `dbcan`, `diamond`, `cluster`, `ncbifam`,
   `kofam`, `interpro`, `signalp`, `tmbed`, `esmfold`, `foldseek`,
   `integrate`, `finalise` and `join`, and both the report and the R object
   were built from it. A 3-plex subset of an 8-plex FragPipe TMT run has also
   completed end to end including the report. What has still never run on real
   data is `smorf`, `context`, `hhblits`, `jackhmmer`, `unipept` and
   `taxonomy` — all off by default — so a run that enables one of those is
   validation rather than production: sanity-check its counts against the
   expectations in `TUTORIAL.md` and say plainly when something looks wrong.
   `3p_profile_only` is fed only by `hhblits` and `jackhmmer`, so no run has
   ever put a protein in it.

9. **A `parsed 0 ... from a non-empty file` warning is not noise.** It means a
   tool's output is truncated or in an unexpected format, so every protein
   silently loses that evidence and the bins shift. Stop and investigate.

10. **Report numbers, do not summarise them away.** If 71% of proteins land in
   the KO-less bins, say 71%. If a stage warns, surface the warning.

11. **To watch a run, use the console; do not write something that touches the
   results directory.** `python3 console/console.py --root <dir>` on the
   machine doing the run serves one read-only page over a mode-0700 UNIX
   socket, reached with `ssh -N -L 8080:<socket> <host>` and a browser. It is
   safe to point at a job that is already running precisely because it writes
   nothing anywhere near it. `tail -f results/metaannot.log` and reading
   `.metaannot_state.json` are fine too. Writing into a results directory in
   order to monitor one is not — see rule 4: one writer per results directory,
   and a watcher that writes is a second writer.

## Where things live

| | |
|---|---|
| tool | `metaannot.py` (this directory) |
| console | `console/console.py` — a read-only watcher, its own program |
| config | `config.yaml` |
| databases | the large storage array, **not** the OS disk |
| results | one directory per project, on the large array |
| inputs | FragPipe output directory + `.fp-manifest` |

## Machines

- **Server (Fedora)**: everything CPU — hmmsearch, DIAMOND, MMseqs2, Foldseek,
  InterProScan, the roll-up, the join.
- **RTX laptop**: the two GPU stages only, `tmbed` and `esmfold`. Results are
  rsynced back and metaannot *adopts* them rather than recomputing.
- **MacBook (R stack)**: the report and the R object, if R is not on the server.

## Things that are deliberate, do not "fix" them

- Global KEGG maps (01100, 01110, …) are excluded from the pathway test.
- `-dp` on InterProScan disables the precalculated lookup, which only holds
  UniParc matches and is useless for novel metagenome ORFs.
- Foldseek searches its targets serially. Two 100 GB+ indices at once thrash.
- `MANUAL` items in `doctor` are not oversights; they are licence-gated or
  version-specific and must not be automated.
- The taxon reference is a median of ratios, not a sum. A sum is biased by any
  strongly changing member of the taxon.
- Isobaric input reached through the **wrong** `quant_format` is refused, not
  read. Returning the pooled MS1 column as the sole "sample", or a
  `tmt-report/` matrix as though its values were linear, was silently wrong; a
  loud death that names `quant_format: fragpipe_tmt` is the intended
  behaviour. Do not weaken a refusal to get a TMT run through — point the
  supported reader at the run directory instead.
- A FragPipe intensity of `0` means "not quantified", not "measured as zero",
  so it is read as missing (`zero_intensity_is_missing: true`, and the run
  logs how many cells that was). Summing zeros as real values turns
  missingness into fold change. Set it false only to reproduce someone else's
  numbers.
- A **FIFO** at an input works where the run reads that input **once**, and is
  refused immediately where it reads it more than once. That is the design, and
  the boundary is measured rather than assumed: a pipe can be drained once, so
  `mkfifo p; zcat big.faa.gz > p &` really does feed this tool from a disk with
  no room on it — at `quant_table`, at `manifest` and at each
  `emapper_precomputed` entry, which a default run opens once each — and cannot
  work at `proteins_faa`, which it opens three times. The counts come from
  `INPUT_READ_SITES` and are computed for the config in hand, so the answer
  changes with the config; `doctor` prints them. **Do not quote the default
  list at somebody**: with `run.taxonomy` on the manifest is read twice and a
  pipe there is refused, and the database paths are in the plan too —
  `unipept.result` once, and each `.dmp` under `db.ncbi_taxonomy` once per
  stage that builds a taxonomy, which is two when `run.taxonomy` and
  `taxon_rank` are both set. On `quant_format: fragpipe_tmt` the plan counts
  the files inside the run directory — each plex's level file, its annotation,
  its `psm.tsv` where `tmt.min_purity` reads one — because `quant_table` names
  a directory there and a directory is not what gets opened. One read in the
  plan is CONDITIONAL and marked: `peptide_features()` re-reads the quant
  table with the peptide-only reader when the full one refuses it, which turns
  on the table and not on the config, so the plan counts it (a pipe cannot be
  promised) and a measurement of a completed run is not held to it. Read the
  row. Where a pipe is allowed, `run`
  opens it non-blocking, says on the log that it is waiting, and waits
  `fifo_wait_s` (6 hours) **for each next byte** before dying with the path in
  the message and with whether anything was holding the write end. A stream
  that ends early is an error, not a short file. `doctor` never waits at all
  and refuses a FIFO on the row, because a command whose job is to answer
  before the run is worthless if it can block. Do not "fix" any of it: a `run`
  that refused every pipe would remove a working workflow, a `run` that waited
  at a multi-read input would sell six hours for a failure it could have
  reported in the first second, and a `doctor` that waited would remove the
  only thing that can tell you about the pipe first.
- A repeated sample name in an `.fp-manifest` is a **fraction**, not a
  duplicate. Fraction rows are collapsed to one sample; the manifest is only
  refused when the rows genuinely disagree about the design.
- The console's three refusals are the design, not gaps to fill. It writes
  nothing into a results directory; it does not import metaannot, taking
  everything it knows from `metaannot.py describe --json`; and it never asserts
  that a run is dead, because a stopped heartbeat is not a stopped process and
  the engine's rule is that unprovable means alive. So there is no `do_POST`,
  no `--force-unlock` button and no death verdict to add. `CONSOLE_VERSION` is
  its own number and does not track `__version__`.
