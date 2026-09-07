# CLAUDE.md — metaannot on the lab server

Project context and standing rules. Read `TUTORIAL.md` for the runbook.

## What this is

`metaannot.py` is a single-file metaproteomics annotation pipeline. It takes
FragPipe (or DIA-NN, or MSstats) **label-free** quantification plus a protein
FASTA, annotates
every protein with whatever evidence exists for it, and bins proteins by that
evidence so the fraction KEGG-based analysis discards stops being invisible.
It then writes an R Markdown report and a QFeatures object.

One file. No package to install. `python metaannot.py --help`.

**FragPipe TMT output is not supported.** Reporter-ion channels are not read.
The per-plex `TMTn/*.tsv` tables are now **refused**: the loader spots the
`Intensity <sample>` reporter columns and dies rather than collapsing to the
single MS1 precursor column and reporting one "sample". The other TMT routes
still fail on their own terms — the `tmt-report/` matrices are log2 and
median-centred and would be log-transformed a second time, and the TMT
`msstats.csv` breaks the parser outright — so neither is safe either. Do not
point this tool at a TMT run and do not read a number out of one; say so
instead. Isobaric support is a separate piece of work.

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

8. **Only the pre-search half has been run on real data.** One label-free
   dataset has gone end to end (FragPipe `combined_peptide.tsv` + manifest +
   precomputed eggNOG, 38,204 proteins): manifest parsing, the id join, the
   shared-peptide rule, binning, roll-up and design recovery all work. Every
   external search stage was off in that run, and the report and R object have
   only ever run on synthetic data. So treat a run that enables search stages,
   or that reaches the report, as validation rather than production:
   sanity-check counts at every phase against the expectations in
   `TUTORIAL.md` and say plainly when something looks wrong.

9. **A `parsed 0 ... from a non-empty file` warning is not noise.** It means a
   tool's output is truncated or in an unexpected format, so every protein
   silently loses that evidence and the bins shift. Stop and investigate.

10. **Report numbers, do not summarise them away.** If 71% of proteins land in
   the KO-less bins, say 71%. If a stage warns, surface the warning.

## Where things live

| | |
|---|---|
| tool | `metaannot.py` (this directory) |
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
- Isobaric input is **refused, not read**. Returning the pooled MS1 column as
  the sole "sample" was silently wrong; a loud death is the intended
  behaviour. Do not weaken it to get a TMT run through.
- A FragPipe intensity of `0` means "not quantified", not "measured as zero",
  so it is read as missing (`zero_intensity_is_missing: true`, and the run
  logs how many cells that was). Summing zeros as real values turns
  missingness into fold change. Set it false only to reproduce someone else's
  numbers.
- A repeated sample name in an `.fp-manifest` is a **fraction**, not a
  duplicate. Fraction rows are collapsed to one sample; the manifest is only
  refused when the rows genuinely disagree about the design.
