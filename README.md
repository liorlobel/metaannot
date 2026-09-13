# metaannot

One file. `metaannot.py` contains the pipeline, the R report, and the config
template.

```bash
python metaannot.py init                        # write config.yaml
python metaannot.py doctor                      # what is missing, and how to get it
python metaannot.py doctor --install-plan i.sh  # a script to review, then run
python metaannot.py doctor --fix                # download and install, after confirming
python metaannot.py describe --json             # the config shape and stage graph, as JSON
python metaannot.py all                         # annotate, report, R object
```

`doctor` prints the exact command for anything missing and can run them for
you. It shows the **total to download and the total on disk** as two separate
figures first, and refuses without confirmation. AFDB50 is the big one: about
**123 GB to download and ~200 GB on disk** after extraction (the Cα-only search
then wants ~150 GB of RAM). Full `Alphafold/UniProt` is ~491 GB and PDB ~2 GB.
The per-item sizes are still hand-maintained, so read the total as indicative
and check a big one against the provider before approving it. Items needing a
licence or a click-through (SignalP 6.0, InterProScan, several targeted
databases) are marked `MANUAL` and never attempted; `doctor` says where to get
each.

Downloads use `curl -fL`, which fails on an HTTP error status but *does* follow
redirects, so a URL repointed at a landing page would otherwise write an HTML
page to disk with exit code 0 and pass the existence check forever. The
generated commands for Pfam-A, the dbCAN HMMs, NCBIfam, UniRef50 and the
targeted DIAMOND FASTAs are now each followed by a guard that sniffs the first
512 bytes, deletes the file and exits non-zero with `is an HTML page, not the
database` if it looks like markup. The two that have no guard — the KOfam
tarball and `ko_list`, and the taxdump — are compressed archives whose `tar` or
`gunzip` step fails loudly on a web page anyway. None of this covers a file you
fetched by hand, so for those check the format yourself (`head -1` of an HMM
file should start with `HMMER3/`). HMM libraries are additionally required to be
`hmmpress`ed, which rejects a truncated download on its own. Paths and URLs from the config are shell-quoted, so a directory
containing a space works and one containing shell metacharacters is treated as
a literal name.

Download URLs live in the config under `sources:`, so a moved link is an edit
there rather than a patch to the tool.

Relative paths in the config resolve against the config file, so a project
runs from any directory. The stage cache is keyed to a separate
`SIGNATURE_VERSION`, not the tool version, so a patch release does not discard
days of compute.

An unrecognised key in `config.yaml` is reported with a spelling suggestion
rather than silently ignored — `run: {unipep: true}` used to leave the stage
disabled with nothing said. The check stops at the free-form blocks
(`tool_args`, `db.diamond`, `sources.diamond`, `diamond_weights`, `vfdb_category_weights`,
`diamond_evalues`, `diamond_min_pidents`), whose keys are user-chosen and
cannot be checked, so a typo there is silent and simply has no effect. Proof-read those blocks by hand.
The `analysis:` block is **not** free-form and **is** checked: its keys are
exactly the Rmd's params, so `fdrr: 0.01` is reported rather than written into
the Rmd header as a spurious param while `fdr: 0.05` stays quietly in force.
Read the `== config ==` block of `doctor` before a run you intend to publish.

A second top-level `run:` or `db:` key in a config that already has one is
refused. YAML itself would retain only the later mapping and silently throw the
earlier block away, so `load_config` installs a no-duplicate loader: the run
exits naming the duplicate key and both line numbers instead. Merge the two
into a single block.

Requires python3 + pandas + pyyaml. R is needed only for the report; external
tools only by the stages that use them. `doctor` says which are missing.

## Tests

```bash
pip install pytest && pytest -q          # a few minutes
pytest -q -m slow                        # the rest: resume, parallel vs serial
```

A healthy default run on this tree is **1753 passed, 1 skipped, 6 xfailed, 38
deselected**, in three to five minutes depending on the machine. Those numbers
are the only yardstick you have for deciding whether your checkout is the one
this document describes, so they are counted rather than estimated. The 38
deselected are the `slow` marker, and they are the second command above.
`pytest -q -m R` selects the 71 R tests, which the default run **already
includes**: they skip rather than fail when `Rscript` or one of its packages is
absent, so on a machine with no R the same run reports 1682 passed and 72
skipped. The single skip here is a Windows-only test pinning a refusal that
cannot happen on POSIX.

One more test skips itself where the machine cannot give it what it asks for:
`test_a_full_log_disk_costs_the_log_line_and_not_the_run` drives a REAL
filesystem with no free blocks — a small ram disk, made, filled, unmounted and
detached inside the test — and the recipe for that is macOS's, so it skips
elsewhere and the count above is one lower there. It skips rather than
substituting a mock because the measurement is the point: on a full
filesystem the buffered write succeeds and the `flush()` after it raises,
which is the reason the guard covers both.

Windows is not run from this machine, so what follows is read off the code and
the test markers rather than measured. `doctor --fix` there is not a recipe
that emits `mkdir -p` any more — since v0.4.0 it is an outright refusal, because
every command it generates is POSIX shell handed to `cmd.exe`, where
`mkdir -p C:\db` creates a directory called `-p` and each step can exit 0 having
done nothing. It names the alternative instead: `doctor --install-plan
install.sh` on Windows, then `wsl bash install.sh` where the tools live.
`--install-plan` itself still works there and every other `doctor` check runs;
only `--fix` is refused. The three `--fix` tests carry that as their skip
reason on Windows, and a fourth runs **only** on Windows to pin the refusal —
it is the one skip in the count above.

Eleven test functions that signal a child process are skipped there as well,
eleven collected items in a default run, and the arithmetic is a coincidence
rather than a rule: `test_a_killed_run_releases_the_results_lock` and
`test_a_killed_run_takes_its_tool_and_the_tools_own_child_with_it` are each
parametrized `SIGTERM` and `SIGHUP`, which adds two, while
`test_a_killed_runs_tail_never_lands_on_the_run_that_replaced_it` and
`test_a_kill_9_leaves_an_orphan_and_the_next_run_refuses_to_join_it` carry the
`slow` marker and are deselected, which takes two away. They do not all carry
the same reason. `test_an_interrupted_run_leaves_parseable_state_and_resumes`
is the `SIGINT` one: `send_signal(SIGINT)` is unsupported on Windows and
`CTRL_C_EVENT` goes to the whole console group including the test runner, and
`test_ctrl_c_unwinds_and_stamps_where_a_sigterm_cannot` and
`test_ctrl_c_still_stops_the_tools_now_that_the_tty_no_longer_does_it` need
both halves — the `SIGINT` above and the `SIGTERM` below. The other eight —
`test_a_killed_run_releases_the_results_lock`,
`test_a_run_killed_that_way_resumes_without_force_unlock`,
`test_a_run_killed_with_sigterm_releases_the_results_lock`,
`test_a_sigterm_releases_the_lock_and_leaves_a_resumable_trace`,
`test_a_sigtermed_run_is_exactly_the_case_the_heartbeat_exists_for`,
`test_a_killed_run_takes_its_tool_and_the_tools_own_child_with_it`,
`test_a_kill_9_leaves_an_orphan_and_the_next_run_refuses_to_join_it` and
`test_a_killed_runs_tail_never_lands_on_the_run_that_replaced_it` — are skipped
because `TerminateProcess` runs no handler on Windows, so a `SIGTERM` cannot be
delivered to a child there at all and none of what they assert (a released
lock, a `128 + N` exit status, the WARN line the handler writes, a dead tool
and a dead grandchild) can happen. There are no `POSIX` process groups to kill
there either, which is the other half of the newer three. All eleven are
verified under Linux and on this Mac.

Five of those markers are new in v0.5.0 and three more arrived with the tool-
group kill, and every one of them was **deduced, not measured**: they describe
a machine that is not the machine they were written on, and the deduction rests
on the same `TerminateProcess` fact the markers written before them give. If
you run the suite on Windows and one of them would in fact have passed, the
marker is wrong and worth removing — that is a better failure than the silent
red it replaced.

**One further test signals a child and carries no Windows marker whatsoever**,
so it is not a skip there: it runs.
`test_the_run_heartbeat_advances_while_a_stage_is_running` sends a `SIGTERM`
only to stop the run, and everything it asserts afterwards — that the lock
survives, that a second run refuses it — is what `TerminateProcess` leaves
behind anyway, so it is left alone deliberately rather than overlooked. Read
the count above as "eleven items are skipped on Windows", not as "the signal
tests are handled on Windows".

The old note here also listed a path test asserting forward slashes among the
Windows failures. That one is gone: v0.4.0 rewrote
`test_home_and_env_vars_in_a_config_path_are_expanded` to compare
`os.path.normcase(os.path.abspath(...))` instead of POSIX strings, so it is
separator-agnostic and passes on both. Stale-lock reclamation is no longer a
Windows exception either: `_holder_is_alive` asks `OpenProcess`, which answers
without touching the process, rather than `os.kill(pid, 0)`, which there calls
`TerminateProcess`. The console's tests are new in this release and have not
been run on Windows at all, so nothing here can say how they behave.

Offline, and needs none of the external tools: where a stage shells out to
hmmsearch, DIAMOND or MMseqs2 the binary is a stub on `PATH` that writes a
canned file, so the stage's own plumbing is exercised without it. Fixtures are
generated from fixed seeds, so nothing binary is committed. The R tests skip
cleanly when `Rscript` or one of its packages is absent, and knit the real
report when they are present.

Each test is named for the defect it protects against and carries a one-line
comment stating the symptom, because most of them exist to stop something
coming back rather than to describe an intended feature. A handful are
`xfail(strict)`: those name guards that are still missing, so a fix turns them
green instead of being forgotten — and because they are strict, a fix that
lands without removing the marker fails the suite rather than passing quietly.
Six are open today, over five tests, and they are live defects the tool
documents rather than hides: an `annotation_pass1.tsv` that is not reproducible
across a resume; a Unipept lineage truncated at the first blank rank; a
`pept2lca` file matching nothing dying with a bare `'verdict'`; an
`emapper.annotations` file with no data rows raising a bare `KeyError` instead
of the "no `#query` header" message every other malformed file gets (two of the
six — the same test over an empty file and a header-only one); and a DIAMOND
database disabled with `""` vanishing without trace, because `resolve_paths()`
drops empty `db.diamond` entries before `doctor` ever sees them, so nothing
anywhere records that you turned it off. That last one is a wish rather than a
regression — the behaviour the tool has today is the one the DIAMOND section
below describes — and the marker is what keeps the wish from being forgotten.

## A worked example

`examples/server-run-plan/` is a real plan for running eight metaproteome datasets on one
server: eight configs plus a four-step runbook (`subset` → `doctor` → run in tmux →
collect). Its paths and dataset names are one lab's, so it is a template rather than
something to run as-is, but it is the shape of a multi-dataset run and it records the
decisions such a run has to make — which stages are on and *why each of the others is
off*, why `min_features_per_protein` is left at 1, which three things to check before
anything long starts (the tool version on the server, the mount point every
path assumes, and the database paths that are inferred rather than confirmed),
and what the whole thing is expected to cost.

## Start from a FragPipe manifest

```yaml
manifest: "experiment.fp-manifest"
quant_table: "combined_peptide.tsv"
quant_format: "fragpipe_peptide"
```

For a **label-free** run the `.fp-manifest` holds path, experiment and
bioreplicate per run, and that is the whole experimental design. metaannot renames the quant columns to the
manifest's sample names, writes `results/quant/design_from_input.tsv`, and
derives the pairwise contrasts. Nothing else is hand-written.

It tries `experiment_bioreplicate`, then `experiment`, then the file basename
against the quant columns, and `doctor` reports whether every manifest run
matched a quant column before a long run rather than failing after one, listing
up to five unmatched runs and, separately, up to four quant columns that no
manifest row claims. It does not print the full mapping.

The manifest is also the sample list, not just a rename table. A quant column
that no manifest row claims is **dropped**, so a run left out of the manifest
is a sample left out of the results — which is why `doctor` reports that second
list before anything long starts.

**A manifest is one row per raw file, not per sample.** A fractionated
acquisition therefore repeats a sample name across its fraction rows — which is
exactly how FragPipe denotes fractions, and FragPipe writes one quant column per
group. metaannot **collapses** those rows to one sample, logs how many samples
were split (`manifest: N sample(s) are split across multiple fraction files`),
and refuses only when the rows of one sample genuinely disagree about the design
(different `data_type`), because collapsing those would invent a sample that
never existed. It also refuses if two *different* manifest samples resolve to the
same quant column, which is a real ambiguity rather than a fraction.

TMT is the case where the `experiment` column is the plex rather than the
condition, so a design derived from it would give plex-versus-plex contrasts —
a batch effect, not a hypothesis. `quant_format: fragpipe_tmt` therefore
ignores the manifest and says so, and takes the condition from the sample names
or from `analysis.metadata` instead; see
[FragPipe TMT](#fragpipe-tmt-isobaric-fragpipe_tmt).

### Identifiers

Everything downstream joins on the FASTA id, which is the first whitespace token
of the header. Two things routinely break that join and both are reported rather
than absorbed:

- FragPipe fills `Protein ID` with `<id> <description>` for a metagenome
  database, a value that can never match the FASTA. The peptide/ion readers
  therefore prefer the **`Protein`** column and reduce whatever they find to its
  first token, warning how many values carried a description.
- A search database often prefixes the ids it was built from (`uhgpSM_`,
  `HUMANHOST_`) while the eggNOG table does not. Set
  `emapper_strip_id_prefix: "uhgpSM_"` (a string or a list) and the prefix is
  stripped from the FASTA id when matching. The coverage diagnostic names the
  prefix that would have matched, so you do not have to guess, and the run
  reports how many rows matched only because of it. It applies to
  `emapper_precomputed` reuse only.

For covariates beyond condition and replicate — sex, age, cage — copy
`design_from_input.tsv`, add columns, and point `analysis.metadata` at it with
a matching `analysis.design_formula`.

## What it does

Every identified protein is binned on the evidence that actually exists for it:

| bin | meaning |
|---|---|
| `1_ko_pathway` | KO, mapped to a specific (non-global) KEGG map |
| `2_ko_orphan` | KO, but no specific pathway map |
| `3_annotated_no_ko` | no KO; informative Pfam (hmmsearch **or** eggNOG's own `PFAMs`) / NCBIfam / InterPro / CAZy / dbCAN / VFDB / MEROPS / CARD / TADB / BAGEL |
| `3d_duf_only` | no KO; only domain evidence is a DUF |
| `3s_structure_only` | no KO and no sequence annotation; confident Foldseek hit to a target that is itself described |
| `3p_profile_only` | no KO, no sequence annotation, no DUF and no fold; only a remote HHblits/jackhmmer profile hit, and only where the hit's target is not itself uncharacterised |
| `4_dark` | no evidence |

**A hit to an uncharacterised target is not a rescue.** `3s_structure_only` and
`3p_profile_only` need more than a hit that clears the thresholds — the
target's *description* has to say something. A description matching `^DUF\d`,
`^UPF\d`, `unknown function`, `uncharacteri…`, `hypothetical`, `family not
named`, `predicted protein` or `putative protein` (case-insensitive; the same
test that decides whether a Pfam or NCBIfam hit is informative) is demoted: the
protein keeps its `foldseek_*` / `hh_*` columns but stays in `4_dark`, or in
`3d_duf_only` if it has a DUF. The run logs the count as a single WARN over
both sources, so quote it next to any 3s/3p number. One exemption, which can
only shrink the rescue claim rather than inflate it: an **empty** description
keeps its evidence, so a bare AlphaFold accession is not demoted — on the
default AFDB50 target the demotion therefore fires rarely, for want of any text
to test.

That is the complete set; `annotation_final.tsv`, `bin_summary.tsv`,
`tier_coverage.tsv` and the report's factor levels use exactly these seven
strings. The tests are applied in
the order KO-with-pathway, KO, sequence annotation, DUF, structure, profile, so a
DUF-only protein stays in `3d_duf_only` however good its fold or profile hit is
— deliberate, since a DUF still names a family, and the run counts how many were
held back so the 3s/3p rescue numbers are read with that in mind.

**eggNOG's own `PFAMs` column counts as domain evidence.** Binning on the
`pfam` stage's `hmmsearch` hits alone meant that with `run.pfam: false` a
protein eggNOG had already assigned a domain to was reported as having no
evidence. On an early eggNOG-only pass of the UC metaproteome — before the
search stages were run, so not the completed UC run quoted later in this file —
6,271 of 17,377 dark proteins (36%) carried an eggNOG Pfam: the dark bin was
inflated by a missing join rather than by biology. The fix is in, so that
figure is a record of the defect rather than something a current run
reproduces; in the completed UC run `4_dark` is 3.8%. An eggNOG Pfam whose only
accession is a DUF/UPF lands the protein in `3d_duf_only` rather than
`3_annotated_no_ko`, on the same rule as an `hmmsearch` DUF. So expect a
smaller `4_dark` than earlier versions of this tool produced, for a join reason
rather than a biological one.

Bins 2 to 4 are what KEGG enrichment silently discards. Global KEGG maps
(01100, 01110, …) are excluded from the pathway test: a protein whose only
"pathway" is *Metabolic pathways* is not pathway-annotated in any useful sense.
In real data almost no protein has *only* global maps, so `2_ko_orphan` means
in practice a KO with no pathway map at all — a module-only or unmapped KO.

## Stages

```
emapper → pfam → dbcan → diamond → signalp → tmbed → cluster
        → ncbifam → kofam → interpro → smorf → context
        → integrate → jackhmmer → hhblits → esmfold → foldseek
        → finalise → unipept → taxonomy → join
```

A fresh `init` config turns on **six** `run:` flags: `eggnog` (which runs the
`emapper` stage), `pfam`, `dbcan`, `diamond`, `cluster` and `join`. `integrate`
and `finalise` have no flag and always run. Everything else is **off** —
`topology` (the `signalp` and `tmbed` stages), `structure` (`esmfold` and
`foldseek`), `context`, `ncbifam`, `kofam`, `interpro`, `hhblits`, `jackhmmer`,
`smorf`, `unipept`, `taxonomy` — because switching each on is a decision about
hardware, a licence or a large database that is not ours to make for you.

`--only` and `--from` take **stage** names, not `run:` flag names: `--only
signalp` works, `--only topology` does not. Naming a stage that its flag has
disabled runs it anyway, which is the point of the option.

So `python metaannot.py all` on a fresh config still needs Pfam-A, the dbCAN
HMMs, the DIAMOND databases and either eggNOG-mapper or an
`emapper_precomputed` table. Run `run --dry-run` and read the plan before the
first real run; `doctor` reports which tools and databases are missing for the
stages you left enabled.

### Added evidence

- **`ncbifam`** — NCBIfam/TIGRFAM HMMs. The stage itself is one more
  `hmmsearch` (`thresholds.ncbifam_cutoff`, default `--cut_tc`), so it is the
  cheapest coverage gain here — but it is not only a coverage gain. `integrate`
  reads each family's `DESC` out of `db.ncbifam_hmm` itself and caches it,
  because `hmmsearch --tblout` records the description of the *target* protein
  and never of the query HMM. A protein whose NCBIfam families **all** describe
  nothing — DUF, UPF, hypothetical, uncharacterised — is marked
  `ncbifam_uninformative` and does not leave `4_dark` on that evidence, the
  same test the Pfam DUF rule applies. Set `ncbifam_uninformative_test: false`
  to count every hit as annotation. The accession is kept alongside the family
  name as `ncbifam_accs`, which is what makes the call comparable with
  InterProScan's NCBIfam member database.
- **`kofam`** — KOfamScan. A *control*, not just coverage. eggNOG assigns KOs
  by DIAMOND search; KOfam uses per-family HMMs with adaptive thresholds. The
  run reports how many proteins KOfam rescues from the KO-less bins and how
  many the two disagree on, and `ko_source` records the provenance per protein.
  If KOfam rescues a large slice, those bins were partly an artefact of
  eggNOG's search rather than a fact about the proteins — quote that number
  next to any claim about the size of the non-KEGG fraction.
- **`interpro`** — InterProScan. Gene3D and SUPERFAMILY are structure-derived
  and catch what Pfam misses. Purely structural analyses (MobiDBLite, Coils,
  TMHMM, Phobius, SignalP) are excluded from the "informative" test so they
  cannot promote a protein out of the dark bin on their own.
- **`hhblits` / `jackhmmer`** — profile-profile and iterative profile search,
  run on the unannotated bins only. This is the real answer to "more sensitive
  than BLAST"; a hit creates the `3p_profile_only` bin. BLASTp is not included.
  Note what that leaves: with the default configuration the DIAMOND stage
  searches only the targeted databases listed under `db.diamond` (VFDB, MEROPS,
  CARD, TADB, BAGEL) with `--max-target-seqs 5`, so it is not a general
  homology search. `db.diamond` is free-form — add a tag and the stage searches
  it, and `integrate` picks it up by globbing `results/diamond/*.tsv` with no
  code change. A tag with no matching `diamond_weights` entry still counts as
  annotation, so it can lift a protein out of `4_dark` while contributing
  nothing to the export score, and the run warns when that happens. The only general-reference search here is
  `jackhmmer` against UniRef50, which is off by default and must be enabled
  explicitly.
- **`context`** — genomic neighbourhood. Needs a `gff` whose identifiers match
  `proteins_faa` exactly; without it the stage has nothing to work from.
  `context_window` counts neighbours on **each side, in genes**, while
  `immunity_max_gap` is in **base pairs** — different units, adjacent keys.
- **`smorf`** — smORFinder and Macrel on the assembly, so it needs
  `contigs_fna`. Note the direction: this produces ORFs that must be **added
  to the search database and the MS data re-searched**. Nothing on the annotation side recovers peptides that
  were never in the search space.

Foldseek now searches several targets (`foldseek_extra_targets` — PDB and
Swiss-Prot carry far better annotation than mostly-unreviewed AFDB50) and
clusters the unannotated structures **against each other**. Fifty dark
proteins sharing a fold is a much stronger signal than fifty singletons, and
needs no reference database.

Both halves run on a **pLDDT-gated subset** of those structures.
`thresholds.esmfold_min_plddt` (default 70) is applied *before* the search: a
45-pLDDT model of a short dark ORF matching a fold at TM 0.5 is noise, and a
hit here is what promotes a protein out of `4_dark` into `3s_structure_only`,
so low-confidence models are never searched rather than filtered afterwards.
Survivors are linked into `{results}/foldseek/query_hq` and the log says how
many of how many passed. On the real UC run that was 993 of 1,821.

The search asks for `qtmscore` and `qlen`, not only `alntmscore`. That matters:
`alntmscore` is normalised by the *alignment*, so a 40-residue local match
inside a 300-residue query can score 0.6 while the two proteins share almost no
fold. A Foldseek build too old to report the query-normalised score falls back
to the legacy columns and says so — the gate then loosens, and the run tells
you it has.

The context stage also detects polysaccharide utilisation loci: several
CAZymes plus a SusC/SusD-like importer in one neighbourhood. That is what
dbCAN-PUL looks for, computed from evidence already in hand rather than adding
a dependency.

Each records a signature over its inputs and the config keys it depends on, so
a rerun skips what has not changed. Outputs produced elsewhere are **adopted**,
not recomputed — fold on the GPU box, rsync `results/structures/` back, rerun.
`--only`, `--from`, `--force`, `--no-adopt`, `--dry-run` control the rest.

`--force` discards **only the state records of the stages it is running**: with
`--only` or `--from` that is the selection, and every unselected stage keeps its
record and is still signature-checked next time. So `--force --only <stage>` is
the supported way to redo one stage — use it instead of deleting outputs, which
triggers the same recomputation with nothing written down. A bare `--force`
still discards everything, which on a real dataset is days of compute.

## Quantification input

| `quant_format` | file |
|---|---|
| `diann` | DIA-NN `report.pg_matrix.tsv` |
| `fragpipe` | `combined_protein.tsv` |
| `fragpipe_peptide` | `combined_peptide.tsv` |
| `fragpipe_ion` | `combined_ion.tsv` |
| `msstats_csv` | FragPipe **label-free** `MSstats.csv` |
| `msstats_feature` | `dataProcess()$FeatureLevelData` as TSV |
| `msstats_protein` | `dataProcess()$ProteinLevelData` as TSV |
| `fragpipe_tmt` | FragPipe **TMT**: the per-plex `TMTn/` directories (`quant_table` is the run directory, not a file) |

Every route above except `fragpipe_tmt` is **label-free (MS1) quantification**,
and everything downstream — the roll-up, the taxon sums, the size factor, the
log2 step — works on **linear** intensity. Most routes supply it directly. The
two `dataProcess()` routes need a word: `msstats_protein` reads
`LogIntensities` (or `ABUNDANCE`), which is log2, and de-logs it silently —
MSstats' own normalisation is baked into those values and survives the
transform, so your numbers are normalised whether or not you wanted that.
`msstats_feature` prefers `INTENSITY` when the file has it; that column is raw,
linear and *not* normalised, so `dataProcess()`'s normalisation is absent from
the run. Only when `INTENSITY` is missing does it fall back to `ABUNDANCE`,
de-log it and warn. Which column your export carries therefore decides whether
MSstats' normalisation is in your numbers, and there is no key to override the
choice. Reporter intensities are linear too, which is why the per-plex tables
and not `tmt-report/` are what `fragpipe_tmt` reads.

FragPipe writes `0` for "not quantified", not for "measured as zero", so zeros
in a `fragpipe`/`fragpipe_peptide`/`fragpipe_ion` table are read as **missing**
by default and the run warns with the cell count. Summed as real zeros they turn
missingness into fold change. `zero_intensity_is_missing: false` restores the
old behaviour if you need to reproduce someone else's numbers.

### FragPipe TMT (isobaric): `fragpipe_tmt`

Isobaric runs are read by exactly one route, and it is the **per-plex** one:
each plex's own `ion.tsv` or `peptide.tsv`, mapped to samples through that
plex's own annotation file and joined here. `quant_table` is then the FragPipe
run directory — the one holding `TMT1/`, `TMT2/`, … — not a file.

- **Read:** `TMTn/ion.tsv` or `TMTn/peptide.tsv`, `TMTn/<PLEX>_annotation.txt`,
  and `TMTn/psm.tsv` when `tmt.min_purity` is set.
- **Not read, and refused by name under every `quant_format`:** the eight
  `tmt-report/` matrices, the TMT flavour of `MSstats.csv`, and the per-plex
  `protein.tsv`. Those are the three files a user reaches for first, and each
  is a worse input than it looks; *The TMT files that are still refused*, at
  the end of this section, says why for each. Nothing here falls back to any of
  them — a missing per-plex table is a failure, not a reason to read a
  different file.

The configuration, in full:

```yaml
quant_format: fragpipe_tmt
quant_table: /data/run2              # the directory holding TMT1/ .. TMT8/
tmt:
  plex_glob: "TMT*"                  # matched case-sensitively on every OS
  level: ion                         # ion.tsv | peptide.tsv
  annotation: "{plex}_annotation.txt"   # or a {plex: path} map
  reference_name: "Pool*"            # or reference_channel: "131C"
  use_reference_ratios: false        # false = the covariate treatment
  condition_from_name: auto          # auto | "" | a regex with one group
  within_plex_normalise: median      # median | none
  min_plexes: 1                      # FEATURE level, before the roll-up
  drop_empty_channels: true
  min_purity: 0                      # 0 = off; needs psm.tsv
analysis:
  min_valid_per_group: 3             # counts SAMPLES
  min_plexes: 1                      # counts PLEXES, protein level
```

Every key of the `tmt:` block, with its default:

| key | default | what it does |
|---|---|---|
| `tmt.plex_glob` | `"TMT*"` | which subdirectories of `quant_table` are plexes. Matched case-sensitively on every OS, so it cannot also pick up `tmt-report/`. Sorted naturally: `TMT10` after `TMT9`. |
| `tmt.level` | `"ion"` | `ion` reads `ion.tsv` (sequence + modified sequence + charge); `peptide` reads `peptide.tsv` (one row per sequence). |
| `tmt.annotation` | `"{plex}_annotation.txt"` | the annotation file of each plex, relative to the plex directory. A `{plex: path}` map is also accepted, and a plex missing from that map is an error rather than a fallback to the pattern. |
| `tmt.reference_name` | `""` (none) | glob on the annotated **sample name** that marks the reference/bridge channel, e.g. `"Pool*"`. |
| `tmt.reference_channel` | `""` (none) | glob on the **channel** instead, e.g. `"131C"`. Setting both is refused: they can disagree per plex. |
| `tmt.use_reference_ratios` | `false` | `false` is the covariate treatment (drop the reference, keep `plex` in the model); `true` divides every channel of a plex by that plex's reference. Refused with no reference named. |
| `tmt.condition_from_name` | `"auto"` | where the condition may come from: `auto` accepts an unambiguous split of the sample names, `""` never derives, anything else is a regex with exactly one capture group. Never the plex. |
| `tmt.within_plex_normalise` | `"median"` | `median` centres each channel of a plex on the plex's median channel before the roll-up; `none` keeps FragPipe's numbers and warns. |
| `tmt.min_plexes` | `1` | keep only **features** identified in at least this many plexes, in the reader, before the roll-up. |
| `tmt.drop_empty_channels` | `true` | drop channels the annotation names `<PLEX>_<CHANNEL>`, which is how FragPipe writes an unassigned one. |
| `tmt.min_purity` | `0` (off) | drop features whose **median** PSM purity is below this. Reads `psm.tsv`, the only table that has purity at all. |

Two keys outside that block behave differently for `fragpipe_tmt`:

| key | default | for `fragpipe_tmt` |
|---|---|---|
| `analysis.min_plexes` | `1` | the **protein**-level companion to `tmt.min_plexes`, applied in the report beside `min_valid_per_group`. Inert without a per-sample plex, so label-free is untouched. |
| `analysis.design_formula` / `analysis.factor_cols` | `"~ 0 + group"` / `"group"` | with two or more plexes the **defaults** become `"~ 0 + group + plex"` and `"group,plex"`. Only the literal defaults are replaced; a formula you wrote is left exactly as written — and `plex` is then added to `factor_cols` only if your formula actually models it, since naming a factor the model never uses stops the knit. |

Each plex is read on its own and the plexes are joined on the feature id — the
peptide sequence, the modified sequence and the charge at `level: ion`, the
peptide at `level: peptide` — so the id is comparable across plexes. **A
feature not identified in a plex is `NA` for every sample of that plex, never
`0`**: FragPipe's `0` is a real value here (5-16% of reporter cells in a real
run) and is itself read as missing, so filling one in for "not identified"
would turn plex-shaped missingness into fold change. Cross-plex overlap is
low — about 45% of ion keys are shared between two plexes — so the run always
logs how many features were seen in 1, 2, … n plexes; `min_plexes` filters on
that count.

Reporter columns are named `Intensity <sample>` after the **annotated sample
name**, which is what each plex's `<PLEX>_annotation.txt` supplies (FragPipe
does not write a plain `annotation.txt`). Two things in that file are easy to
get wrong and are handled explicitly:

- **The reference channel does not sit at a fixed position.** In a real 8-plex
  design the pool is at `131C` in six plexes and at `131N` in the other two, so
  `reference_channel` alone cannot describe the run. `reference_name` globs the
  sample name (`Pool*`), which is the stable signal; set one or the other, and
  the log says which was used and what it resolved to per plex.
- **An unassigned channel is named `<PLEX>_<CHANNEL>`** (`TMT7_131C`) in the
  annotation. Its signal is isotope carry-over, not a sample, so it is dropped
  by default and the log names it. `drop_empty_channels: false` keeps it.

#### Which designs are read

Every plex is described by its own annotation, and nothing assumes a common
channel count, a common channel set, or a bridge, so all of these read:

- **Reference-free** — no pool anywhere. Nothing is divided by anything and no
  channel is held back as a denominator; the `plex` term carries the batch. This is what you get
  when neither `tmt.reference_name` nor `tmt.reference_channel` is set. A lone
  channel named `Pool*` is pointed out in the log rather than being treated as
  a reference behind your back.
- **Multi-plex without a bridge** — several plexes with no channel in common.
  Read, and honest as long as each condition appears in more than one plex.
  What links the plexes is then the plex coefficient and the report's median
  normalisation: a per-plex mean shift, not the per-protein correction a shared
  channel would give. A bridge is the better design; its absence is not a
  reason to refuse the data.
- **A bridge in every plex** — name it (`reference_name: "Pool*"`) and pick the
  covariate or the ratio treatment below. Under both it stops being a sample.
- **Mixed plex sizes** — 16 channels in one plex, 11 in the next, 6 in the
  third read as one experiment; that exact run was read to check it. Sample
  names must be unique across the whole run, since they are the columns of the
  joined matrix, and two plexes claiming one name is refused.
- **TMTpro 16- and 18-plex, and iTRAQ** — channel labels are strings out of the
  annotation file and no channel set is hard-coded anywhere, so a TMTpro
  annotation reads exactly as a TMT-11 one does. Only TMT-11 has been run on
  real data (8 plexes, 88 channels); a 16-channel annotation was read to check
  that nothing counts channels, and nothing wider than 16 has been tried.
- **A single plex** — read. `plex` is then not added to `design_formula`, because
  one level is not a batch effect, and the cross-plex filters have nothing to do.

#### The designs whose statistics cannot be made honest

These are not supported-with-a-caveat. Two of them stop the run:

- **A condition that does not cross plexes** — one condition per plex. **The run
  stops**, in the reader when the condition came from the sample names and again
  in the report's `plex_confounding()` before the model is fitted, with the
  cross-tabulation printed. The batch and the biology are the same vector, so
  every fold change would be both; dropping `+ plex` does not fix it, it reports
  the batch as biology. There is no flag to override this and there should not be.
- **A model that is rank deficient once `plex` is in it** — the same problem
  arriving through a hand-written formula. Refused with the non-estimable
  coefficients named, and the message says that in a TMT run the plex is the
  batch a condition has to cross.
- **A condition that crosses only some plexes.** This one fits. The report
  prints the condition-by-plex table and gates on the empty cells, because the
  contrast is then carried by whichever plexes hold both levels, and a level
  present in one plex only is estimated from that plex's batch as much as from
  its biology. Nothing in the numbers separates the two.
- **A plex holding a single sample.** Its plex coefficient fits that one channel
  exactly, so the channel contributes nothing to the condition. The report gates
  it as a factor level with one sample. Adding such a plex adds no power.
- **A protein quantified in one plex only.** Its fold change is that plex's batch
  as much as the condition, and `min_valid_per_group` cannot see it because it
  counts samples. That is what the two `min_plexes` keys are for, and the report
  states the exposure whether or not either is set.
- **Effect sizes compared with label-free ones.** Reporter ratios are compressed
  toward 1 by co-isolation, metaannot does not correct that, and the compression
  is worst exactly where this database is weakest — near-identical paralogues in
  one isolation window. Direction and ranking survive the compression; magnitude
  does not. A TMT log2 fold change and a DIA one are not the same quantity and
  must not be pooled or compared.

#### The reference channel is not a sample

A pooled bridge channel is not a biological sample and never enters the design
as one. Which of the two treatments is used is yours to choose, and both are
written into `design_record.txt` beside the numbers they produced:

- **covariate** (the default, `use_reference_ratios: false`). The reference is
  dropped from the sample columns; the plex stays in the model as a batch term
  and absorbs the plex effect. Nothing is divided, so no value is lost.
- **ratios** (`use_reference_ratios: true`). Every channel of a plex is divided
  by that plex's reference and the reference column is dropped. This is the
  classic bridge design and it removes the plex effect directly, but it assumes
  the same pool went into every plex, it discards the reference's own variance,
  and it **propagates the reference's missingness**: a feature with no
  reference value in a plex becomes `NA` for every channel of that plex — and
  a reference FragPipe wrote as `0` counts as no reference whatever
  `zero_intensity_is_missing` says, since dividing by it is not a ratio. The
  run reports what that cost, per plex and in total — on the real 8-plex run,
  6,576 of 1,237,468 values, 0.53% — and warns when it is large.

The result of either is still linear, so log2, the roll-up and the
median-of-ratios size factor are unchanged.

#### The condition, and the plex in the model

The **condition is not in these files and is never inferred from the plex** —
a plex is a batch, and a plex-versus-plex contrast is a batch effect presented
as a hypothesis. `design_from_input.tsv` carries `sample`, `plex` and
`channel`, plus `group` when a condition could be had honestly. A `manifest` is
ignored for TMT and says so: a FragPipe TMT manifest names LC-MS runs, and its
experiment column is the plex.

Two places the condition can come from:

- **The annotated sample names**, when they carry it unambiguously.
  `condition_from_name: auto` (the default) accepts the part of the name before
  a separator only when every name has one, every level has at least two
  samples, and no other separator groups the samples differently; it logs the
  result as the guess it is. `condition_from_name: ""` never derives; a regular
  expression with one capture group states the rule explicitly.
- **`analysis.metadata`**, keyed on sample. This is required whenever the names
  do not carry the condition — the real dataset's `MF####` codes do not — and
  the run then stops and prints the exact file to write, starting from
  `design_from_input.tsv`, which already lists every sample and its plex.

For TMT the default model is **`~ 0 + group + plex`** and `factor_cols` is
`group,plex`: the condition is the hypothesis and the plex is a nuisance term
fitted alongside it, so contrasts come out over the condition. A
`design_formula` you write yourself is left exactly as written, and a
single-plex run keeps `~ 0 + group` because one level is not a batch effect.
**A plex perfectly confounded with the condition stops the run**, naming the
plex and printing the cross-tabulation: with one condition per plex the batch
and the biology are the same vector and no model separates them.

#### Within-plex normalisation

The channels of one plex are the same LC-MS run, so what differs between them
is how much peptide was loaded and how completely it was labelled: a
per-channel constant with no biology in it, which the roll-up would otherwise
sum straight into the protein. `within_plex_normalise: median` (**the
default**) divides each channel by its own median and multiplies by the plex's
median channel, and logs the log2 scale factors it applied. Because the median
commutes with log2, that is exactly a per-channel median centring, applied one
plex at a time and before the roll-up rather than after it; because it centres
on the plex's own median rather than on 1, the values stay linear, which is
what the size factors need.

The **between-plex** difference is deliberately left alone. That one is the
batch, and the plex term in the model — or the report's own `normalise:
median` — is what removes it; taking it out here would hide it from both.
`within_plex_normalise: none` keeps FragPipe's numbers exactly and says so.
Either way the choice is written into `design_notes.txt` and from there into
`design_record.txt`, because a matrix that has been median-centred per channel
and one that has not are different data and nothing downstream can tell them
apart by looking.

One case this handles worst is named in the log rather than hidden: the median
is taken over **observed** values, so a channel whose low end went missing has
a median above its true centre and is scaled up too little. A channel that is
both far off the plex scale and much emptier than its neighbours is therefore
reported as under-corrected — on the real 8-plex run that is three channels of
eighty-eight, one of them 10.7x off with 42% of its cells missing, which is a
loading failure rather than an imbalance. An unequal but complete load is
exactly what this step is for and raises no alarm.

#### `min_plexes`: what `min_valid_per_group` cannot see

`min_valid_per_group` counts **samples**. An isobaric run's missingness is
shaped by the plex: a protein identified in one plex only is all-`NA` in every
other, so "3 valid values in every group" can be satisfied entirely inside one
batch, and the difference the model then reports is that batch. `min_plexes`
counts **plexes** instead, and both filters are applied:

- `tmt.min_plexes` filters **features**, in the reader, before the roll-up.
- `analysis.min_plexes` filters **protein groups**, in the report, beside
  `min_valid_per_group`.

The report counts each filter against the same starting set and prints them
separately, so it is visible which one bit, along with the histogram of
proteins by number of plexes. Both default to 1, so no run loses proteins to a
filter it did not ask for — and with the protein-level filter at 1 the report
still says how many of the proteins it kept are quantified in a single plex,
which is the number to set it on, **and what setting it would leave**.

Above 1 the report says the other half, because at that point the setting is
what decided which proteins the document is about: how many protein groups
that is, out of how many `min_valid_per_group` passed, and that at the default
every one of the rest would be in the model. It is one comparison — the value
in force against the default — and not a ladder, because the value in force is
the only one the document can be exact about and 1 is not a point on a ladder
but the absence of a choice. The sentence is silent where the filter took
nothing `min_valid_per_group` would have kept, since the line above it already
carries that as its `0 of which min_valid_per_group would have kept` clause —
the removal count on that line can be any number, and it is the clause rather
than the count that says the filter took nothing worth pricing.

**That sentence is not the histogram cumulated, and the histogram cannot be
made into it.** `n_plex` is counted over every quantified row, *before*
`min_valid_per_group`, and the retained set is both filters at once, so a
protein quantified in every plex can still be one `min_valid_per_group` drops.
The caption says so where the table is printed.

The counterfactual is only ever offered for `analysis.min_plexes`.
`tmt.min_plexes` filtered features in the reader, before the roll-up, so every
count in the report is already after it and what it cost in protein groups is
not derivable from the document at all — only from a re-run. The report says
that, in those words, on the runs where the feature-level filter was in force
AND a protein-level sentence fired: the caveat is attached to that sentence
rather than standing alone, so a run where `analysis.min_plexes` had nothing
to say prints neither. That is deliberate — a caveat about a number the
document did not print is a sentence with no referent — but it does mean the
feature-level filter can be in force with the document silent about it, and
the run log is where its cost in features is recorded either way. `analysis.min_plexes` is inert without a plex column, so label-free
runs are unaffected; setting it above 1 where no per-sample plex exists stops
the report rather than passing everything.

#### `min_purity`, and where purity actually lives

Precursor purity is written **only into `psm.tsv`** — not into `ion.tsv`,
`peptide.tsv` or `protein.tsv`. `tmt.min_purity` is therefore implemented as a
join: psm.tsv is read per plex, keyed on the same columns the feature id was
built from, and a feature is judged by the **median** purity of the PSMs that
produced it, because its reporter intensities are a sum over those PSMs and no
single one describes it. `0` (the default) disables the filter and psm.tsv is
not read at all — and is not a stage input either, so a rewritten psm.tsv does
not invalidate a cached join that never opened it.

Two consequences are logged rather than assumed. A feature that matches no PSM
row is **kept** and counted: an unmatched key is a join failure (FragPipe
leaves `Modified Peptide` empty on rows `ion.tsv` writes a modified sequence
for; about 1.5% of ion keys in a real plex), and dropping those would look
exactly like a purity filter working. And this is an approximation of the
per-PSM filter TMT-Integrator would apply *before* summarising — FragPipe has
already summed by the time this reader sees the file.

`min_purity` limits how co-isolated the accepted spectra were. It does not
correct the ratio compression that co-isolation causes, and nothing here does:
the report states that as a limitation (fold changes are lower bounds, ranking
is more trustworthy than magnitude) rather than dividing by an estimate of the
contamination and turning a known bias into an unknown variance.

#### The taxon size factor under a plex effect

The taxon size factor is a median of ratios computed across **all** samples, so
an isobaric run has to be asked whether the batch got into it. Measured on a
fixture with no biology at all — the same random draw run twice, once with
equal plex loading and once with a 7.5x spread between plexes — the plex effect
lands entirely in the part every taxon shares (recovered 1.575 log2 against a
true 1.585, and −1.289 against −1.322), which is exactly what a size factor is
for, and the taxon-by-taxon part is **identical to 9.4e-05 log2**, four orders
of magnitude below the effect. So a plex effect does not reach a taxon's ratio
model.

What does reach it is plex-shaped **missingness**. With fewer than
`taxon_min_proteins_for_factor` members observed in every plex, a taxon loses
the complete-case reference and falls back to the poscounts variant, whose
median then mixes proteins referenced inside one plex with proteins referenced
across all of them; in the same fixture a taxon with 7 of 10 members confined
to one plex had no size factor at all in the other two and was displaced 1.40
log2 in the one it lived in, while the unaffected taxa stayed within 0.03. The
join stage counts those taxa and warns, naming `analysis.min_plexes` and
`tmt.min_plexes` as the two ways to remove them.

#### The TMT files that are still refused

A TMT run also writes files this tool must not quantify, and each is refused by
name rather than read:

- The `tmt-report/` matrices are **refused** on their `ReferenceIntensity`
  column, under every `quant_format`. All eight of them: `abundance_` and
  `ratio_`, at `gene`, `protein`, `peptide` and `modified-peptide` level — the
  peptide-level ones carry `Peptide` and `Mapped Proteins` and so look more
  like readable feature input than the protein ones do, and they are refused
  just the same. They are already log2 and median-centred, so the log2 step
  downstream would log them twice; they are already rolled up, so the
  shared-peptide rule, `peptide_assignment` and `peptide_evidence.tsv` have
  nothing to work on; and their inference is TMT-Integrator's, which is the
  inference a strain-redundant metagenome database makes least trustworthy.
  `fragpipe_tmt` never falls back to them — a missing per-plex table is a
  failure, not a reason to read a different file.
- The TMT flavour of `msstats.csv` is **refused** on its `Channel <mass>`
  columns, recognised from the header before pandas parses the body (unquoted
  commas in `Protein.Description` used to kill it with a raw tokenising error).
- A per-plex `ion.tsv`/`peptide.tsv`/`psm.tsv` handed to `fragpipe_peptide` or
  `fragpipe_ion` is **still refused**, naming the reporter channels it found,
  and now points at `fragpipe_tmt`. Those formats keep only columns *ending* in
  `Intensity`, so the run would otherwise reduce to the single MS1 precursor
  intensity, pooled over all channels, and report one "sample". Do not work
  around that refusal.
- A per-plex `protein.tsv` handed to `fragpipe` (the protein-level format) is
  **refused** on its `Intensity <sample>` columns. It is the one isobaric file
  that carries neither marker above — no `ReferenceIntensity`, no
  `Channel <mass>` — and the protein-level column detector takes every numeric
  column that is not declared metadata, so it would have quantified a SINGLE
  plex as the whole experiment with `Length`, `Protein Qvalue` and
  `Razor Intensity` in the matrix beside the channels. The discriminator is the
  prefix: FragPipe writes `Intensity <sample>` for a reporter channel and
  `<sample> Intensity` for a label-free run, so a label-free
  `combined_protein.tsv` is untouched.
- `fragpipe_tmt` pointed at a single FILE — any of the four per-plex tables, or
  a `tmt-report/` matrix — is refused naming the file, because `quant_table` is
  the run **directory** for this format. Pointed at one plex directory it
  refuses too, naming the `tmt.plex_glob` that matched nothing and listing what
  is there instead.

`analysis.msstats_comparison` additionally accepts a
`groupComparison()$ComparisonResult` export, which replaces the abundance model
in the report.

### The shared-peptide rule

Shared peptides are the central problem in a strain-redundant metagenome
database, so `peptide_assignment` is explicit:

- `protein_unique` — only features matching exactly one protein.
- `taxon_unique` (default) — features whose candidates all resolve to one
  taxon, assigned to the razor protein. A peptide shared between two strains
  of one organism still quantifies that organism; one shared across taxa
  quantifies neither. A candidate with **no** taxonomy makes a feature
  `shared_unknown_taxon`, never taxon-unique.
- `taxon_or_family_unique` — as `taxon_unique`, plus features whose candidates
  cannot be compared by taxon at all because at least one has none, but which
  all share one MMseqs `family_id`. A family is a sequence cluster at
  `cluster_min_seq_id`, **not** an organism, so the intensity is attributed to
  one representative of the cluster. Opt-in only: it needs `run.cluster` on —
  metaannot refuses the mode if `family_id` is missing, and refuses again if no
  family has more than one member, since with the cluster stage off every
  family is a singleton and the rule could never fire. Those features are
  recorded as their own `family_unique` class and the run logs the count, as a
  WARN when any feature was kept that way and INFO when none. It does not
  rescue features whose candidates resolve to two *different* taxa; those stay
  `shared` and are dropped.
- `razor` — FragPipe's own behaviour. Arbitrary here; useful for measuring how
  much it changes the answer.

#### `quant_funnel.tsv` — where the proteins went

The first full real run narrowed from 455,571 proteins to 1,282 rows, and said
why across a dozen lines of an 8,990-line log, in three stages, some of them
counting features and some counting proteins. The join stage now writes those
steps down as it takes them, beside the table they produce:

| step | unit | before | after | dropped |
|---|---|---|---|---|
| proteins in the annotation table | protein | | 455,571 | |
| proteins named by at least one quantified feature | protein | 455,571 | … | … |
| features read from the quant table | feature | | … | |
| `peptide_assignment=taxon_unique` | feature | … | … | … |
| proteins carrying at least one assigned feature | protein | … | … | … |
| `min_features_per_protein=1` | protein | … | … | 0 |
| rows written to `annotated_quant.tsv` | protein | … | 1,282 | … |

**The `unit` column is the load-bearing one.** A funnel that chained a feature
count straight into a protein count would read as one number shrinking while
being a different claim at every step. So a row's `before` is the last `after`
recorded *for its own unit* — the protein chain steps over the feature rows in
the middle of the stage rather than restarting after them — and the first row
of each unit has no `before` at all. A `why` column carries the reason in
words; a filter that removes nothing is still listed, priced at zero, because a
funnel silent about a filter reads as a funnel with no such filter.

A **negative** `dropped` is not an arithmetic bug: it means a step's
population is not a subset of the one above it, which here has exactly one
cause — protein ids the quant table names that `annotation_final.tsv` does
not. A run where that is true of most ids dies; a run where it is true of some
warns and carries on, and this is where the consequence becomes visible.

Two things it does not cover, and says so rather than leaving them to be
inferred. Rows the quant reader refused before this stage was handed anything
— decoys, contaminants, unusable columns — are already gone from its first
count, which says so; the reader logs them. And every filter the **report**
applies afterwards (`min_valid_per_group`, `analysis.min_plexes`) runs in R
over the file this funnel ends at, and is counted there. The join stage says
that out loud on the log beside the path, because a funnel that stopped at
1,282 without saying it would be read as the end of the narrowing when it is
the middle of it.

`peptide_evidence.tsv` has one row per protein and nine columns:
`protein_id`, `n_features_used`, `n_unique`, `n_taxon_unique`,
`n_family_unique`, `n_features_dropped`, `taxon_unique_dominated`,
`rollup_method` and `peptide_assignment`. The last two record how the numbers
were made, because a `sum` run and a `median_polish` run are otherwise
indistinguishable once the log is gone. `taxon_unique_dominated` flags proteins
resting more on shared-but-taxon-unique features than on their own unique ones
— the ones whose intensity is most sensitive to the assignment rule you chose.

The file has a row for every razor protein the rule considered, **including
those whose features were all dropped**. Those rows read `n_features_used = 0`,
carry no number in `annotated_quant.tsv` and appear in no report table, and
they can never be `taxon_unique_dominated` either, because a row of zeros
cannot rest more on one kind of feature than on another. So the rate the join
stage logs is over the proteins with at least one assigned feature — the ones
the rule actually decided something about — and the run prints how many rows it
left out, so a rate computed straight off this file will not match it. Where
`min_features_per_protein` is set above its default of 1 the retention line
beside it reports over that same population, which is the comparison to make;
at the default that line does not print at all, because the filter is inert.

The report states the rate again over ITS protein set, which is the one a claim
about your results rests on. Whether that number differs from the stage's
depends entirely on your config: `min_features_per_protein` defaults to 1 and
`analysis.min_features` to 0, and both filters are skipped when they are at
those values, so on a default run the two populations are the same and the two
rates agree. Raise either and they part company — on the 455,571-protein run,
at `min_features_per_protein: 2`, the stage logged 54.7% and the report 61.2%,
because the single-feature proteins that filter removed were less often
dominated than the ones it kept. That direction is a property of that dataset
and not a rule.

MSstats-format inputs carry only the razor protein per feature, so
shared-peptide filtering is unavailable from them; metaannot says so.

`ProteinLevelData` inherits whatever `dataProcess()` did. Check `MBimpute` and
`censoredInt`: imputed low-abundance values are exactly where the KO-less
fraction sits.

## Taxonomy: Unipept vs eggNOG

eggNOG's `seed_ortholog` taxid is the taxon of the best-matching *reference*
protein, not of the organism in your sample. Unipept computes a peptide LCA
against UniProt. Independent estimates, different failure modes — and both the
`taxon_unique` roll-up and the ratio model stand on whichever you use.

```bash
npm install -g unipept-cli     # Node 22+; the Ruby gem is the legacy client
unipept pept2lca --equate --all -i results/unipept/peptides.txt -o pept2lca.csv
```

The current CLI writes `domain_id`/`domain_name` where the old gem wrote
`superkingdom_id`/`superkingdom_name`, which is what the parser reads. The
unipept.ugent.be web export is **not** an accepted input: it is name-based and
lacks the taxid columns entirely.

Set these inside the blocks your config already has — do not paste a second
top-level `run:` or `db:` key (see above; a duplicate key is refused):

```yaml
run:                              # in the existing run: block
  unipept: true
  taxonomy: true
unipept: {result: "pept2lca.csv", split_missed_cleavages: true}
db:                               # in the existing db: block
  ncbi_taxonomy: "/data/db/taxdump"          # nodes/names/merged/delnodes.dmp
taxonomy_source: "concordant"                # eggnog | unipept | concordant
```

Per protein the consensus is the deepest rank where a majority of its peptides
agree — not the LCA of the LCAs, which one spurious peptide would drag to the
root. `concordant` keeps a taxon only where **both** methods produced a taxid
and either the two taxids are the same (verdict `identical`, which is bare
taxid equality — no rank test is applied to it) or their lineages agree down to
genus or species (`concordant`). Everything else is blanked, not only the
disagreements: `unipept_missing`, `eggnog_missing`, a seed taxid that resolves
to nothing even after `merged.dmp` (`eggnog_unresolved` — a taxid that merely
*merged* is followed to its current node and compared normally),
`concordant_above_genus`, `no_common_rank`, no peptide LCA, or fewer than
`consensus_min_peptides` peptides. So the taxon-based steps can cover far fewer
proteins than you expect.

`db.ncbi_taxonomy` is not optional for this, even though it defaults to empty.
Without a taxdump there is no lineage to compare and every verdict collapses to
`identical` or `differ_no_lineage`, so `concordant` degrades to exact taxid
identity and blanks everything else — including the cases it exists to keep,
such as eggNOG *E. faecalis* against a Unipept LCA of genus *Enterococcus*. The
run only WARNs. Check the verdict counts in `taxonomy_comparison.tsv`.

## Databases

```bash
hmmpress /data/db/Pfam-A.hmm ; hmmpress /data/db/dbCAN-HMMdb-V13.txt
diamond makedb --in VFDB_setA_pro.fas -d /data/db/vfdb_core     # and MEROPS, CARD, TADB, BAGEL
foldseek databases Alphafold/UniProt50 /data/db/foldseek/afdb50 tmp
```

`db.diamond` and `sources.diamond` **replace** the built-in defaults rather
than merging into them. Whatever you write under `db.diamond` is the complete
list for the run, and the run names what you dropped:

```
WARN config: db.diamond lists ['vfdb'], so the default entries
     ['bagel', 'card', 'merops', 'tadb'] are NOT used
```

So to drop one database, **list the ones you keep**. Setting an entry to `""`
or `null` also drops it, silently, before any stage or `doctor` sees it — but
only within the block you supply. Writing a block whose only key is
`merops: ""` therefore leaves **no** DIAMOND databases at all: merops is
dropped as falsy and the other four were never merged back in. The run does not
fail; the stage logs `no diamond databases configured, nothing to do` and
`annotation_final.tsv` comes out with no virulence, protease, AMR,
toxin-antitoxin or bacteriocin evidence — indistinguishable from a real
absence. A database dropped either way does not appear in `doctor` at all; it
is **not** reported as `MANUAL: no path configured`.

Give each database you do keep an entry in `diamond_weights` or it scores 0 and
the run warns.

### The identity a hit has to reach

Not every DIAMOND database is read the same way. A MEROPS or TADB hit is a
family assignment — "this looks like a protease" — and 30% identity over good
coverage supports that. A CARD, VFDB or BAGEL hit is read as a claim about a
**particular** protein: this confers resistance, this is a virulence factor.
A 32%-identity match to a beta-lactamase over half a query is a hit against
the fold, and reporting it as "carries an AMR gene" is the claim a reviewer
would check first.

```yaml
thresholds:
  diamond_min_pident: 30     # merops, tadb, and anything else
diamond_min_pidents:
  card: 50                   # defaults; override or empty to change them
  vfdb: 50
  bagel: 50
```

On the 455,571-protein run this is the difference between 39,138 CARD hits and
4,661, and between 107,219 VFDB hits and 21,623.

`diamond_min_pidents` mirrors `diamond_evalues` and is applied **twice on
purpose**: as DIAMOND's own `--id` during the search, so a floored database
writes thousands of rows instead of hundreds of thousands, and again when the
table is read, because a `<tag>.tsv` adopted from another machine or written
before the floor existed never saw `--id`. It is part of the `diamond`,
`integrate` and `finalise` signatures.

One consequence worth knowing: `thresholds.diamond_strong_pident` (50) halves
the weight of a hit below it, and a database floored at or above that number
can no longer produce one. Every VFDB hit that reaches the scoring takes the
full weight, and the run says so once rather than leaving the config implying
a grading that cannot happen.

### A database that cannot hit

Three things a DIAMOND database can do that look exactly like "no virulence
factors here", and that both `doctor` and the `diamond` stage now refuse or
warn about before the search starts:

* **A file too small to be a database.** A `diamond makedb` that failed leaves
  a zero-byte `.dmnd` behind. Searching it reports 0 hits, which is
  indistinguishable in `annotation_final.tsv` from a real absence. Anything
  under 128 bytes — smaller than DIAMOND's own header — or that
  `diamond dbinfo` reports as holding no sequences is now a refusal, naming the
  `diamond makedb` line that rebuilds it. `doctor` reports it as `MISS`.

* **A motif seed set built as though it were a sequence database.** BAGEL4
  ships two different things: its bacteriocin sequence files, and the motif
  seed set its HMM/regex step is built from — headers like `LE-nisin`,
  `ggmotif`, `lasso`, median length 15 residues. The seed set is what got
  built here, and it returned exactly 0 hits against 38,204 proteins. A
  `blastp` against 15-residue seeds searches for those fifteen residues, not
  for the molecules they mark, so **no e-value makes it a bacteriocin
  search** — lowering the threshold would have produced meaningless hits
  instead of meaningless silence. A database whose typical sequence is under
  25 residues, or whose source FASTA headers carry those markers, is now
  called out as a seed set, and the e-value advice below is deliberately
  *replaced* rather than added to. (25, not 40: mature nisin is 34 residues,
  so a genuinely short bacteriocin database must not be accused of this.)

* **A database whose sequences are too short for the e-value.** Distinct from
  the above, and now the narrower case: real short peptides that still cannot
  reach the threshold set for them. The estimate is a perfect self-match of
  the database's typical sequence, so it fires only when a hit is essentially
  impossible rather than merely unlikely. The run says so with the number, and
  names the weight the database is holding while it cannot hit.

The fix for the third one is a per-database e-value:

```yaml
thresholds:
  diamond_evalue: 1e-10     # vfdb, merops, card, tadb
diamond_evalues:
  bagel: 1e-3               # short peptides cannot reach 1e-10
```

Check *what* you built before reaching for this knob. On this pipeline's own
run the 0 hits were the seed-set mistake above, not the threshold, and the
per-database e-value would have papered over it.

`diamond_evalues` overrides `thresholds.diamond_evalue` for that tag alone, in
the search *and* in the filter `integrate` applies to the hit table, and the
run logs each database that is searched at a threshold other than the headline
one. The estimate behind the warning uses DIAMOND's own BLOSUM62 constants
(Lambda 0.267, K 0.041) against a perfect self-match, so it fires only when a
hit is essentially impossible, not merely unlikely. The typical sequence length
comes from `diamond dbinfo`, which reports Sequences and Letters, so the length
it yields is the **mean**. Whenever dbinfo cannot answer — diamond is not
installed, the call fails, or the `.dmnd` is too old to report those fields —
the check falls back to a source FASTA beside the database or named in
`sources.diamond`, and when neither is available it says the length is unknown
rather than guessing one. A mean hides the distribution, so a database of
mostly-long sequences with a short tail can pass this check while its short
entries remain unreachable at the configured threshold.

## The R object

```bash
python metaannot.py object          # results/metaannot.rds
python metaannot.py all             # run, report, then object
```

A **QFeatures** with two linked assays when the input was peptide- or
ion-level:

```r
qf <- readRDS("results/metaannot.rds")
qf                                     # peptides -> proteins
assay(qf, "proteins")                  # intensities, proteins x samples
rowData(qf[["proteins"]])$bin          # annotation evidence bin
rowData(qf[["proteins"]])$effector_score
rowData(qf[["peptides"]])$assignment_class  # five values, see below
colData(qf)                            # the design in quant/design_from_input.tsv
metadata(qf)$taxon_size_factors
metadata(qf)$normalisation_risk
```

`assignment_class` takes five values: `unique`, `taxon_unique`,
`family_unique`, `shared` and `shared_unknown_taxon`. `family_unique` appears
only under `peptide_assignment: taxon_or_family_unique`, where a feature whose
candidates could not be compared by taxon was kept because they all share one
MMseqs `family_id` — a sequence cluster, not an organism.

`colData` is read from `results/quant/design_from_input.tsv` and from nothing
else: the manifest for a label-free FragPipe or DIA-NN run, `Condition` and
`BioReplicate` off the long table for the MSstats formats, and `sample`, `plex`
and `channel` for `fragpipe_tmt` — plus `group` only when it could be derived,
which for TMT usually means it comes from `analysis.metadata` at report time
rather than from the input.

The link between the assays records the assignment metaannot actually made,
not one re-derived in R, so the shared-peptide decisions stay inspectable
rather than becoming implicit.

Falls back to a **SummarizedExperiment** when QFeatures is not installed or the
input was protein-level (the peptide assay, if any, goes in
`metadata()$peptides`), and to a plain list if neither Bioconductor package is
present — the data is never trapped behind a missing dependency.

If the report has already run, its differential-abundance results are folded
into `rowData` as `DE.<contrast>.<column>`, so identification, quantification,
annotation and statistics live in one object. `doctor` reports which R packages
are installed.

`export_feature_quant: false` skips the feature-level matrix if its size is a
problem; the object then carries a protein assay only.

## The report

`metaannot.py report` writes the Rmd into `results/analysis/` with its params
filled in from the config, then renders it if `Rscript` is available. Two
models are fitted: protein abundance, and abundance relative to the source
organism. The taxon reference is a median of ratios, not a sum, and it is
DESeq2's: the per-protein reference is the mean across samples computed over
the taxon's proteins observed in **every** sample, and the median is taken over
those same complete-case proteins, so the reference does not move with the
sample set. The protein being tested is not left out, so a taxon carried by
only a few complete proteins partly regresses against itself. When too few of a
taxon's proteins are complete, the poscounts variant is used instead — each
protein referenced against its own observed samples — and a sample whose median
would then rest on fewer than `taxon_min_proteins_for_factor` ratios is left
blank rather than guessed. The factor degrades to a plain sum when the taxon
has fewer than `taxon_min_proteins_for_factor` proteins at all (4 by default),
when neither variant is supported, or when poscounts leaves no usable sample;
the join log and the report each print how many taxa landed there.

Just above the threshold is not safe either — the median then rests on as few
as four ratios, so the protein being tested is a quarter of its own reference.
Treat calls resting on a taxon near the threshold as provisional, and read the
number out of your own run rather than a fixed protein count. This Python
threshold and the R usability threshold `analysis.taxon_min_proteins` are
separate keys, but the defaults deliberately hold them equal at 4 and nothing
enforces it. Keep them equal: the report's "reference is the plain sum" flag is
only evaluated over taxa the R threshold already admitted, so setting
`analysis.taxon_min_proteins` higher drops exactly the taxa that flag exists to
mark.
Proteins whose significance disappears under the second model were tracking
their organism, not being regulated.

The shortlist is not the top of `effector_score`: significant, no KO, and
predicted to reach the host. `surface_or_secreted` gates the list; the score
only orders it. That gate is the OR of a signal peptide, a beta-barrel, an
LPxTG motif and an anchor domain — so with `topology` off it narrows to the
last two rather than closing, and an empty shortlist is as likely to mean
nothing was significant as it is to mean the tool was missing. Check the
differential-abundance table before concluding either. See `docs/signalp-6.md`.

**The two KEGG-invisibility rates are now subtracted rather than left on two
pages.** The bin composition chunk states what fraction of the quantified
protein groups is invisible to KEGG pathway enrichment; the finalise stage has
always logged the same fraction over the whole search database. Both numbers
have existed for as long as there have been bins, in different places, and
nothing ever put them side by side — so on the first full real run, 43.6% of
what was quantified against 29.4% of the database it was searched against went
unremarked. That gap is the point: the KO-less fraction is not merely large,
it is **over-represented among the proteins that were actually expressed and
measured**. The report now prints the database rate on the line under the
quantified one, with the ratio between them, and the join stage logs the same
pair in a single line rather than half of it here and half of it two stages
ago. Both are scored off `bin` on both sides, because `bin == "1_ko_pathway"`
is what `kegg_enrichment_visible` is defined as, and both count only
quantified groups that HAVE an annotation row — a group with none is not a
KEGG-invisible protein, it is one the pipeline knows nothing about, and it is
named separately. A NOTE carries the caveat the ratio needs: the two
populations are selected very differently, the database being every predicted
ORF and the quantified set being what was identified and survived every filter
above, so the ratio is an observation about that run and not a general rate.

The report also says when the statistics do not reach the bins this tool exists
for, rather than leaving it in a table. `retention by bin` prints each bin's
quantified count beside the count that survives `min_valid_per_group` and
`min_plexes` and is actually fitted; a bin that arrives with proteins and
leaves with none raises a GATE naming it and saying that nothing below — no
table, no enrichment, no shortlist — is about it. The KO-less bins raise a
second GATE stating what fraction of that whole population the statistics
cover, and it is conditional in both directions: only when the run quantified
at least `COVERAGE_MIN_N` KO-less groups, so that the fraction has something
behind it, and only when the fraction is under `COVERAGE_MIN_PCT` — one per
cent — because above that the coverage is a number the reader can work with
and not a failure. The percentage is of the KO-less groups the run QUANTIFIED,
and the model it is set against is the whole model, the unbinned groups
included, since those are fitted too. A NOTE names the knobs behind it:
`analysis.min_valid_per_group` with its value, `analysis.min_plexes` where an
isobaric run applies it as well, and `min_features_per_protein` — a top-level
key, not part of any stage block — which decided which proteins were written
to `annotated_quant.tsv` before the report saw anything. The escalation
is by DENOMINATOR, so a bin that is empty because nothing was ever put in it
stays silent — `3p_profile_only` is fed only by `hhblits` and `jackhmmer` and
no run has put a protein in it, and a gate that fires on that every time is how
a reader learns to skip the line that matters. Under `COVERAGE_MIN_N`
quantified groups the same zero is a NOTE rather than a GATE, because "none of
them" over a handful of proteins is an anecdote. The ratio model and the
effector shortlist repeat it where it changes how they are read, under the same
floor and distinguishing the same two causes: a ratio model whose usable set
holds no KO-less protein says whether none reached the model at all — the
coverage failure above, and nothing to do with taxonomy — or whether the ones
that did have no usable taxon, which is a limit of the taxonomy and gated only
over a population; and an empty shortlist prints the population it was drawn
from, so "nothing was significant" and "there was nobody to rank" stop looking
the same on the page.

## Parallelism

Stages form a dependency graph and independent ones run concurrently. Almost
every stage before `integrate` consumes no other stage's output — just the
protein FASTA and its own reference database — so eggNOG-mapper, hmmsearch,
DIAMOND, InterProScan, KOfamScan, SignalP, TMbed and MMseqs2 all dispatch in
the first wave. `context` is the exception that actually waits: it reads the
annotations produced by `emapper`, `pfam`, `signalp` and `dbcan` (plus the
`gff`), so it cannot start until all four have finished. `jackhmmer`, `hhblits`
and `esmfold` then run after `integrate`, because their input is the dark set
it writes.

Concurrency is capped by `stage_workers`, so the first wave is a queue rather
than a stampede: with the default 4, twelve ready stages compete for four
slots. The scheduler orders that queue longest-first, by two numbers per
stage. The first is a coarse cost rank — hours, minutes, or seconds, measured
on real runs — and it is also the ratio the machine is divided by, which is
why it stays coarse. The second is what that stage took on the release's
reference run, and it exists to break the ties in the first: eleven stages
are hours-class, so the rank on its own left all of them in table order, and
`interpro` — the longest stage in the pipeline by an order of magnitude, and
the stage the ordering was written for — was dispatched seventh of them.
Dispatching in pure table order instead, as it did before v0.4.0, gave the
first wave to `dbcan` (10 min) and `diamond` (5 min) while `signalp` and
`tmbed` (about an hour each) queued.

Both numbers are published per stage in `describe --json`, as `cost` and
`order_s`, so the dispatch order is reproducible from the release alone —
without reading anybody's results directory, and identically on a first run, a
fresh clone and a `--force`. Only their ORDER is used: list scheduling never
reads the magnitudes, so the same table orders a faster machine and a dataset
an order of magnitude larger the same way, which is what makes one run's
measurements a fair table for every run. The run's log names the order it
took, once, before the first dispatch.

The order is a starting order and not a schedule — nothing in it decides when a
stage finishes — but the sentence that used to stand here, "it makes nothing
faster", was wrong, and it was wrong in the same way the code was: on every
selection this repository has published durations for, dispatching `interpro`
first instead of seventh takes hours off the run, and at the default
`stage_workers: 4` it lands the whole pipeline exactly on the dependency
graph's critical path — the shortest that any ordering of the same work can
make it. THE SLOT COUNT IS PART OF THAT CLAIM AND NOT A DETAIL: at
`stage_workers: 3`, which is what the 455,571-protein run used, the same
reordering saves 13.7 h and is still 13.4 h above the critical path, because
there the run is bounded by the work rather than by the graph. What is left
after that is `stage_workers` and InterProScan itself.
See [Sizing your run](TUTORIAL.md#sizing-your-run) for both.

```yaml
threads: 32
stage_workers: 4      # 4 stages at once, 8 cpu each
ram_gb: 128           # 32 GB each; 0 = 80% of detected RAM
```

`--threads` and `--ram` override both from the command line on `run` and `all`
(`doctor` reads only the config). The memory budget takes the same spellings in
both places, because the flag and the config key go through the same parser: a
bare number of gigabytes (`64`), or a suffixed size (`64G`, `64GB`, `512M`,
`1T`). So `ram_gb: 64G` in config.yaml is valid. A size below 1 GB rounds up
rather than down, since `0` is reserved to mean "auto-detect 80% of detected
RAM", and a value that is not a size at all exits with a message rather than a
traceback. Note that `init` writes the config with `yaml.safe_dump`, which
drops the comments that explain these units, so the generated file shows a bare
`ram_gb: 0`.

### Where the memory budget goes

The budget is split across concurrent stages exactly as threads are, then
translated into each tool's own flag:

| tool | flag | note |
|---|---|---|
| DIAMOND | `-b` | peak memory is roughly 6 GB per unit of block size |
| MMseqs2, Foldseek | `--split-memory-limit` | |
| hhblits | `-maxmem` | default is 3 GB and silently truncates alignments |
| InterProScan | `_JAVA_OPTIONS=-Xmx…` | read by the JVM whichever launcher is used |
| eggNOG-mapper | `--dbmem` | only above `emapper_dbmem_min_gb`; below it, loading the annotation database just swaps |
| hmmsearch, KOfamScan | — | no memory flag; footprint tracks `--cpu`, so the CPU split is the control |

`doctor` prints the resulting per-stage allocation before you commit to a long
run.

### When a flag is wrong

This code is pinned to the CLI of a dozen programs whose options move between
versions. `tool_args` appends flags verbatim, so a wrong or missing one is a
config change rather than an edit:

```yaml
tool_args:
  diamond: ["--comp-based-stats", "0"]
  hmmsearch: "--nonull2"
```

A results directory takes a lock for the run, so two processes cannot
interleave their writes. Reclamation is deliberately conservative: a lock is
removed automatically only when it names a pid **on this host that is provably
gone**. A lock written on another host, one owned by another user, a garbled
lock file, and *every* lock on Windows — where `os.kill(pid, 0)` calls
`TerminateProcess`, so asking whether the holder is alive would kill it — all
read as live. Trampling a live run is silent corruption; refusing is a
message. So an ordinary crashed run on Windows exits with `another metaannot is
already running here` and needs `--force-unlock`, which is the escape hatch for
a holder you are certain is gone.

`--force-unlock` is refused in the one case where this host can prove you are
wrong about that: the lock names a pid **here**, and the process table says it
is running. The refusal prints what you would otherwise go and assemble — pid,
host, when it started, when it last stamped `_run`, the stages it has recorded
running, its command line, and the `ps -p` to run — and `--force-unlock-live`
takes the directory anyway if that is really what you mean. Note the asymmetry
with the paragraph above, which is deliberate: reclaiming needs proof of
DEATH, refusing needs proof of LIFE, and everything unprovable — another node,
another user's garbled file, every lock on Windows — is left exactly where it
was. A refusal that fired on "not provably dead" would take `--force-unlock`
away on the cluster, which is the machine this tool runs on.

The `_run.last_seen` heartbeat below does **not** change that. It is advisory:
it tells you *how long* a lock has been silent, and nothing reclaims a lock on
the strength of it. A heartbeat that stopped is not proof that a process
stopped — one failed write ends the timer, not the run — and from another host
the two are indistinguishable, so acting on it would trade a stale lock for two
runs writing one directory. Deciding that a silent holder is dead is
`--force-unlock`, which is a person deciding; the heartbeat is there to give
that person the number.

The one exception is a **zero-byte** lock — what a power loss, an OOM kill or a
host bugcheck leaves behind between the `O_EXCL` create and the write. A real
writer closes that gap in microseconds, so an empty lock that has sat unchanged
for more than a minute cannot belong to a live run, and it is removed on sight
with a WARN naming when it was left. A freshly created empty lock still blocks,
so the genuine race stays safe.

A stage whose dependencies were excluded by `--only`/`--from` refuses rather
than producing a confident answer from inputs that do not exist.

`--serial` forces one at a time. Parallel and serial runs are verified to
produce the same output: `tests/test_outputs.py` digests a parallel run against
a `--serial` run for every `quant_format`, and for the `protein_unique`,
`taxon_unique` and `razor` assignment modes; `tests/test_scheduler.py` compares
every file under `results/` byte-for-byte between the two.
`taxon_or_family_unique` is not in that sweep. The full
`quant_format` × `peptide_assignment` cross-product carries the `slow` mark and
is excluded from the default run.

Within stages: DIAMOND runs its databases concurrently (`diamond_workers`,
sublinear thread scaling makes 4×N/4 faster than 4 sequential N), and hhblits
parallelises its per-query searches (`hhblits_workers`). Foldseek's multiple
targets stay **serial on purpose** — a target index is tens to hundreds of GB
and two at once will thrash. ESMFold is GPU-bound and serial by nature; see
"The GPU is leased, not shared" below for why it does not run beside `tmbed`.

The scheduler's own overhead is negligible next to the search tools. No
benchmark script or speed-up measurement ships with this file, so no figure is
quoted here.

**Peak memory scales with `stage_workers`.** Four hmmsearch jobs against
Pfam-A alongside InterProScan is the usual squeeze. `doctor` says so.

### The GPU is leased, not shared

CPU and RAM are split between the stages running at a time; the GPU is not
divisible in the same way. `tmbed` held 15.5 GB of a 16 GB card and ESMFold
peaked at 13.3 GB on a single short sequence, so with `run.topology` and
`run.structure` both on they do not fit together — and the one that loses dies
of a CUDA OOM that names no cause. Stages marked as needing the GPU (`tmbed`,
`esmfold`) therefore take an exclusive lease, and only `gpu_workers` of them
run at once:

```yaml
gpu_workers: 1     # GPU stages at a time; CPU-only stages keep running
```

The rest of the pipeline is *not* serialised: every CPU-only stage is
dispatched in the same round as before. A stage that is waiting says so, so an
enabled stage that has not started is explained rather than mysterious:

```
[  312.0s] INFO  --- esmfold: waiting for the GPU — tmbed is using it and
                     gpu_workers is 1. It starts when that stage finishes;
                     everything else carries on meanwhile.
```

`gpu_workers` is a **lease count, not a device map.** Every GPU stage is pinned
to the single `gpu_device`, so on a two-card machine raising it to 2 runs two
stages on the *same* card rather than one per card. One stage per device is not
implemented; the limit here is policy, not hardware.

**The lease follows the stage, not the device.** `gpu=True` is a fixed
property of the stage table, not a probe of what the stage does at run time.
`tmbed` takes the lease whether or not it ends up on the GPU — with
`tmbed_use_gpu: false`, or on a machine with no CUDA device, it still holds the
slot that would otherwise let `esmfold` start. If you are running tmbed on CPU
deliberately, raise `gpu_workers` so the lease stops serialising two stages
that are no longer competing.

### When a model is missing

`structure_requested` says a protein was on the fold work-list;
`structure_attempted` says a `.pdb` for it exists. When the second falls short
of the first the run says which of three things happened, because they are not
the same event:

| what the log says | what it means |
| --- | --- |
| `INFO … none of the N proteins in dark.faa have been folded yet` | `esmfold` has not run: no `.done` marker, no `plddt.tsv`, no models. Nothing has been lost; the pass carries no structural evidence. |
| `WARN … esmfold has not finished` | Models exist but the stage did not complete. The rest are pending — rerun `esmfold`, which resumes. |
| `WARN … requested structures exist … although esmfold has finished` | The only case where a model can be missing for a bad reason, and the message names its causes rather than guessing one. Proteins longer than `max_len_structure` are counted separately as never submitted — that count uses the **static** limit only, so a protein excluded by the tighter VRAM cap falls into the remainder instead. The remainder is then split by `results/structures/esmfold_failed.tsv`: proteins listed there were attempted and failed twice, with the error in the table; proteins absent from it were added to `dark.faa` after the last fold and were never attempted at all. |

The first structure run of the real dataset printed *"0/1913 requested
structures exist …; the rest were skipped (OOM) or never folded"* while the
directory was empty because nothing had been folded **yet** — which reads as
1,913 models lost to the OOM killer.

**A skipped protein is now a record, not an inference.** A protein whose two
fold attempts both fail is written to `results/structures/esmfold_failed.tsv`
with its length and the error, so it can be told apart from one added to
`dark.faa` after the last fold. `esmfold` also logs `WARN skipping <id>` as it
happens. The count itself is taken over the proteins this pass requested
rather than over every `.pdb` in the directory, so models left from an
earlier, larger `dark.faa` cannot mask a shortfall — but when `finalise` takes
its "no structure or profile evidence, reusing the first pass" path, no
shortfall message is printed at all.

### Did the parts of a merged database annotate alike?

The search database for a metaproteomics run is normally a merge, and one
headline coverage over it is an average of things that are not alike. The run
this tool was built on is four tiers by identifier prefix: `uhgpL_` (395,467
UHGP proteins), `OIDECCNN_` (43,139 Prokka calls off the matched metagenome),
`uhgpSM_` (14,797) and `ampS_` (2,168). One of those arrives with precomputed
eggNOG annotations and another is ORFs nobody has ever seen — which is the
reason the proteome was merged in the first place, and exactly what "94% have
an eggNOG hit" hides.

`finalise` writes `results/tier_coverage.tsv`: one row per prefix with the
count, the share of the proteome, the percentage carrying each kind of
evidence, the percentage dark, and the median export score. `emapper` logs and
records the same split for its own coverage, which is where the difference is
starkest.

**The tag is not the key.** A tier tag says which *source* a protein came
from; the identifier key says what the id actually *is*, and it lives behind
the tag:

| id | tier tag | identifier key |
|---|---|---|
| `uhgpL_MGYG000004906_01237` | `uhgpL_` | `MGYG#_#` |
| `uhgpSM_MGYG000009567_01280` | `uhgpSM_` | `MGYG#_#` |
| `OIDECCNN_00158` | `OIDECCNN_` | `#` |
| `ampS_AMP10.000_478` | `ampS_` | `AMP#.#_#` |

On the real database that is **four tiers over three key spaces**: `uhgpL_`
and `uhgpSM_` (and the `ent_` entrapment set) all wrap the same MGnify
`MGYG…` namespace, 31.8M of the search database's 36.6M records. Two tiers
sharing a key are one namespace under two labels — the same protein appears
once per tag, one row of a precomputed annotation table annotates all of
them, and **every** such tag must be in `emapper_strip_id_prefix` or its tier
loses that table entirely and reports as unannotated. `tier_coverage.tsv`
carries `key_shape` and `key_shape_pct` columns, the run names any key shared
by more than one tier, and `emapper` warns while it can still be fixed if one
sharing tier is listed and another is not.

Digit *runs* are masked rather than digits, because widths vary inside one
namespace: the Prokka tier runs `OIDECCNN_00001` to `OIDECCNN_1712297` — 5-,
6- and 7-digit accessions, all one key space. Masking per digit would report
three.

Tiers are detected, not configured: the prefix is the text up to the first
`_`, `|`, `:` or `.`. Two cases are declined rather than guessed at — one
prefix over everything (contig ids all share `k141_`, so splitting says
nothing) and more than twelve (one Prokka locus tag per MAG is not a source
label). When that happens no file is written **and the run says why**, so an
absent `tier_coverage.tsv` is never ambiguous between "not applicable" and "a
stage failed". A UniProt-style merge splits into `sp|` and `tr|` for free.

A column absent from the annotation frame is skipped rather than reported as
0%: "the stage did not run" and "the stage found nothing" are different
answers and must not share a cell.

### Do the sources agree, or merely overlap?

Two sources reaching the same protein is not the same as the two of them
agreeing about it, and only the second says the evidence is corroborated.
`finalise` writes `results/source_agreement.tsv`, comparing every pair of
columns that shares an identifier namespace:

| comparison | what it tests |
| --- | --- |
| `pfam_accs` vs InterProScan's `Pfam:` | one Pfam-A library, two implementations |
| `ncbifam_accs` vs InterProScan's `NCBIfam:` | one NCBIfam library, two implementations |
| `ko` vs `kofam_ko` | orthology by DIAMOND against orthology by per-family HMM |
| `pfam_hits` vs `pfams_emapper` | domains found directly against domains carried by the ortholog |

Each row counts `identical`, `overlapping` (they share calls but one has more)
and `disjoint` (both called something and they share nothing), plus which side
is the superset when they differ. On the UC run:

```
comparison                          both  identical  overlapping  disjoint  pct_agree
pfam: hmmsearch vs interproscan    34106      24834         9272         0      100.0
ncbifam: hmmsearch vs interproscan 12677      12364          311         2      100.0
ko: eggnog vs kofamscan            13797      12414         1021       362       97.4
pfam names: hmmsearch vs eggnog    25393      16300         8081      1012       96.0
```

Read `disjoint` first. A few are ordinary — a paralogue boundary, a threshold
near the edge of a family. **A rate near 100% is almost never real
disagreement**, and the run says so rather than reporting a number that invites
the wrong conclusion: it means the two columns hold different kinds of
identifier and nothing is being compared at all.

That is not hypothetical. `ncbifam_hits` holds family *names* (`PorV_fam`),
while InterProScan reports *accessions* (`NF033709`). Comparing them scored
97.2% "conflict" between two searches of one library. Adding `ncbifam_accs`
— the accession was in the `hmmsearch` output all along and was simply being
discarded — turned that into 100% agreement over the same 12,677 proteins.

The direction matters too. Where hmmsearch and InterProScan differ on Pfam,
hmmsearch is the superset in **every** case and InterProScan in none: the two
never contradict each other, one is simply more sensitive. A pair like that is
called out in the report, because "they disagree 27% of the time" and "one
finds more than the other" are very different claims.

### The length TMbed can actually embed

TMbed embeds with ProtT5, whose attention score matrix is length-squared ×
heads, so one sequence costs roughly `len² × 32 × 4` bytes: 1.1 GB at 3,000
residues and 141 GB at titin's 34,350. That is a **single** allocation, so no
device setting saves it — observed twice on real data, 8.79 GiB refused on a
16 GB card and then 151 GB refused on a 94 GB host under `--cpu-fallback`, each
time hours in and with nothing written, because TMbed writes its predictions
only at the end. So the input is capped rather than the failure retried:

```yaml
tmbed_max_len: 3000     # residues; 0 disables the cap
```

Proteins longer than this are excluded from the search and listed with their
lengths in `results/topology/tmbed_excluded.tsv`, and the run says how many
were cut and how long the longest was. They get **no topology evidence** — no
`n_tmh`, no `n_tmb` — so for those proteins the export score's beta-barrel term
and the TMbed half of the `surface_or_secreted` gate are missing evidence
rather than evidence against, exactly as they are on a run with `topology` off.
It cuts both ways: `n_tmh` also drives `multi_tm_helix_penalty`, so an excluded
protein escapes that penalty as well as forgoing `tm_beta_barrel`. Its
`export_score` is not simply lower than it should be — it is uninformed, and
can land either side of the score it would have had.
The key is part of the `tmbed` stage's signature, so changing it re-runs that
stage, and the predictions it writes are an input to `integrate`, so that and
`finalise` follow. No other stage recomputes.

### TMbed writes nothing until it finishes

Not "buffers a bit" — nothing. On the 455,571-protein run it produced no
progress bar and no partial file for 31 hours, and the two attempts before it
died at 2 h 36 min with nothing recoverable. A single invocation over a whole
proteome is therefore an all-or-nothing bet measured in days.

```yaml
tmbed_chunk_residues: 5000000   # 0 = one invocation, the old behaviour
tmbed_allow_partial: false      # true: finish without the failed chunks
tmbed_max_consecutive_failures: 2
```

The input is grouped into chunks of about that many residues and each is
committed as it lands, so an interrupted run resumes from the last finished
chunk. Under 5M residues — roughly 17k average proteins — there is exactly one
chunk and the behaviour is what it always was.

Chunks are **length-sorted, longest first**, for two reasons: ProtT5 pads every
sequence in a batch out to the longest one in it, so a chunk of similar lengths
wastes less work than one mixing 30-residue peptides with 3,000-residue
proteins; and whatever is going to exhaust the device is then in the *first*
chunk, where it costs one chunk to discover instead of the whole stage. The
count is bounded at both ends — the budget is a floor, 256 parts a ceiling —
because ProtT5 is loaded once per chunk.

A chunk that dies keeps whatever TMbed wrote. What ends up in the committed
file is reconciled against what was handed over, from the files rather than
from the loop's bookkeeping, because a chunk can exit 0 and still come back
short; anything missing is named in `results/topology/tmbed_failed.tsv`. By
default that shortfall is fatal, so nothing downstream reads a partial
topology set by accident, and `tmbed_max_consecutive_failures` stops a wedged
card from failing every remaining chunk the same way, slowly.

`tmbed_chunk_residues` is deliberately **not** part of the stage signature. It
changes how the work is divided, and therefore the order of the records, but
not one prediction in them; listing it would discard a 30-hour stage because
someone tuned a checkpoint size. `tmbed_allow_partial` and
`tmbed_max_consecutive_failures` *are* listed, because they change what is in
the file.

### The length a card can actually fold

ESMFold's cost does not rise smoothly with length. It rises smoothly until the
working set stops fitting in VRAM, and then falls off a cliff. Measured over
1,819 folds on one 16 GB card:

| sequence length | median fold time |
| --- | --- |
| 450-470 aa | 20.7 s |
| 470-478 aa | 22.1 s |
| **481 aa** | **140 s** |
| 486-491 aa | 949 s, then 2,053 s |
| 495-510 aa | ~305 s |

A 0.6% increase in length cost 6x, and shortly after that 90x. **Nothing
reports an out-of-memory error**, because the driver pages device memory to
host RAM rather than failing — so the run does not stop, it stops being
finishable, and the only outward sign is that a stage which was going to take
an hour is now going to take a week.

So the fold work-list is capped by the memory that is actually free, measured
once the weights are resident:

```
max_len = sqrt((free_vram - esmfold_vram_reserve_gb) / esmfold_bytes_per_residue_pair)
```

Peak footprint above the weights is dominated by terms quadratic in length —
the pair representation and the triangular attention over it — hence the
square root. `esmfold_bytes_per_residue_pair` is empirical, calibrated so the
measurement above (478 aa at 4.8 GB free) comes out at exactly 478. It is a
config key because one card is not a law:

```yaml
esmfold_vram_cap: true                  # false: use max_len_structure alone
esmfold_bytes_per_residue_pair: 20200   # raise = more conservative
esmfold_vram_reserve_gb: 0.5            # left for the driver and the display
```

Whichever of this and `max_len_structure` is tighter wins, and the log says
which. Sequences above the cap are reported as **never attempted**, not as
failures — and they are now written down, because "fold them elsewhere" is
advice a reader cannot act on when the log gives only a count. Every run of
the stage writes both files, empty ones included, so that their absence means
the stage did not run rather than meaning nothing was skipped:

| | |
|---|---|
| `results/structures/not_folded.tsv` | `protein_id`, `length`, `limit_aa`, `limit_from` |
| `results/structures/not_folded.faa` | the same sequences, ready to fold |

Fold the `.faa` on a card with more memory, in the cloud, or on CPU, and drop
the models into `results/structures/` before rerunning `foldseek`.

**A low cap is a statement about free VRAM, not about the coefficient.** On
the first full real run the cap came out at 193 aa and excluded 1,462 of 2,000
sequences, which reads like a badly-guessed constant and is not one: 20200
reproduces both measured points exactly, and running the estimate backwards,
a 193 aa cap means about **1.2 GB** was free with the weights resident, where
the 478 aa measurement had 4.8 GB. On a 16 GB card holding an 11.2 GB trunk,
4.8 GB is what should be left — so roughly 3.6 GB was held by something else,
and that, not `esmfold_bytes_per_residue_pair`, is what to go and look at. The
warning now prints the free VRAM the longest skipped sequence would have
needed, one line under how much was actually free, so the comparison is on the
page rather than left to be worked out. Lower the coefficient only where a
card **demonstrably** folds longer sequences at a smooth rate; the guard is
there because crossing the cliff bugchecked a host twice.

Whether the skipped tail was worth much is a separate question and this run
answered it for one dataset: 538 completed folds produced only 14 proteins
whose *sole* evidence was structural. That is an argument about how much to
spend chasing the cap, not an argument for raising it blind.

It can get worse than slow. On the machine these numbers came from, folding
above the cliff also produced repeated `CUDA driver error: device not ready`
faults and then took the whole host down twice inside eleven minutes with a
hypervisor bugcheck — the GPU there is reached through a virtualisation layer
(WSL2), and sustained paging across it is what broke. That is a defect in
somebody else's code and not something this tool can fix, but staying under
the cliff avoids it. For the same reason, a stage with nothing left to fold
now returns **without loading the weights at all**, rather than uploading
~11 GB to the card and then discovering it had no work.

### When the card fails rather than the protein

Not every fold failure is an out-of-memory. A long sequence at a large
`esmfold_chunk_size` can run a single attention kernel for long enough that
the display driver resets the device under it, which arrives as
`CUDA driver error: device not ready` — a plain `RuntimeError`, not an
`OutOfMemoryError`. Both answer to the same remedy, so both are retried once
at half the chunk size, and a sequence that fails twice is skipped rather than
taken as a reason to abandon the stage.

Every `.pdb` is written as it is folded, so a card that dies at protein 1,800
of 1,900 has cost the tail and nothing else; rerunning `esmfold` resumes.
What a rerun cannot fix is a card that has stopped responding altogether, and
walking the remaining list to fail every one takes hours to produce nothing.
So after `esmfold_max_consecutive_failures` failures in a row (default 5) the
stage stops and says how many structures it has:

```yaml
esmfold_max_consecutive_failures: 5     # in a row = a wedged card, not hard proteins
esmfold_allow_partial: false            # true: go on with what folded
```

With `esmfold_allow_partial: true` the stage finishes with the structures it
has and Foldseek searches those. That is a real reduction in evidence, not a
neutral setting: a protein with no structure hit may simply never have been
folded. Leave it off unless the alternative is no structural evidence at all.

### Telling a live stage from a hung one

Every minute, a running tool's newest line of stderr is echoed under its stage
tag with how long it has been going:

```
[ 9421.3s] INFO       tmbed | tmbed running 2h36m | 61%|######    | 23310/38204 [2:36:04<1:39:41]
```

Set `progress_interval_s` to change the interval, or to `0` to turn it off:

```yaml
progress_interval_s: 60   # seconds between progress lines; 0 = silent
```

Only the newest line is kept, not the output: stdout still goes to
`/dev/null`, because InterProScan and friends emit tens of MB of chatter, and
stderr is held in a small ring whose only other use is the tail quoted when a
tool fails. That tail is unchanged. Before this, tmbed could run for 2 h 36 min
and InterProScan for 2.9 h with nothing between the command and its failure,
and the only way to tell either apart from a hang was to watch its CPU ticks
accumulate in `/proc`. Tools that write a `tqdm` bar (tmbed, InterProScan,
ESMFold's own loop) are the ones this shows; a tool that writes nothing still
gets `no output yet on stderr` on the same schedule, which is the heartbeat.

#### How far through InterProScan is

InterProScan is the exception that paragraph does not cover. It is the longest
stage in the pipeline — 57 h of the 455,571-protein run — and what reaches
stderr from it says that it is alive without ever saying how far through it
is. What it does leave is a trail under the `-T` directory this tool hands it:
a `.fasta` per chunk it splits out, and a `.raw` beside that chunk once the
chunk has been analysed. The heartbeat counts them:

```
[31840.7s] INFO    interpro | interproscan.sh running 8h50m | chunk 312/380, ~2h08m left | ...
```

What that line does not claim, deliberately:

- The unit is **chunks**, not sequences and not a percentage of the stage.
  InterProScan chose the slices, they are not equal, and the merge and write
  that follow the last one are not counted at all — so it reaches `380/380`
  with real work still to do. Do not turn it into a percentage.
- **No remaining time is offered while the denominator is still moving.**
  InterProScan goes on splitting while it analyses, so a rate taken then is
  measured against a number that is about to grow — and it is perfectly
  computable and perfectly wrong, because chunks really are completing. The
  count has to hold still for `_IPS_ETA_STABLE_TICKS` consecutive heartbeats
  before a rate is anchored, and the rate is then a long-run average from that
  anchor rather than an instantaneous one: it lags a machine that slows down
  later, which on a stage measured in days is the trade worth taking.
- **`.raw` beside `.fasta` is the layout of one observed run.** No
  InterProScan was available where this was written to check it against, so a
  build that writes a different tree gets no progress on the line at all
  rather than a wrong one — and the stage says so once when it finishes:
  `the heartbeat never found a chunk file under …`. Silence from a probe that
  never engaged must not be read as a stage that made no progress. The same
  goes for a census too large to walk: it reports nothing rather than a
  partial count.

Any stage can pass `run_cmd` a probe of its own the same way; InterProScan is
simply the one where the tool's own output leaves the most unsaid.

Log lines cannot themselves kill a run. Tool output is decoded with
`errors="replace"`, so a description can carry U+FFFD, and printing one to a
Windows console in cp1252 raises `UnicodeEncodeError` — on the log line rather
than on the work, hours in. `stdout` and `stderr` are reconfigured to UTF-8
with `errors="replace"` at startup, and each write falls back to replacing
what the stream cannot encode.

**What the progress line is not.** It is the tool's own newest line of stderr,
verbatim — not a parsed percentage — so nothing here estimates time remaining,
and a tool that prints a spinner or a fixed banner repeats it every interval.
Three limits are worth knowing before relying on it:

* **Only stderr is watched.** stdout still goes to `/dev/null`, so a tool that
  reports progress there gets `no output yet on stderr` and nothing more. That
  heartbeat proves the process is alive; it says nothing about how far it has
  got.
* **Only external commands are watched.** `integrate`, `finalise`, `join`, the
  quant readers and the emapper join do their work in-process and emit no
  progress lines however long they take, so a silent stretch during those is
  not evidence of a hang either way.
* **Progress is not a checkpoint.** Watching a stage does not make it
  resumable. `esmfold` and `hhblits` resume, because they work one protein at a
  time and skip what is already on disk (`<id>.pdb`, `<id>.hhr`). Every stage
  that shells out to one long command — `hmmsearch` for `pfam`, `dbcan` and
  `ncbifam`, InterProScan, KOfamScan, `tmbed`, DIAMOND, `jackhmmer`, Foldseek,
  MMseqs2 — starts again from the beginning if it is killed at 90%, and a stage
  recorded as `running` when the process died is always recomputed rather than
  adopted.

`progress_interval_s` is read once, at the start of `run` (and `all`), and
nothing else reads it, because nothing else calls the wrapper that watches a
running tool. `doctor --fix` does not use this machinery at all: it runs each
recipe through the shell with the terminal attached, so you see curl's or
conda's own progress rather than a metaannot heartbeat. Standalone `report` and
`object` are outside it too — they shell out to Rscript through their own
wrapper, which collects stderr and quotes it only on failure.

### A FIFO at an input: only where the run reads that input once

Every parser in `metaannot.py` reads through one function, `opener()`, and a
named pipe is the one input state where what it should do is genuinely
contested. The rule is short and the reason it is this rule was measured, not
argued:

> **A FIFO can be drained exactly once, so it works at an input this run opens
> exactly once — and is refused, immediately, at an input it opens more than
> once.**

`mkfifo p; zcat big.faa.gz > p &` is how you feed this tool from a disk you
have no room on, and it is a supported input where the rule allows it, for a
plain stream and a gzipped one alike. Where it does not, the refusal is at the
first open and says why:

```
input/proteins.faa: this is a FIFO, and this run reads that input 3 times:
  - stage_emapper -> prepare_emapper -> read_fasta
  - stage_integrate_pass1 -> build_annotation -> read_fasta (ids and lengths)
  - stage_integrate_pass1 -> build_annotation -> read_fasta (the dark work-lists)
A FIFO can be drained EXACTLY ONCE. ... Write the stream to a real file first
```

**Which inputs those are is worked out from your config, not from a list in a
document.** `INPUT_READ_SITES` in `metaannot.py` names, per input, every place
a run opens it and the condition under which that read happens; `run` and
`doctor` both compute the plan for the config in front of them, because the
answer changes with it — a project with `run.unipept` on reads its quant table
once more than one without, so a pipe that works in the first config does not
work in the second. `doctor`'s row for a FIFO prints the count and the reads.

On a default label-free project the plan is: `proteins_faa` **3** (the emapper
stage, then the integrate stage twice), `quant_table` **1**, `manifest` **1**,
each `emapper_precomputed` entry **1**. So a pipe works at all of those but
the FASTA. That is a change: `quant_table` and `manifest` were read three
times and twice before, in a way nothing could see — `header_columns()` opened
the table and then `read_delim_table()` opened it for the delimiter and pandas
opened it for the body, and the join read the manifest once for the column
mapping and again for the design. Those reads are now one open each.

**Turn a stage on and the plan changes, which is the point of computing it.**
With `run.taxonomy` on and a feature-level table, `manifest` is **2** — the
join reads it through `read_feature_table()`, and so does the taxonomy stage
through `peptide_features()`, which reaches the same reader. With
`unipept.result` set, `stage_unipept` ingests that export and returns before
it ever opens the quant table, so `quant_table` does **not** gain a read from
it. With `run.eggnog` off, `emapper_precomputed` is **0**. The plan is derived
from the same functions `doctor` uses to decide which stages open which input,
so a row and a refusal cannot disagree about it.

**The database paths are in the plan too.** `db.ncbi_taxonomy` names a
directory, and `nodes.dmp`, `names.dmp`, `merged.dmp` and `delnodes.dmp`
inside it are read by `NCBITaxonomy()` — once with `run.taxonomy` on, once
more when `taxon_rank` is set, because `stage_join` builds its own. So on a
config with both, a pipe cannot work at any of them and the refusal says so;
on a config with one, it can. `unipept.result` is read once. All of them used
to be read through bare `open()` calls outside the choke point, which is to
say a FIFO at any of them hung the run with no exit status and the results
lock still held.

**A TMT quant tree is in the plan file by file.** With `quant_format:
fragpipe_tmt`, `quant_table` names the run DIRECTORY and not a table, and a
directory is not something any reader opens — so what the plan counts is what
the readers really open: each plex's `ion.tsv` (or `peptide.tsv`), each plex's
annotation file, and `psm.tsv` where `tmt.min_purity` reads one. Every stage
that consumes the tree opens each of them once, so with the taxonomy stage and
`run.join` both on a pipe cannot work at any of them and the refusal says so;
with one consumer it can. The annotation and `psm.tsv` are counted separately
from the level file because the peptide-only reader never opens them.

**One read in the plan is conditional, and it is marked as one.** Under
`peptide_only_reader: auto`, which is the default, a taxonomy or unipept stage
that meets a quant table the full reader REFUSES re-reads it with the
peptide-only reader — a second open of the same path, and whether it happens
is a property of the table rather than of your config. The plan cannot promise
it will not, so it counts it: a pipe at a quant table a peptide stage opens is
refused rather than blessed, and the refusal lists that read as the
conditional one it is. `peptide_only_reader: always` skips the full reader
altogether, which is one open again.

**`report` reads your metadata through the same open.** On a TMT project the
contrasts are derived from `analysis.metadata`, which is a path you wrote in
the config, and it goes through the same gate every reader here uses. `report`
sets no read plan — it is not a run — so the honest single-read default
applies: a FIFO there is announced and streamed, and a second open of the same
drained pipe is refused by name instead of blocking with no exit status.

The **signature digest** is deliberately not in that plan, and it is the
subtle one. Every input listed in a stage's signature is also read by
`_content_digest()`, so on a regular file each of these paths is opened once
more than the numbers above. It is not counted because it never touches a
stream: `_open_regular_binary()` opens the path, `fstat`s the **descriptor**
and refuses anything that is not a regular file, so the digest cannot drain a
pipe and cannot be raced into doing so by a path that changes under a stat.

```yaml
fifo_wait_s: 21600   # seconds to wait for the NEXT BYTE on a FIFO; 0 refuses one
```

At a single-read input, the read is announced before it waits:

```
[   12.4s] INFO     eggnog | input/cat.emapper.annotations.gz: this is a FIFO, and
                            this run reads it exactly once, so it is read as a live
                            stream. Each read waits up to 21600s for the next byte
                            (`fifo_wait_s`); a live writer streams straight through
                            and nothing else will be logged until one appears.
```

(one line on stderr and in `metaannot.log`; wrapped here to fit the page.)

The path is **opened once, non-blocking, and that same descriptor is what the
parser reads.** That is the whole design, and it is not an implementation
detail: opening the read end of a FIFO and closing it again sends the writer
on the other side `EPIPE`, so a check that probed the path and then let the
reader open it again would destroy the stream it was meant to protect. There
is no probe. `O_NONBLOCK` makes the open return instead of wait, and it also
releases a writer that is already blocked in its own `open()`; `fstat()` on
the descriptor — never a second stat on the path — says what was really
opened; and the flag comes off again before the read, because a non-blocking
descriptor raises `EAGAIN` the moment a healthy writer pauses, which would
turn a working pipe into an intermittent failure.

**`fifo_wait_s` bounds every read, not only the first.** A writer that
attaches, sends half a table and then holds the write end open without writing
again is the same hang one read further in, so the setting is an *idle*
timeout: each read waits that long for the next byte. A writer that is
streaming never comes near it.

**When the wait runs out, the message says which silence this is.** These are
two different problems with two different remedies, and telling the operator
the wrong one is worse than telling them nothing:

```
input/x.faa: this is a FIFO and nothing is holding the write end - nothing has
opened it in 21600s (`fifo_wait_s`) ... Start the writer before the run

input/x.faa: this is a FIFO, something IS holding the write end, and it has not
written a byte in 21600s (`fifo_wait_s`). The writer is attached, so starting
another one is not the fix ... raise `fifo_wait_s`
```

**A stream that stops early is a failure, not a shorter result.** A writer
that wrote half a FASTA and died used to yield the half and no error at all,
which is the worst outcome available here, because the number at the end of it
is wrong and looks right. A pipe now fails if the writer closed without
writing a byte, and — for a plain stream — if it stopped in the middle of a
line. A **gzipped** stream is the stronger form and is what to pipe when you
have the choice: gzip carries an end-of-stream marker, so a truncation is
caught wherever it happens rather than only when it happens mid-line. A plain
writer killed exactly on a line boundary is a case no amount of looking at the
bytes can catch, and that limit is stated here rather than left to be found.

**Six hours is the default because the two costs are not symmetric.** Waiting
too long costs only the tail of a mistake already made, and it costs it
visibly — the wait is announced before it starts. Waiting too briefly costs a
workflow that works: the writer this has to survive is not the shell one-liner
above but a producer that is itself a queued job on a shared cluster, where
waits run to hours. The upper bound is a run holding an exclusive lock on its
results directory the whole time it waits, so a job started at the end of a
working day should have failed with a message by the next morning rather than
still be sitting on the open. Set `fifo_wait_s: 0` if you never pipe anything
in and would rather a FIFO were refused on the spot. A **socket**, a **device
node** or a **directory** is refused at once whatever this is set to: nothing
there is ever going to become a table, and `doctor` names which of the three
it found rather than offering you all of them.

`doctor` is unaffected and never waits at all — it reads every path through
its own non-blocking probe, and that asymmetry is deliberate: a command whose
job is to answer *before* the run is worthless if it can block, while a `run`
that refused every pipe would be removing something that works.

### What a results directory says about itself

A results directory now says what produced it and whether that is still
happening. Neither file below is ever read back by metaannot, and neither is in
any stage's signature, so neither can make a stage recompute.

**`config.effective.yaml`** is the merged configuration the run actually used:
the built-in defaults, then your config file, then the command line. It is not
the same as the config file beside it, which carries none of the defaults, none
of `--threads`/`--ram`/`--faa`, and says whatever it says today rather than
what it said in March.

**`_run`** is a record at the head of `results/.metaannot_state.json`, above the
per-stage entries and deliberately not one of them — the key starts with an
underscore, no stage is named that, and `--force` never clears it:

```json
"_run": {
 "run_id": "20260901T144530-31284",
 "version": "0.3.0",
 "config_path": "/data/projects/gut2/config.yaml",
 "argv": ["<script>", "run", "--config", "config.yaml"],
 "host": "server", "pid": 31284,
 "started": "2026-09-01T14:45:30",
 "last_seen": "2026-09-01T19:12:00", "last_seen_epoch": 1788289920.0,
 "heartbeat_s": 30, "finished": null, "final_status": "running"
}
```

`argv` is `sys.argv` verbatim, so the `<script>` above is really the path of
the `metaannot.py` that ran — useful when more than one copy is installed.
`config_path` is the absolute path of the `--config` file, or **`null`** when
the run was given none and took the built-in defaults — a parser that types it
as a string meets a real file it cannot read the first time someone runs
`metaannot.py run` without `--config`. It is the only nullable field here
besides `finished`, which is `null` until the run ends.
`final_status` is `ok`, `failed`, `interrupted` (Ctrl-C) or `running`. A record
still saying `running` with a `last_seen` from hours ago is most often the
**ordinary** trace of a `kill`: `SIGTERM` releases the lock and exits without
unwinding, so nothing stamps a verdict on the way out and `running` is simply
the last thing the record was ever told (the next section spells out why that
is the right trade). It is also what a `SIGKILL`, a lost machine or a wedged
process leaves behind, and what a run whose **lock** was removed under it
leaves behind — such a run cannot prove it still owns the directory, so it
stops writing `_run` altogether and keeps recording its stages, and it says
that once in the log. From the record alone those cases are indistinguishable —
which is the point of writing it down rather than interpreting it, and the log
is where they are told apart. `last_seen` tells you how long ago the process
last wrote a `_run` it could prove it was entitled to write. It is for reading, not for deciding: nothing in metaannot reclaims a
lock because a heartbeat went quiet (see the lock section above). Both timestamps are the same
instant: the string is local time for reading, the epoch is for arithmetic,
because two hosts sharing one filesystem cannot subtract each other's local
clocks. `heartbeat_s` appears in the record exactly as you set it, so an
integer there stays an integer. Set the interval, or turn it off, with:

```yaml
heartbeat_s: 30   # seconds between last_seen stamps; 0 = no heartbeat
```

Both files carry more about your setup than the state file used to. `_run`
records the host name, the pid, the absolute config path and the full command
line, and `config.effective.yaml` materialises every merged setting including
absolute input and database paths — and `.metaannot_state.json` is embedded
whole in the `.rds` object people publish as supplementary data. Nothing in the
defaults is credential-shaped, but `sources.*` and `tool_args` are free-form:
if you have put a presigned URL or a token there, it now travels with the
results.

**The lock is released on `SIGTERM`**, not only on Ctrl-C. `kill`,
`systemctl stop` and `wsl --terminate` used to end the process where it stood,
leaving `.metaannot.lock` behind. `run` and `all` now install a handler for
`SIGTERM`, `SIGHUP` and (on Windows) `SIGBREAK` that does exactly four things,
in this order: it `SIGKILL`s the process **group** of every tool the run
started, it removes the lock file, it writes one line to stderr, and it calls
`os._exit(128 + N)`. The kill is first on purpose — `ResultsLock.__exit__`
reads and parses the lock file before it unlinks it, which on a wedged NFS
mount blocks indefinitely, and a wedged mount is exactly when a run gets
killed. Losing the lock release costs one `--force-unlock`; losing the kill
costs core-hours.

**That is a different path from Ctrl-C, deliberately, and the difference is
worth knowing before you `kill` a run.** The handler does not unwind, so there
is no `interrupted` message, no `_run` stamp, and no waiting for the stage that
is running. Ctrl-C is the opposite trade: `SIGINT` is left as Python's default,
raises `KeyboardInterrupt`, unwinds through the stage pool's `with` — which
**waits for its workers** — and only then prints `interrupted` and stamps the
record. That wait used to be the length of the stage, because nothing in the
program stopped the tool: a keyboard Ctrl-C reached it only because the tty
broadcasts `SIGINT` to the whole foreground process group. It is now short,
because the interrupt is caught **inside** the executor's `with` and kills the
tool groups there before re-raising, and a latch stops the workers from
starting the next chunk or the next query while it unwinds. `kill -INT` from a
script, which has no terminal and so never got that broadcast at all, now
behaves the same way. Unwinding on `SIGTERM` was tried and taken out again, because waiting
is precisely what a supervisor cannot afford: a tmbed chunk or an InterProScan
stage is an hour, `systemd` hits `TimeoutStopSec` long before that and sends
`SIGKILL`, and the lock release that is the whole point of handling the signal
is then lost. So `SIGTERM` buys the lock at the cost of the trace, and Ctrl-C
buys the trace at the cost of the wait. The one line the handler does write is
pre-formatted and pre-encoded at registration and goes out through
`os.write(2, ...)` — a handler runs between two bytecodes of the main thread,
so touching `sys.stderr`'s buffer lock can deadlock the process it was meant to
release — and it names the signal, the lock file, the fact that every tool the run
started was killed group and all, and the dot-prefixed `.part` file a
part-written output is left behind as. It reaches stderr but not
`results/metaannot.log`, whose buffer cannot be flushed from a handler.

The exit status is **128 + the signal**: `130` for Ctrl-C, `143` for `SIGTERM`,
`129` for `SIGHUP`. That is what a shell and `systemd` both expect — units carry
`SuccessExitStatus=143` precisely so a `systemctl stop` is not recorded as a
failed unit — and it is the channel in which a supervisor tells an operator's
stop from a person at a keyboard. Note that only `run` and `all` install the
handler, and only once they have taken the lock: every other subcommand, and
`run --dry-run`, still takes the unwinding path on `SIGTERM`, which costs
nothing because none of them holds a lock.

**So what a `kill`ed directory looks like, and what to do with it.** The lock is
gone, so the next run starts without `--force-unlock`. `_run` still says
`"final_status": "running"` with `"finished": null` — see above; that is the
signature of this path, not evidence of anything worse. The stage that was
mid-flight is still recorded `"running"` because `mark_running()` stamped it
before it started, and the next run reads that, says `the previous run was
interrupted while this stage was writing, or could not write the stage's
record, so its output may be truncated; recomputing`, and redoes it. Both
readings are in that sentence because from the next run's position they are the
same bytes: `mark_running()` writes before the stage starts and `finish()`
writes when it ends, so a state file that goes unwritable in between leaves
this record for a stage that finished perfectly well (the I/O part below). The
verdict is the same either way and it is the safe one. You do not have to
`--force` anything, and you should not delete anything.

One thing a `kill` does **not** leave is an account of what the run could not
record. A run that reaches its own exit says which records never got there;
`SIGKILL`, and the `SIGTERM` handler that releases the lock and `_exit`s
without unwinding, say nothing, because that handler may not format a string or
take a lock. On Windows `TerminateProcess` runs nothing at all. That is why
each such loss is also said in the log at the moment it happens.

**Under a supervisor the question does not arise, and that is worth knowing
before reaching for anything cleverer.** `systemd`'s default is
`KillMode=control-group`, so a run started as a unit or under
`systemd-run --scope` has its whole cgroup torn down on `systemctl stop`
whatever this program does or fails to do — `setsid` does not change cgroup
membership — and the same holds for a Slurm step. The leak this section is
about was always specific to a bare `kill` from a shell or a tmux pane, which
is how the run that burned the core-hours was started. If you have a machine
where runs get killed and you want the guarantee to hold even for `kill -9`,
the answer is a unit (or `systemd-run --scope -p KillMode=control-group -- \
python metaannot.py all --config config.yaml`), not code in here.

**The deaths no handler runs for are the ones left, and the next run is
what covers them.** A `kill` now takes the tools with it, but `SIGKILL`, the
OOM reaper, a host reset and Windows' `TerminateProcess` run no handler at all
— and putting each tool in a session of its own has made that case slightly
worse, because a tool is no longer in the terminal's foreground group and so
no longer dies with a closing tmux pane either. What such a death leaves is an
orphan: an `hmmsearch` whose parent is gone, still holding its cores, still
writing into this results directory, and the lock that would have kept a
second writer out is reclaimable as stale by anyone, because the process it
names really is dead.

So every run, immediately after it takes the lock and **before** it dispatches
anything, lists the in-progress `.part` files under its declared output
directories, says one aggregated `WARN` naming each directory, the minting pid,
the size, the mtime and which stage writes those names — and then stats them
again two seconds later. A file whose **size changes** is positive proof that
something which is not this run is writing here, and the run is refused with
the `lsof`/`fuser` to find the writer — with the one exception in the next
paragraph. A file that does not change is a corpse:
one line, no refusal, and **nothing is deleted, renamed or touched** (rule 5,
and because it is the only artefact showing what happened — an `unlink` would
not even free the space while the writer holds the inode, and it would destroy
the one handle that gets from the file back to the process). The two seconds
are spent only when a candidate exists. `find results -name '.*.part.*'` is the
same list by hand.

**`--force-unlock-live` is the one exemption, and it is per file.** That flag
means "take the directory even though this host can see the holder is still
running", and a live holder's tool is exactly what leaves a `.part` file
growing — so refusing there would refuse the flag's only case, every time,
with no escape, and would make the handover this section goes on to describe
unreachable. Where the pid that MINTED the growing file is the pid of the lock
this run has just taken from a holder **this host proved alive**, the run is
not refused: it says what it found, at `WARN`, names the pid and the growth,
says that you now have two writers in one results directory deliberately, and
goes on. Any other growing file still refuses the run, in the same census, on
the same evidence — including one lying beside it, because the flag is an
assertion about one process and not about the directory. And the exemption is
keyed on what the lock really displaced rather than on the flag being typed: a
`--force-unlock-live` over a lock that was vacant, garbled, on another node, or
held by a process this host proved **dead** displaced no live run at all, so
the orphan case above is refused exactly as before. The one thing it cannot
tell apart is a file minted by a dead run whose pid was later recycled onto the
live holder; the message prints the number so you can see that, and `lsof` is
what settles it.

The number in `.pfam.336.140234.part.tblout` is **metaannot's** pid, not the
tool's: `atomic_out` mints the name with `os.getpid()`, and a tool's own pid is
recorded nowhere. So `ps -p 336` may show something unrelated, or nothing — a
pid can be recycled by anything once its owner is gone — and nothing in the
program ever reads that number as evidence about a process. One surface
disagrees and will until a separate change lands: the console picks the newest
`.part` file beside a stage's output and shows it as that stage's live
progress, so it will attribute an orphan's bytes to a healthy run, and
eventually report that run as stalled.

**A run that is unwinding stops writing when it is superseded.** That is the
Ctrl-C case, and the reason it matters is that a run can still be inside a stage
when you decide it has hung and `--force-unlock` the directory for a
replacement. The lock read has three answers and they are not one rule:

* The lock **now holds somebody else's token** — a replacement took the
  directory and is still running. The old run writes nothing further into
  `.metaannot_state.json` at all, removes no lock, says
  `this run no longer holds ...` once in the log, and exits. The latch is
  permanent. So the `_run` record and the lock you see afterwards belong to the
  run you started.
* The lock is **gone** — the replacement took the directory and has since
  exited, or something removed the lock outright. That is not proof of either,
  so the old run stops writing `_run` (which is an ownership claim) and carries
  on recording its own stage records (which are not). It says
  `is no longer there, so this run cannot prove it still owns this directory`
  once. Its `_run` then stays at `"final_status": "running"` for good.
* The lock is **unreadable** — an NFS `EIO`, an `EACCES`. That is no evidence
  that anything changed hands, so nothing changes; a transient error must not
  be able to strand a lock or silence a run.

A `SIGTERM`ed run cannot reach any of that — it is gone before it could write
anything — which protects the replacement just as effectively and says nothing
about it.

**What that protects, and what used to be left out of it.** The lock was never
the hole here: `__exit__` removes on `is_still_ours() is True`, and that answer
covers a lock it could not read at all, deliberately — a lock we cannot read is
one we have no evidence has changed hands, and leaving it behind was the
regression that rule was written against. `_run`
was *believed* safe and was not — the gate on it declined only the write that
would create the state file, and a run's own stage record creates that file one
call earlier, so a superseded run's `_run` went in behind its own `pfam` record
with nothing latched and nothing warned. That is fixed here; the paragraph on
`_run` below says what the real gate costs. Two more things were never covered
at all, and both are addressed now.

The first was every OTHER stage's record. A state write used to hand the whole
document to `save_state` from one process's in-memory dict, so a superseded run
recording the single stage it had just finished replaced the file with a
snapshot that had never heard of the stages the replacement completed
meanwhile. A reviewer reproduced it: run B finished `pfam`, `dbcan`, `diamond`
and `cluster`, run A wrote its own older one-stage view over the top, the
other records were gone, and the next run printed `adopting output this run did
not produce` for every one of them — the warning that says outright it cannot tell a
finished file from an interrupted one. Writes name the keys they change now:
the file is re-read immediately before each one, the changed key is merged into
what is there, and the file is written back and read back to check it still
says what we wrote. Where that re-read comes back **missing or unparseable** —
an operator's `rm`, a remount, a page of NULs where a crashed writer's payload
should be — the write creates a document holding *only the keys it names*, and
every other record stays gone. It does **not** rebuild the document from the
writing run's own snapshot. Three earlier revisions of this fix kept that
rebuild and tried to gate it, and there is no gate that works: the losing
ordering needs no race at all, because a replacement that has *finished*
leaves a vacant lock and no `_run` to read, and a vacant lock reads exactly
the same as a replacement that has not started yet.

Read what that buys the way the outputs half below is read — as a guarantee, a
clock and an uncovered part — because it does not all fall in one. An earlier
version of this section said **"That is a guarantee, and no part of it is a
clock"** about the whole of it, and that sentence was wrong:

* **The guarantee.** A run never writes a key it did not name. Into a readable
  document every other record is read, kept and written back; where the
  document is missing or unparseable the write creates one holding only the
  keys it names and rebuilds nothing from this run's snapshot. That is what
  stops the reviewer's failure, in which a whole document was rebuilt from one
  run's snapshot. It is not a promise that a record can never go backwards: a
  write that lands in the read-to-rename gap below is not in what this run
  merged, so this run's older copy of it goes back over the top and nothing is
  said. `--force`'s discard obeys the same rule: it carries
  the record it was decided about, not just the stage name, and deletes only
  while the document still holds that record, so a discard decided from one
  read at the start of the run cannot delete a record another run wrote
  afterwards. Driven.
* **The clock, and it is not the thirty-second one.** A record another process
  writes **between this run's pre-write read and its `os.replace`** is dropped
  rather than merged: a plain lost update, one read-to-rename gap wide, not
  configurable, and nothing detects it. The read-back after the write catches
  only a writer that lands *after* the rename; one that landed before it is
  simply not in what we merged, and the file we rename is a complete document
  without it. This is what `update_state()`'s own docstring means by shrinking
  the lost-update window rather than closing it. The mechanism is demonstrated
  by injecting a write into that gap; the race has not been won at shipped
  speeds, so read the window as *unmeasured*, not as *small*. For a
  **superseded** run the succession check closes this — but only from the
  moment that run has READ the replacement's `_run`, and a read that came back
  missing or unparseable is not that moment.
* **Not covered at all.** Whatever destroyed the file. If the state file is
  destroyed under a run, the records in it are destroyed with it, by whatever
  destroyed the file; the run logs `is no longer there` once, naming the path,
  and carries on writing a document that holds only what it wrote after that
  moment. The price is a **recomputation, not a wrong answer** — the outputs
  are still on disk and only their provenance is gone — and over-invalidating
  rather than under-invalidating is the trade this tool takes everywhere (see
  the signature section). The alternative, a dying run putting a live run's
  records back to its own older view of them, costs a results table that looks
  fine and is not.

**And what happens when the write cannot be made at all**, which is a
different failure from every one above: not another writer, but this run's own
filesystem refusing it — a full disk, an NFS `EIO` or `ESTALE`, a permissions
flip, a Windows sharing violation from an indexer. Both halves of a state write
can meet it: the read before it, and the write itself, which is a temp file,
its bytes, its close and a rename. ENOSPC can only meet the second of those —
a read needs no blocks — which is why the arm that declines an unreadable
document could never have covered the errno this section leads with.

* **It is declined, and it never raises.** The pre-write read gets 2 attempts
  (`STATE_READ_TRIES`) before it is called unreadable — one immediate retry,
  no sleep, which is the whole remedy for an `ESTALE` the next `open()`
  revalidates and for a sharing violation that clears when the other holder
  closes. A write that fails is
  retried through the same loop that redoes a lost merge, and each attempt
  re-reads and re-judges ownership, so a write that lands on the last attempt
  lands under the answer derived from the read before it. When the attempts run
  out the write is declined: nothing is written blind, and nothing is held for
  a later write either.
* **A stage that has already succeeded is never failed by it.** That is the
  rule this part exists for, and it was being broken two ways: an `OSError` out
  of the state write left `finish()` and killed the run — at hour thirty, with
  the stage's output already complete on disk — and in the drain loop, where
  `finish()` runs inside a broad `except`, the stage that had *succeeded* was
  reported as the stage that failed and the run exited `1`. A run whose records
  are all refused now does its work and exits on its stages' verdict.
* **The run says what it lost, twice.** Once per key when it happens, naming
  the stage and quoting the error, because an interrupt exits by a door that
  reaches no end-of-run report; and once as an account before the run exits,
  which names the records that never landed and says which of the two residues
  the next run will act on — a stage already recorded `running` is
  **recomputed**, a stage with no record at all is **adopted**, with the
  warning that says outright it cannot tell a finished file from an interrupted
  one. That second one is the one to look at: those outputs were produced while
  the filesystem was refusing writes. `--force --only <stage>` redoes one,
  `--no-adopt` refuses the adoptions wholesale, and there is no supported way
  to write a record by hand.
* **A declined `_run` leaves the same trace a `kill` does.** `final_status`
  stays `running`, which a console reads as a run that never ended, and
  `last_seen` stops advancing. The heartbeat says so out loud rather than going
  quiet — that warning had never been able to fire on a read error, because it
  counted exceptions and a declined write raises none.
* **The log file is part of this and not beside it.** `results/metaannot.log`
  is on the same filesystem as the state file, so the same ENOSPC stops the
  line that explains the declined write. A log file that will not take a line
  now costs the line: it is said once that the rest of the run's log is on
  stderr only, later lines are still attempted so the log resumes if space is
  freed, and **stderr itself is left raising** — a log nobody is watching live
  and the channel `tmux` keeps are different facts.

What none of that recovers is the record itself. A write that cannot be made is
not made, so the cost is the same one the rest of this section trades in: a
recomputation, or an adoption with a warning, never a wrong answer.

`_run` is the one key held to a stricter rule, because it is not a record of
work but a **claim about who owns the directory**. A run writes it only on
positive proof of ownership: a lock it reads and finds is still its own, or one
it could not read, which is no evidence that anything changed hands. A run
whose lock has gone **vacant** writes no `_run` at all from that moment — no
further `last_seen`, no `final_status`, no `finished`. That vacancy is exactly
what a replacement that took the directory and then finished leaves behind, and
a state file naming the old run the owner is read by the next run, which adopts
it as its predecessor and admits its writes for the rest of its life.

**And that gate costs the run its own last word, which is the price of having
it hold at all.** The version this replaced declined only the write that would
*create* the state file, so it never fired for a run that had recorded a stage
— and every real run records a stage, which creates the document one call
earlier. A run whose lock and document were both removed exited `0` leaving its
own `_run` in a file it had just made. Now it does not: it records the stage,
logs one line saying it can no longer prove it owns the directory and is
writing no more `_run`, and its verdict is lost. **No ordinary run pays that.**
The lock is released by an `atexit` hook that runs *after* `stamp_run()`, and
`kill` releases the lock and exits without stamping at all; driven end to end,
the ownership read answered *ours* at the final stamp of a clean run and of a
Ctrl-C'd one. The run that does pay it is one whose lock an operator, a
tmp-reaper or a remount removed under it, and on a console that run reads
afterwards as one that never ended.

Each stage record also carries the `run_id` that wrote it, so a document two
runs have both written into can be read rather than guessed at.

The second was the stage output itself. The executor waits for a stage that is
already running, and `atomic_out` renamed its result into place at the end; if
the replacement had meanwhile finished that same stage and recorded it `"ok"`,
the old run's output would land **under the new run's valid signature**, and
the run after that would report `cached` and read it. The rename now asks
whether this run still owns the directory, and on proof that it does not a
**declared** output is *parked* beside its target as
`.superseded.<stem>.<run_id><ext>` rather than renamed over it — nothing is
deleted, and the log says where the work went. Declared is the whole of it: the
check lives in `atomic_out`, so a file a stage writes some other way is not
covered, and for several stages the declared entry is a `.done` sentinel rather
than the table beside it. That check is an in-memory flag plus, at most, one rate-limited read
of the state file; it renames on every doubt, so a transient error can never
abort a stage that is legitimately finishing. A stage that *starts* after the
handover cannot rename at all; one already running when the directory changed
hands is only **noticed**, and the difference between those two words is the
whole of what this buys. Read them separately:

* **The guarantee.** A stage that *starts* after the handover cannot rename.
  `mark_running` is a merged write, the succession check reads the document
  before the stage begins, and the run is stood down there. No timing is
  involved. (A document that is missing or unreadable at that instant gives
  the check nothing to read, and that stage falls into the next bullet.)
* **The clock.** A stage *already running* when the directory changed hands is
  noticed at whichever comes first of the next heartbeat tick and the fallback
  probe — `min(heartbeat_s, STATE_PROBE_S)`, thirty seconds as shipped. The
  fallback probe is not what does it on a default run: `seen` is refreshed by
  every merge read, a heartbeat tick *is* a merge read, and the two intervals
  are the same, so the probe essentially never opens the file and the
  heartbeat's own write is the detector — essentially, because a tick that
  lands a few milliseconds late does let one probe through, which costs a read
  and never a wrong answer. The probe covers the runs the heartbeat does not —
  `heartbeat_s: 0`, or an interval longer than the probe's. Either way,
  **inside that window nothing detects the handover at all**: a takeover that
  completes in under a second is caught by nothing, and that stage's output is
  renamed over the live run's. Narrowing the window means lowering
  `heartbeat_s`, at one small file rewrite per interval.
* **What the check never sees.** It lives in `atomic_out`, so it reaches only
  outputs written through it. `hhblits` renames its per-query `<id>.hhr`
  itself, `esmfold` appends `plddt.tsv` and writes `esmfold_failed.tsv` in
  place, and `diamond`, `hhblits` and `esmfold` each create a `.done` sentinel
  with a plain `open()`. A superseded run still puts all of those into the
  live run's directory, with no check and no warning.

And one case is beyond any ownership check whatever — though not the one it
looks like. A run `SIGKILL`ed between its rename and its record leaves the
`"running"` record `mark_running` wrote *before* the stage started, and the next
run reads that and recomputes: the solo kill is covered, by over-invalidating.
What is left is the case where the record that **survives** belongs to a
different run than the **bytes** do. A marks `pfam` running; the operator
`--force-unlock`s; B recomputes `pfam` and records it `ok` over A's `running`;
A's rename then lands on top of B's file inside the clock above; A is killed and
writes nothing at all. B's record, A's bytes, a signature that agrees, and
`cached` over an output nothing was ever told about. The same shape needs no
kill whatever: a table copied in by hand, a parked `.superseded.*` moved back,
an over-broad copy back from the GPU box. `signature()` hashes inputs and config
and never an output, so nothing else here is capable of noticing any of it.

**So `--force-unlock` on a run that is still alive remains unsupported, and
none of this changes that** — what changed is that it now says so. Where the
lock names a pid **on this host that this host can see running**,
`--force-unlock` refuses with a message naming the pid, the host, when the run
started, when it last stamped `_run`, which stages it has recorded running, and
the `ps -p` to run; taking the directory anyway needs `--force-unlock-live` as
well. The refusal fires on proof and nothing else, so the case the flag exists
for — a stale lock from another node of the array, which this host can never
disprove — is taken over with the single flag exactly as before. It is what
`CLAUDE.md` rule 4 and the troubleshooting table already say: find out what the
other process is first. A run that is genuinely gone writes nothing, and none
of the above applies to it.

**So the file is dated against its own record, and the run says so when the two
disagree.** On the `cached` branch, and for an `ok` record only, `decide()`
compares each declared output's mtime against the `finished` stamp of the record
that describes it. **Nothing else happens.** The verdict is still `cached`, the
stage is still reused, no record is rewritten and nothing is deleted —
recomputing would rename over the one artefact that shows anything happened, and
would spend hours of InterProScan or ESMFold on evidence that has legitimate
ways to be wrong. `--force --only <stage>` is how to act on it, and the message
says so.

**It is one `WARN` per run, not one per stage.** Whatever the directory holds,
the run says this once: a heading, then a line per stage naming the file, its
size, when it was last written and the run whose record it is, then the reading.
*"cached, and the record is no longer about the file that is there."* The size of
that report grows with the number of late stages and its **count does not**,
which is the whole of why it is collected rather than said where it is found: a
report whose line count grows with the directory is the report that gets turned
off, and it would be turned off by *following the advice in it* — one
`--force --only <stage>` on a copied directory re-records one stage and leaves
every other one to be re-reported on every run for ever. The line arrives at the
end of the run rather than at the moment the stage was reused, and that is the
price: this check changes nothing about the run either way, so it is something to
read afterwards.

Read it in one direction, because that is the only direction it holds in:

* **A disagreement is evidence, not a verdict.** It says the record is not about
  that file. It cannot say what wrote it, or whether what is there is any good.
* **Agreement proves nothing at all.** An overwrite that landed before the
  record, or one that kept the file's timestamp, leaves no trace here — nor does
  anything inside `OUTPUT_STAMP_SLACK_S` of the record, which is the hole the
  slack buys. The slack is derived rather than chosen: `strftime` floors a stamp
  to the whole second and a filesystem's timestamp granularity can be coarser
  still, so with no slack at all every healthy stage of every resume would be
  reported. It is deliberately far smaller than `min(heartbeat_s,
  STATE_PROBE_S)`, the window a superseded run's rename has to land inside.
* **Where several stages are late at once, the cause is not knowable and the
  line does not claim one.** A results directory copied with `cp -r`, unpacked
  from an archive or restored from a backup has every timestamp in it rewritten
  at one moment and looks exactly like this — and so do that many separate
  replacements. There is nothing in a timestamp that separates the two: measured
  on this project's own results directory, a `cp -r` puts every output inside
  twenty milliseconds of every other, and so does overwriting the two outputs
  that a `--only pfam dbcan` run leaves. So the line names both readings, says
  that rebuilding stage by stage would spend hours of compute for nothing if it
  is the first, and leaves the choice with whoever knows the directory's history.
  (`rsync -a`, which the two-machine workflow in `TUTORIAL.md` prescribes
  everywhere, preserves mtimes and reaches none of this.)
* **Declared outputs only, which is narrower than it sounds.** `diamond`,
  `hhblits` and `esmfold` declare a `.done` sentinel, so for those stages a
  sentinel that dates cleanly says nothing whatever about the per-database
  tables, the per-query `<id>.hhr` files, `plddt.tsv`, `esmfold_failed.tsv` or
  the per-protein PDBs beside it. (The tables and the PDBs *are* covered by the
  ownership gate above — they go through `atomic_out`. They are outside *this*
  check because they are not declared.) The sentinel is still the only witness
  those files have anywhere: a superseded run that completes one of those stages
  re-touches its `.done` with a plain `open()` that no gate covers at all.
* **`adopted` records are never dated.** `finished` on one of those is when this
  box *noticed* the file, not when the box that made it wrote it, and `rsync -a`
  preserves the source mtime. Dating the GPU hand-off against it would report
  that hand-off on every run, for ever.
* **`--force` compares nothing**, having already decided to recompute.
* **A stamp that names no single instant is declined rather than dated, and the
  run says which record.** A `finished` stamp is a naive local time, so the hour
  a zone repeats when it leaves summer time is two moments carrying one text and
  the hour it skips entering summer time is none at all. Dating one of those
  wrong by an hour is far outside the slack, and a stage that went unjudged in
  silence could not be told from one that dated cleanly, so it is named.
* **One case reports itself on every run, and nothing silences it.** A stage
  this box ran and recorded `ok`, over which the operator later `rsync`s a
  newer output from the GPU box, is a record that really is no longer about the
  file that is there — so it is reported, correctly, and nothing re-records it.
  The remedy the message names, `--force --only <stage>`, would recompute a GPU
  stage on the wrong machine. Re-recording a stage *without* recomputing it is
  a change of its own and is not in this one.
* **`emapper` is the one declared output a kill can genuinely truncate**: its
  live branch lets `emapper.py` write `eggnog/emapper.emapper.annotations`
  itself, with no temp file and no rename. A timestamp says nothing about
  truncation, and this check does not claim to.

Where the timestamps themselves cannot be read, the run says *that* instead,
once, and compares nothing further — and there is exactly one such case, because
it is the only one that is a **measurement** rather than a reading. A filesystem
whose clock leads this machine's — measured against a file the run has just
written — makes every output look newer than its record for as long as the mount
is skewed, and two clocks compared against each other say nothing at all.

**A `--dry-run` cannot take that measurement, and says so rather than leaving it
to be inferred.** The skew is measured on a file the run has just written; a dry
run writes nothing, on purpose, and a plan check that created files in order to
measure them would stop being one. So a dry run — which calls `decide()` for
every stage, and is therefore the one command that reports a whole directory
without running anything — reports it and adds that the clock question was not
asked at all.

None of this needs a migration, and nothing in the state file moved for it: the
check is built out of `finished`, which every `ok` record has always carried, so
a results directory written by an earlier build is read exactly as one written
by this one, and a record this build cannot parse at all is compared against
nothing and reported as nothing.

### `describe`: what this build is, as JSON

```bash
python metaannot.py describe                    # a human summary
python metaannot.py describe --json             # the machine-readable contract
python metaannot.py describe --json --config config.yaml
```

`--json` emits one object: `default_config` (the whole config with its
defaults), `stages` (each stage's `enabled` flag, `deps`, `outputs` and the
config `keys` that decide whether it recomputes), `requirements` (the
tool-and-database half of `doctor`, as data), the config vocabulary —
`path_keys`, `db_path_keys`, `replace_blocks`, `freeform_keys`, `retired_keys`
— and `paths`, the files a watcher polls. Those paths are absolute: as
configured when you pass `--config`, and against the current directory
otherwise.

`requirements` is **not** everything `doctor` checks. It is the tools and
databases, with how to obtain each. `doctor` additionally checks the inputs
(`proteins_faa`, `quant_table`, `gff`, `contigs_fna`), the
`emapper_precomputed` files,
unrecognised and retired config keys, whether each DIAMOND database is usable
rather than merely present, the CUDA probe for `structure`/`topology`, the
`tmt:` block against the plexes actually present, the `manifest` and whether
its runs map to the quant columns, the `taxonomy` inputs, the `resources`
split (threads and RAM per concurrent stage, and whether eggNOG gets
`--dbmem`), and the `R` packages the report and the object need. Those are
`doctor`'s own section headings — `== inputs ==`, `== precomputed emapper ==`,
`== config ==`, `== gpu ==`, `== tmt ==`, `== manifest ==`, `== taxonomy ==`,
`== resources ==`, `== R ==` — and the only two that `requirements` does feed
are `== tools ==` and `== databases ==`. A preflight screen built on
`requirements` alone can be green for a config `doctor` fails, so run `doctor`
too.

It is versioned: `describe_version` changes when a key is removed or its
meaning changes, so a reader can refuse a shape it does not understand instead
of guessing. The exact key set is pinned by a test, which is what makes that
promise true rather than aspirational.

It exists so that anything driving metaannot from outside reads a contract
rather than importing private functions or scraping `--help`, which is how a
front end drifts from the engine and starts lying. The stakes are concrete:
most stages hash their database path **by value**, so a wrapper that rewrites
one cosmetically — `D:\db\Pfam-A.hmm` into `/mnt/d/db/Pfam-A.hmm`, the same
file — invalidates those stages and restarts InterProScan. Each stage's `keys`
is the list that says which of its settings do that.

Half the answer is static and half is a probe of the machine it ran on:
`default_config`, `stages` and the versions are the same everywhere, while
`requirements[].ok` is `shutil.which` and `os.path.exists` on `host` at
`generated`. Read `ok` as a fact about that machine, not about metaannot.

### `doctor --json`: what this machine can do with this config

```bash
python metaannot.py doctor --json --config config.yaml
python metaannot.py doctor --json --config config.yaml --install-plan i.sh
```

The sibling of `describe --json`, and deliberately shaped like it: the header
is spelled key for key the same way (`metaannot_version`, `signature_version`,
`generated`, `host`, `config_path`, plus `doctor_version` and
`describe_version`, so one call tells you the level of both contracts), and
`requirements` is `describe --json`'s `requirements` array element for element
— the same function, not a re-shaping of it. `cmds`, `size_gb`, `disk_gb` and
`note` live there and nowhere else; every check points back into it by
`requirement_id`.

What `doctor` adds is `checks`: **one entry per line the printed report
prints**, in the same order, under the same section headings. `detail` is the
printed sentence rather than a second wording of it, so the terminal and a
front end cannot come to disagree about the same config. `sections` carries the
`== ... ==` headings so a table reproduces the grouping instead of inventing
one.

Each check names **which enabled stages it kills**:

| field | means |
|---|---|
| `status` | `ok` / `warn` / `fail` / `skip`. The printed report has a fifth mark: `doctor_mark()` prints `MANUAL` rather than `MISS` for a `fail` whose remedy is one doctor will not fetch, so `status: "fail"` maps to `MISS` **or** `MANUAL` and the other three are one for one with `OK` / `WARN` / `-` |
| `blocks` | enabled stage names that die without this — joins to `describe --json`'s `stage_names` |
| `blocks_commands` | `run`, `report` or `object`: a whole command refuses or dies, and no stage does. A `fail` can name **only** these — a missing `limma` is `blocks: []`, `blocks_commands: ["report", "object"]` |
| `degrades` | enabled stages that run anyway and produce less |
| `expect` / `found` | the right kind of thing, and what is actually there |
| `depth` | how far the check looked: `config`, `existence`, `kind`, `header`, `parsed`, `probe`. A ratchet enforced in both directions — `header`/`parsed` must carry a `caveat`; `existence`/`kind` must carry a `found` that is evidence of a stat (an `unset` one is not); and `config`, which is a claim that **nothing** outside the config was consulted, forbids a real `found` — so a row can neither hide a read nor claim one it did not make |
| `remedy` | `auto` / `manual` / `config` / `input` / `none` — which affordance to offer |
| `fails_reason` | `stage_or_command_dies` or `setting_ignored` |

**The exit status is derivable from the document and nothing else.**
`verdict.rule` says it in a sentence: 1 if and only if some check has
`status: "fail"`, otherwise 0. `ok`, `exit_status` and `verdict.fails` all come
off that one expression, and `$?` is that expression too, so a caller reading
the JSON and a supervisor reading the status cannot disagree. A check fails
exactly when something the config asks for dies on it, and the row names what:
an enabled **stage** in `blocks`, or a whole **command** in `blocks_commands`,
which is a failure with no stage in it at all. `status == "fail"` if and only
if `fails_reason` is set, and `fails_reason == "stage_or_command_dies"` if and
only if `blocks` or `blocks_commands` is non-empty — the value is named for
both halves because it was named for one, and three of the four places the
rule was written down then copied that name's sentence and dropped the command
class entirely. There is **exactly one** declared
exception and it names itself: an unrecognised config key is ignored by `run`,
so no stage dies, but the setting the config asks for is silently not in
effect and the run answers a different question than the config asked; that
check fails with `fails_reason: "setting_ignored"` and an empty `blocks`. It
stays a failure deliberately — `doctor` has always exited 1 on a misspelled
key, and demoting it to a warning would change the exit status for a real
class of config, which is not a call to make inside a contract change. Every
other row in the document obeys the rule without exception: if it fails, an
enabled stage or a whole command dies on it.

`scope.statement` is a literal sentence the engine authors and a front end
renders verbatim:

> A check that passes says the thing exists, is the right kind of thing, and is
> not empty. It does not say the thing is correct inside: doctor does not parse
> FASTA files, quant tables or sequence databases FOR THEIR CONTENT. Four
> checks do read into a file to answer their own question - the DIAMOND
> usability check, the TMT plex annotations, the manifest and the quant table's
> header line - and each says which on its own row, in depth and caveat. An
> input that passes every check here can still fail the stage that reads it.

That is the boundary, and `depth` is what keeps it honest per row rather than
as a promise in a docstring: four checks do read past a file's existence — the
DIAMOND usability check, each plex's TMT annotation file, the manifest itself
(`read_manifest`, in full, because the mapping check has no runs to map
without it) and `pd.read_csv(nrows=0)` on the quant table's header line — and
each of those carries a `caveat` naming exactly what it opened. The DIAMOND one is the
reason the caveat is **derived and not declared**: what it reads depends on
the host. With `diamond` on PATH it reads the database header
(`diamond dbinfo`, `depth: "header"`); without it — the ordinary state the
first time anyone runs `doctor` — it falls back to the source FASTA beside the
database and reads up to 200,000 records for a length profile, plus 200
deflines for the motif-seed signal, which is `depth: "parsed"`. A zero-byte
`.dmnd` is refused on `os.path.getsize` alone and reads nothing at all
(`depth: "kind"`). Nothing in the shape can carry a parse result; there is no
`rows`, no `columns`, no record count and no id list, which is what stops this
growing into a different program.

**Every one of those reads goes through one gate, and `doctor` opens nothing
the gate has not already opened.** `regular_readable()` proves a path with an
`os.open(O_RDONLY | O_NONBLOCK)` and an `fstat()` on the descriptor, and
`_deep_readable()` wraps it for `doctor`, handing back the row's sentence
alongside the verdict. That is not fastidiousness about stats. `os.path.exists()`
is **true** for a FIFO and a read-only `open()` of a FIFO with no writer does
not return until somebody writes: a FIFO at a plex's annotation file used to
make `doctor` wait for ever — no document and no exit, which is worse than any
traceback, because a caller waiting on the process has nothing to time out
against either. `O_NONBLOCK` is
what makes that open return instead of wait; `fstat()` on the descriptor is
the only "is this a regular file" test that cannot be raced by the path
changing underneath; and opening at all is the only test a permission cannot
lie to, since `os.path.getsize()` answers happily for a file nobody may read.
Every deep check keeps a wide `except` behind the gate as well, because the
whole history of this class is that the next path state was one nobody had
thought of. **No state of any path `doctor` is pointed at can stop it printing
a document or stop the process exiting** — the six configured input keys and
the TMT tree's plex directories, level files and annotation files, through
absent, empty, directory, dangling symlink, FIFO, socket, device node,
unreadable directory and unreadable file — and a test walks all of them, with
a timeout, so a future hang fails a test rather than wedging a suite.

"The right kind of thing" is **derived, never assumed**. `expect.kind` is
`file` for every quant format but one, and `dir` for `fragpipe_tmt`, which
reads the run directory holding the per-plex folders;
`expect.derived_from: ["quant_format"]` names the key that decided, so a ninth
format changes the value and not the schema. Its full vocabulary is `file`,
`dir`, `on_path`, `probe` (the host was asked — CUDA, an import,
`requireNamespace`), `setting` (a claim about the config and nothing else) and
`r_package`. `expect.members` carries the siblings the engine's **own** test
requires beside the thing named — the four `hmmpress` files for an HMM
library, `.dbtype` and `.index` for a Foldseek target, the `nodes.dmp` inside
a taxdump directory — and is empty where that test names none, as for
`hhblits_db`, whose `_prefix_exists()` accepts any non-empty sibling of the
stem. `found.kind` names what is really there — `file`, `dir`, `empty_file`,
`empty_dir`, `symlink_broken`, `absent`, `unset`, `other` (present, and
neither a regular file nor a directory: a FIFO, a socket, a device node) and
`unreadable` (there, and this process may not read it — a directory with no
read bit, or a file whose mode or ACL refuses an `open()`) — so a consumer
never has to do arithmetic on `bytes`, and a dangling symlink is reported as
one instead of as a plain absence. `unreadable` is reached for a **file** as
well as a directory, which it was not until now: `os.path.getsize()` is a stat
and answers happily for a mode-000 file, so four rows that branch on
`kind == "file"` used to report `ok` for one while `run` died on "Permission
denied" in the stage that opened it. `other` is not a curiosity: `os.path.exists()`
is **true** for a FIFO, so `run` does not refuse one and every row that read
`absent` for it told an operator the opposite of what happens. What `run` does
with one is its own section below; what matters here is that the two commands
answer differently on purpose, and `doctor` is the one that may never wait.

**`other` is three things, and `found.other_kind` says which.** It is `fifo`,
`socket`, `char_device` or `block_device`, and `null` for every other `kind`.
It is a new KEY and not a new `kind`, so `DOCTOR_VERSION` does not move and a
consumer that has never heard of it keeps every answer it had. It exists
because a sentence derived from `kind` alone is written for one of the three
and printed over all of them: a UNIX socket at `proteins_faa` published "dies,
but not at once ... A FIFO in particular is not refused on sight ... waits
`fifo_wait_s`", and a socket cannot be opened as a file at all — `os.open()`
itself fails on it. (The **errno** is deliberately not published any more: this
was written down as `ENXIO` and measured as `EOPNOTSUPP`, and which one you get
is the OS's business. What decides who dies is that it is an `OSError`, which
nothing downstream catches.) The verb,
the FIFO paragraph and the phrase naming what is there are all keyed on it
now.

**Two questions, not one: WHEN the open ends and WHAT it ends with.** They
part company on a character device, which answers at once — so the "when"
predicate is true of it — and raises a `StageError` from the `die()` at the
bottom of `_open_for_read()`, which `peptide_features()` catches and re-reads
around. `input:manifest` is the row where that matters, because it is the row
that names an exception and computes `blocks` from it: on a char device it
published "it dies with an `OSError` … every stage that opens a quant table on
this format goes with it" beside its own `blocks: ["join"]`, and `blocks` was
the right half — driven, the taxonomy stage fell back to the peptide-only
reader and finished. Both are computed from one predicate now, which also
corrects the opposite error: a **socket** and a **block device** share the
`other` bucket with the FIFO and end in an `OSError`, so they cost every
reader and `blocks` used to under-claim them.

**`found.kind` has a tenth value and it is `null`.** A row can carry a
non-null `found` whose `kind` is `null`, and it means *this build did not
compute it*, never *nothing is there*: every `requirements` row, the CUDA probe
and every row in the `R` block answer with one boolean — `have()`,
`_exists()`, `_pyhas()`, `torch.cuda`, `requireNamespace()` — which cannot
tell an absent database from a directory or an unpressed library. In the
default document those are the **majority** of the rows, so a consumer
switching on `found.kind` needs a `null` arm, and `found.present` is what
carries the answer there. The `db:diamond:<tag>:usable` rows are **not** in
that list and were wrongly named in it in three places: they are built from
`_found()` and carry a real kind, which makes them the document's own
counterexample to the sentence they were named in. The ordinal above is
derived from `DOCTOR_FOUND_KINDS`, and the families that really carry a null
are `DOCTOR_NULL_KIND_ROWS`; a test recomputes the first and drives a document
to check the second.

`--json` **replaces** the printed report on stdout, exactly as
`describe --json` does, and `log()` already writes to stderr, so the stream
stays clean even when the config has typos in it. `--install-plan FILE` still
writes its file — it is a side effect the operator asked for, not output — and
`install_plan` in the document records the path, what went in and what was left
out. **`--fix` with `--json` is refused by argparse, which exits 2**, not by
`die()`, which exits 1 and so could not be told apart from "problems found":
`--fix` reads stdin for its confirmation and streams progress to stdout, and a
document written while a 123 GB download is half done describes nothing. The
honest loop is `doctor --json`, decide, `doctor --install-plan`, run it,
`doctor --json` again.

`totals` carries both sizes as raw floats, plus the two things that make them
readable: `counted` / `not_counted` (MANUAL items are in neither total, so the
figures are a floor rather than the whole job) and `unsized` — the ids where
the `0.0` in `requirements[]` means "nobody estimated this", not "free". The
printed report hides anything under 0.05 GB, which is how a large Java
distribution came to look like no download at all.

It is versioned the way `describe --json` is: `doctor_version` moves when a key
is removed or its meaning changes, never when one is added. One clause is
added, because this document has enums a consumer switches on — **adding a
value to a closed enum is a meaning change and is a bump**. The closed sets are
`status`, `remedy`, `fails_reason`, `depth`, the `blocks_commands` vocabulary,
`found.kind` and `expect.kind` — seven, and the last two joined the list after
they were found outside it, which is the difference between a contract and a
suggestion: this page tells you to switch on `found.kind` rather than compute
`bytes > 0`, and a `fifo` or an `unreadable` appearing there would not have
been a version bump. Both of those states now have a name **in** the set, and
`fails_reason`'s first value was renamed in the same pass, neither of which is
a bump only because `doctor_version` 1 has not shipped: it and this page are
in the same unreleased change set, so there is no consumer to break.
`finding`, the stage names in `blocks`, and every `detail` and `caveat` are
open: a consumer that meets an unfamiliar value there keeps `status` and
renders `detail`.

## The console: watching a run without touching it

`console/console.py` is a read-only watcher for results directories. It is one
stdlib-only file — Python 3.9 or newer, nothing to install, deployable by `scp`
— and it serves one page: a list of the directories it watches, and per
directory a stage table read out of `.metaannot_state.json`, a log tail, the
lock, and what the run record says about itself. It is `CONSOLE_VERSION`
`0.1.2`, and that number is deliberately not `__version__`: the console and the
engine ship in one repository but they are two programs with two audiences, and
tying their versions together would mean either lying about one of them or
bumping a number nobody asked about.

**The NEXT rows carry the engine's dispatch order, taken from the engine.**
A stage with no record and no unmet dependency is ready, and several are ready
at once; which one actually starts is `stage_priority`, the cost rank and then
the stage's seconds on the release's reference run, both of which
`describe --json` publishes for exactly this purpose. The console reads both
and says which ready stages are ranked ahead of which — it does not reorder the
table, because a dependency graph read out of order is harder rather than
easier. It says nothing about how many workers are free or whether a GPU stage
will be deferred, because it cannot know either; the note under the table
carries those two unknowns once. Reading only the rank was not enough and the
page was wrong for a while because of it: eleven stages declare the same rank,
so the page annotated the hours class in its own table order and told an
operator that several ready stages were ahead of InterProScan when none were.

**It never writes a byte into a results directory.** Not a lockfile, not a
cache, not a temp file, not a log line. That is not politeness, it is the
entire reason it is safe to point at a three-day job that is already running:
every reader opens `O_RDONLY`, the console never `chdir()`s into a watched
directory so not even a core dump can land there, there is no `do_POST`, and
the one file it does create — a `flock` that stops two consoles fighting over
one socket — lives in the console's own runtime directory and is refused
outright if you try to put it anywhere near a watched tree. The rule is
enforced rather than intended: `tests/test_console_contract.py` proves it by
AST scan and by snapshotting a results directory around every route.

**It never imports metaannot either.** Stage order, dependency edges, the
`_run` key and the names of the files it polls all come from
`describe --json`, shelled out once at startup and cached. This is the
difference between a front end that stays true and one that starts lying: most
stages hash their database path by value, a console that hard-codes what a
stage is named or where a state file lives drifts the day either changes, and
the drift is invisible until somebody reads a stale answer off a page that
looks authoritative. `docs/gui-design.md` is the design record for that choice.

Run it on the host that writes the results directories — it reads local files,
so it has to be where they are:

```bash
# on the machine that runs the pipeline
python3 console/console.py --root /data/projects

# from your workstation, once
ssh -N -L 8080:/run/user/1000/metaannot.sock lab-fedora
# then open http://localhost:8080/
```

It binds a **mode-0700 UNIX socket**, never a TCP port. File permissions are
the whole access control, by decision: a `127.0.0.1` listener on a shared lab
server is reachable by every other account on that box, and a token would leak
into `ps`, shell history and the URL bar. A UNIX-socket forward target needs
OpenSSH 6.7 or newer on your side. The console prints the exact `ssh -L` line
for the socket it actually bound, so you can copy it rather than reconstruct
it, and any free local port will do.

The flags are few. Results directories are named as bare arguments or with
`--project` (repeatable, and identical); `--root` (repeatable) scans a tree
three levels deep for them instead — deep enough for
`<root>/<dataset>/results`, shallow enough not to wander into a 200 GB Foldseek
tmp directory — and re-scans every 30 seconds, so a directory created later is
picked up and keeps the number it was given. `--socket PATH` overrides where to
bind. Leave it off **and `$METAANNOT_CONSOLE_SOCK` unset** and the console picks
the first private directory it can find from `$XDG_RUNTIME_DIR`,
`/run/user/<uid>`, `~/.cache/metaannot-console` and
`/tmp/metaannot-console-<uid>`, and puts `metaannot.sock` in it. If
`$METAANNOT_CONSOLE_SOCK` **is** set, that is the flag's default value and the
search never happens: the variable is treated exactly as if you had typed
`--socket` yourself, refusals included. That matters because the two paths are
not equivalent — the console creates the last component of a directory of its
own choosing, and creates nothing at all for a path you named, so a
`$METAANNOT_CONSOLE_SOCK` pointing into a directory that is not there yet is a
refusal rather than a `mkdir`. Either way the directory has to be one no other
account can reach, and a socket path inside or under anything the console is
watching is refused outright — but a *candidate* that fails the privacy test is
only passed over, with a line on stderr, for the next one in the list.
`--metaannot PATH` is the engine to ask for the contract, defaulting to the
`metaannot.py` next to the `console/` directory; `--python PATH` is the
interpreter that runs it, defaulting to the one running the console.
`--interval S` is the fastest the page will poll, default 3 seconds, floor
0.5 — `0` used to be accepted and turned every open tab into a fetch loop
against an NFS mount.

**It will not tell you a run is dead.** That is a decision, not an omission.
The heartbeat in `_run.last_seen` is advisory — one failed write ends the
heartbeat thread while the run carries on — and the engine's own rule is that
unprovable means alive. So the console says how long it has been since the last
sign of work, names what it is reading, hands you the `ps -p <pid> -o
pid,etime,stat,args` line for the pid in the lock file, and stops there. Where
the evidence really is informative it says so and no more: a `FATAL` as the
newest log line while the record still says `running` is the killed-mid-write
shape, and it is called that rather than called death. The verdict is yours,
and the action that follows from it is `--force-unlock` on the engine, which is
a person deciding.

Which is the other half of the same decision: **there is no button that acts.**
No `--force-unlock`, no launch, no config edit, no `do_POST` — the HTTP handler
implements `GET` and `HEAD` and nothing else. A watcher that cannot act cannot
act wrongly on a directory you care about, and that is what makes "point it at
the running job and see" a reasonable first thing to do rather than a decision.

The preflight checklist, the server-side directory picker, a `doctor --json` to
feed them and the config authoring described in `docs/gui-design.md` are later
milestones and are **not** in this release. This is M1, a watcher, and it is
deliberately the whole of it: if the pane is not useful, the loss is one file.

## What has actually been run

**On real data:** one label-free dataset, end to end — a FragPipe
`combined_peptide.tsv` plus its `.fp-manifest` plus a precomputed eggNOG table,
38,204 proteins, 122,278 features, 36 samples, 6 groups. Fifteen of the
twenty-one stages ran on that input: `emapper` (reuse), `pfam`, `dbcan`,
`diamond` over five databases (VFDB, MEROPS, CARD, TADB3, BAGEL4), `cluster`
(MMseqs2), `ncbifam`, `kofam`, `interpro`, `signalp` (SignalP 6), `tmbed`,
`esmfold`, `foldseek`, `integrate`, `finalise` and `join`.

Final bins: `1_ko_pathway` 47.6%, `2_ko_orphan` 28.2%, `3_annotated_no_ko`
19.0%, `3d_duf_only` 1.0%, `3s_structure_only` 0.4% (141 proteins), `4_dark`
3.8%. Dark rescue took 11,106 proteins dark on eggNOG alone down to 1,462 —
9,644 rescued, of which KOfam alone supplied a KO to 8,129 that eggNOG missed,
and 6,372 of the 11,106 got their KO from KOfam. The eggNOG-only figure
reconciles with the join described above: 17,377 proteins had no eggNOG KO, and
6,271 of them carried an eggNOG `PFAMs` entry, leaving 11,106.
ESMFold built 1,821 models, 993 of them passed the pLDDT gate, and Foldseek
returned 21,791 hits over PDB and AlphaFold Swiss-Prot plus 664 self-clustered
fold groups.

A FragPipe TMT run is in progress on a second dataset — 8 plexes, 88 channels,
75 biological samples, 74,051 features, 455,571 proteins, eggNOG coverage
98.4% — and a 3-plex subset of it has already completed end to end including
the report. So `quant_format: fragpipe_tmt` is not a paper path.

**Not on real data:** `smorf`, `context`, `hhblits`, `jackhmmer`, `unipept` and
`taxonomy` — all off by default — have still never run. So the profile
searches, the Unipept API and the peptide-LCA taxonomy comparison remain
unexercised outside their output parsers, and `3p_profile_only` has never been
produced by any run: it is fed only by `hh_hit` and `jackhmmer_hit`, so it
stays empty until one of those two stages runs. `3s_structure_only` **has**
been populated — 141 proteins, from the Foldseek search above.

**Report and R object: built from the real run.** Both were produced from the
38,204-protein dataset above on R 4.6.1 — the report knits to HTML with figures
and tables, and `object` returns a QFeatures of about 12.5 MB with linked
`peptides` and `proteins` assays. Nothing in the tool pins an R version; R 4.3
with pandoc 3.1 was simply the earlier validation environment. The object
script can still degrade in two independent places and each names the real
error rather than guessing: if `addAssayLink` fails you still get a QFeatures,
both assays present but unlinked, logged as `assay link not added (<the real
condition>)`; only a failure of the QFeatures constructor itself falls back to
a SummarizedExperiment, with the peptide assay in `metadata()$peptides`. Check
the class of the object you get back.

A test suite does ship, in `tests/` — `conftest.py`, `fixtures.py` and sixteen
test modules — so the plumbing described here is reproducible offline from what
you have; see **Tests** above. What does not ship is a benchmark script, and
the datasets themselves are not redistributable, so the run figures above
cannot be re-derived.

## Scale

Wall times are measured, and the tool measures them itself. Every stage that
actually **runs** records its duration as `seconds` in
`results/.metaannot_state.json` — also embedded in the R object as
`metadata(obj)$metaannot$state` — so your own run reports its own numbers. A
stage that is *adopted* or *skipped* on a resume records none.

One label-free run, 38,204 proteins, on a single workstation, with stages
running concurrently, so these are shares of a parallel run rather than
single-stage benchmarks. This is also the run the scheduler's dispatch order
comes from: each figure below is that stage's `order_s` in `describe --json`,
where it is used for its ORDER alone and never as a prediction of your run.
The stages this run did not enable carry a placement rather than a measurement
and the comment at `STAGE_ORDER_S` says which is which:

| stage | wall time |
| --- | --- |
| `interpro` | 2.84 h |
| `signalp` | 0.93 h |
| `tmbed` | 0.87 h |
| `kofam` | 0.49 h |
| `pfam` | 0.45 h |
| `ncbifam` | 0.40 h |
| `foldseek` | 0.02 h |
| `dbcan` | 0.01 h |
| `emapper`, `cluster`, `diamond`, `integrate`, `finalise`, `join` | under a minute each |

ESMFold is quoted separately because its cost depends on length rather than on
protein count: 1,805 folds at or under 478 aa took 1.6 h in total, while the 91
sequences above that machine's VRAM cliff were projected at 7.4 h on their own.
See **The length a card can actually fold**.

**They do not scale linearly.** On a 455,571-protein run — 11.9× the size —
the stages that finished came in at 14–30×, not 12×: kofam 14.3 h (29×), pfam
7.2 h (16×), ncbifam 5.6 h (14×), dbcan 617 s (15×), diamond 291 s (30×). The
cause is contention rather than size: at 38k a long stage rarely overlaps
another long stage, and at 455k every one of them overlaps every other for its
whole life. SignalP measured 11.4 sequences/s with the machine mostly to
itself and 2.0–3.2 sequences/s alongside tmbed and kofam — the same work,
three to five times slower. Doubling the linear estimate past ~100k proteins
is a fair planning rule and an optimistic one for the worst stage.

MMseqs2 is the exception: 455,571 proteins clustered into 49,347 families in
109 seconds, 7.6× for 11.9× the proteins, because clustering scales with
redundancy rather than with count.

Three of that run's stages had not finished when this was written, and are
deliberately not quoted rather than rounded: SignalP was 58% through after
30.7 h, InterProScan had been going 3.7 h, and TMbed had written nothing at
all in 30.7 h. See [TUTORIAL.md](TUTORIAL.md#resource-guide--and-how-to-size-a-run).

The hot paths are vectorised: protein-group explosion, bin assignment and
export scoring are array operations, not row-wise `apply`, and sequences are
not retained after the single FASTA pass, so memory scales with protein count
rather than total residues.

Quant tables are read by `read_delim_table`, which takes the delimiter from the
header line alone and hands the body to pandas' C parser. The old
delimiter-sniffing path (`sep=None, engine="python"`) is gone from every quant
read — roughly 18× slower and 4–6× the memory for no benefit, since the first
line already says which delimiter this is. It survives only for the Unipept
result file and a header-only peek in `doctor`. Reading is no longer the
dominant cost.

Scratch is **partly** cleaned up. The two trees that can reach hundreds of GB
are removed by the run itself: each Foldseek target search deletes its
`{results}/foldseek/tmp{i}` when the search returns, and self-clustering
deletes `{results}/foldseek/tmpc`, so scratch peaks at one target's tree rather
than the sum over targets — but that peak is real: against AFDB50 budget tens
to hundreds of GB of free space for the duration of the stage. A Foldseek that
is killed, or that exits non-zero, still leaves the tree it died in.

What survives a clean run is smaller and fixed in kind: `{results}/kofam/tmp`
is the largest leftover at 229 MB on a 38k-protein run;
`{results}/interpro/query.faa` is the sanitised copy of the whole proteome
InterProScan is actually searched against, ~14 MB at 38k proteins and scaling
with it; `{results}/cluster/tmp` and `{results}/interpro/tmp` are created and
never removed but MMseqs2 and InterProScan empty them themselves, so they
survive as empty directories; `{results}/foldseek/query_hq` is a hardlink farm
of the pLDDT-passing models, so it costs no real space. Budget a few hundred
MB, not hundreds of GB.

## Three things to check before believing any of it

1. **Gene calling.** Prodigal in meta mode has a hard floor of 90 nt and calls
   genes below about 100 aa with reduced sensitivity rather than discarding
   them, so without a smORF-augmented database bacteriocins, TA toxins and
   RiPPs are under-sampled rather than absent — and their absence from your
   results is not evidence of absence. The report does not print this as a
   blanket: it counts quantified groups at or under 100 aa and uses the
   "effectively never in the search space" wording only when that count is
   zero, printing the count and a softer warning otherwise. That count is a
   proxy for what the database contained, not a check of it, so still confirm
   whether yours has a smORF tier. Going further needs a **re-search, not a
   re-annotation**: `run.smorf: true` with `contigs_fna` set calls small ORFs
   into `results/smorf/smorf_proteins.faa`, which you must append to your MS
   search database and search the raw data against again yourself. No
   downstream stage consumes that file.
2. **Group and taxonomy conflicts.** `group_conflicts.tsv` and
   `taxonomy_comparison.tsv`. KO-less proteins are strain-specific, so theirs
   are the least trustworthy.
3. **Taxon confounding.** `taxon_intensity.tsv` as a covariate, or a protein
   that merely tracks its source organism reads as a regulatory finding.
