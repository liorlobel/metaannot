# metaannot

One file. `metaannot.py` contains the pipeline, the R report, and the config
template.

```bash
python metaannot.py init                        # write config.yaml
python metaannot.py doctor                      # what is missing, and how to get it
python metaannot.py doctor --install-plan i.sh  # a script to review, then run
python metaannot.py doctor --fix                # download and install, after confirming
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
(`tool_args`, `db.diamond`, `sources.diamond`, `diamond_weights`,
`effector_predictions`), whose keys are user-chosen and cannot be checked, so a
typo there is silent and simply has no effect. Proof-read those blocks by hand.
The `analysis:` block is **not** free-form and **is** checked: its keys are
exactly the Rmd's params, so `fdrr: 0.01` is reported rather than written into
the Rmd header as a spurious param while `fdr: 0.05` stays quietly in force.
Read the `== config ==` block of `doctor` before a run you intend to publish.

Do not paste a second top-level `run:` or `db:` key into a config that already
has one. `yaml.safe_load` keeps only the last mapping with a given key, so the
earlier block is discarded, every setting in it reverts to the default, and the
key check cannot see it because the dict has already collapsed. Edit the
existing block instead.

Requires python3 + pandas + pyyaml. R is needed only for the report; external
tools only by the stages that use them. `doctor` says which are missing.

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
up to five unmatched runs and four candidate columns. It does not print the
full mapping.

**A manifest is one row per raw file, not per sample.** A fractionated
acquisition therefore repeats a sample name across its fraction rows — which is
exactly how FragPipe denotes fractions, and FragPipe writes one quant column per
group. metaannot **collapses** those rows to one sample, logs how many samples
were split (`manifest: N sample(s) are split across multiple fraction files`),
and refuses only when the rows of one sample genuinely disagree about the design
(different `data_type`), because collapsing those would invent a sample that
never existed. It also refuses if two *different* manifest samples resolve to the
same quant column, which is a real ambiguity rather than a fraction.

TMT is the case that still needs a hand-written design: there the `experiment`
column is the plex rather than the condition, so the derived contrasts would be
plex-versus-plex — a batch effect, not a hypothesis. Skip the manifest-derived
design and write `analysis.metadata` by hand (sample, condition, and any batch
column), with `analysis.design_formula` to match. Note that the quant side
refuses isobaric input outright (below), so this only arises if you are
quantifying the channels elsewhere.

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
| `3p_profile_only` | rescued by profile-profile or iterative profile search |
| `3s_structure_only` | no sequence annotation; confident Foldseek hit |
| `4_dark` | no evidence |

That is the complete set; `annotation_final.tsv`, `bin_summary.tsv` and the
report's factor levels use exactly these seven strings. The tests are applied in
the order KO-with-pathway, KO, sequence annotation, DUF, structure, profile, so a
DUF-only protein stays in `3d_duf_only` however good its fold or profile hit is
— deliberate, since a DUF still names a family, and the run counts how many were
held back so the 3s/3p rescue numbers are read with that in mind.

**eggNOG's own `PFAMs` column counts as domain evidence.** Binning on the
`pfam` stage's `hmmsearch` hits alone meant that with `run.pfam: false` a
protein eggNOG had already assigned a domain to was reported as having no
evidence; on the UC metaproteome that was 36% of the dark bin. An eggNOG Pfam
whose only accession is a DUF/UPF lands the protein in `3d_duf_only` rather than
`3_annotated_no_ko`, on the same rule as an `hmmsearch` DUF. So expect a smaller
`4_dark` than earlier versions of this tool produced, for a join reason rather
than a biological one.

Bins 2 to 4 are what KEGG enrichment silently discards. Global KEGG maps
(01100, 01110, …) are excluded from the pathway test: a protein whose only
"pathway" is *Metabolic pathways* is not pathway-annotated in any useful sense.
In real data almost no protein has *only* global maps, so `2_ko_orphan` means
in practice a KO with no pathway map at all — a module-only or unmapped KO.

## Stages

```
emapper → pfam → dbcan → diamond → signalp → tmbed → cluster
        → ncbifam → kofam → interpro → smorf → effectors → context
        → integrate → jackhmmer → hhblits → esmfold → foldseek
        → finalise → unipept → taxonomy → join
```

A fresh `init` config turns on **six** of them: `eggnog`, `pfam`, `dbcan`,
`diamond`, `cluster` and `join`. `integrate` and `finalise` have no flag and
always run. Everything else is **off** — `topology` (SignalP 6 + TMbed),
`structure` (ESMFold + Foldseek), `context`, `ncbifam`, `kofam`, `interpro`,
`hhblits`, `jackhmmer`, `smorf`, `effectors`, `unipept`, `taxonomy` — because
switching each on is a decision about hardware, a licence or a large database
that is not ours to make for you.

So `python metaannot.py all` on a fresh config still needs Pfam-A, the dbCAN
HMMs, the DIAMOND databases and either eggNOG-mapper or an
`emapper_precomputed` table. Run `run --dry-run` and read the plan before the
first real run; `doctor` reports which tools and databases are missing for the
stages you left enabled.

### Added evidence

- **`ncbifam`** — NCBIfam/TIGRFAM HMMs. Cheapest coverage gain: one more
  `hmmsearch`.
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
  Note what that leaves: the DIAMOND stage searches targeted databases only
  (VFDB, MEROPS, CARD, TADB, BAGEL) with `--max-target-seqs 5`, so it is not a
  general homology search. The only general-reference search here is
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
- **`effectors`** — ingest for Bastion3/4/6, EffectiveDB, T4SEpp, SecretomeP.
  One generic reader (`{file, id_col, score_col, threshold}` — there is no
  per-predictor `weight`; every predictor over threshold contributes a flat
  `effector_prediction_weight`) rather than five bespoke parsers for tools
  that are web services. Matters because many
  bacterial effectors have no signal peptide and the `surface_or_secreted`
  filter misses them.

Foldseek now searches several targets (`foldseek_extra_targets` — PDB and
Swiss-Prot carry far better annotation than mostly-unreviewed AFDB50) and
clusters the unannotated structures **against each other**. Fifty dark
proteins sharing a fold is a much stronger signal than fifty singletons, and
needs no reference database.

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

Every route above is **label-free (MS1) quantification**, and every quant input
is treated as linear intensity.

FragPipe writes `0` for "not quantified", not for "measured as zero", so zeros
in a `fragpipe`/`fragpipe_peptide`/`fragpipe_ion` table are read as **missing**
by default and the run warns with the cell count. Summed as real zeros they turn
missingness into fold change. `zero_intensity_is_missing: false` restores the
old behaviour if you need to reproduce someone else's numbers.

### Limitation: FragPipe TMT output is not supported

Reporter-ion channels are not read. A TMT run does not write the `combined_*`
files listed above at all: it writes per-plex `TMTn/{psm,ion,peptide,protein}.tsv`,
the `tmt-report/` matrices (`abundance_*_MD.tsv`, `ratio_*_MD.tsv`) and a
TMT-flavoured `msstats.csv`. There is no row in that table a TMT user can
honestly pick, and pointing the tool at any of these files quantifies the wrong
column rather than failing:

- The per-plex tables are **refused**. Their reporter channels are named
  `Intensity <sample>` while the column selector keeps only columns *ending* in
  `Intensity`, so the run would reduce to the single MS1 precursor intensity,
  pooled over all channels, exit 0, and report **one sample**. The loader now
  detects that shape — reporter-style columns present, and the only match a bare
  `Intensity` — and dies naming the channels it found. This refusal is the
  intended behaviour; do not work around it.
- `quant_format: fragpipe` on `abundance_protein_MD.tsv` turns the metadata
  columns (`NumberPSM`, `MaxPepProb`, `ReferenceIntensity`, …) into sample
  channels.
- The `tmt-report/` matrices are already log2 and median-centred; they are read
  as linear and log-transformed a second time.
- `msstats_csv` on the TMT `msstats.csv` raises a raw pandas `ParserError` on
  unquoted commas in `Protein.Description`.

Do not use this tool for reporter-ion quantification until a TMT reader exists.
Isobaric support is out of scope for this pass.

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
- `razor` — FragPipe's own behaviour. Arbitrary here; useful for measuring how
  much it changes the answer.

`peptide_evidence.tsv` records per protein how many features were unique,
taxon-unique and dropped.

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
top-level `run:` or `db:` key (see above; YAML keeps only the last one):

```yaml
run:                              # in the existing run: block
  unipept: true
  taxonomy: true
unipept: {result: "pept2lca.csv", split_missed_cleavages: true}
db:                               # in the existing db: block
  ncbi_taxonomy: "/data/db/taxdump"          # nodes.dmp + names.dmp
taxonomy_source: "concordant"                # eggnog | unipept | concordant
```

Per protein the consensus is the deepest rank where a majority of its peptides
agree — not the LCA of the LCAs, which one spurious peptide would drag to the
root. `concordant` keeps a taxon only where **both** methods produced a taxid
and agreed at genus or below. Proteins where either method simply had no answer
(`unipept_missing`, `eggnog_missing`, an unresolved or merged taxid, no peptide
LCA, or fewer than `consensus_min_peptides` peptides) are blanked as well, not
only those that disagree, so the taxon-based steps can cover far fewer proteins
than you expect. Check the verdict counts in `taxonomy_comparison.tsv`.

## Databases

```bash
hmmpress /data/db/Pfam-A.hmm ; hmmpress /data/db/dbCAN-HMMdb-V13.txt
diamond makedb --in VFDB_setA_pro.fas -d /data/db/vfdb_core     # and MEROPS, CARD, TADB, BAGEL
foldseek databases Alphafold/UniProt50 /data/db/foldseek/afdb50 tmp
```

To disable a default DIAMOND database, set its value to an **empty string**
(`merops: ""`). Deleting the key has no effect — the defaults are merged back
in and keep pointing at `/data/db/*.dmnd` — and `merops: null` reaches
`os.path.exists(None)` in the diamond stage and raises a TypeError. `doctor`
still lists a disabled database as `MANUAL: no path configured`. Give each
database you do keep an entry in `diamond_weights` or it scores 0 and the run
warns.

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
rowData(qf[["peptides"]])$assignment_class  # unique / taxon_unique / shared
colData(qf)                            # the design, from your manifest
metadata(qf)$taxon_size_factors
metadata(qf)$normalisation_risk
```

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
organism. The taxon reference is the median of ratios across **all** proteins
assigned to that taxon, the protein being tested included, so a taxon carried
by only a few proteins partly regresses against itself. Below
`taxon_min_proteins_for_factor` proteins the factor degrades to a plain sum and
is not trustworthy — treat calls for taxa with fewer than about six proteins as
provisional. Note that this Python threshold and the R usability threshold
`analysis.taxon_min_proteins` are separate keys and are not kept in step.
Proteins whose significance disappears under the second model were tracking
their organism, not being regulated.

The shortlist is not the top of `effector_score`: significant, no KO, and
predicted to reach the host. `surface_or_secreted` gates the list; the score
only orders it. That gate is the OR of a signal peptide, a beta-barrel, an
LPxTG motif and an anchor domain — so with `topology` off it narrows to the
last two rather than closing, and an empty shortlist is as likely to mean
nothing was significant as it is to mean the tool was missing. Check the
differential-abundance table before concluding either. See `docs/signalp-6.md`.

## Parallelism

Stages form a dependency graph and independent ones run concurrently. Every
stage up to `integrate` depends only on the protein FASTA, so hmmsearch,
DIAMOND, InterProScan, KOfamScan, SignalP and the rest run together instead of
one after another; `jackhmmer`, `hhblits` and `esmfold` then run together after
it.

```yaml
threads: 32
stage_workers: 4      # 4 stages at once, 8 cpu each
ram_gb: 128           # 32 GB each; 0 = 80% of detected RAM
```

`--threads` and `--ram` override both from the command line. The suffixed
spellings (`64G`, `64GB`, `512M`) are accepted by the `--ram` **flag only**;
`ram_gb` in config.yaml must be a plain integer number of gigabytes
(`ram_gb: 64`), and a suffixed value there raises a ValueError traceback. Note
that `init` writes the config with `yaml.safe_dump`, which drops the comments
that explain these units, so the generated file shows a bare `ram_gb: 0`.

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
interleave their writes; a lock from a dead process is reclaimed
automatically. A stage whose dependencies were excluded by `--only`/`--from`
refuses rather than producing a confident answer from inputs that do not exist.

`--serial` forces one at a time. Parallel and serial execution are *intended*
to produce identical output; no shipped test verifies this.

Within stages: DIAMOND runs its databases concurrently (`diamond_workers`,
sublinear thread scaling makes 4×N/4 faster than 4 sequential N), and hhblits
parallelises its per-query searches (`hhblits_workers`). Foldseek's multiple
targets stay **serial on purpose** — a target index is tens to hundreds of GB
and two at once will thrash. ESMFold is GPU-bound and serial by nature.

The scheduler's own overhead is negligible next to the search tools. No
benchmark script or speed-up measurement ships with this file, so no figure is
quoted here.

**Peak memory scales with `stage_workers`.** Four hmmsearch jobs against
Pfam-A alongside InterProScan is the usual squeeze. `doctor` says so.

## What has actually been run

**On real data:** one label-free dataset, end to end — a FragPipe
`combined_peptide.tsv` plus its `.fp-manifest` plus a precomputed eggNOG table,
38,204 proteins, 122,278 features. The stages that ran were `emapper` (reuse),
`integrate`, `finalise` and `join`: manifest parsing, the id join, the
shared-peptide rule, binning, the roll-up and the design recovery. That is the
path most runs take before any search tool starts, and it works.

**Not on real data:** every external search stage was disabled in that run, so
hmmsearch, DIAMOND, InterProScan, KOfamScan, HHblits, jackhmmer, ESMFold,
Foldseek, SignalP/TMbed and the Unipept API remain unexercised outside their
output parsers. `3s_structure_only` and `3p_profile_only` have never been
populated from a real search.

**Report and R object: synthetic data only.** The report **knits** to HTML with
figures and tables and the object builds, under R 4.3, pandoc 3.1, rmarkdown,
limma, SummarizedExperiment and the tidyverse. The QFeatures assay-link branch
— the headline deliverable of the object script — has not been exercised; it
sits inside two nested `try(..., silent = TRUE)` calls and on failure prints a
hard-coded "QFeatures version differs" regardless of the real error and degrades
to a SummarizedExperiment. Check the class of the object you get back.

No test file, benchmark script, synthetic generator or fixture ships with this
file, so none of the above is reproducible from what you have.

## Scale

No benchmark data accompanies this file, so no wall-time or peak-memory figure
is quoted. The hot paths are vectorised — protein-group explosion, bin
assignment and effector scoring are array operations, not row-wise `apply` —
and sequences are not retained after the single FASTA pass, so memory scales
with protein count rather than total residues. Reading a large quant table on
the delimiter-sniffing path (`sep=None, engine="python"`) is the expensive step
and a full run makes three such reads; budget on the order of a gigabyte and
tens of seconds per read at 400k features.

Scratch is **not cleaned up**: each Foldseek target search leaves
`{results}/foldseek/tmp{i}`, self-clustering leaves `{results}/foldseek/tmpc`,
and `cluster` leaves `{results}/cluster/tmp`. Against AFDB50 that is tens to
hundreds of GB. Delete them yourself when the run is done.

## Three things to check before believing any of it

1. **Gene calling.** Prodigal in meta mode has a hard floor of 90 nt and
   reduced sensitivity below about 100 codons, so unless the search database
   was augmented with a dedicated smORF call, bacteriocins, TA toxins and RiPPs
   were largely absent from the search space. Check whether your database has a
   smORF tier before quoting this caveat — the report prints it unconditionally.
2. **Group and taxonomy conflicts.** `group_conflicts.tsv` and
   `taxonomy_comparison.tsv`. KO-less proteins are strain-specific, so theirs
   are the least trustworthy.
3. **Taxon confounding.** `taxon_intensity.tsv` as a covariate, or a protein
   that merely tracks its source organism reads as a regulatory finding.
