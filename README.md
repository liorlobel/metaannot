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
(`tool_args`, `db.diamond`, `sources.diamond`, `diamond_weights`, `vfdb_category_weights`,
`diamond_evalues`), whose keys are user-chosen and cannot be checked, so a
typo there is silent and simply has no effect. Proof-read those blocks by hand.
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

A healthy default run is about **589 passed, 35 skipped, 8 xfailed, 33
deselected**, in two to four minutes depending on the machine. On Windows four
of those come back as failures instead: a path test asserting forward slashes,
a `doctor --fix` recipe emitting `mkdir -p`, a SIGINT test, and stale-lock
reclamation, which `_holder_is_alive` deliberately disables on Windows because
`os.kill(pid, 0)` there calls `TerminateProcess` — asking whether the holder is
alive would kill it. All four are POSIX assumptions in the tests, not defects
in the tool.

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
The ones open today are live defects the tool documents rather than hides: an
`annotation_pass1.tsv` that is not reproducible across a resume, a Unipept
lineage truncated at the first blank rank, and a `pept2lca` file matching
nothing dying with a bare `'verdict'`.

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

That is the complete set; `annotation_final.tsv`, `bin_summary.tsv` and the
report's factor levels use exactly these seven strings. The tests are applied in
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
which is the number to set it on. `analysis.min_plexes` is inert without a
plex column, so label-free runs are unaffected; setting it above 1 where no
per-sample plex exists stops the report rather than passing everything.

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

`peptide_evidence.tsv` has one row per protein and nine columns:
`protein_id`, `n_features_used`, `n_unique`, `n_taxon_unique`,
`n_family_unique`, `n_features_dropped`, `taxon_unique_dominated`,
`rollup_method` and `peptide_assignment`. The last two record how the numbers
were made, because a `sum` run and a `median_polish` run are otherwise
indistinguishable once the log is gone. `taxon_unique_dominated` flags proteins
resting more on shared-but-taxon-unique features than on their own unique ones
— the ones whose intensity is most sensitive to the assignment rule you chose.

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
slots. The scheduler orders that queue longest-first. Every stage carries a
coarse cost rank — hours, minutes, or seconds, measured on real runs — and the
hours-class stages claim the workers while the seconds-class ones fill in
behind them as slots free. Dispatching in table order instead, as it did
before v0.4.0, gave the first wave to `dbcan` (10 min) and `diamond` (5 min)
while `signalp` and `tmbed` (about an hour each) queued.

The order is a starting order, not a schedule: it makes nothing faster, and
InterProScan still paces a large run. See
[Sizing your run](TUTORIAL.md#sizing-your-run) for what to do about that.

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
read as live. Trampling a live run is silent corruption; refusing is a message.
So an ordinary crashed run on Windows exits with `another metaannot is already
running here` and needs `--force-unlock`, which is the escape hatch for a
holder you are certain is gone.

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
failures — fold them on a card with more memory, in the cloud, or on CPU, and
drop the models into `results/structures/` before rerunning `foldseek`.

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
3.8%. Dark rescue took 9,731 proteins dark on eggNOG alone down to 1,462 —
8,686 rescued, of which KOfam alone supplied a KO to 8,129 that eggNOG missed.
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

A test suite does ship, in `tests/` — `conftest.py`, `fixtures.py` and thirteen
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
single-stage benchmarks:

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

Those figures scale roughly with protein count. On a 455,571-protein run —
twelve times the size — expect InterProScan alone to pace the whole thing at
well over a day.

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
