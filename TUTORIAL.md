# TUTORIAL — running metaannot on the lab server

A runbook. Each phase ends with a **check** whose result decides whether to
continue. Do not skip the checks; the whole point of them is that a wrong
assumption at phase 3 is cheap and at phase 7 is not.

Estimated wall time for a first full run on a real dataset: **1–3 days**, most
of it InterProScan and Foldseek. Phases 0–4 take under an hour and catch
almost everything that goes wrong.

---

## Phase 0 — Survey the machine

Nothing is configured yet. Find out what there is.

```bash
nproc                                  # cores
free -g | awk '/^Mem:/ {print $2" GB"}'
df -h --output=target,size,avail -x tmpfs -x devtmpfs | sort -k2 -h
lsblk -o NAME,SIZE,FSTYPE,MOUNTPOINT
nvidia-smi 2>/dev/null || echo "no GPU on this host"
```

**Decide two things and write them down:**

- `THREADS` = cores, minus 2 if others share the machine.
- `DBROOT` = a directory on the **large array**, not the OS disk. AFDB50 is
  about **123 GB to download and ~200 GB on disk**; full `Alphafold/UniProt` is
  ~491 GB and PDB ~2 GB. Everything else together is well under 100 GB.
  Leave room for scratch: Foldseek and MMseqs2 temp directories under
  `results/` are never cleaned up and can run to tens of GB.

```bash
export DBROOT=/path/to/large/array/db          # edit
export PROJ=/path/to/large/array/metaannot     # edit
mkdir -p "$DBROOT" "$PROJ" && cd "$PROJ"
df -h "$DBROOT" | tail -1
```

**CHECK.** `$DBROOT` has ≥ 400 GB free if Foldseek/AFDB50 is wanted (download
plus extraction plus scratch), or ≥ 100 GB without it. If not, stop and decide which databases to skip — that is
a scientific choice, not a technical one, and it belongs to Lior.

---

## Phase 1 — Environment

```bash
mamba create -y -n metaannot -c conda-forge -c bioconda \
  python=3.11 pandas pyyaml numpy \
  hmmer diamond mmseqs2 foldseek
mamba activate metaannot
python -c "import pandas, yaml, numpy; print('python deps ok')"
```

Only `pandas`, `pyyaml` and `numpy` are hard requirements. Everything else is
needed solely by the stage that uses it.

Optional, add as needed:

```bash
mamba install -y -c bioconda kofamscan hmmer          # kofam stage
mamba install -y -c bioconda hhsuite                  # hhblits stage
pip install tmbed && tmbed download                   # GPU host only
pip install "fair-esm[esmfold]"                       # GPU host only
```

SignalP 6.0 and InterProScan are separate installs (SignalP needs an academic
licence; InterProScan is a large Java distribution). Skip both on the first
pass — `run.topology: false` and `run.interpro: false`.

**CHECK.**

```bash
python metaannot.py --version && python metaannot.py --help | head -20
```

---

## Phase 2 — Config skeleton

```bash
python metaannot.py init --out config.yaml
```

Edit these first; leave the rest at defaults for now.

```yaml
results_dir: "results"
threads: 32            # your THREADS
ram_gb: 128            # leave 20% for the OS; 0 auto-detects
stage_workers: 4       # 4 stages at once, threads/4 cpu each

proteins_faa: ""       # filled in at phase 3
quant_table: ""
quant_format: "fragpipe_peptide"
manifest: ""

emapper_precomputed: []   # your existing eggNOG output — see phase 3

db:
  pfam_hmm: "/DBROOT/Pfam-A.hmm"
  ncbi_taxonomy: "/DBROOT/taxdump"
  diamond: {}          # empty for now
run:
  pfam: true
  dbcan: false
  diamond: false
  topology: false      # needs SignalP + tmbed
  cluster: true
  structure: false     # needs GPU + Foldseek DB
  interpro: false
  kofam: false
```

Substitute the real `$DBROOT` path — YAML does not expand shell variables.

A fresh `init` config already has `eggnog`, `pfam`, `dbcan`, `diamond`,
`cluster` and `join` **on** and everything else off. So of the block above, the
two lines that actually change something are `dbcan: false` and
`diamond: false` — they save you two database downloads on the first pass;
`pfam: true` and `cluster: true` restate the default, and the rest are already
off and are written out only so the file says what you decided. `integrate` and
`finalise` have no flag and always run. Edit this block in place; do not append
a second top-level `run:` key later in the file. A duplicate key is refused —
the run exits naming the key and both line numbers.

**Relative paths in the config are resolved against the config file**, not the
current directory, so `metaannot.py run --config /path/to/project/config.yaml`
works from anywhere. Paths given on the command line (`--results-dir`) are
relative to where you are standing.

---

## Phase 3 — Inputs, and the identifier check

This phase is where runs are saved or lost.

### 3a. Point at the FragPipe output

```bash
FP=/path/to/fragpipe/output          # edit
ls "$FP"/*.fp-manifest "$FP"/combined_peptide.tsv
head -3 "$FP"/*.fp-manifest
head -1 "$FP"/combined_peptide.tsv | tr '\t' '\n' | grep -c Intensity
```

The manifest carries path, experiment and bioreplicate per run — that is the
whole experimental design, so no metadata file needs writing.

A manifest lists one row per **raw file**, not per sample. A fractionated
acquisition therefore repeats the experiment name across its fraction rows,
which is how FragPipe denotes fractions; it writes one quant column per group.
metaannot collapses those rows and logs it:

```
manifest: 4 sample(s) are split across multiple fraction files; each group is
one quant column in FragPipe's output and is treated here as one sample
```

Read that line and confirm the count is the number of samples you expect. It
refuses only when rows of one sample disagree about the design (different
`data_type`), or when two different samples resolve to the same quant column.

**TMT does not come through this path at all.** In a `.fp-manifest`
`experiment` is the plex, so a design derived from it would contrast plex
against plex. Set `quant_format: fragpipe_tmt` and point `quant_table` at the
FragPipe run directory holding the `TMTn/` folders: that reader builds the
design itself with `plex` as a covariate, and derives the condition from the
annotated sample names (`tmt.condition_from_name`, `auto` by default). It logs
the derivation as a **guess** — check it, and set `tmt.condition_from_name` or
write `analysis.metadata` by hand (sample, condition, batch) where the names
cannot carry it. See README, "FragPipe TMT (isobaric)".

Set in `config.yaml`:

```yaml
manifest: "/path/to/experiment.fp-manifest"
quant_table: "/path/to/combined_peptide.tsv"
quant_format: "fragpipe_peptide"
```

### 3b. Build the identified-protein FASTA

Annotate what was identified, not the whole catalogue — it is far faster and
the bin percentages then describe the actual result.

```bash
python metaannot.py subset \
  --db /path/to/search_database.fasta \
  --quant "$FP"/combined_peptide.tsv \
  --format fragpipe_peptide \
  --out input/proteins.faa
grep -c '^>' input/proteins.faa
```

A warning about ids absent from the database usually means decoy or
contaminant prefixes (`rev_`, `sp|`, `CON__`). A handful is normal; thousands
means the wrong FASTA.

### 3c. Reuse the existing eggNOG-mapper output

```yaml
proteins_faa: "input/proteins.faa"
emapper_precomputed:
  - "/path/to/catalogue.emapper.annotations.gz"
```

### 3d. **The identifier check — do this before anything long**

```bash
python metaannot.py run --config config.yaml --only emapper 2>&1 | tail -20
```

Read the coverage line. Three outcomes:

| coverage | meaning | action |
|---|---|---|
| ≥ 90% | fine | continue |
| < 90%, a transform is named | identifiers are formatted differently | set `emapper_id_transform` to the named value and rerun this step |
| < 90%, a prefix is named | the search database prefixed its ids (`uhgpSM_`, `HUMANHOST_`) and the eggNOG table did not | set `emapper_strip_id_prefix` to the named prefix (string or list) and rerun this step |
| < 50% | aborts | **stop.** Wrong eggNOG file, or the FASTA is from a different assembly. Do not raise `emapper_min_coverage`. |

`emapper_strip_id_prefix` applies to `emapper_precomputed` reuse only; the log
says so if you set it without one. When it fires, the run reports how many rows
matched only because of it — quote that number, it is the size of the bridge.

**CHECK.** Coverage ≥ 90%, or a clear reason why not. Report the number.

---

## Phase 4 — Doctor, dry run, and a subsample

### 4a. Doctor — and let it install what is missing

```bash
python metaannot.py doctor --config config.yaml
```

Every enabled stage must show `OK` for its tool and database, every manifest
run must map to a quant column, and the resources block must show a sane
per-stage CPU/RAM split. A `== config ==` block means a key is misspelled and
that setting is silently not in effect — fix it before anything else. The
`analysis:` block **is** covered by that check (its keys are the Rmd's params,
so a misspelled `fdr` is reported rather than leaving the default quietly in
force). What the check still cannot see inside are the free-form blocks
(`tool_args`, `db.diamond`, `sources.diamond`, `diamond_weights`,
`vfdb_category_weights`, `diamond_evalues`), whose keys you choose — proof-read
those by hand.

For anything missing, `doctor` prints the exact commands. Three ways to act on
them, in increasing order of trust:

```bash
python metaannot.py doctor --config config.yaml --install-plan install.sh
less install.sh && bash install.sh          # review, then run

python metaannot.py doctor --config config.yaml --fix       # asks first
python metaannot.py doctor --config config.yaml --fix --yes # no prompt
```

`--fix` prints every item and two totals — **to download** and **on disk** —
before doing anything, and refuses without confirmation. Read them, but treat
them as indicative: the per-item sizes are hand-maintained. Check any item over
about 10 GB against the provider before approving it.

Two things it will not do:

- **`MANUAL` items are never attempted.** SignalP 6.0 needs an academic
  licence, InterProScan is a version-specific Java distribution, and several
  targeted databases are behind click-throughs. `doctor` tells you where to get
  each and what to set afterwards.
- **It does not trust exit codes.** Downloads use `curl -fL`, which fails on an
  HTTP error status, and after installing it re-runs the checks. An item whose
  commands succeeded but whose file is still absent is reported `UNVERIFIED`
  rather than passing. `-L` also follows redirects, so a link repointed at a
  landing page used to write an HTML page to disk with exit code 0 and pass
  forever. The generated commands for Pfam-A, dbCAN, NCBIfam, UniRef50 and the
  DIAMOND FASTAs now each sniff the first 512 bytes and abort with `is an HTML
  page, not the database`, deleting the file. If you see that, fix the URL under
  `sources:` — do not retry. The KOfam and taxdump recipes have no such guard,
  but both go through `tar`/`gunzip`, which fail loudly on a web page. A file
  you fetched by hand is covered by neither, so check it: `head -1` of an HMM
  file should start with `HMMER3/`.
- **A database path that is not configured is reported, not guessed.** An empty
  `db.<name>` shows as `MANUAL: no path configured` instead of producing
  commands that download into the working directory.

URLs live in the config under `sources:`. If a provider moves a link, edit it
there rather than patching the tool — and check the release numbers (`dbCAN
V13`, the NCBIfam archive layout) against the provider before a big download.

**Prefer `--install-plan` on a shared server.** It gives you a script to read
before anything touches the filesystem.

### 4b. Dry run

```bash
python metaannot.py run --config config.yaml --dry-run
```

Confirms the stage plan. `disabled` is expected for anything switched off.

### 4c. Subsample first — do not skip this

A full run is a day or more. A 5,000-protein subsample takes minutes and
catches format problems that no amount of static checking will.

```bash
mkdir -p sub/input
awk '/^>/{n++} n<=5000' input/proteins.faa > sub/input/proteins.faa
grep -c '^>' sub/input/proteins.faa
python metaannot.py run --config config.yaml \
    --faa sub/input/proteins.faa --results-dir sub/results 2>&1 | tail -30
cat sub/results/bin_summary.tsv
```

Override the two paths on the command line rather than copying the config.
Relative paths inside a config resolve against **that config file's own
directory**, so a `sub/config.yaml` containing `sub/input/proteins.faa`
resolves to `sub/sub/input/proteins.faa` and the run dies with
`proteins_faa not found`. If you do want a separate config, write it into
`sub/` with paths relative to `sub/` (`proteins_faa: input/proteins.faa`,
`results_dir: results`) — and never edit `results_dir` with a quote-sensitive
`sed`: `init` emits it unquoted, so the substitution silently misses and the
subsample overwrites the production `results/`.

The bin table and the `% of proteins are invisible to KEGG pathway enrichment`
line both go through the logger, so they land in `results/metaannot.log` as well
as on the console — you can read them back after the fact with
`grep -A 10 'invisible to KEGG' results/metaannot.log`.

**CHECK — what a sane bin summary looks like.**

- `1_ko_pathway` is typically **25–50%** of proteins in a gut metaproteome.
- The KO-less bins together are typically **50–75%**. A high number here is the
  expected result, not an error — it is the whole reason the tool exists.
- **`4_dark` at 95%+ means something is broken**, almost always the identifier
  mismatch from phase 3d rather than genuinely unannotated proteins.
- **`1_ko_pathway` at 0%** means the eggNOG table has no `KEGG_Pathway` column
  or the wrong column set.

Only bins with at least one protein appear as rows, so a missing
`3s_structure_only` / `3p_profile_only` means the structure and profile stages
were off, not that the table is broken.

A real reference point, from a label-free UC metaproteome run with **every
search stage off** (eggNOG reuse only, 38,204 proteins):

```
1_ko_pathway       15591  40.8%
2_ko_orphan         5236  13.7%
3_annotated_no_ko   5484  14.4%
3d_duf_only          787   2.1%
4_dark             11106  29.1%
```

`3_annotated_no_ko` and `3d_duf_only` are non-zero there only because eggNOG's
own `PFAMs` column counts as domain evidence; turning `pfam`, `ncbifam`,
`interpro` and `dbcan` on moves proteins out of `4_dark` and into those two.

Report the actual percentages. If they fall outside these ranges, stop and say
so rather than continuing to a full run.

---

## Phase 5 — Full CPU run

### 5a. Turn on what you want, then let doctor fetch it

Enable the stages in `run:`, point `db:` at paths under `$DBROOT`, then:

```bash
python metaannot.py doctor --config config.yaml --install-plan install.sh
less install.sh          # read it: this is where the AFDB50 decision is made
bash install.sh
python metaannot.py doctor --config config.yaml    # confirm
```

Sensible order if you are adding things gradually: **Pfam-A** (~2 GB, the
largest single gain), then **NCBI taxdump** (~0.1 GB to download, ~0.5 GB once
extracted; needed for the eggNOG/Unipept comparison), then the targeted DIAMOND
databases (small), and
only then decide about Foldseek. `foldseek databases PDB <path> tmp` is a few
GB instead of AFDB50's ~123 GB download and still finds most classical toxin
folds.

Give every DIAMOND database an entry in `diamond_weights` or it contributes 0
to the effector score — `doctor` warns about this.

### 5b. Run

```bash
tmux new -s metaannot
mamba activate metaannot
cd "$PROJ"
python metaannot.py run --config config.yaml --threads 32 --ram 128 \
  2>&1 | tee run_$(date +%F_%H%M).log
```

Detach with `Ctrl-b d`. Reattach with `tmux attach -t metaannot`.

A results directory takes a lock for the duration. A second run against the
same directory refuses rather than interleaving its writes. A lock is reclaimed
automatically only when it names a pid on this host that is provably gone —
which on Windows is never, because asking there would kill the process. A
crashed run may therefore need `--force-unlock`. The exception is a zero-byte
lock, which is what a power loss or a hard crash leaves behind: those are
removed on sight once they are more than a minute old.

Monitor from another shell:

```bash
tail -f "$PROJ"/results/metaannot.log
python -c "import json;print(json.dumps(json.load(open('results/.metaannot_state.json')),indent=1))"
```

The state file shows which stages finished, how long each took, and any that
failed with the reason.

**If it dies:** rerun the same command. Completed stages are cached and skipped;
only the failed one and its dependents rerun. Do **not** add a bare `--force` —
it discards every cached stage.

To redo one stage on purpose — you fixed a database, or a tool wrote garbage —
use the scoped form, not `rm`:

```bash
python metaannot.py run --config config.yaml --force --only pfam
```

`--force` discards the state records of only the stages the selection names, so
every other stage keeps its record and is still signature-checked next time.
`--only` on its own reports the stage as `cached` and skips it, and says as
much. Deleting the stage's output files does trigger a recomputation, but
leaves nothing written down about why — use the flag.

**CHECK.** `done: N run, M adopted, K skipped` and `annotation_final.tsv`
exists. Compare `bin_summary.tsv` against the subsample — the percentages
should be similar. A large divergence means the subsample was not
representative, which is worth knowing.

---

## Phase 6 — GPU stages

`tmbed` and `esmfold` want CUDA. Run them on a GPU host and bring the results
back; metaannot **adopts** outputs it finds rather than recomputing.

### No CUDA GPU? Read this first

`doctor` now tells you before a run starts, rather than letting you find out
when `esmfold` finally runs. With `run.structure` or `run.topology` on it
prints a `== gpu ==` block:

```
== gpu ==
  WARN   a card is present (NVIDIA GeForce RTX 5080) but torch reports CUDA
         unavailable - usually a CPU-only torch build
  MISS   run.structure needs CUDA: stage_esmfold exits rather than fold on CPU
  WARN   run.topology: SignalP 6 is CPU-only and unaffected. tmbed will fall
         back to CPU, where it is one to two orders of magnitude slower
```

It separates the two cases that look identical and are not: **no card**, and
**a perfectly good card with a CPU-only `torch` wheel**. The second is the more
confusing failure, because the hardware is right there. Fix it by reinstalling
torch from the CUDA index matching your driver.

What to do, by stage:

| stage | without CUDA |
|---|---|
| `signalp` | **unaffected** — SignalP 6 is CPU-only anyway |
| `tmbed` | runs, but 1–2 orders of magnitude slower. Fine for a few thousand proteins; not for a few hundred thousand |
| `esmfold` | **refuses.** Folding on CPU is impractical at any real scale, so the stage exits rather than pretend |
| `foldseek` | **unaffected** — CPU-only, and it searches whatever models are present |

So a machine with no GPU can still run 19 of the 21 stages, and can still get
structural evidence: fold elsewhere, copy `results/structures/` across, and
`foldseek` will search what it finds. That split is deliberate — the expensive
GPU step is separable from the search that uses it.

Three settings decide what happens:

```yaml
tmbed_use_gpu: auto     # GPU if present, CPU if not (default)
tmbed_use_gpu: true     # a missing GPU is FATAL - use this when slow is worse than absent
tmbed_use_gpu: false    # CPU deliberately, no warning
run.structure: false    # skip esmfold and foldseek entirely
```

`auto` is the default because losing topology evidence to a missing card is
worse than being slow — but at a few hundred thousand proteins that judgement
inverts, which is why the run says the protein count out loud before starting
the CPU path.

**On Apple Silicon**, `torch.cuda.is_available()` is False regardless of how
good the chip is: neither `tmbed` nor `esmfold` has an MPS path today. Treat a
Mac as a no-GPU host for these two stages. See `docs/gui-design.md` for the
full per-platform picture.

On the laptop:

```bash
rsync -av server:$PROJ/input/ ./input/
rsync -av server:$PROJ/config.yaml ./
rsync -av server:$PROJ/results/dark.faa \
         server:$PROJ/results/annotation_pass1.tsv ./results/   # after phase 5
pip install tmbed "fair-esm[esmfold]" && tmbed download

python metaannot.py run --config config.yaml --only tmbed --serial
python metaannot.py run --config config.yaml --only esmfold --serial --ram 16
```

`annotation_pass1.tsv` has to come across as well. `esmfold` and `tmbed` read
only `dark.faa`, but the dependency check requires **every** output of the
`integrate` stage to be present before it lets the stage run, and refuses
otherwise.

Phase 2 turned `topology` and `structure` off. Naming a stage with `--only`
overrides its run flag — the two commands above therefore work as written and
warn `run.topology is false, but the stage was named with --only, so it is
running anyway`. That override is per-invocation and does **not** help adoption:
on the server the dry run below plans without `--only`, so `topology` and
`structure` must be `true` in the server's config or both stages read
`disabled` and nothing is adopted.

Enabling `topology` also arms the `signalp` stage, which shares that flag and
needs the licence-gated `signalp6` binary; if you have no licence, keep
`topology: false` on the server for normal runs and set it true only for the
adoption dry run.

Back on the server:

```bash
rsync -av laptop:results/topology/tmbed.pred  results/topology/
rsync -av laptop:results/structures/          results/structures/
python metaannot.py run --config config.yaml --dry-run | grep -E "tmbed|esmfold"
```

Both should read `adopt`. Note that `--dry-run` does not run the dependency
check, so it can show `RUN` for a stage a real invocation would refuse.
If either reads `RUN`, the file is missing or empty —
adoption deliberately refuses a zero-byte output, because that is usually an
interrupted writer.

Then Foldseek on the server, which is CPU-only and wants the big array:

```bash
foldseek databases Alphafold/UniProt50 "$DBROOT"/foldseek/afdb50 tmp   # ~123 GB, hours
python metaannot.py run --config config.yaml --only foldseek finalise --ram 128
```

---

## Phase 7 — Unipept taxonomy

Run Unipept yourself; the HTTP route depends on a reachable server and on the
API version matching.

**Turn the stages on first.** `run.unipept` is false by default, so `--only
unipept` before this edit reports `disabled (run.unipept)`, writes nothing, and
the rest of the phase has no file to work on. Inside the **existing** `run:`
block set `unipept: true` and `taxonomy: true`; inside the **existing** `db:`
block set `ncbi_taxonomy: "/DBROOT/taxdump"`; and at the top level set
`unipept: {result: "pept2lca.csv"}` and `taxonomy_source: "concordant"`
(`eggnog | unipept | concordant`). Do not paste a second top-level `run:` or
`db:` key: a duplicate key is refused, so rather than quietly discarding your
phase 2 settings (including `structure: false` and `topology: false`), the run
exits naming the key and both line numbers.

Then:

```bash
python metaannot.py run --config config.yaml --only unipept 2>&1 | tail -5
# writes results/unipept/peptides.txt, then aborts asking for the result
wc -l results/unipept/peptides.txt

npm install -g unipept-cli   # Node 22+; the Ruby gem is the legacy client
unipept pept2lca --equate --all -i results/unipept/peptides.txt -o pept2lca.csv
# the current CLI writes domain_id/domain_name where the gem wrote
# superkingdom_id/superkingdom_name, which is what the parser reads
```

The unipept.ugent.be web export is **not** an accepted input: it is name-based
and carries no taxid columns, so the parser cannot read it. Use the CLI.

```bash
python metaannot.py run --config config.yaml --only unipept taxonomy join
column -t results/unipept/taxonomy_comparison.tsv | head
```

**CHECK.** The verdict counts. Expect a mix of `identical`, `concordant` and
`concordant_above_genus`. A large `conflict` fraction is a real finding worth
reporting — eggNOG has no close reference for divergent gut organisms — but
first confirm `ncbi_taxonomy` is set. Without it the comparison degrades to
exact taxid identity, emits the verdict `differ (no lineage)`, and
*Enterococcus* vs *E. faecalis* reads as a conflict. That verdict is not one of
the report's factor levels, so those rows show as NA in the verdict table and
the plot, and the "disagree below family" gate ignores them. Run the comparison
with a taxdump.

`taxonomy_source: concordant` makes the taxon-based steps use only proteins
where both methods produced a taxid and agreed at genus or below. Proteins
where either method had **no answer** are dropped as well, not only those that
disagree, so far fewer proteins may survive than you expect — read the verdict
counts, not just the conflict fraction. That is the conservative choice.

---

## Phase 8 — Report and R object

Needs R. If R is on the MacBook rather than the server, rsync the results
directory there **and regenerate the Rmd locally** with
`python metaannot.py report --no-render --config config.yaml`: the Rmd written
on the server embeds absolute server paths in its YAML header, so an rsynced
one points at directories that do not exist there. The laptop therefore needs
python3, pandas and pyyaml as well as R.

```r
install.packages(c("tidyverse","rmarkdown","knitr","patchwork"))
BiocManager::install(c("limma","clusterProfiler","QFeatures","SummarizedExperiment"))
```

```bash
python metaannot.py doctor --config config.yaml   # the "== R ==" block
```

That block probes four packages (SummarizedExperiment, QFeatures, limma,
rmarkdown) and checks neither knitr nor pandoc, while the Rmd also loads readr,
dplyr, tidyr, tibble, stringr, ggplot2 and purrr. A clean machine can pass
`doctor` and then fail at the first `library()` call; install what it names by
hand and rerun.

Set the design in the config's `analysis:` block. Contrasts are derived from
the manifest automatically; override only for a non-pairwise comparison or to
add covariates:

```yaml
analysis:
  design_formula: "~ 0 + group"
  factor_cols: "group"
  block_col: ""              # "cage" for co-housed animals
  contrasts: ""              # "" = every pairwise contrast
  min_valid_per_group: 3
  fdr: 0.05
```

For covariates beyond condition and replicate, copy
`results/quant/design_from_input.tsv`, add columns, point `analysis.metadata`
at it, and extend `design_formula` to match.

```bash
python metaannot.py report --config config.yaml     # writes and knits the Rmd
python metaannot.py object --config config.yaml     # writes results/metaannot.rds
```

If the knit fails, `report --no-render` writes the Rmd without running it —
open it in RStudio and step through, which gives a far better error. The
generated document carries absolute paths, because `rmarkdown::render()`
evaluates a document with the working directory set to its own folder.

**CHECK the report, in this order.**

1. **Design diagnostics.** Rank deficiency stops the run. Empty cells in a
   factor cross-tabulation mean partial confounding.
2. **Feature support by bin.** KO-less proteins carry fewer peptides. If the
   median for `4_dark` is 1, those quantifications are not evidence.
3. **`% invisible to KEGG enrichment`.** This is the headline number.
4. **Taxonomy agreement by bin.**
5. **Normalisation risk.** Any taxon marked `AT RISK` needs both models
   reported, not just the ratio model.
6. **Effector shortlist** — significant, no KO, predicted secreted or surface.

---

## Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `manifest: N sample(s) are split across multiple fraction files` | a fractionated acquisition, one manifest row per fraction | not an error — the fractions are collapsed into one sample. Check N is the number of samples you expect |
| `manifest: several samples map to the same quant column` | two different experiment/bioreplicate pairs resolve to one column | rename them so each run resolves to its own column; this is not fractionation and is not collapsed |
| `have rows with different data_type` from the manifest | rows sharing a sample name are not fractions of one sample | fix the manifest, or write `analysis.metadata` by hand |
| `this looks like isobaric (TMT/iTRAQ) output` | a per-plex FragPipe TMT table read under a label-free `quant_format` | the format is wrong, not the file. Set `quant_format: fragpipe_tmt` and point `quant_table` at the run directory holding the `TMTn/` folders, which the message names; that reader takes the reporter channels of every plex and joins them |
| `is an HTML page, not the database` from an install command | the download URL now redirects to a landing page | fix that entry under `sources:` in `config.yaml` and rerun the install; do not retry the same URL |
| `only N% of proteins matched` | identifier format differs | apply what metaannot names — an `emapper_id_transform`, or an `emapper_strip_id_prefix` when the search database prefixed its ids. Never raise `emapper_min_coverage` |
| `emapper_strip_id_prefix 'X' matches no fasta id` | the prefix is wrong, or the FASTA never carried it | check the prefix against `head -1 input/proteins.faa`; the coverage diagnostic names one that would have matched |
| `N fasta ids collide once emapper_strip_id_prefix is removed` | the prefix was the only thing distinguishing those ids | the first is kept and the rest cannot be annotated. Do not strip that prefix; fix the FASTA ids instead |
| `4_dark` is 95%+ | same as above, or the wrong eggNOG file | recheck phase 3d |
| a stage produced none of the evidence you expected | its dependency was disabled | a disabled dependency does **not** block a dependent stage; the dependent runs with that evidence missing. Check `--dry-run` for `disabled` lines |
| `not adopting X — the file is empty` | interrupted writer | correct; let it rerun |
| `stage '<name>' failed: <tool> not found` | tool not installed | the message carries the install command; run it, or disable the stage. The stage is recorded `failed`, so the next run retries just it |
| `N manifest run(s) match no column in <path>` | FragPipe named columns differently | `doctor` lists up to five unmatched runs and four candidate columns, not the full mapping; check experiment vs bioreplicate naming |
| OOM, or the machine swaps | too many concurrent stages | lower `stage_workers`, or `--serial` |
| InterProScan runs out of Java heap | heap below its need | raise `--ram`; it is passed as `_JAVA_OPTIONS` |
| a tool rejects a flag | version drift | add the correct flag via `tool_args:` rather than editing the tool |
| report: `contrast term(s) not among the design coefficients:` | coefficient names differ from the guess | the error prints the available names; `make.names()` turns `high fiber` into `high.fiber` |
| everything reruns unexpectedly | a config key in that stage's signature changed | expected; check `git diff config.yaml` |
| `unrecognised key 'x' — did you mean 'y'?` | typo in config.yaml | the setting is **not** in effect; fix the spelling |
| `parsed 0 <records> from a non-empty file` | a tool wrote a truncated or wrong-format output | never noise: every protein silently loses that evidence. Find out why the output is wrong, then redo the stage with `--force --only <stage>` rather than deleting the file |
| report: `0 residual degrees of freedom` | no replication left after the model | you need replicates, or a simpler `design_formula` |
| report: `factor(s) with only one level` | the manifest does not distinguish conditions | check the `experiment` column of the manifest |
| report: `N protein(s) have zero variance within every group` | identical in every replicate | listed in `zero_variance.tsv`; usually one shared peptide or an imputed constant. `drop_zero_variance: true` removes them |
| `another metaannot is already running here` | a second run on the same results directory | wait for it; `--force-unlock` only if you are certain the other is gone |
| `cannot run: [...] produced nothing and were not selected` | `--only`/`--from` skipped a stage this one needs | rerun without the selection, or add the named stages to it |
| `deadlock: [...] can never become ready` | the scheduler has stages left but none whose dependencies are all satisfied, and nothing is still running | the message names the stuck stages and, per stage, which dependencies are unsatisfied. Normally that means those deps were excluded by `--only`/`--from`: add them to the selection, or drop the selection. If they *were* selected this is a bug in the stage graph, not a config error — keep the message |
| a stage has said nothing for hours | it is working, or it is hung — you cannot tell from silence alone | every `progress_interval_s` seconds (60 by default) the newest line the tool wrote to stderr is echoed with its elapsed time; `no output yet on stderr` is the heartbeat from a tool that prints nothing. Set `progress_interval_s: 0` to switch it off |
| `have been folded yet` after `integrate` | `esmfold` has not run | **not** a loss: nothing has been folded, so this pass simply carries no structural evidence. It is INFO, not WARN |
| `requested structures exist` … `esmfold has not finished` | the stage was interrupted, or has not been rerun since `dark.faa` grew | the missing models are **pending, not lost**. Every PDB already written is kept, so rerun the `esmfold` stage and it resumes from there |
| `requested structures exist` … `although esmfold has finished` | `esmfold` finished and models are still missing | the line names each cause separately rather than guessing between them: proteins over `max_len_structure` `were never submitted`; proteins ESMFold `attempted and failed twice` are listed with the real error in `results/structures/esmfold_failed.tsv` (`skipping <id>` in the esmfold log is the same event); proteins `absent from that list` were added to `dark.faa` after the last fold and were never attempted at all |
| a stage runs even though its `run.<stage>` flag is false | it was named with `--only` | intended — `--only` is a clearer statement of intent than a flag left off for another machine, and the run warns. This is how the GPU box runs `tmbed`/`esmfold` against the server's config |
| `waiting for the GPU — tmbed is using it` | `esmfold` and `tmbed` both need the card | expected, not stuck. Each wants essentially a whole GPU (15.5 GB and 13.3 GB of a 16 GB device were observed), so they are leased one at a time; every CPU-only stage keeps running. `gpu_workers` raises the count, but they all still share the single `gpu_device` |
| `refusing to search a DIAMOND database that cannot answer` | a `.dmnd` that is empty, truncated, or holds no sequences | a failed `diamond makedb` leaves one behind, and searching it reports 0 hits — the same output as a real absence. Rebuild it with the `diamond makedb` line the message prints, or drop the tag with `<tag>: ""` |
| `incapable of a hit before it starts` | the database's sequences are too short for the configured e-value | that search's 0 hits say nothing about the biology. Give the tag its own threshold under `diamond_evalues:`, or record that this database needs different settings. The warning also names the `diamond_weights` entry that claims the database can score |
| `looks like a motif or seed set rather than a protein sequence database` | the `.dmnd` was built from motifs, not sequences | check **what** you built before touching any threshold. BAGEL4 ships bacteriocin sequence files *and* the motif seed set its HMM step uses (`LE-`, `MA-`, `ggmotif`, `lasso`, median 15 residues); building the seed set gives 0 hits and no e-value fixes that, because a `blastp` against 15-residue seeds searches for those residues and not for the molecules they mark. Rebuild from the sequence files |
| `doctor --fix cannot run on Windows` | `--fix` on a Windows host | its install commands are POSIX shell and cmd.exe mis-executes them (`mkdir -p C:\db` makes a directory called `-p`). Use `--install-plan install.sh` here and `wsl bash install.sh` where the tools live |
| `got no prediction; listed in tmbed_failed.tsv` | one or more tmbed chunks failed | the finished chunks are kept and a rerun resumes from them. If the card is wedged, restart the machine or the WSL session rather than rerunning. `tmbed_allow_partial: true` accepts the shortfall deliberately |

---

## Resource guide — and how to size a run

These are **measured**, not estimated. Every stage that runs records its own
duration as `seconds` in `results/.metaannot_state.json`, so your run reports
its own numbers; quote those rather than these.

Two real runs on one workstation (22 cores, 94 GB to WSL2, one 16 GB card),
three stages at a time:

| stage | 38,204 proteins | 455,571 proteins | ratio |
|---|---|---|---|
| `emapper` (reuse) | 3.0 s | 76 s | 25× |
| `diamond` (5 databases) | 9.7 s | 291 s | 30× |
| `cluster` (MMseqs2) | 14 s | **109 s** | **7.6×** |
| `dbcan` | 41 s | 617 s | 15× |
| `ncbifam` | 0.40 h | **5.6 h** | 14× |
| `pfam` | 0.45 h | **7.2 h** | 16× |
| `kofam` | 0.49 h | **14.3 h** | 29× |
| `signalp` | 0.93 h | *still running* | — |
| `tmbed` | 0.87 h | *still running* | — |
| `interpro` | 2.84 h | *still running* | — |
| `foldseek` | 57 s | — | — |
| `esmfold` | 1.6 h for 1,805 models under 478 aa | — | — |

The three *still running* cells are not estimates withheld for tidiness —
that run was in its 31st hour when this table was written and those stages had
not finished. What was known at that moment: SignalP was 58% through
(262,200 of 455,571 sequences after 30.7 h, so on the order of 53 h);
InterProScan had been going 3.7 h; TMbed had written **nothing at all** in
30.7 h, which is what that tool does — see `tmbed_chunk_residues` below.
Rather than round those into the table, read your own: every stage records
its `seconds` in `results/.metaannot_state.json`.

### The thing that surprises people: it is not linear

The protein count went up **11.9×**. The stages that finished went up
**14–30×**, and the spread between them is as informative as the numbers.

The reason is contention, not size. At 38k every stage finishes quickly, so
three-at-a-time rarely means three long stages overlapping. At 455k every long
stage runs concurrently with every other long stage for its entire life, and
they share 22 cores. SignalP measured **11.4 sequences/s** with the machine
mostly to itself and **2.0–3.2 sequences/s** with tmbed and KOfam alongside —
the same work, three to five times slower.

So: **scale by observed contention, not by protein count.** The measured
overshoot against a linear estimate ran from 1.2× (ncbifam) to 2.5× (kofam),
so doubling the linear estimate past about 100k proteins is a fair planning
rule and an optimistic one for the worst stage.

MMseqs2 is the exception worth noting — 455,571 proteins clustered into 49,347
families in 109 seconds, **7.6×** for 11.9× the proteins, because clustering
scales with redundancy rather than with count.

### The other thing: when a stage starts is not when it is reached

On that same run, with `stage_workers: 3`, the stages started here:

| start | stage | why then |
|---|---|---|
| 0 s | `emapper`, `pfam`, `dbcan` | first three **in table order** |
| 76 s | `diamond` | emapper finished |
| 291 s | `signalp` | diamond finished |
| 617 s | `tmbed` | dbcan finished |
| 7.2 h | `cluster`, `ncbifam` | pfam finished |
| 12.8 h | `kofam` | ncbifam finished |
| **27.1 h** | **`interpro`** | kofam finished |

InterProScan is the longest stage in the pipeline and it started **27 hours
in**, because it is tenth in a table whose first three entries include one
that takes ten minutes. Since v0.4.0 the scheduler sorts each round's ready
stages longest-first, so the hours-class stages claim the workers and dbcan
waits instead. It is still worth giving InterProScan its own pass at this
scale — see below — but it no longer queues behind a ten-minute stage.

### Sizing your run

**Measure the database first.** `grep -c '^>' proteins.faa` decides everything
below.

| database | what to do |
|---|---|
| **under ~50k** | run everything; a full pass is hours |
| **50k–150k** | run everything, but expect a day; the default `stage_workers: 4` is fine |
| **over ~150k** | read the rest of this section before starting |

Past roughly 150k proteins, three decisions matter more than any tuning:

1. **Give InterProScan its own pass.** It is the longest stage by a wide
   margin. Since v0.4.0 the scheduler dispatches the hours-class stages before
   the short ones, so it no longer starts *behind* `dbcan` and `diamond` — but
   with only `stage_workers` running at once it still may not get a slot in
   the first round, and when it shares it gets `threads / stage_workers` CPU.
   Running it alone gives it the whole machine:
   ```bash
   python metaannot.py run --config config.yaml --only interpro
   python metaannot.py run --config config.yaml          # everything else
   ```
   On a 455k run, letting it start last cost most of a day.

2. **Ask whether the whole database needs annotating.** The expensive stages
   run over `proteins_faa`, but only the proteins that survive to the quant
   table can reach the report. On one real TMT run that was **8,668 razor
   proteins out of 455,571 — 1.9%**. The other 98% were searched only to
   populate a database-wide bin table.

   `subset` reduces the FASTA to what a quant table references, and is worth
   trying — but check what it gives you before relying on it. On a smORF
   database it saved nothing (455,571 of 457,611), because short peptides map
   almost everywhere. And restricting the set is **not free**:
   `peptide_assignment: taxon_unique` needs a taxid for every shared-peptide
   *candidate*, so a smaller universe turns features into
   `shared_unknown_taxon` and drops them. eggNOG is the cheap stage that has to
   stay broad; the search stages are the expensive ones that need not.

3. **Consider fewer concurrent stages, not more.** The runs measured above used
   `stage_workers: 3` — not the default of 4 — on 22 cores, so three tools took
   7 each. If one of them takes 10 anyway (tmbed does), everything else is
   squeezed. Lowering to 2 does not reduce total CPU
   work, but it does make each stage finish sooner, which matters because
   several stages write nothing until they are done.

### Stages that write nothing until they finish

The `hmmsearch` stages buffer their output and write once at the end. On a
large database that means hours with an empty output file, which cannot be
told from a hang — check CPU time (`ps -o time -C hmmsearch`) rather than file
size. A crash loses the whole stage.

**TMbed is the extreme case**, and since v0.4.0 the tool works around it. On
the 455,571-protein run it wrote nothing at all for the first 31 hours: not a
progress bar, not a partial file. Two earlier attempts died at 2 h 36 min with
nothing recoverable. So the input is now split into chunks of about
`tmbed_chunk_residues` (5,000,000 by default, roughly 17k average proteins),
each committed as it lands:

```yaml
tmbed_chunk_residues: 5000000   # 0 = one invocation, the old behaviour
tmbed_allow_partial: false      # true: finish without the failed chunks
tmbed_max_consecutive_failures: 2
```

An interrupted run resumes from the last finished chunk instead of starting
over, and the log says which chunk it is on and roughly how long is left. The
chunks are length-sorted longest-first, so the sequences most likely to
exhaust the card are in the *first* chunk, where that costs minutes to
discover rather than the whole stage. The price is one ProtT5 load per chunk,
which is why the budget is millions of residues and not thousands — under 5M
there is exactly one chunk and nothing changes.

A chunk that dies keeps whatever TMbed managed to write, and what came back is
reconciled against what was handed over — from the files, because a chunk can
exit 0 and still be short. Anything missing is named in
`results/topology/tmbed_failed.tsv`, and by default that is fatal, so nothing
downstream reads a short topology set by accident.

`esmfold` checkpoints every model as it is written, so a machine that dies at
protein 1,800 of 1,900 costs the tail and nothing else. Each structure is
**renamed** into place rather than written in place — a run killed between the
open and the flush used to leave a 0-byte `.pdb` that every later run counted
as folded, and a file with no coordinate line in it is now refolded rather
than trusted.

---

## What has actually been run

**Two real datasets, end to end.**

A label-free FragPipe run — 38,204 proteins, 122,278 features, 36 samples in six
groups — went through fifteen of the twenty-one stages: `emapper` (reuse),
`pfam`, `dbcan`, `diamond` over five databases, `cluster`, `ncbifam`, `kofam`,
`interpro`, `signalp`, `tmbed`, `esmfold`, `foldseek`, `integrate`, `finalise`
and `join`. The report knitted and the QFeatures object built from that run on
R 4.6.1, not from synthetic data. Final bins were 47.6 / 28.2 / 19.0 / 1.0 /
0.4 / 3.8%, and the dark bin fell from 9,731 proteins on eggNOG alone to 1,462.

A FragPipe **TMT** run followed on a second dataset — 8 plexes, 88 channels, 75
biological samples, 455,571 proteins — so `quant_format: fragpipe_tmt` is not a
paper path either.

**What still has not run on real data:** `smorf`, `context`, `hhblits`,
`jackhmmer`, `unipept` and `taxonomy`, all off by default. So the remote-profile
searches and the peptide-LCA taxonomy comparison remain unexercised outside
their output parsers, and `3p_profile_only` has never been populated by any run
— it is fed only by `hh_hit` and `jackhmmer_hit`.

A test suite ships in `tests/`, so the plumbing described in this tutorial is
reproducible offline. What does not ship is a benchmark script, and the datasets
themselves are not redistributable, so the timings above cannot be re-derived
from what you have.
