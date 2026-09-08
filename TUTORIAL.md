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

**TMT is different and still needs a hand-written design**: there `experiment`
is the plex, so the derived contrasts would be plex-versus-plex. Write
`analysis.metadata` by hand (sample, condition, batch) — and note that
reporter-ion quantification is not supported at all: a per-plex FragPipe table
is now refused outright by the loader. See README, "FragPipe TMT output is not
supported".

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
`effector_predictions`), whose keys you choose — proof-read those by hand.

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

## Phase 6 — GPU stages, on the RTX laptop

`tmbed` and `esmfold` need CUDA. Run them there and bring the results back;
metaannot **adopts** outputs it finds rather than recomputing.

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
| `this looks like isobaric (TMT/iTRAQ) output` | a per-plex FragPipe TMT table | correct, and deliberate. Reporter-ion channels are not read; use a label-free or DIA-NN input |
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
| `requested structures exist` … `were skipped (OOM) or never folded` | `esmfold` finished and models are still missing | proteins over `max_len_structure` are counted separately in the same line; the remainder are OOM skips (`skipping <id>` in the esmfold log) or proteins added to `dark.faa` after the last fold |
| a stage runs even though its `run.<stage>` flag is false | it was named with `--only` | intended — `--only` is a clearer statement of intent than a flag left off for another machine, and the run warns. This is how the GPU box runs `tmbed`/`esmfold` against the server's config |
| `waiting for the GPU — tmbed is using it` | `esmfold` and `tmbed` both need the card | expected, not stuck. Each wants essentially a whole GPU (15.5 GB and 13.3 GB of a 16 GB device were observed), so they are leased one at a time; every CPU-only stage keeps running. `gpu_workers` raises the count, but they all still share the single `gpu_device` |
| `refusing to search a DIAMOND database that cannot answer` | a `.dmnd` that is empty, truncated, or holds no sequences | a failed `diamond makedb` leaves one behind, and searching it reports 0 hits — the same output as a real absence. Rebuild it with the `diamond makedb` line the message prints, or drop the tag with `<tag>: ""` |
| `incapable of a hit before it starts` | the database's sequences are too short for the configured e-value | a 15-residue peptide cannot reach `1e-10`, so that search's 0 hits say nothing about the biology. Give the tag its own threshold under `diamond_evalues:`, or record that this database needs different settings. The warning also names the `diamond_weights` entry that claims the database can score |

---

## Resource guide

| stage | cost | memory | notes |
|---|---|---|---|
| emapper (reuse) | minutes | low | streams; a 50 M-row table is fine |
| pfam / ncbifam | hours | moderate | scales with `--cpu` |
| diamond | hours | `-b` × 6 GB per job | databases run concurrently |
| interpro | **longest** | large JVM heap | disable on the first pass |
| cluster | minutes–hours | `--split-memory-limit` | |
| tmbed / esmfold | GPU-bound | 16 GB VRAM | laptop only; one at a time (`gpu_workers`) |
| foldseek | hours | large | serial across targets on purpose |
| integrate / finalise / join | minutes | ~1 GB per 100k proteins | |

No benchmark ships with the tool, so the non-tool costs above are estimates,
not measurements. Reading a large quant table is the one expensive non-tool
step (order of a gigabyte and tens of seconds per read, and a full run reads
three times). The wall time of a real run is essentially the search tools.

---

## First-run caveat

Phases 0–5 with the search stages **off** are the part that has been run on real
data: one label-free FragPipe dataset, 38,204 proteins and 122,278 features,
through `emapper` (reuse), `integrate`, `finalise` and `join`. Manifest parsing,
the id join, the shared-peptide rule, binning, the roll-up and the recovered
design all worked there.

What has still never run on real data is the search tools themselves —
hmmsearch, DIAMOND, InterProScan, KOfamScan, HHblits, jackhmmer, ESMFold,
Foldseek, SignalP/TMbed — and the Unipept API. Only their output parsers are
tested, and `3s_structure_only` / `3p_profile_only` have never been populated
from a real search.

The **report knits and the R object builds** on synthetic data under R 4.3 with
pandoc 3.1, rmarkdown, limma, SummarizedExperiment and the tidyverse — synthetic
only, never on a real result. None of that is reproducible from what you have:
no test file, benchmark or synthetic generator ships with `metaannot.py`. The
QFeatures assay-link branch of the object script has not been exercised; when it
fails it prints a generic "QFeatures version differs" and falls back to a
SummarizedExperiment, so check the class of what you get back.

So anything past phase 5, and any run that enables a search stage, is still
validation. Check counts at every phase, and when something looks wrong, say so
and stop rather than pressing on.
