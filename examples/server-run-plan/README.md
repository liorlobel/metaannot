# Server run plan — eight metaproteome datasets, full stage set

Generated 2026-09-07 on the MacBook, from the eggNOG-only baseline runs. Nothing here has
been executed on the server; every path is a claim to be checked by `doctor` before
anything long starts.

## Why

The eggNOG-only baselines showed `4_dark` is small in the catalogue datasets (2–4%), but
that is not the same as annotated: **14–28% of proteins have no KO at all** and **27–44%
are invisible to KEGG pathway enrichment**, sitting in `3_annotated_no_ko` and
`2_ko_orphan` rather than in the dark bin. Those are proteins eggNOG-mapper *did*
annotate — a COG category, a description, a Pfam — but never gave a KO.

That population is what the other search programs address, and the rescue is measured,
not assumed. On the one dataset where all eight sources have run, KOfamScan assigned KOs
to 21,926 proteins, **8,129 of which eggNOG missed entirely**, 2,598 of them reaching a
specific KEGG map. eggNOG assigns KOs by DIAMOND search; KOfam uses per-family HMMs with
adaptive thresholds.

Expect these datasets to respond differently from UC. There, most of the movement was
`4_dark` collapsing 29.1 → 4.2%. Here there is almost no dark bin to collapse, so the
movement has to come from `3_annotated_no_ko` → `1_ko_pathway`, which is exactly the
KOfam-versus-eggNOG axis.

| dataset | identified | no KO | KEGG-invisible | rescuable |
|---|---|---|---|---|
| Results_Combined_UC_metaprotoemics | 674,327 | 23.6% | 41.7% | 159,187 |
| Acute_colitis_cohort | 571,550 | 23.9% | 41.6% | 136,685 |
| UC_nitsan_2023_1st_batch | 418,384 | 19.0% | 36.8% | 79,466 |
| DDA_mouse_microbiome_acetylation2 | 282,450 | 19.8% | 36.2% | 55,968 |
| Heyer_IBD_metaproteomics | 242,724 | 17.8% | 33.4% | 43,263 |
| DDA_mouse_microbiome_acetylation | 233,668 | 17.9% | 33.8% | 41,747 |
| mouse_dietary_intervnetion | 146,244 | 27.7% | 44.2% | 40,456 |
| DDA_mouse_microbiome_nitrosylation | 160,590 | 14.3% | 27.4% | 23,010 |

## Check these three things first

**1. The tool.** The UC run's log reported `metaannot 1.0.0`, which matches no release —
v0.1.0 and v0.2.0 are the only ones that exist. Whatever is on the server is of unknown
provenance. Put v0.2.0 there before running anything:

```bash
cd ~/metaannot && git fetch --tags && git checkout v0.2.0
python metaannot.py --version        # must print 0.2.0
```

v0.2.0 matters for this run specifically: v0.1.0 dies on a single non-UTF-8 byte in a
tool's output, which is what killed `integrate` on the UC run after InterProScan had
already spent three hours.

**2. The mount.** Every path assumes the Extreme Pro drive is at `/mnt/d` (`D:\` in
Windows), which is how the UC run addressed it. If it is elsewhere:

```bash
sed -i 's#/mnt/d#/your/mount#g' */config.yaml */FASTA_PATH
```

**3. Four database paths are inferred.** Only `ncbifam_hmm` and `diamond.vfdb` were
confirmed from the UC log. `pfam_hmm`, `dbcan_hmm`, `kofam_profiles`, `kofam_ko_list` and
`interproscan_sh` are marked `# INFERRED` in each config. If `doctor` says MISS for one,
the database is probably present under another name — **fix the path in the config
first**, and only download if it is genuinely absent.

## Running it

```bash
rsync -av ~/metaannot-server-plan/ server:~/metaannot-runs/     # from the MacBook
```

```bash
cd ~/metaannot-runs
./00_subset.sh          # identified-protein FASTAs; reads ~46 GB, minutes
./01_doctor.sh          # checks tools, databases, manifest→column mapping. Downloads NOTHING.
                        # Read every install.sh it writes before running one.
tmux new -s metaannot
./02_run.sh 2>&1 | tee run_$(date +%F_%H%M).log
./03_collect.sh         # one table across all eight, plus the KOfam control
```

`00_subset.sh` skips a dataset whose FASTA already exists, so it is safe to re-run.
`02_run.sh` is sequential on purpose: each dataset already takes 22 cpu / 80 GB across 3
concurrent stages, and four hmmsearch jobs against Pfam-A alongside InterProScan is the
usual squeeze. It runs the small datasets first, so a path mistake surfaces in minutes
rather than after the 674k-protein set.

If a dataset dies, just re-run `02_run.sh`: finished stages are cached and skipped, and
only the failed stage and its dependents rerun. **Do not add a bare `--force`** — it
discards every cached stage, which here is days. To redo one stage on purpose:
`metaannot.py run --config <d>/config.yaml --force --only <stage>`.

## Cost

InterProScan is the long pole: 2.8 hours for 38,204 proteins on this server. These sets
are 146k–674k proteins, so it dominates everything else and the total is **days, not
hours**. Two ways to cut it, both a real choice rather than a default:

- Set `run.interpro: false`. You keep the KOfam control, which is the point of the run,
  and lose Gene3D/SUPERFAMILY — the structure-derived signatures that catch what Pfam
  misses.
- Annotate only what gets quantified. `subset` currently takes every candidate protein
  per peptide, which in a redundant catalogue is 2–8 proteins per feature; only 5k–37k of
  them ever reach a protein group. Annotating the quantified set instead would cut the
  work roughly tenfold, at the cost of bin percentages that describe the quantified
  subset rather than the identification.

## What is off, and why

Each is a decision recorded in the configs rather than an oversight: `topology` (SignalP
6.0 is licence-gated and not installed; tmbed needs a GPU), `structure` (ESMFold needs a
GPU, Foldseek AFDB50 is ~123 GB down / ~200 GB on disk), `hhblits` and `jackhmmer` (they
query `dark_all.faa`, which is 2–8% here, and UniRef50 is 8.8 GB down / 27 GB on disk),
`context`, `smorf`, `effectors`, `unipept`, `taxonomy` (each needs an input that does not
exist yet).

`min_features_per_protein` is left at the documented default of **1**. The earlier UC run
set it to 2, which dropped 8,201 of 24,271 proteins before the report saw them —
disproportionately the small, KO-less ones this tool exists to study, and the reason no
dark protein reached its statistics. See issue #5.

Only `vfdb` is configured under `db.diamond`, matching the UC run. Listing it alone means
the four stock defaults (merops, card, tadb, bagel) are **not** used — `deep_merge`
replaces that block rather than merging it, and warns. Add them if you have them, and
give each an entry in `diamond_weights` or it scores 0.

## Afterwards

The report and the R object need R, which is on the MacBook rather than the server. Bring
the results back and generate the report locally — do **not** rsync the Rmd, because the
one written on the server embeds absolute server paths in its YAML header:

```bash
rsync -av server:~/metaannot-runs/<dataset>/results/ ./results/
python metaannot.py report --config <dataset>/config.yaml
python metaannot.py object --config <dataset>/config.yaml
```
