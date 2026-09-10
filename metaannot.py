#!/usr/bin/env python3
"""
metaannot — single-file metaproteome functional annotation pipeline.

Rescues the fraction of an identified metaproteome that KEGG-based analysis
discards, by binning every protein on what evidence actually exists for it.
The bins are listed here in the order the classifier tests them — the first
match wins, so a protein with both a fold and a profile hit is 3s, never 3p:

  1_ko_pathway       KO, mapped to a specific (non-global) KEGG map
  2_ko_orphan        KO, but no specific pathway map
  3_annotated_no_ko  no KO; informative Pfam / CAZy / VFDB / MEROPS / TADB / ...
  3d_duf_only        no KO; only domain evidence is a DUF
  3s_structure_only  no sequence annotation; confident Foldseek hit
  3p_profile_only    no sequence annotation and no fold; only a remote
                     profile hit (HHblits / jackhmmer)
  4_dark             no evidence

This list, BIN_ORDER below, and the report's BIN_LEVELS/BIN_COLS are one
vocabulary; a mismatch is checked for at import.

Usage
  metaannot.py init                      write a config template
  metaannot.py doctor                    check tools, databases and inputs
  metaannot.py subset --db D --quant Q   build the identified-protein fasta
  metaannot.py run                       run every enabled stage, resumably
  metaannot.py run --dry-run             show the plan without executing
  metaannot.py run --only pfam diamond   run named stages
  metaannot.py run --from integrate      run from a stage onwards
  metaannot.py run --force               ignore cached stage state

Requires: python3, pandas, pyyaml. External tools are needed only by the
stages that use them; `doctor` reports which are missing.
"""

__version__ = "0.4.0"

# Bumped only when the MEANING of a stage's output changes, so that existing
# results become genuinely invalid. It is deliberately not __version__: tying
# the cache to the tool version means a patch release that fixes a log message
# discards days of InterProScan and Foldseek compute.
SIGNATURE_VERSION = 1

# The evidence bins, in the order build_annotation() tests them. This is the
# single source of the bin vocabulary: the module docstring above, the report's
# BIN_LEVELS/BIN_COLS and bin_summary.tsv all take their order from it, because
# a table ordered 3p-before-3s next to a classifier that decides 3s-before-3p
# reads as if the two disagree about what a protein is.
BIN_ORDER = ("1_ko_pathway", "2_ko_orphan", "3_annotated_no_ko",
             "3d_duf_only", "3s_structure_only", "3p_profile_only", "4_dark")

import argparse
import atexit
import concurrent.futures
import contextlib
import fnmatch
import glob
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
from collections import defaultdict, deque

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("metaannot requires numpy and pandas:  pip install numpy pandas")
try:
    import yaml
except ImportError:
    yaml = None


# ======================================================================
# configuration
# ======================================================================
DEFAULT_CONFIG = {
    "proteins_faa": "input/proteins.faa",
    "quant_table": "input/report.pg_matrix.tsv",
    # diann | fragpipe            protein-level tables
    # fragpipe_peptide | fragpipe_ion   combined_peptide.tsv / combined_ion.tsv
    # msstats_csv                 FragPipe MSstats.csv (long, feature level)
    # msstats_feature             MSstats dataProcess()$FeatureLevelData
    # msstats_protein             MSstats dataProcess()$ProteinLevelData
    # fragpipe_tmt                FragPipe TMT: the per-plex TMTn/ directories
    #                             (quant_table is then the RUN DIRECTORY that
    #                             holds them, not a file). See the tmt: block.
    "quant_format": "diann",
    # FragPipe .fp-manifest. When set, sample names, conditions and
    # bioreplicates all come from it: quant columns are renamed to the
    # manifest's sample names and results/quant/design_from_input.tsv is
    # written for the report.
    "manifest": "",
    # Feature-level input only. Shared peptides are the central problem in a
    # strain-redundant metagenome database, so the rule is explicit.
    #   protein_unique         only features matching exactly one protein
    #   taxon_unique           (default) also features whose candidates all
    #                          share one taxon
    #   taxon_or_family_unique taxon_unique, plus a fallback for features whose
    #                          candidates cannot be shown to share a taxon
    #                          because at least one of them has none: if they
    #                          all share one MMseqs family_id from the cluster
    #                          stage, they are treated as one unit. A family is
    #                          a SEQUENCE CLUSTER at cluster_min_seq_id, NOT an
    #                          organism — two strains' orthologues and two
    #                          paralogues in one genome both land in one family.
    #                          Opt in only when the alternative (dropping the
    #                          shared peptides of exactly the unannotated
    #                          proteins this tool studies) is worse for the
    #                          question being asked. Needs run.cluster on.
    #   razor                  keep everything, assign to the razor protein
    "peptide_assignment": "taxon_unique",
    # How assigned features become a protein number.
    #   sum            (default) plain sum over observed features. Simple and
    #                  what every published metaannot number was computed with.
    #                  A peptide seen in only some samples still contributes
    #                  its whole intensity to those samples, so missingness
    #                  becomes fold change.
    #   median_polish  Tukey median polish in log2 space (the MaxLFQ argument):
    #                  each feature's own response level is removed before the
    #                  samples are compared, so a peptide present in a subset
    #                  of samples cannot by itself create a ratio. The result
    #                  is rescaled to preserve each protein's summed intensity,
    #                  so magnitudes stay comparable with "sum" but ratios do
    #                  not. Changes every quantified number: opt in, and say so
    #                  in any methods section.
    "rollup_method": "sum",                 # sum | median_polish
    # 1, not 2: most proteins <= 100 aa yield a single peptide, so a default
    # of 2 removed the population this tool exists to study before anything
    # could report it. Filter in the report instead (analysis.min_features),
    # where the drop is visible per bin.
    "min_features_per_protein": 1,
    # The unipept and taxonomy stages need only peptide sequences and the
    # proteins each peptide could come from — never an intensity. Reading
    # them through the full quant reader makes them fail for reasons that
    # have nothing to do with peptides (no intensity column, an unmappable
    # manifest, isobaric input this tool refuses to quantify).
    #   auto   full reader first, peptide-only reader if it refuses (default)
    #   always peptide-only reader, never touching the quant columns
    #   never  full reader only; a quant problem aborts taxonomy as before
    "peptide_only_reader": "auto",
    # FragPipe writes 0 for "not quantified", not for "measured as zero".
    # Summed as a real zero it turns missingness into fold change.
    "zero_intensity_is_missing": True,
    # Decoys and contaminants travel in the same FragPipe table. They must not
    # become quantified "proteins" and must not veto a shared peptide.
    # Decoy, entrapment, contaminant and host ids travel in the same FragPipe
    # table. They must not become quantified "proteins", must not veto a
    # shared peptide, and are unannotated by construction so they do not
    # belong on the dark work-list either.
    "exclude_id_prefixes": ["rev_", "decoy_", "ent_", "contam_", "Cont_",
                            "CON__", "HUMANHOST_"],
    # Write the feature-level matrix as well as the roll-up. It is the raw
    # material for the peptide assay of the R object; set false if the file
    # size is a problem.
    "export_feature_quant": True,
    "feature_intensity_suffix": "Intensity",
    "feature_exclude_suffixes": ["MaxLFQ Intensity", "Spectral Count"],
    # quant_format 'fragpipe_tmt' only. Reporter intensities are read from the
    # PER-PLEX tables, never from tmt-report/ (see read_fragpipe_tmt for why).
    "tmt": {
        # Plex directories inside quant_table. Sorted naturally, so TMT10
        # follows TMT9 rather than TMT1.
        "plex_glob": "TMT*",
        # Which per-plex table: ion.tsv keeps the modified sequence and the
        # charge state apart, peptide.tsv is one row per sequence.
        "level": "ion",                     # ion | peptide
        # FragPipe writes '<PLEX>_annotation.txt' (TMT1/TMT1_annotation.txt),
        # never a plain 'annotation.txt'. A {plex: path} map overrides the
        # pattern when the files live somewhere else.
        "annotation": "{plex}_annotation.txt",
        # The reference/bridge channel, resolved PER PLEX because it does not
        # sit at a fixed position: in a real 8-plex design the pool is at 131C
        # in six plexes and at 131N in the other two, so no single channel
        # names it. reference_name globs the ANNOTATED SAMPLE NAME ('Pool*'),
        # which is the stable signal; reference_channel globs the channel
        # ('131C'). Set one, not both.
        "reference_name": "",
        "reference_channel": "",
        # How the reference is treated. False (the default) is the COVARIATE
        # treatment: the reference is dropped from the sample columns, because
        # a pooled bridge is not a biological sample and would otherwise get a
        # condition and a row in the design; the plex stays in the model to
        # absorb the batch. True is the RATIOS treatment: every channel of a
        # plex is divided by that plex's reference before the roll-up, which
        # removes the plex effect directly (the classic bridge design) but
        # assumes the same pool went into every plex and throws away the
        # reference's own variance. It also propagates the reference's
        # missingness: a feature with no reference value in a plex becomes NA
        # for that whole plex, and the reader reports how many values that
        # cost. Either way the reference is never a modelled sample.
        "use_reference_ratios": False,
        # Where the CONDITION comes from. The plex is a batch and is never
        # used as one, so without this the condition has to be written by
        # hand into analysis.metadata.
        #   "auto" - try to split the annotated sample names ("resp_01" ->
        #            "resp"), and accept the split only when it is
        #            unambiguous; otherwise say so and require the metadata.
        #   ""     - never derive; always require analysis.metadata.
        #   regex  - a regular expression with one capture group, applied to
        #            each sample name; the capture is the condition.
        "condition_from_name": "auto",
        # Keep only features identified in at least this many plexes. 1 keeps
        # everything, which is the honest default: cross-plex overlap is low
        # (~45% of ion keys between two plexes), so raising this trades
        # features for a fuller matrix.
        # This is the FEATURE-level filter, applied here before the roll-up.
        # analysis.min_plexes is the PROTEIN-level one, applied in the report
        # beside min_valid_per_group; they are different questions and both
        # are reported.
        "min_plexes": 1,
        # Bring the channels of one plex to a common scale before the join.
        # Channels differ by how much peptide was loaded and how completely it
        # was labelled, which is a per-channel constant carrying no biology,
        # and the roll-up sums features across it. "median" divides each
        # channel by its own median and multiplies by the plex's median
        # channel, so the values stay linear (the size factors need that), the
        # step is exactly a per-channel median centring after the report's
        # log2, and the BETWEEN-plex difference is deliberately left alone for
        # the plex term to absorb. "none" leaves FragPipe's numbers untouched.
        "within_plex_normalise": "median",   # median | none
        # FragPipe names an unassigned channel '<PLEX>_<CHANNEL>' in the
        # annotation. That is not a sample: its signal is isotope-impurity
        # carry-over from the neighbouring channels.
        "drop_empty_channels": True,
        # Drop features whose precursor was too co-isolated to trust. Purity
        # is written ONLY into psm.tsv, so this is implemented by reading that
        # file and aggregating it onto the feature key: a feature is judged by
        # the MEDIAN purity of the PSMs that produced it, because its reporter
        # intensities are a sum over those PSMs and no single one describes
        # it. 0 disables the filter and psm.tsv is not read at all. A feature
        # that matches no PSM row is KEPT and counted, never dropped — an
        # unmatched key is a join failure, not a dirty spectrum. This limits
        # how much co-isolation contributed; it does not correct for it.
        "min_purity": 0,
    },
    # Intensity columns are auto-detected (numeric, not a known metadata
    # column) and always logged. Override when auto-detection is wrong —
    # an unrecognised numeric metadata column would otherwise be summed as
    # if it were sample signal.
    "intensity_columns": [],   # explicit list wins
    "intensity_regex": "",     # else a regex on the column name
    "results_dir": "results",

    # How a stage decides that one of its inputs changed. False (the default,
    # and what every result so far was produced with) digests a file over 64
    # MB at its two ends plus its size: reading a 36 GB Foldseek or UniRef
    # index end to end on every signature would cost more than the stage it
    # protects, and a rebuilt database changes both ends. It cannot see an
    # edit that keeps a huge file's size and both 8 MB edges — an in-place
    # patch to the middle of a database. Set true to hash every byte and
    # catch that; expect the first run afterwards to rehash everything, and
    # every stage whose digest changes to recompute.
    "full_content_digest": False,

    "emapper_precomputed": "",
    "emapper_id_transform": "exact",
    # A search DB may prefix the ids it was built from
    # (uhgpSM_, HUMANHOST_) while the eggNOG table has not;
    # stripped off the fasta id when matching. str or list.
    "emapper_strip_id_prefix": "",
    "emapper_min_coverage": 0.50,
    "emapper_warn_coverage": 0.90,
    "emapper_diagnose_rows": 2_000_000,

    # Inputs of the context and smorf stages. Their ids must match the ids in
    # proteins_faa, or every neighbourhood comes back empty.
    "gff": "",
    "contigs_fna": "",
    "context_window": 5,        # neighbours on EACH side, counted in GENES
    "pul_min_cazymes": 2,       # CAZymes inside that window to call a PUL
    # aa, the same small-protein cap the rest of the tool uses; 150 admitted
    # so many ordinary proteins that the +2 immunity weight meant little.
    "immunity_max_len": 100,
    "immunity_max_gap": 60,   # bp between the two CDS

    "threads": 8,
    # How often run_cmd echoes the newest line a running tool has written to
    # stderr, in seconds. 0 turns it off. tmbed ran for 2 h 36 min and
    # InterProScan for 2.9 h with nothing on the log, because stderr was
    # captured and only quoted on failure: the only way to tell either apart
    # from a hang was to watch its CPU ticks accumulate in /proc. Both write a
    # tqdm bar to stderr the whole time, so one line a minute is the
    # difference between a silent process and a progress bar. Nothing is
    # hoarded - see run_cmd.
    "progress_interval_s": 60,
    # How many independent stages may run at once. Each gets threads //
    # stage_workers CPUs. Parallel stages multiply peak memory.
    "stage_workers": 4,
    # Total memory budget in GB, split across concurrent stages exactly as
    # threads is. 0 means auto-detect. Propagated to the flags that actually
    # control memory: diamond -b, mmseqs/foldseek --split-memory-limit,
    # hhblits -maxmem, InterProScan's JVM heap, and whether eggNOG-mapper may
    # use --dbmem.
    "ram_gb": 0,
    # eggNOG --dbmem loads the whole annotation database into RAM. Only used
    # when the stage's own budget reaches this, because on a smaller machine
    # it swaps instead of speeding anything up.
    "emapper_dbmem_min_gb": 64,
    # Appended verbatim to the named tool's command line.
    "tool_args": {},

    # Download URLs, in the config so a moved link is an edit here rather than
    # a patch to the tool. `doctor --fix` uses these; check them against the
    # provider's current release before a big download.
    "sources": {
        "pfam": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
        "taxdump": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        "uniref50": "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz",
        "kofam_profiles": "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz",
        "kofam_ko_list": "https://www.genome.jp/ftp/db/kofam/ko_list.gz",
        # hmm_PGAP.LIB, not hmm_PGAP.HMM.tgz: the .tgz unpacks to ~19,000
        # single-model files that hmmpress cannot take.
        "ncbifam": "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB",
        # Both of these were dead: the old dbCAN link 302s to an HTML landing
        # page, which `curl -fL` accepts with exit 0 and writes to disk as if
        # it were the database, and the old InterProScan directory 404s.
        # bcb.unl.edu itself now redirects to pro.unl.edu, so the download host
        # is named directly rather than trusting a chain of redirects to end
        # somewhere that is still the database and not a web page.
        "dbcan": "https://pro.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt",
        "interproscan": "https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/5/",
        "diamond": {
            "vfdb": "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz",
        },
    },
    "gpu_device": 0,
    # How many GPU stages may run at once. The CPU and RAM budgets are split
    # between concurrent stages; the GPU was not modelled at all, and tmbed and
    # esmfold each want essentially a whole card — tmbed held 15.5 GB of a
    # 16 GB device and ESMFold peaked at 13.3 GB on a single short sequence.
    # With run.topology and run.structure both on they cannot share it, and
    # the loser dies of an OOM that names no cause.
    # This is a LEASE COUNT, not a device map: every GPU stage is pinned to the
    # single `gpu_device` above, so raising this on a two-card machine runs two
    # stages on the SAME card rather than one per card. Per-device assignment
    # is not implemented; the limit here is a policy, not the hardware.
    "gpu_workers": 1,
    "max_len_structure": 700,
    "esmfold_chunk_size": 64,
    "max_dark_structures": 2000,
    # A CUDA fault partway through folding is not a reason to throw away the
    # structures already on disk. Every PDB is checkpointed as it is written,
    # so a failure at protein 1800 of 1900 has cost nothing but the tail.
    # With this on, ESMFold skips the sequences it cannot fold, records them,
    # and lets Foldseek search what did fold. Off, it stops and says how many
    # it has, so a rerun resumes rather than a partial result passing silently.
    "esmfold_allow_partial": False,
    # Cap the fold work-list by the card's ACTUAL free VRAM, measured once the
    # weights are resident. ESMFold does not raise an OOM when it stops
    # fitting - the driver pages device memory to host RAM and the fold simply
    # becomes one to two orders of magnitude slower, so a run does not fail,
    # it stops being finishable. See vram_fit_length for the measurements.
    "esmfold_vram_cap": True,
    # Empirical, and calibrated to reproduce one measurement exactly: on a
    # 16 GB card holding the 11.2 GB fp32 trunk, 4.8 GB was free and 478 aa
    # was the longest length that still folded at the smooth rate. With the
    # 0.5 GB reserve below, this coefficient puts the cap at that 478. It is a
    # config key because that is one card, not a law: raise it to be more
    # conservative, lower it if a card demonstrably folds longer sequences at
    # a smooth rate.
    "esmfold_bytes_per_residue_pair": 20200,
    # Left free for the driver, the display and fragmentation.
    "esmfold_vram_reserve_gb": 0.5,
    # Consecutive failures that mean the device is wedged rather than the
    # sequence being hard. Folding does not recover from that on its own, so
    # stop instead of walking the rest of the list failing every one.
    "esmfold_max_consecutive_failures": 5,

    "run": {
        "eggnog": True, "pfam": True, "dbcan": True, "diamond": True,
        "cluster": True, "join": True,
        # Off by default because switching them on is a decision about someone
        # else's hardware: structure is ESMFold on a GPU plus a ~1 TB Foldseek
        # target, topology needs SignalP 6 and tmbed installed. With them on,
        # `run all` on a CPU box died in these two after every cheap stage had
        # already burned its hours.
        "topology": False, "structure": False,
        "context": False,
        "unipept": False, "taxonomy": False,
        # Added evidence, all off by default: each needs its own tool or
        # database, and doctor reports exactly which.
        "ncbifam": False,     # NCBIfam/TIGRFAM HMMs — cheapest coverage gain
        "kofam": False,       # KOfamScan: a CONTROL on eggNOG's KO calls
        "interpro": False,    # InterProScan: Gene3D, SUPERFAMILY, PANTHER, ...
        "hhblits": False,     # profile-profile, the real answer to "past BLAST"
        "jackhmmer": False,   # iterative profile search, cheaper alternative
        "smorf": False,       # small-ORF calling on contigs (see the warning)
    },

    "db": {
        "eggnog_data": "/data/db/eggnog",
        "pfam_hmm": "/data/db/Pfam-A.hmm",
        "dbcan_hmm": "/data/db/dbCAN-HMMdb-V13.txt",
        # One or more Foldseek targets. AFDB50 is mostly unreviewed; PDB and
        # Swiss-Prot carry far better functional annotation per hit.
        "foldseek_target": "/data/db/foldseek/afdb50",
        "foldseek_extra_targets": [],
        # Directory holding NCBI nodes.dmp and names.dmp, for resolving eggNOG
        # seed taxids to a lineage offline. Without it the eggNOG/Unipept
        # comparison degrades to exact taxid identity.
        "ncbi_taxonomy": "",
        # The full NCBIfam/TIGRFAM library. NCBIfam-AMRFinder.LIB is the
        # AMR-only ~1k-model subset and was never what sources.ncbifam fetched.
        "ncbifam_hmm": "/data/db/hmm_PGAP.LIB",
        "kofam_profiles": "/data/db/kofam/profiles",
        "kofam_ko_list": "/data/db/kofam/ko_list",
        "interproscan_sh": "/opt/interproscan/interproscan.sh",
        "hhblits_db": "/data/db/hhsuite/pfam",
        "jackhmmer_db": "/data/db/uniref50.fasta",
        "diamond": {
            "vfdb": "/data/db/vfdb_core.dmnd",
            "merops": "/data/db/merops_scan.dmnd",
            "card": "/data/db/card_protein_homolog.dmnd",
            "tadb": "/data/db/TADB3_toxin.dmnd",
            "bagel": "/data/db/bagel4.dmnd",
        },
    },

    "unipept": {
        # Offline route (recommended): run Unipept yourself and point at the
        # output.  unipept pept2lca --equate --all -i peptides.txt -o out.csv
        "result": "",
        # HTTP route. Off by default: the API version, its field names and its
        # rate limits are outside this tool's control.
        "allow_http": False,
        "api_url": "https://api.unipept.ugent.be/api/v2/pept2lca.json",
        "equate_il": True,          # MS cannot distinguish I from L
        "batch_size": 100,
        "sleep": 0.3,
        "retries": 3,
        "timeout": 180,
        # Off: the suffix-array API behind pept2lca matches whole peptides,
        # missed cleavages included. Splitting therefore trades a species-level
        # answer for two short fragments that resolve to root or domain, and
        # those votes then outvote the protein's real consensus. Turn it on
        # only against an older Unipept that indexed fully tryptic peptides.
        "split_missed_cleavages": False,
        "consensus_min_fraction": 0.5,
        "consensus_min_peptides": 2,
    },
    # Which taxonomy drives taxon_unique roll-up and the taxon sums.
    #   eggnog     seed_ortholog taxid (default; always available)
    #   unipept    peptide LCA consensus (needs the unipept stage)
    #   concordant eggNOG taxid, but only for proteins where the two agree —
    #              proteins in conflict get no taxon and drop out of the
    #              taxon-based steps rather than being silently wrong
    "taxonomy_source": "eggnog",
    # Below this many proteins a within-taxon median is not robust either, so
    # the size factor falls back to the plain sum.
    "taxon_min_proteins_for_factor": 4,
    # A seed_ortholog taxid is a reference *strain*, so two ORFs of one gut
    # organism can carry different ones. Set species/genus/family to collapse
    # them through the taxdump first; "" keeps the raw seed taxid.
    # The default "" is kept because collapsing needs a taxdump that
    # db.ncbi_taxonomy often does not point at, and because changing it would
    # redefine "taxon" for every existing run. Read the consequence out of the
    # log: with "" the taxa in taxon_intensity.tsv / taxon_size_factors.tsv
    # are eggNOG REFERENCE-GENOME taxids, not organisms, so one gut species
    # can appear as several of them and each gets its own size factor.
    "taxon_rank": "",

    # Everything the R report needs. These become the Rmd's params, so the
    # whole tool is driven by this one config file.
    "analysis": {
        "metadata": "",                  # "" = use the design from the manifest
        "sample_col": "sample",
        "design_formula": "~ 0 + group",
        "factor_cols": "group",
        "numeric_cols": "",
        "block_col": "",
        "contrasts": "",                 # "" = every pairwise group contrast
        "primary_contrast": "",
        "msstats_comparison": "",
        "msstats_label": "",
        "min_features": 0,
        "min_valid_per_group": 3,
        # Protein-level companion to min_valid_per_group, for isobaric input.
        # min_valid_per_group counts SAMPLES, and an isobaric run's
        # missingness is shaped by the plex: a protein identified in one plex
        # only is all-NA in every other, so "3 valid values in every group"
        # can be satisfied entirely inside one batch and the difference the
        # model then reports is that batch. This counts PLEXES instead. 1 is
        # the default so that no run loses proteins to a filter it never
        # asked for; the report always says how many of the retained proteins
        # live in a single plex, so the number to set this on is in the
        # report whether or not the filter bites. Inert without a plex
        # column, so label-free runs are unaffected.
        "min_plexes": 1,
        "drop_zero_variance": False,
        "group_col_for_filtering": "group",
        "normalise": "median",
        "fdr": 0.05,
        "min_lfc": 0.585,
        # Kept equal to taxon_min_proteins_for_factor: with 5 here the report
        # dropped every taxon that could ever have used the sum fallback, so
        # the "fragile" flag it documents could never appear.
        "taxon_min_proteins": 4,
        "deviation_lfc": "",
        "assumption_frac_deviating": 0.5,
        "assumption_opposing": 0.6,
        # 3, not 1: a taxon flagged on one opposing protein is noise.
        "assumption_min_deviating": 3,
        "assumption_opposing_frac": 0.3,
        "family_min_size": 3,
        "require_taxonomy_concordance": False,
        "out_subdir": "analysis",
    },

    # InterProScan member databases; "" runs all installed ones (slow).
    "signalp_mode": "fast",       # fast | slow | slow-sequential
    "signalp_batch_size": 0,      # 0 = tool default
    "tmbed_batch_size": 0,        # residues per batch; 0 = tool default
    "tmbed_use_gpu": "auto",      # auto = GPU with CPU fallback | true | false
    # Proteins longer than this are excluded from tmbed and listed in
    # results/topology/tmbed_excluded.tsv. ProtT5's attention is
    # length-squared, so cost per sequence is roughly (len^2 x 32 x 4) bytes:
    # 1.1 GB at 3000 residues, 141 GB at titin's 34,350. One such sequence
    # aborts the stage on any device, hours in, with nothing written. 0
    # disables the cap and restores the old behaviour.
    "tmbed_max_len": 3000,
    # TMbed writes nothing until it finishes, so one invocation over a large
    # proteome risks days of work on a single process. The input is split into
    # chunks of about this many residues and each is committed as it lands, so
    # an interrupted run resumes instead of starting over. It costs one ProtT5
    # load per chunk, which is why the default is millions of residues and not
    # thousands: under 5M (roughly 17k average proteins) there is exactly one
    # chunk and the behaviour is what it always was. 0 disables the split.
    "tmbed_chunk_residues": 5000000,
    # A chunk that fails leaves its finished predictions behind and its
    # proteins in results/topology/tmbed_failed.tsv. False stops the run there
    # so nothing downstream reads a partial topology set by accident; True
    # accepts the shortfall deliberately.
    "tmbed_allow_partial": False,
    # A device that has stopped responding fails every remaining chunk the
    # same way, slowly. Stop after this many in a row rather than working
    # through all of them. 0 = never stop early.
    "tmbed_max_consecutive_failures": 2,
    "interpro_applications": "Pfam,NCBIfam,Gene3D,SUPERFAMILY,PANTHER,SMART,CDD,PIRSF",
    "hhblits_iterations": 2,
    "hhblits_workers": 4,
    "diamond_workers": 4,
    "jackhmmer_iterations": 3,
    "smorf_mode": "meta",         # smorf meta = assembly | single = isolate


    "diamond_weights": {"vfdb": 4, "tadb": 3, "bagel": 3, "merops": 2, "card": 0},
    # VFDB's own category, by its stable numeric code. A hit whose category is
    # not listed falls back to diamond_weights.vfdb. Empty disables the split
    # and every VFDB hit scores the flat weight, which is the old behaviour.
    # This map cannot resurrect a database the user switched off: with
    # diamond_weights.vfdb at 0 or absent, VFDB contributes nothing and the
    # categories are not applied, because the WARN, doctor and the README all
    # promise that a zero weight means the database does not score.
    #   VFC0086 effector delivery system   VFC0235 exotoxin
    #   VFC0001 adherence                  VFC0204 motility
    #   VFC0258 immune modulation          VFC0272 nutritional/metabolic
    #   VFC0282 stress survival            VFC0301 regulation
    "vfdb_category_weights": {
        "VFC0086": 4, "VFC0235": 4,          # exported effectors and toxins
        "VFC0001": 3,                        # adherence: surface, host-facing
        "VFC0204": 2,                        # motility
        "VFC0258": 2,                        # immune modulation: broad
        "VFC0282": 1, "VFC0301": 1,          # stress, regulation: mostly cytoplasmic
        "VFC0272": 1,                        # housekeeping in a virulence coat
    },
    # Per-database e-value, overriding thresholds.diamond_evalue for that tag
    # alone. One threshold cannot fit every database: a bacteriocin database
    # is built from peptides an order of magnitude shorter than VFDB's
    # 350-residue proteins, and no short alignment can reach 1e-10, so at the
    # pipeline default that search is close to incapable of a hit before it
    # starts. Empty = one threshold for all.
    #   diamond_evalues: {bagel: 1e-3}
    #
    # A caution the 0-hit BAGEL run on this pipeline earned: check WHAT was
    # built before tuning the threshold. Those 262 "sequences" of median
    # length 15 were BAGEL4's motif SEED set, not its bacteriocins, and no
    # e-value turns a blastp against 15-residue seeds into a bacteriocin
    # search. diamond_db_check now says so instead of offering this knob.
    "diamond_evalues": {},

    # Per-database minimum percent identity, overriding
    # thresholds.diamond_min_pident. Applied twice on purpose: as DIAMOND's
    # own --id during the search, so the table stays small, and again when the
    # table is read, so a file adopted from elsewhere or produced before the
    # floor existed keeps the same promise.
    #
    # 50 for these three because their hits are read as claims about a
    # PARTICULAR protein, not as a family assignment. On a 455,571-protein
    # run CARD returned 39,138 hits and VFDB 107,219 at the 30% default; above
    # 50% they are 4,661 and 21,623. The rest are the usual metagenome background
    # -- a 32%-identity match to a beta-lactamase over half a query is a hit
    # against the fold, not evidence that this protein confers resistance, and
    # reporting it as "carries an AMR gene" is the error a reviewer would find
    # first. BAGEL's peptides are short enough that a low-identity alignment
    # over them is close to meaningless.
    "diamond_min_pidents": {"card": 50, "vfdb": 50, "bagel": 50},

    "thresholds": {
        "diamond_evalue": 1e-10,
        "diamond_min_qcov": 50,
        "diamond_min_pident": 30,
        "diamond_strong_pident": 50,   # below this a hit scores half weight
        "dbcan_evalue": 1e-15,
        "dbcan_min_cov": 0.35,
        "foldseek_evalue": 1e-3,
        "foldseek_min_prob": 0.90,
        "foldseek_min_tmscore": 0.50,
        # Alignment coverage of the query. Only applied when the hit table
        # carries qlen; see FOLDSEEK_COLS.
        "foldseek_min_qcov": 0.50,
        # A short dark ORF is exactly where ESMFold confidence collapses, so
        # models below this mean pLDDT are not searched at all. 0 disables.
        "esmfold_min_plddt": 70.0,
        "foldseek_cluster_evalue": 1e-2,
        "foldseek_cluster_tmscore": 0.50,
        "foldseek_cluster_coverage": 0.80,
        "hhblits_min_prob": 90.0,
        "jackhmmer_evalue": 1e-5,
        "ncbifam_cutoff": "--cut_tc",
        "interpro_ignore_analyses": "MobiDBLite,Coils,Phobius,TMHMM,SignalP",
        "cluster_min_seq_id": 0.50,
        "cluster_coverage": 0.80,
        "smorf_max_len": 100,
        "lpxtg_min_offset": 15,
        "lpxtg_max_offset": 60,
    },

    # Effector priority weights. Every entry rewards positive evidence except
    # "no_ko", which rewards the ABSENCE of a KO. That looks wrong in
    # isolation and is kept deliberately: both populations the score is used
    # on are KO-less by construction — the report's shortlist filters
    # has_ko == FALSE before it sorts, and dark.faa is drawn from
    # 4_dark / 3d_duf_only / 3p_profile_only, none of which can carry a KO —
    # so no_ko is a constant inside each list and changes no ranking. It only
    # sets the offset between an effector_score printed here and one computed
    # over the whole table. Set it to 0 if you want scores that are comparable
    # across annotated and unannotated proteins; that shifts every KO-less
    # protein by -1 and will not match previously published scores.
    "weights": {
        "signal_sec_spi": 3, "signal_lipo_spii": 2, "signal_tat": 2,
        "signal_pilin_spiii": 2,
        "lpxtg_anchor": 2, "slh_or_anchor_domain": 2, "tm_beta_barrel": 3,
        "multi_tm_helix_penalty": -1, "cazy_hit": 1, "small_protein": 1,
        "foldseek_toxin_fold": 3, "context_mge_or_secretion": 2,
        "context_immunity_pair": 2, "context_pul": 1, "no_ko": 1,
    },

    "foldseek_self_cluster": True,
    # Which Foldseek target database wins when the same protein hits several.
    # hits.tsv records the target database per row, so the merged table can be
    # ranked by provenance rather than by raw bitscore alone: an unannotated
    # AFDB50 model routinely outscores a described Swiss-Prot or PDB hit, and
    # the description is what the report, toxin_fold and the uninformative
    # test all read. Names are matched against the basename of each target
    # path, most-preferred first. Empty (the default) means "the order you
    # listed your targets in", i.e. db.foldseek_target then
    # db.foldseek_extra_targets, which is the order a user already expresses a
    # preference in. Targets not named here rank last and fall back to
    # bitscore among themselves.
    "foldseek_target_priority": [],
    # NCBIfam/TIGRFAM families whose DESC line says "hypothetical protein" or
    # "DUF1234 domain-containing protein" name a family but no function, and
    # counting them as annotation promotes a protein out of the dark bin on
    # nothing. hmmsearch --tblout carries the description of the TARGET (the
    # protein), never of the query HMM, so the family descriptions are read
    # from db.ncbifam_hmm itself and cached beside the tblout. Set false to
    # restore the old behaviour, where every NCBIfam hit counted as
    # informative.
    "ncbifam_uninformative_test": True,
    # Matched as whole words against the Foldseek target description. Bare
    # "deaminase", "phospholipase", "patatin" and "hemolysin" hit cytidine
    # deaminases, patatin-like housekeeping lipases and hemolysin-III, each
    # scoring the largest single effector weight, so they are spelled out.
    # This list scored 0 of 38,204 proteins on a real gut metaproteome while
    # carrying the joint-largest weight in the table, so the additions below
    # are the families that were actually missing rather than a widening of
    # the ones already here. "holotoxin" is listed explicitly because the
    # whole-word rule means "Tc toxin" cannot match inside it, which is how
    # PDB 2vse - a genuine Tc-family holotoxin - was missed. The contact-
    # dependent and T6SS families matter most for gut commensals, which carry
    # CDI and LXG systems far more often than they carry classical exotoxins.
    # `Ntox\d*` and not `Ntox`: every family in that set is Ntox followed by a
    # number (Ntox15, Ntox28, Ntox47), and a digit is a word character, so the
    # trailing \b of the whole-word rule means a bare "Ntox" cannot match any
    # of them - the same failure that kept "Tc toxin" out of "holotoxin".
    # These are joined into one alternation, so a pattern is a REGEX, not a
    # literal; anything added here that contains a metacharacter must be
    # written as one.
    "toxin_fold_patterns": [
        "aerolysin", "MACPF", "cholesterol-dependent cytolysin",
        "perfringolysin", "ADP-ribosyltransferase", "ADP-ribosylating",
        "RTX", "alpha-hemolysin", "alpha-haemolysin", "hemolysin BL",
        "leukocidin", "Tc toxin", "holotoxin", "pore-forming", "colicin",
        "pyocin", "VgrG", "Rhs", "MARTX", "delta-endotoxin", "cytolysin",
        "insecticidal toxin", "nuclease toxin", "contact-dependent",
        "CdiA", "LXG", r"Ntox\d*", "zeta toxin", "pierisin",
    ],

    "anchor_pfams": ["PF00395", "PF01473", "PF00746", "PF13715",
                     "PF05738", "PF17998"],
}


# Blocks that are a list of the things the user actually has, so the config
# REPLACES the defaults there instead of merging into them. Merging made the
# five stock DIAMOND databases unremovable: listing only vfdb still left four
# /data/db/*.dmnd phantoms in doctor, in the diamond signature and in the
# stage's "database missing" warnings.
REPLACE_BLOCKS = {"db.diamond", "sources.diamond"}


def deep_merge(base, override, prefix=""):
    """Recursive merge. A shallow update would silently drop every default
    inside a nested block the user partially overrides (e.g. supplying only
    `thresholds.smorf_max_len` would delete all other thresholds)."""
    out = dict(base)
    for k, v in (override or {}).items():
        path = f"{prefix}{k}"
        if v is None and isinstance(out.get(k), dict):
            # `db:` with every child commented out parses as None. Letting that
            # replace the block only turns an ordinary config edit into a
            # NoneType traceback three functions later, so keep the defaults
            # and say out loud that the block did nothing.
            log(f"config: '{path}:' is empty - keeping the built-in defaults "
                f"for that block", "WARN")
            continue
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            if path in REPLACE_BLOCKS:
                dropped = sorted(x for x in out[k] if x not in v)
                if dropped:
                    log(f"config: {path} lists {sorted(v) or 'nothing'}, so "
                        f"the default entries {dropped} are NOT used", "WARN")
                out[k] = dict(v)
            else:
                out[k] = deep_merge(out[k], v, path + ".")
        else:
            out[k] = v
    return out


# Blocks whose keys the user defines: a database tag, a tool name, a predictor
# name. Everything else has a fixed key set, so an unrecognised key there is a
# typo rather than an addition. `analysis` is deliberately NOT here: its keys
# are exactly the Rmd's params, so `min_lfC` or `design_formla` is a typo that
# used to run the wrong statistics with nothing said.
FREEFORM = {"tool_args", "diamond_weights", "vfdb_category_weights",
            "diamond_evalues", "db.diamond", "sources.diamond"}


def unknown_keys(user, default, prefix=""):
    """Keys in the user's config that no default recognises.

    A mistyped key used to be silently ignored, so `run: {unipep: true}` left
    the stage disabled with nothing said, and the run simply did something
    other than what was asked.
    """
    out = []
    if not isinstance(user, dict) or not isinstance(default, dict):
        return out
    for k, v in user.items():
        path = f"{prefix}{k}"
        if k not in default:
            out.append(path)
            continue
        if path in FREEFORM or prefix.rstrip(".") in FREEFORM:
            continue
        out += unknown_keys(v, default[k], path + ".")
    return out


def nearest_config_key(bad):
    """A ready-to-append `did you mean ...` clause, or "" if nothing is close.

    Shared so that `doctor` and a real run say the same thing about the same
    typo: doctor used to print the key without the suggestion, which is the
    half of the message that actually fixes it.
    """
    import difflib

    def flat(d, pre=""):
        for k, v in d.items():
            yield pre + k
            if isinstance(v, dict) and (pre + k) not in FREEFORM:
                yield from flat(v, pre + k + ".")

    known = list(flat(DEFAULT_CONFIG))
    near = difflib.get_close_matches(bad, known, n=1, cutoff=0.7) or \
        difflib.get_close_matches(bad.split(".")[-1],
                                  [k.split(".")[-1] for k in known],
                                  n=1, cutoff=0.7)
    return f" — did you mean '{near[0]}'?" if near else ""


# Keys this tool used to accept. Without this, a config carrying one reads as a
# TYPO - unknown_keys reports it and nearest_config_key helpfully suggests the
# closest surviving name - and the user goes looking for their own mistake
# instead of learning the key was deliberately removed. Delete an entry once
# the release that removed it is old news; this is not a deprecation framework.
RETIRED_KEYS = {
    "run.effectors":
        "the effectors stage was removed: it ingested predictions from "
        "Bastion/EffectiveDB/T4SEpp, web services this tool has no way to "
        "invoke, and every one of them is now unreachable. Nothing is lost - "
        "the stage contributed nothing unless you configured it, and the "
        "score never depended on it. See CHANGELOG.",
    "effector_predictions":
        "removed with the effectors stage; there is nothing left to ingest. "
        "A pred_* column you join in yourself still reaches the shortlist as "
        "a column. See CHANGELOG.",
    "effector_prediction_weight":
        "removed with the effectors stage. See CHANGELOG.",
}


def report_unknown_keys(user, path):
    bad = unknown_keys(user, DEFAULT_CONFIG)
    if not bad:
        return []
    for b in bad:
        if b in RETIRED_KEYS:
            # Named as removed, and deliberately WITHOUT a "did you mean"
            # suggestion: there is no key to mean instead.
            log(f"{path}: '{b}' is no longer a setting - {RETIRED_KEYS[b]}",
                "WARN")
            continue
        # `or "."` so the sentence ends once, whether or not there is a
        # suggestion: "...key 'x' — did you mean 'y'? It is being ignored."
        log(f"{path}: unrecognised key '{b}'"
            + (nearest_config_key(b) or ".")
            + " It is being ignored, so this setting is NOT in effect.", "WARN")
    return bad


# Path-valued config keys, resolved relative to the config file rather than to
# the current working directory. A config that only works when you happen to
# be standing in the right directory is a trap, especially for an agent.
PATH_KEYS = ["proteins_faa", "quant_table", "manifest", "gff", "contigs_fna",
             "results_dir"]
DB_PATH_KEYS = ["eggnog_data", "pfam_hmm", "dbcan_hmm", "foldseek_target",
                "ncbi_taxonomy", "ncbifam_hmm", "kofam_profiles",
                "kofam_ko_list", "interproscan_sh", "hhblits_db",
                "jackhmmer_db"]


def resolve_paths(cfg, base):
    """Make every path in `cfg` absolute against `base`."""
    def R(v):
        if not isinstance(v, str) or not v:
            return v
        # `~` is not absolute, so without this `~/db/Pfam-A.hmm` became
        # <configdir>/~/db/Pfam-A.hmm and doctor reported a database the user
        # does have as missing. $VARS are expanded for the same reason.
        ex = os.path.expanduser(os.path.expandvars(v))
        if os.path.isabs(ex):
            # normalise only what expansion touched, so a path the user wrote
            # absolute still reaches the tools exactly as they wrote it.
            return os.path.normpath(ex) if ex != v else ex
        return os.path.normpath(os.path.join(base, ex))

    for k in PATH_KEYS:
        if cfg.get(k):
            cfg[k] = R(cfg[k])
    pre = cfg.get("emapper_precomputed")
    if isinstance(pre, str):
        cfg["emapper_precomputed"] = R(pre)
    elif isinstance(pre, list):
        cfg["emapper_precomputed"] = [R(x) for x in pre]
    db = cfg.get("db") or {}
    for k in DB_PATH_KEYS:
        if db.get(k):
            db[k] = R(db[k])
    if isinstance(db.get("diamond"), dict):
        # An entry set to null or "" means "I do not have this one". Dropping
        # it here keeps None out of os.path.exists() in the diamond stage and
        # out of doctor's MANUAL list.
        db["diamond"] = {k: R(v) for k, v in db["diamond"].items() if v}
    if isinstance(db.get("foldseek_extra_targets"), list):
        db["foldseek_extra_targets"] = [R(x) for x in db["foldseek_extra_targets"]]
    u = cfg.get("unipept") or {}
    if u.get("result"):
        u["result"] = R(u["result"])
    a = cfg.get("analysis") or {}
    for k in ("metadata", "msstats_comparison"):
        if a.get(k):
            a[k] = R(a[k])
    return cfg


def load_config(path):
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))  # deep copy
    if path:
        if yaml is None:
            sys.exit("pyyaml is needed to read a config file:  pip install pyyaml")
        if not os.path.exists(path):
            sys.exit(f"config not found: {path}")
        if os.path.isdir(path):
            sys.exit(f"{path} is a directory; --config takes a YAML file")
        # yaml.safe_load keeps the LAST of two identical keys, so a config
        # that grew a second `run:` block - the shape every doc snippet has -
        # silently threw the first one away and re-enabled the stages the
        # operator had turned off. Nobody means that, so it is fatal.
        class _NoDupLoader(yaml.SafeLoader):
            def construct_mapping(self, node, deep=False):
                seen = {}
                for kn, _ in node.value:
                    key = self.construct_object(kn, deep=deep)
                    try:
                        dup = key in seen
                    except TypeError:      # unhashable: SafeLoader's problem
                        continue
                    if dup:
                        sys.exit(
                            f"{path}: duplicate key '{key}' at line "
                            f"{kn.start_mark.line + 1}, already set at line "
                            f"{seen[key]}. YAML keeps only the last one, so "
                            "the earlier block would be discarded - merge the "
                            "two into a single block.")
                    seen[key] = kn.start_mark.line + 1
                return super().construct_mapping(node, deep=deep)

        try:
            with open(path, encoding="utf-8") as fh:
                user = yaml.load(fh, Loader=_NoDupLoader) or {}
        except yaml.YAMLError as e:
            mark = getattr(e, "problem_mark", None)
            where = (f" at line {mark.line + 1}, column {mark.column + 1}"
                     if mark else "")
            sys.exit(f"{path}: not valid YAML{where}\n  "
                     f"{getattr(e, 'problem', e)}\n  "
                     "A flow mapping written over two lines needs the "
                     "continuation indented inside the braces.")
        if not isinstance(user, dict):
            sys.exit(f"{path}: top level must be a mapping")
        report_unknown_keys(user, path)
        cfg = deep_merge(cfg, user)
        cfg = resolve_paths(cfg, os.path.dirname(os.path.abspath(path)))
    return cfg


# ======================================================================
# small utilities
# ======================================================================
_START = time.time()
_LOGFH = None
_LOGLOCK = threading.Lock()
_CTX = threading.local()


def set_log_context(name):
    """Tag this thread's log lines, so output from stages running concurrently
    stays attributable."""
    _CTX.stage = name


def configure_console_streams():
    """Make stdout/stderr survive a character the console cannot encode.

    Tool output is decoded with errors="replace" (see opener), which puts
    U+FFFD into descriptions. Python on Windows writes stdout in the console
    code page - cp1252 here - and printing U+FFFD to a cp1252 stream raises
    UnicodeEncodeError. That is not hypothetical: it happened while printing a
    parsed table. A stage that logs a protein description containing one
    replacement character could therefore abort a run that had been going for
    hours, on the LOG LINE rather than on the work.

    UTF-8 first, because that is what the log FILE already is and what a
    modern Windows console (PEP 528) and every POSIX terminal want; falling
    back to leaving the encoding alone and only relaxing the error handler,
    which still cannot raise. Both are best effort: a stream that pytest or a
    caller replaced may not be reconfigurable at all, which is why log()
    additionally writes defensively.
    """
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:
            continue
        for kwargs in ({"encoding": "utf-8", "errors": "replace"},
                       {"errors": "replace"}):
            try:
                reconfigure(**kwargs)
                break
            except (ValueError, OSError, AttributeError, TypeError):
                continue


def _write_safely(stream, text):
    """write(), but a character the stream cannot encode costs that character
    rather than the run. configure_console_streams normally makes this
    unreachable; it stays because the streams are not always ours to
    reconfigure and the logging path must never be what kills a stage."""
    try:
        stream.write(text)
    except UnicodeEncodeError:
        enc = getattr(stream, "encoding", None) or "ascii"
        stream.write(text.encode(enc, "replace").decode(enc, "replace"))


def log(msg, level="INFO"):
    tag = getattr(_CTX, "stage", "")
    prefix = f"[{time.time()-_START:7.1f}s] {level:5s} " + (f"{tag:>10s} | " if tag else "")
    with _LOGLOCK:
        for i, part in enumerate(str(msg).split("\n")):
            line = prefix + part if i == 0 else " " * len(prefix) + part
            _write_safely(sys.stderr, line + "\n")
            if _LOGFH:
                _write_safely(_LOGFH, line + "\n")
        sys.stderr.flush()
        if _LOGFH:
            _LOGFH.flush()


class StageError(RuntimeError):
    """Raised by die(). A real exception, not SystemExit: inside a worker
    thread SystemExit ends only that thread and its str() is just the exit
    code, so the diagnosis was lost before it reached the summary."""


def die(msg):
    raise StageError(msg)


def opener(path):
    """Text reader for a plain or gzipped file, tolerant of bad bytes.

    Every parser in this file reads through here, so a strict decode makes one
    stray byte anywhere in a tool's output fatal to the stage that reads it.
    That is not hypothetical: a 0xa0 (latin-1 non-breaking space) in a search
    result killed `integrate` on a real run after InterProScan had already
    spent three hours, with a message — "'utf-8' codec can't decode byte 0xa0"
    — that named neither the file nor the stage. VFDB subject titles, InterPro
    signature descriptions, HMM DESC lines and FASTA headers all carry latin-1
    in the wild, and none of them is ours to re-encode.

    errors="replace", the same as read_delim_table already uses, so a bad byte
    costs one character of a description instead of the stage. The blast
    radius is small and already covered: these bytes turn up in free text, and
    if one ever landed in an identifier the id-overlap diagnostics would say
    so — an identifier that is not valid UTF-8 is broken whatever we do.

    encoding is stated for the gzip branch too, where it was missing
    altogether. TextIOWrapper falls back to the LOCALE's codec, so that branch
    was never UTF-8 by construction — it was whatever the machine happened to
    be set to, which is UTF-8 on a modern desktop and ASCII under a bare C
    locale. emapper_precomputed is routinely a .gz, so this is the one input
    most likely to be read differently on the server than on the laptop.
    """
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_fasta(path):
    """Yield (id, seq). Header id = first whitespace-delimited token."""
    pid, buf = None, []
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(buf)
                tok = line[1:].split()
                if not tok:
                    raise StageError(
                        f"{path}: a record has an empty identifier (a bare "
                        "'>'). Every protein needs a unique id or the join to "
                        "the quant table is meaningless.")
                pid, buf = tok[0], []
            elif line:
                buf.append(line)
    if pid is not None:
        yield pid, "".join(buf)


def have(tool):
    return shutil.which(tool) is not None


def resolve_tool(name):
    """Turn a tool name into the absolute path PATH says it is.

    On POSIX this changes nothing: execvp would have found the same file.
    On Windows it is load-bearing. CreateProcess, which is what subprocess
    uses without a shell, does its own PATH search and only ever appends
    ".exe" - it ignores PATHEXT. So a ".cmd" or ".bat" earlier on PATH loses
    to an ".exe" later on it, and the tool that runs is not the tool
    shutil.which reported. Resolving here means one answer to "which binary is
    this", used by the availability check, by the logged command line and by
    the process that actually starts.

    An absolute path, or a name PATH cannot resolve, is handed back unchanged
    so the caller still fails with its own not-found message.
    """
    if not isinstance(name, str) or os.sep in name or (os.altsep and os.altsep in name):
        return name
    return shutil.which(name) or name


# Seconds between progress lines from a running tool; 0 disables. Set from
# config.progress_interval_s at the start of a run, because run_cmd is called
# from every stage and threading a cfg through all of them would touch code
# that has nothing to do with logging.
_PROGRESS_INTERVAL = float(DEFAULT_CONFIG["progress_interval_s"])
# Stderr lines kept while a command runs. Bounded on purpose: the failure tail
# has only ever quoted the last 15, and the reason stdout goes to devnull -
# not buffering tens of MB of chatter - applies here too.
_STDERR_KEEP = 200
# CSI escapes, which tqdm uses to colour the bar and to erase the line.
_ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")


def set_progress_interval(seconds):
    global _PROGRESS_INTERVAL
    _PROGRESS_INTERVAL = max(0.0, float(seconds or 0))


def _elapsed_str(seconds):
    s = int(seconds)
    return f"{s // 3600}h{s % 3600 // 60:02d}m" if s >= 3600 else \
        f"{s // 60}m{s % 60:02d}s"


def _progress_line(text, width=160):
    """One printable line out of whatever a tool last wrote.

    tqdm redraws in place with carriage returns and erases with CSI codes, so
    the raw text is a wall of partial redraws rather than a message. Universal
    newlines have already turned each \\r into its own line by the time we get
    here, so what is left is to drop the escape codes, the stray control
    characters and any excess width.
    """
    text = _ANSI_RE.sub("", text)
    # Two passes, not one: a tab is not printable, so filtering before
    # translating deleted it and glued the columns of a table together.
    text = "".join(" " if c in "\t\x08\x0b\x0c" else c for c in text)
    text = "".join(c for c in text if c == " " or c.isprintable()).strip()
    return text[:width - 1] + "…" if len(text) > width else text


def run_cmd(cmd, cwd=None, env=None):
    """Run a command, raising with the tail of stderr on failure.

    stdout goes to /dev/null: InterProScan and friends emit tens of MB of
    progress chatter, and buffering it in memory bought nothing.

    stderr is read as it arrives into a bounded ring instead of being
    collected in one go at the end, and every progress_interval_s seconds the
    newest line is echoed with the command and how long it has been running.
    Nothing is hoarded - the ring holds _STDERR_KEEP lines, far more than the
    15 the failure tail quotes - so the reasoning above still holds.

    The silence this fixes was real: tmbed ran 2 h 36 min and then died,
    twice, and InterProScan 2.9 h, with the log saying nothing between the
    command and the failure. Both write a tqdm bar to stderr the entire time;
    none of it reached the operator, so the only way to tell a live stage from
    a hung one was to watch its CPU ticks in /proc. The tutorial tells people
    to expect 1-3 DAY runs.

    Return value ("") and failure behaviour (RuntimeError quoting the tail)
    are deliberately unchanged: every stage depends on both.
    """
    log("$ " + " ".join(str(c) for c in cmd))
    name = os.path.basename(str(cmd[0]))
    started = time.time()
    kept = deque(maxlen=_STDERR_KEEP)
    # A deque raises if it is mutated while being iterated, and the reader
    # thread appends to this one continuously, so every read of it is a
    # snapshot taken under the lock.
    kept_lock = threading.Lock()
    stage = getattr(_CTX, "stage", "")

    def pump(stream):
        # The reader thread inherits nothing from threading.local, so the
        # stage tag is carried over by hand or these lines lose their owner
        # exactly when several stages are running at once.
        set_log_context(stage)
        try:
            for line in iter(stream.readline, ""):
                with kept_lock:
                    kept.append(line.rstrip("\n"))
        except (ValueError, OSError):
            pass          # pipe closed under us; the exit status still speaks
        finally:
            try:
                stream.close()
            except OSError:
                pass

    def newest_line():
        with kept_lock:
            snapshot = list(kept)
        for raw in reversed(snapshot):
            clean = _progress_line(raw)
            if clean:
                return clean
        return ""

    with open(os.devnull, "w") as null:
        # errors="replace" for the same reason opener() uses it: a tool that
        # writes latin-1 to stderr must not turn into a UnicodeDecodeError
        # that hides its actual failure.
        argv = [str(c) for c in cmd]
        argv[0] = resolve_tool(argv[0])
        proc = subprocess.Popen(argv, cwd=cwd, env=env,
                                stdout=null, stderr=subprocess.PIPE,
                                text=True, errors="replace")
        reader = threading.Thread(target=pump, args=(proc.stderr,),
                                  daemon=True)
        reader.start()
        try:
            while True:
                try:
                    proc.wait(timeout=_PROGRESS_INTERVAL or None)
                    break
                except subprocess.TimeoutExpired:
                    newest = newest_line()
                    log(f"{name} running {_elapsed_str(time.time() - started)}"
                        + (f" | {newest}" if newest else
                           " | no output yet on stderr"))
        except BaseException:
            # What subprocess.run did on Ctrl-C: do not leave a GPU job or an
            # InterProScan behind, still running, after the pipeline exits.
            proc.kill()
            proc.wait()
            raise
        reader.join(timeout=5)

    if proc.returncode != 0:
        with kept_lock:
            blob = "\n".join(kept)
        tail = "\n".join(blob.strip().splitlines()[-15:])
        raise RuntimeError(
            f"{cmd[0]} exited {proc.returncode}\n--- stderr tail ---\n{tail}")
    return ""


def nonempty(path):
    try:
        return os.path.getsize(path) > 0
    except OSError:
        return False


# Suffix for every in-progress output. Fixed, not random, so an operator
# looking at a results directory after a crash can tell a leftover from a
# result at a glance: `find results -name '.*.part.*'`.
ATOMIC_SUFFIX = ".part"


@contextlib.contextmanager
def atomic_out(path):
    """Yield a temp path to write, renamed onto `path` only on success.

    Every stage output used to be written in place, so a writer killed by the
    OOM reaper, a SIGKILL or a dropped ssh left a truncated file that the next
    run happily adopted as a finished stage — a silently short annotation
    table, hmmsearch tblout or quant matrix. The temp file lives in the SAME
    directory as its target because os.replace is only atomic within one
    filesystem; a temp under /tmp would degrade to a cross-device copy that
    can itself be interrupted half way.

    The extension is preserved (foo.tsv -> .foo.<pid>.<tid>.part.tsv) because
    several of the tools we hand this path to (InterProScan, foldseek) decide
    what to write from the name they are given. The leading dot matters just
    as much: results/diamond/*.tsv is globbed as "every DIAMOND database we
    searched", and glob skips dot-files, so a temp left behind by a killed
    writer cannot come back as a database tag called "vfdb.9134.7.part".

    Usage:
        with atomic_out(p.pfam) as tmp:
            run_cmd([... "--tblout", tmp ...])
    or, for our own writers:
        with atomic_out(p.context) as tmp:
            df.to_csv(tmp, sep="\t", index=False)
    """
    d = os.path.dirname(path) or "."
    os.makedirs(d, exist_ok=True)
    stem, ext = os.path.splitext(os.path.basename(path))
    tmp = os.path.join(
        d, f".{stem}.{os.getpid()}.{threading.get_ident()}{ATOMIC_SUFFIX}{ext}")
    _atomic_rm(tmp)
    try:
        yield tmp
    except BaseException:
        # Leave nothing a later run could mistake for a result. Best effort:
        # the original exception is what the operator needs to see.
        _atomic_rm(tmp)
        raise
    if not os.path.exists(tmp):
        # A tool that exited 0 without writing anything is a failure we can
        # name here, rather than an empty output adopted three stages later.
        die(f"nothing was written to {tmp}, so {path} was left alone. "
            "The tool exited successfully but produced no output — check the "
            "log above for what it was asked to do.")
    os.replace(tmp, path)


def _atomic_rm(path):
    try:
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        elif os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


def exists_all(paths):
    for path in paths:
        if not os.path.exists(path):
            return False
        if os.path.isdir(path) and not os.listdir(path):
            return False
    return True


class Paths:
    def __init__(self, cfg):
        R = cfg["results_dir"]
        self.R = R
        self.emapper = f"{R}/eggnog/emapper.emapper.annotations"
        self.emapper_report = f"{R}/eggnog/reuse_report.tsv"
        self.pfam = f"{R}/hmm/pfam.tblout"
        self.dbcan = f"{R}/hmm/dbcan.domtblout"
        self.diamond_dir = f"{R}/diamond"
        self.diamond_done = f"{R}/diamond/.done"
        self.feature_quant = f"{R}/quant/feature_quant.tsv"
        self.robject = f"{R}/metaannot.rds"
        self.signalp = f"{R}/topology/signalp/prediction_results.txt"
        self.tmbed = f"{R}/topology/tmbed.pred"
        self.cluster = f"{R}/cluster/fam_cluster.tsv"
        self.context = f"{R}/context.tsv"
        self.pass1 = f"{R}/annotation_pass1.tsv"
        # Two files on purpose. dark.faa is the ESMFold WORK-LIST: capped at
        # max_dark_structures and filtered to max_len_structure because GPU
        # time is finite. dark_all.faa is the DEFINITION of "unannotated":
        # every protein in the dark/DUF-only/profile-only bins, uncapped.
        # hhblits and jackhmmer must query the second, or a structure budget
        # silently truncates profile search and the "rescued by profile
        # search" fraction is computed on an effector_score-ranked subset.
        self.dark = f"{R}/dark.faa"
        self.dark_all = f"{R}/dark_all.faa"
        self.structures = f"{R}/structures"
        self.struct_done = f"{R}/structures/.done"
        self.foldseek = f"{R}/foldseek/hits.tsv"
        self.fold_clusters = f"{R}/foldseek/fold_clusters.tsv"
        self.final = f"{R}/annotation_final.tsv"
        self.summary = f"{R}/bin_summary.tsv"
        self.tier_coverage = f"{R}/tier_coverage.tsv"
        # Whether the evidence sources that reach the same protein agree
        # about it, which coverage numbers alone cannot say.
        self.agreement = f"{R}/source_agreement.tsv"
        self.quant_dir = f"{R}/quant"
        self.ncbifam = f"{R}/hmm/ncbifam.tblout"
        # Cache of NAME/ACC -> DESC read out of the NCBIfam HMM library. The
        # library is gigabytes; the descriptions are a few MB.
        self.ncbifam_desc = f"{R}/hmm/ncbifam_family_desc.tsv"
        self.kofam = f"{R}/kofam/kofam.tsv"
        self.interpro = f"{R}/interpro/interproscan.tsv"
        self.hhr_dir = f"{R}/hhblits"
        self.hhr_done = f"{R}/hhblits/.done"
        self.jackhmmer = f"{R}/hmm/jackhmmer.tblout"
        self.smorf_faa = f"{R}/smorf/smorf_proteins.faa"
        self.unipept_peptides = f"{R}/unipept/peptides.txt"
        self.unipept_cache = f"{R}/unipept/pept2lca_cache.tsv"
        self.unipept_lca = f"{R}/unipept/pept2lca.tsv"
        self.protein_taxonomy = f"{R}/unipept/protein_taxonomy.tsv"
        self.taxonomy_comparison = f"{R}/unipept/taxonomy_comparison.tsv"
        self.state = f"{R}/.metaannot_state.json"
        self.lock = f"{R}/.metaannot.lock"
        self.logfile = f"{R}/metaannot.log"

    def mkdirs(self):
        if os.path.exists(self.R) and not os.path.isdir(self.R):
            die(f"results_dir {self.R} exists and is not a directory. Choose "
                "another path, or move that file out of the way.")
        for d in [self.R, f"{self.R}/eggnog", f"{self.R}/hmm", self.diamond_dir,
                  f"{self.R}/topology/signalp", f"{self.R}/cluster",
                  f"{self.R}/foldseek", self.structures, self.quant_dir,
                  f"{self.R}/unipept", f"{self.R}/kofam",
                  f"{self.R}/interpro", f"{self.R}/hhblits", f"{self.R}/smorf"]:
            os.makedirs(d, exist_ok=True)


# ======================================================================
# parsers
# ======================================================================
GLOBAL_MAPS = {
    "map01100", "map01110", "map01120", "map01200", "map01210", "map01212",
    "map01230", "map01232", "map01250", "map01240", "map01220", "map01310",
    # map01320 (Sulfur cycle) is a current global map; map01130 (Biosynthesis
    # of antibiotics) was retired by KEGG but eggnog-mapper 2.1.x still emits
    # it from its eggNOG 5 snapshot, so it has to be listed to be excluded.
    "map01320", "map01130",
}
MAP_RE = re.compile(r"\bmap\d{5}\b")
KO_RE = re.compile(r"K\d{5}")
LPXTG_RE = re.compile(r"LP.TG")
DUF_RE = re.compile(r"^DUF\d|unknown function|uncharacteri", re.IGNORECASE)

CANONICAL_V2 = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]

ID_TRANSFORMS = {
    "exact": lambda s: s,
    "strip_after_first_pipe": lambda s: s.split("|")[0],
    "strip_before_last_pipe": lambda s: s.rsplit("|", 1)[-1],
    "strip_version_suffix": lambda s: re.sub(r"\.\d+$", "", s),
    "lowercase": lambda s: s.lower(),
    "basename_after_last_slash": lambda s: s.rsplit("/", 1)[-1],
}


# eggnog-mapper 2.0.x (the version behind the 2020-2021 catalogues) names its
# first column '#query_name' and spells three annotation columns differently
# from 2.1.x. Without the aliases set_index('protein_id') raises a bare
# KeyError and col() would silently return '' for Description/COG_category.
EMAPPER_V20_ALIASES = {
    "query_name": "protein_id",
    "seed_eggNOG_ortholog": "seed_ortholog",
    "seed_ortholog_evalue": "evalue",
    "seed_ortholog_score": "score",
    "best_tax_level": "max_annot_lvl",
    "eggNOG OGs": "eggNOG_OGs",
    "COG Functional cat.": "COG_category",
    "eggNOG free text desc.": "Description",
}


def emapper_header(line):
    """Canonicalise an emapper header line onto the 2.1.x column names.

    The first column is renamed unconditionally: it is the query id whatever
    the version calls it, and guessing wrong there costs the whole join.
    """
    header = line.lstrip("#").rstrip("\n").split("\t")
    if not header:
        return header
    header = [EMAPPER_V20_ALIASES.get(h.strip(), h.strip()) for h in header]
    if header[0] in ("query", "protein_id"):
        header[0] = "protein_id"
    else:
        log(f"emapper header starts with '{header[0]}', not 'query'; treating "
            "it as the protein id column", "WARN")
        header[0] = "protein_id"
    return header


def parse_emapper(path):
    header, rows = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = emapper_header(line)
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                raise ValueError(
                    f"{path}: no '#query' header line; run the emapper stage, "
                    "which reconstructs it for headerless files")
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                f += [""] * (len(header) - len(f))
            rows.append(f[:len(header)])
    df = pd.DataFrame(rows, columns=header)
    # 'query' is already renamed by emapper_header; the rename stays for
    # tables reconstructed by an older run of this script.
    df = df.rename(columns={"query": "protein_id"}).set_index("protein_id")
    return df.replace("-", "")


def from_dict(d, index, fill=""):
    """Aligned column from a sparse dict. Building the Series once and
    reindexing is vectorised; a list comprehension over the index does a
    scalar pandas lookup per protein, which dominated the profile."""
    if not d:
        return pd.Series([fill] * len(index), index=index)
    return pd.Series(d).reindex(index).fillna(fill)


def col(df, name, index):
    """Safe column access: always a str Series aligned to `index`.

    DataFrame.get() returns the bare default when a column is missing, which
    then breaks every downstream .fillna()/.str call. eggnog-mapper column
    sets differ between versions, so this guard is load-bearing.
    """
    if df is None or name not in df.columns:
        return pd.Series("", index=index, dtype=object)
    return df[name].reindex(index).fillna("").astype(str)


def warn_if_no_records(path, n, kind):
    """A tolerant parser turns a corrupted file into zero records, not an error.

    Zero records from a non-empty file means every protein silently loses this
    evidence and the bins shift wholesale — the same damage as a truncated
    file, without the truncation to notice.
    """
    if n == 0 and nonempty(path):
        log(f"{path}: parsed 0 {kind} from a non-empty file. The format is "
            "probably not what was expected; every protein will look "
            "unannotated by this evidence.", "WARN")
    return n


# hmmsearch --tblout: 18 fixed, whitespace-delimited columns and then a
# free-text "description of target" that runs to end of line. Verified against
# HMMER 3.4: the description is the TARGET's FASTA description (for hmmsearch
# the target is the protein and the query is the HMM, so this is NOT the HMM's
# DESC line), and it is a bare "-" when the target header carries none.
HMM_TBLOUT_FIXED_COLS = 18


def parse_hmm_tblout(path):
    """(query name, query accession, full E-value, full score, target desc).

    The description used to be dropped. It is the only free text the tblout
    carries, and for jackhmmer (where the target is a UniRef entry) it is what
    tells "Uncharacterized protein" apart from a named homolog, so it is kept
    rather than reconstructed later from a file that is gigabytes wide.
    """
    out = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 6:
                continue
            # Only split off the description once the fixed columns are all
            # present; a truncated line must not have its 7th field read as
            # free text.
            desc = ""
            if len(f) > HMM_TBLOUT_FIXED_COLS:
                desc = " ".join(f[HMM_TBLOUT_FIXED_COLS:])
            if desc == "-":
                desc = ""         # hmmsearch's "no description" placeholder
            try:
                rec = (f[2], f[3], float(f[4]), float(f[5]), desc)
            except ValueError:
                continue          # build first: out[f[0]] would create an entry
            out[f[0]].append(rec)
    warn_if_no_records(path, len(out), "HMM hits")
    return out


def parse_hmm_lib_desc(path, limit_names=None):
    """NAME/ACC -> DESC for every model in an HMMER library.

    hmmsearch --tblout cannot report the query HMM's description (see above),
    so the only place a family's DESC line exists is the library itself. Keyed
    on both NAME and ACC, with and without the version suffix, because the
    tblout's query name and query accession columns may carry either.
    """
    out = {}
    name = acc = desc = None

    def flush():
        if desc:
            for k in (name, acc):
                if not k or k == "-":
                    continue
                if limit_names is None or k in limit_names:
                    out[k] = desc
                base = k.rsplit(".", 1)[0]
                if base != k and (limit_names is None or base in limit_names):
                    out.setdefault(base, desc)

    with opener(path) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line[5:].strip()
            elif line.startswith("ACC "):
                acc = line[4:].strip()
            elif line.startswith("DESC "):
                desc = line[5:].strip()
            elif line.startswith("//"):
                flush()
                name = acc = desc = None
    flush()                       # a library truncated before its last "//"
    return out


def parse_hmm_domtblout(path, min_cov, max_evalue):
    out = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 22:
                continue
            try:
                qlen = float(f[5])
                i_evalue = float(f[12])
                hmm_from, hmm_to = float(f[15]), float(f[16])
            except ValueError:
                continue
            cov = (hmm_to - hmm_from + 1) / qlen if qlen else 0.0
            if i_evalue <= max_evalue and cov >= min_cov:
                # dbCAN-HMMdb NAME fields end in '.hmm' (run_dbcan strips them
                # too); keeping the suffix makes dbcan_hits string-incomparable
                # to the eggNOG CAZy column and to any family grouping.
                name = re.sub(r"\.hmm$", "", f[3])
                out[f[0]].append((name, i_evalue, round(cov, 3)))
    warn_if_no_records(path, len(out), "domain hits")
    return out


def parse_diamond(path, min_evalue, min_qcov, min_pident):
    best = {}
    if not os.path.exists(path):
        return best
    cols = ["qseqid", "sseqid", "pident", "length", "evalue",
            "bitscore", "qcovhsp", "scovhsp", "stitle"]
    n_rows, n_shape = 0, 0
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(cols):
                continue
            n_rows += 1
            r = dict(zip(cols, f))
            # stitle is free text. A number there means the file is not in the
            # 9-field custom format stage_diamond asks for — the DIAMOND
            # default --outfmt 6 has 12 numeric columns — and every field would
            # be read from the wrong position, so refuse the row rather than
            # filtering it on nonsense.
            try:
                float(r["stitle"])
            except ValueError:
                pass
            else:
                n_shape += 1
                continue
            try:
                if (float(r["evalue"]) > min_evalue
                        or float(r["qcovhsp"]) < min_qcov
                        or float(r["pident"]) < min_pident):
                    continue
                bs = float(r["bitscore"])
            except ValueError:
                continue
            q = r["qseqid"]
            if q not in best or bs > best[q][0]:
                best[q] = (bs, r["sseqid"], float(r["pident"]),
                           float(r["evalue"]), r["stitle"])
    if n_shape:
        log(f"{path}: {n_shape}/{n_rows} rows are not in the column format "
            f"this pipeline writes ('--outfmt 6 {' '.join(cols)}') and were "
            "skipped. A file copied in from a hand-run DIAMOND with the "
            "default 12 columns would be read from the wrong positions.",
            "WARN")
    warn_if_no_records(path, len(best), "alignments")
    return {q: v[1:] for q, v in best.items()}


def parse_signalp6(path):
    out, header = {}, None
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if "Prediction" in line:
                    header = line.lstrip("#").strip().split("\t")
                continue
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            cs = ""
            if header and "CS Position" in header:
                i = header.index("CS Position")
                cs = f[i] if i < len(f) else ""
            elif f[-1].startswith("CS pos"):
                cs = f[-1]
            out[f[0].split()[0]] = (f[1].strip(), cs)
    warn_if_no_records(path, len(out), "predictions")
    return out


def _count_segments(labels, chars):
    """TMbed encodes orientation in the CASE of the label (H = in->out,
    h = out->in), so a hairpin's two helices can abut with no loop residue
    between them. Counting runs of the H/h class merges them into one segment;
    a change of the exact character starts a new one."""
    n, prev = 0, ""
    for c in labels:
        if c in chars and c != prev:
            n += 1
        prev = c
    return n


def parse_tmbed(path):
    """tmbed --out-format 0 (3-line): header, sequence, per-residue labels.
    H/h transmembrane helix, B/b transmembrane beta strand, S signal."""
    # Streamed through iter_tmbed_records, which is also what stage_tmbed
    # counts chunks with: one definition of a complete record, so a file
    # truncated mid-write cannot be read as short here and full there.
    out = {}
    for hdr, _seq, lab in iter_tmbed_records(path):
        out[hdr[1:].split()[0]] = (_count_segments(lab, "Hh"),
                                   _count_segments(lab, "Bb"))
    warn_if_no_records(path, len(out), "topologies")
    return out


# Foldseek --format-output field lists, told apart by field count.
# FOLDSEEK_COLS is what stage_foldseek should ask for; FOLDSEEK_COLS_LEGACY is
# what older result files and hand-run searches carry. The distinction matters
# scientifically: alntmscore is normalised by the ALIGNMENT, so a 40-residue
# local match inside a 300-residue query can score 0.6 while the two proteins
# share almost no fold. qtmscore is normalised by the query and is the number
# the "TM >= 0.5" promise in the report actually needs; qlen gives the
# alignment-coverage backstop that a local 3Di+AA search otherwise lacks.
FOLDSEEK_COLS = ["query", "target", "fident", "alnlen", "qlen", "tlen",
                 "evalue", "bits", "prob", "alntmscore", "qtmscore",
                 "ttmscore", "lddt", "theader"]
FOLDSEEK_COLS_LEGACY = ["query", "target", "fident", "alnlen", "evalue",
                        "bits", "prob", "alntmscore", "lddt", "theader"]
# Appended by stage_foldseek when it concatenates the per-target result
# files. Not a Foldseek --format-output field: Foldseek never reports which
# target database a row came from, and the accession cannot stand in for it.
FOLDSEEK_DB_COL = "target_db"


def parse_foldseek(path, min_evalue, min_prob, min_tm, min_qcov=0.0,
                   target_priority=None):
    """Best structural hit per query, with the target database taken seriously.

    stage_foldseek concatenates one result file per target database, so a
    protein can hold hits from AFDB50, PDB and Swiss-Prot at once. Ranking
    those by bitscore alone lets an unannotated AFDB50 model beat a described
    Swiss-Prot or PDB hit, and the description is what toxin_fold, the
    uninformative test and the report all read. The accession cannot be used
    to tell the databases apart (AFDB50 and AFDB-SwissProt both look like
    AF-<acc>-F1-model_v4), so hits.tsv now carries the target database per row
    and `target_priority` (most-preferred first, matched against those names)
    decides between databases before bitscore does.

    Files written before that column existed still parse: every row then
    carries an empty database name, no priority can apply, and selection is by
    bitscore exactly as it was.
    """
    best = {}
    if not os.path.exists(path):
        return best
    n_legacy = 0
    hdr_cols = None
    prio = [str(x).strip().lower() for x in (target_priority or [])
            if str(x).strip()]

    def rank(db):
        db = (db or "").strip().lower()
        return prio.index(db) if db in prio else len(prio)

    with opener(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                # stage_foldseek names its columns, so target_db is found by
                # name rather than by counting fields: with the column
                # appended, a 10-field legacy row and an 11-field new one are
                # not distinguishable by width.
                if hdr_cols is None:
                    hdr_cols = line.lstrip("#").strip().split("\t")
                continue
            f = line.split("\t")
            if hdr_cols and len(f) >= len(hdr_cols):
                cols = hdr_cols
            elif len(f) >= len(FOLDSEEK_COLS):
                cols = FOLDSEEK_COLS
            elif len(f) >= len(FOLDSEEK_COLS_LEGACY):
                cols = FOLDSEEK_COLS_LEGACY
            else:
                continue
            if "qtmscore" not in cols:
                n_legacy += 1
            r = dict(zip(cols, f))
            try:
                ev, prob = float(r["evalue"]), float(r["prob"])
                b = float(r["bits"])
                tm = float(r.get("qtmscore") or r["alntmscore"])
                qlen = float(r.get("qlen") or 0)
                qcov = float(r["alnlen"]) / qlen if qlen else 1.0
            except ValueError:
                continue
            if ev > min_evalue or prob < min_prob or tm < min_tm:
                continue
            if qcov < min_qcov:
                continue
            q = re.sub(r"\.(pdb|cif)(_[A-Za-z0-9]+)?$", "", r["query"])
            db = (r.get(FOLDSEEK_DB_COL) or "").strip()
            # Best hit per (query, database) first, so the choice between
            # databases is made between each one's own best hit.
            k = (q, db)
            if k not in best or b > best[k][0]:
                best[k] = (b, r["target"], prob, tm, r["theader"], db)
    if n_legacy:
        log(f"{path}: {n_legacy} rows are in the 10-column Foldseek format, so "
            "the TM gate falls back to alntmscore (normalised by the "
            "alignment, not the query) and no coverage filter can be applied; "
            "a short local match can pass it. Re-run the foldseek stage to get "
            "qtmscore and qlen.", "WARN")

    per_query = defaultdict(list)
    for (q, _db), v in best.items():
        per_query[q].append(v)
    out, n_prio, dbs = {}, 0, set()
    for q, cands in per_query.items():
        dbs.update(c[5] for c in cands)
        by_bits = max(cands, key=lambda c: c[0])
        pick = min(cands, key=lambda c: (rank(c[5]), -c[0]))
        if pick[1] != by_bits[1]:
            n_prio += 1
        out[q] = pick[1:]
    named = {d for d in dbs if d}
    if len(named) > 1 and prio:
        log(f"foldseek: {len(named)} target databases in {path}; preference "
            f"order {', '.join(prio)}. {n_prio} protein(s) kept a hit from a "
            "preferred database over a higher-scoring one elsewhere.")
        unranked = sorted(d for d in named if d.lower() not in prio)
        if unranked:
            log("foldseek target database(s) not named in "
                "foldseek_target_priority, so ranked last: "
                + ", ".join(unranked), "WARN")
    elif out and not named and prio:
        log(f"{path} predates the target_db column, so which database each hit "
            "came from is unknown and foldseek_target_priority cannot be "
            "applied; the best bitscore wins. Re-run the foldseek stage to "
            "record provenance.", "WARN")
    warn_if_no_records(path, len(out), "structure hits")
    return out


def parse_cluster(path):
    out = {}
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[1]] = f[0]
    warn_if_no_records(path, len(out), "cluster members")
    return out


def parse_context(path):
    if not os.path.exists(path) or not nonempty(path):
        return {}
    df = pd.read_csv(path, sep="\t", encoding="utf-8", encoding_errors="replace")
    if "protein_id" not in df.columns:
        return {}
    df = df.set_index("protein_id")
    df = df[~df.index.duplicated(keep="first")]
    return df.to_dict("index")


# ======================================================================
# additional annotation evidence
# ======================================================================
def parse_kofam(path):
    """KOfamScan detail-tsv. A leading '*' marks a hit above the family's own
    adaptive threshold; anything else is below it and is not a KO call."""
    out = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                f = re.split(r"\s{2,}", line.rstrip("\n"))
            if len(f) < 6 or f[0].strip() != "*":
                continue
            gene, ko = f[1].strip(), f[2].strip()
            try:
                score, ev = float(f[4]), float(f[5])
            except ValueError:
                continue
            # KOfamScan's detail-tsv template wraps the definition in double
            # quotes; left in, they are quoted a second time by to_csv.
            desc = f[6].strip().strip('"') if len(f) > 6 else ""
            out[gene].append((ko, score, ev, desc))
    warn_if_no_records(path, len(out), "KO assignments")
    return out


def parse_interproscan(path):
    """InterProScan TSV. Columns 12/13 (InterPro accession and description)
    are only present when -iprlookup was used."""
    out = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rec = {"analysis": f[3], "signature": f[4], "desc": f[5],
                   "ipr": f[11] if len(f) > 11 else "",
                   "ipr_desc": f[12] if len(f) > 12 else "",
                   "go": f[13] if len(f) > 13 else ""}
            out[f[0]].append(rec)
    warn_if_no_records(path, len(out), "InterPro matches")
    return out


# The hit table is printed as '%3d %-30.30s %5.1f ...'. When the truncated
# name fills all 30 columns and Prob is 100.0 there is exactly one space
# between them, so a mandatory description group cannot match and the BEST hit
# is silently skipped in favour of the next one. The group is optional here.
HHR_HIT = re.compile(
    r"^\s*\d+\s+(\S+)(?:\s+(.*?))?\s+(\d+\.\d)\s+(\S+)\s+(\S+)\s+"
    r"(-?[\d.]+)\s+(-?[\d.]+)\s+(\d+)\s")


def parse_hhr_dir(d, min_prob):
    """HH-suite .hhr files, one per query. Probability is the metric that
    matters; E-value alone under-reports remote profile-profile hits."""
    out = {}
    for path in sorted(glob.glob(os.path.join(d, "*.hhr"))):
        q = os.path.splitext(os.path.basename(path))[0]
        best = None
        with opener(path) as fh:
            in_table = False
            for line in fh:
                if line.startswith(" No Hit"):
                    in_table = True
                    continue
                if in_table:
                    if not line.strip():
                        break
                    m = HHR_HIT.match(line)
                    if not m:
                        continue
                    hit = m.group(1)
                    desc = (m.group(2) or "").strip()
                    prob = float(m.group(3))
                    if prob >= min_prob and (best is None or prob > best[1]):
                        best = (hit, prob, desc, m.group(4))
        if best:
            out[q] = best
    return out


# A search database is often built by renaming the fasta it came from
# ("MGYG000001_00023" -> "uhgpSM_MGYG000001_00023", "HUMANHOST_..."), while
# the eggNOG table still carries the original, unprefixed id. No entry in
# ID_TRANSFORMS can add or remove a prefix, so those rows never match and the
# coverage gate ends up blaming biology for what is a naming mismatch. Hence
# both a configurable prefix strip and, in the diagnostic, prefix detection.
PREFIX_DELIM_RE = re.compile(r"[_|:.]")


def id_prefix_candidates(pid, limit=3):
    """(prefix, tail) pairs for the first few delimiters in pid.

    Only used to name a likely prefix in the diagnostic: if the eggNOG table
    carries pid's tail, the search database is what added the prefix.
    """
    out = []
    for m in PREFIX_DELIM_RE.finditer(pid):
        tail = pid[m.end():]
        if tail:
            out.append((pid[:m.end()], tail))
        if len(out) >= limit:
            break
    return out


# A merged search database is the normal case for a metaproteomics run, and
# its tiers do not annotate alike. The Pittsburgh one is four: uhgpL_ (395,467
# UHGP proteins), OIDECCNN_ (43,139 Prokka calls off the matched metagenome),
# uhgpSM_ (14,797) and ampS_ (2,168). A single "94% have an eggNOG hit" hides
# that one of those tiers arrives with precomputed annotations while another
# is ORFs nobody has ever seen, and the difference is the whole reason the
# proteome was merged in the first place.
#
# 12 is where a prefix stops being a source label and starts being part of
# the id: one Prokka locus tag per MAG would give hundreds, and a table with
# a row each answers no question anyone asked.
MAX_ID_TIERS = 12


def id_key_shape(pid):
    """The identifier under a tier tag, with its digit runs masked.

    A merged database has two things worth telling apart and they are easy to
    confuse. The TIER TAG says which source a protein came from and is what
    id_tiers reports. The KEY says what the identifier actually is, and it
    lives BEHIND the tag:

        uhgpL_MGYG000004906_01237   ->  tier uhgpL_    key MGYG#_#
        uhgpSM_MGYG000009567_01280  ->  tier uhgpSM_   key MGYG#_#
        OIDECCNN_00158              ->  tier OIDECCNN_ key #
        ampS_AMP10.000_478          ->  tier ampS_     key AMP#.#_#

    Two tiers sharing a key are the same identifier namespace under two
    labels, which is exactly the case emapper_strip_id_prefix exists for: one
    eggNOG row annotates a protein under every tag it appears with, and a tag
    left out of that list loses its whole tier.

    Digit RUNS are masked, not digits, because widths vary within one
    namespace: the Prokka tier of the database this was written for runs
    OIDECCNN_00001 to OIDECCNN_1712297, and 5-, 6- and 7-digit accessions are
    all one key space. Masking per digit would split it into three.
    """
    m = PREFIX_DELIM_RE.search(pid)
    tail = pid[m.end():] if m else pid
    return re.sub(r"\d+", "#", tail) or "(empty)"


def id_tiers(ids, max_tiers=MAX_ID_TIERS):
    """{prefix: count} when the ids look like a merged database, else {}.

    The prefix is everything up to and including the first _ | : or . - the
    same delimiter set id_prefix_candidates uses. Returns {} when there is
    only one prefix (nothing to split) or more than max_tiers (the prefix is
    part of the id, not a label). Ids with no delimiter at all group under "".
    """
    counts = {}
    for pid in ids:
        m = PREFIX_DELIM_RE.search(pid)
        key = pid[:m.end()] if m else ""
        counts[key] = counts.get(key, 0) + 1
        if len(counts) > max_tiers:
            return {}
    return {} if len(counts) < 2 else counts


def prepare_emapper(sources, faa, out, report, transform, min_cov, warn_cov,
                    diagnose_rows, strip_prefixes=()):
    want = {pid for pid, _ in read_fasta(faa)}
    if not want:
        die(f"no sequences in {faa}")
    log(f"emapper reuse: {len(want)} proteins to annotate")

    want_lc = {i.lower() for i in want}
    tf = ID_TRANSFORMS[transform]
    lc = transform == "lowercase"
    alt = {n: set() for n in ID_TRANSFORMS if n != transform}

    # bridge: eggNOG id (the fasta id with the prefix removed) -> [fasta ids].
    # Everything downstream joins on the fasta id, so a row matched this way
    # has to be written under the fasta id, not the id the table carried.
    # A list, not a single id: a composite search database holds the same bare
    # id under two tier prefixes (uhgpSM_MGYG..._00100 and uhgpL_MGYG..._00100
    # are one UHGP protein in two tiers), and the eggNOG table carries one row
    # for it. Keeping only the first fasta id sent the other one to the dark
    # bin even though its annotation was sitting right there.
    bridge = defaultdict(list)
    for pref in strip_prefixes:
        n_pref = 0
        for pid in want:
            if not pid.startswith(pref) or len(pid) == len(pref):
                continue
            key = pid[len(pref):]
            key = key.lower() if lc else key
            n_pref += 1
            if pid not in bridge[key]:
                bridge[key].append(pid)
        if not n_pref:
            log(f"emapper_strip_id_prefix '{pref}' matches no fasta id; "
                "check the prefix", "WARN")
        else:
            log(f"emapper reuse: prefix '{pref}' stripped from {n_pref} "
                "fasta ids when matching the eggNOG table")
    n_shared_key = sum(1 for v in bridge.values() if len(v) > 1)
    if n_shared_key:
        log(f"emapper reuse: {n_shared_key} bare id(s) occur under more than "
            "one emapper_strip_id_prefix; one eggNOG row is copied out once "
            "per fasta id it annotates")

    # Candidate tails of every fasta id, so an unmatched row can name the
    # prefix that would have matched it. Built only when the diagnostic can
    # actually run, because it costs a few entries per protein.
    probe_map = {}
    if diagnose_rows > 0:
        for pid in want:
            for pref, tail in id_prefix_candidates(pid):
                probe_map.setdefault(tail.lower() if lc else tail, (pref, pid))
    pref_hits = {}

    header, n_scanned, n_probed, n_shared_rows = None, 0, 0, 0
    seen, bridged = set(), set()

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    # Written to a temp path and renamed only once the coverage gate passes,
    # so a killed or rejected run cannot leave a truncated table behind that
    # the next run adopts as a finished stage.
    tmp = out + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh_out:
        fh_out.write("## reconstructed by metaannot\n")
        for path in sources:
            if not os.path.exists(path):
                die(f"emapper_precomputed path not found: {path}")
            with opener(path) as fh:
                for line in fh:
                    if line.startswith("#query"):
                        if header is None:
                            header = line.lstrip("#").rstrip("\n").split("\t")
                            fh_out.write("#" + "\t".join(header) + "\n")
                        continue
                    if line.startswith("#") or not line.strip():
                        continue
                    n_scanned += 1
                    f = line.rstrip("\n").split("\t")
                    if header is None:
                        if len(f) != len(CANONICAL_V2):
                            die(f"{path}: no '#query' header and {len(f)} columns, "
                                f"expected {len(CANONICAL_V2)}")
                        header = list(CANONICAL_V2)
                        fh_out.write("#" + "\t".join(header) + "\n")
                        log("no header found; assumed canonical eggnog-mapper v2 "
                            "column order", "WARN")
                    raw = f[0]
                    q = tf(raw)
                    # Every fasta id this row annotates, not just the first:
                    # the direct match (if any) plus every prefixed id whose
                    # bare form is this row's query.
                    targets = [q] if ((q.lower() in want_lc) if lc
                                      else (q in want)) else []
                    for pid in bridge.get(q.lower() if lc else q, ()):
                        if pid not in targets:
                            targets.append(pid)
                            bridged.add(pid)
                    if targets:
                        for pid in targets:
                            if pid in seen:
                                continue
                            seen.add(pid)
                            # Write the id the fasta uses, or the downstream
                            # join to the fasta silently finds nothing.
                            row = list(f)
                            row[0] = pid
                            fh_out.write("\t".join(row) + "\n")
                        if len(targets) > 1:
                            n_shared_rows += 1
                    elif n_probed < diagnose_rows:
                        n_probed += 1
                        for name in alt:
                            v = ID_TRANSFORMS[name](raw)
                            if (v in want_lc) if name == "lowercase" else (v in want):
                                alt[name].add(v)
                        got = probe_map.get(q.lower() if lc else q)
                        if got:
                            pref_hits.setdefault(got[0], set()).add(got[1])

    cov = len(seen) / len(want)
    log(f"emapper reuse: scanned {n_scanned} rows, matched {len(seen)} "
        f"({100*cov:.1f}% of the protein set)")
    # Per tier, because one headline coverage over a merged database is an
    # average of things that are not alike: a public catalogue tier arrives
    # with precomputed annotations and a tier assembled from this study's own
    # reads does not, and 94% overall can be 99% and 30%.
    tiers = id_tiers(want)
    tier_cov = {}
    if tiers:
        hit = {k: 0 for k in tiers}
        for i in seen:
            m = PREFIX_DELIM_RE.search(i)
            k = i[:m.end()] if m else ""
            if k in hit:
                hit[k] += 1
        shape_of = {}
        for pref, n in sorted(tiers.items(), key=lambda kv: -kv[1]):
            tier_cov[pref] = (n, hit[pref])
            ex = next((i for i in want if i.startswith(pref)), pref)
            shape_of.setdefault(id_key_shape(ex), []).append(pref)
            log(f"emapper reuse:   {pref or '(no prefix)':20s} "
                f"{hit[pref]:>8,}/{n:<8,} {100.0 * hit[pref] / n:5.1f}%")
            if not hit[pref]:
                # Zero is not a small number, it is a different kind of
                # answer. A tier at 5% has a table that mostly misses; a tier
                # at 0% has no table keyed on its identifiers at all, and
                # every one of its proteins will be reported 4_dark by
                # construction rather than by biology. On the run this was
                # written for that is the AMPSphere tier -- 2,168 identified
                # proteins, 30% of the whole dark fraction, and no eggNOG
                # table anywhere on the machine is keyed on AMP/SPHERE ids.
                # The headline coverage was 98.4%, so nothing else said it.
                log(f"emapper reuse: tier {pref or '(no prefix)'} matched "
                    f"NONE of its {n:,} protein(s). Not a low number -- a "
                    "zero, which is what a tier whose identifiers no "
                    "configured table is keyed on looks like. Every one of "
                    "them will bin as 4_dark for want of a join, not for "
                    "want of biology. Either add a table for it to "
                    "emapper_precomputed, or record that this tier is "
                    "unannotated by construction so the dark fraction is "
                    "read with that in mind", "WARN")
        # A tier tag left out of emapper_strip_id_prefix while a tier sharing
        # its key space is in it loses the whole table for that tier, and the
        # symptom is a plausible-looking coverage number rather than an error.
        strip = tuple(strip_prefixes or ())
        for shape, prefs in shape_of.items():
            if len(prefs) < 2:
                continue
            miss = [p for p in prefs if p not in strip]
            if strip and miss and len(miss) < len(prefs):
                log(f"emapper reuse: {', '.join(prefs)} share the identifier "
                    f"key {shape}, but emapper_strip_id_prefix lists only "
                    f"{[p for p in prefs if p in strip]}. "
                    f"{', '.join(miss)} will not be bridged to the same "
                    "table rows, so that tier reports as unannotated when it "
                    "is only unjoined. Add it.", "WARN")
    if bridged:
        log(f"emapper reuse: {len(bridged)} of those proteins matched only "
            "after emapper_strip_id_prefix was removed from the fasta id")
    if n_shared_rows:
        log(f"emapper reuse: {n_shared_rows} eggNOG row(s) annotated more "
            "than one fasta id (the same bare id under several prefixes)")

    recommend, best_pref = "", ""
    if cov < warn_cov:
        log(f"coverage {100*cov:.1f}% is below emapper_warn_coverage "
            f"({100*warn_cov:.0f}%); proteins recoverable under other id "
            f"transforms (probed {n_probed} unmatched rows):", "WARN")
        log(f"    {transform:28s} {100*cov:5.1f}%  (used)", "WARN")
        ranked = sorted(((len(v) / len(want), k) for k, v in alt.items()),
                        reverse=True)
        for frac, name in ranked:
            log(f"    {name:28s} {100*frac:5.1f}%", "WARN")
        pranked = sorted(((len(v) / len(want), k)
                          for k, v in pref_hits.items()), reverse=True)
        for frac, pref in pranked[:5]:
            log(f"    fasta-id prefix {pref!r:19s} {100*frac:5.1f}%", "WARN")
        if ranked and ranked[0][0] > cov + 0.05:
            recommend = ranked[0][1]
            log(f"id mismatch detected — set emapper_id_transform: {recommend}",
                "WARN")
        if pranked and pranked[0][0] > 0.05:
            best_pref = pranked[0][1]
            log("the search database appears to have prefixed the fasta ids — "
                f"set emapper_strip_id_prefix: {best_pref}", "WARN")
        if not recommend and not best_pref:
            # Only claim this once both an id transform and a prefix have been
            # ruled out, and say what was actually looked at: probing stops
            # after emapper_diagnose_rows unmatched rows.
            capped = "" if n_probed < diagnose_rows else \
                ", the emapper_diagnose_rows cap was hit so later rows went " \
                "unprobed"
            log("neither an id transform nor a common fasta-id prefix explains "
                f"the shortfall ({n_probed} unmatched rows probed{capped}); "
                "the unmatched proteins look genuinely unannotated", "WARN")

    if report:
        with open(report, "w", encoding="utf-8") as fh:
            fh.write("metric\tvalue\n")
            fh.write(f"proteins_in_fasta\t{len(want)}\n")
            fh.write(f"emapper_rows_scanned\t{n_scanned}\n")
            fh.write(f"proteins_annotated\t{len(seen)}\n")
            fh.write(f"proteins_bridged_by_prefix\t{len(bridged)}\n")
            fh.write(f"rows_shared_by_several_fasta_ids\t{n_shared_rows}\n")
            fh.write(f"coverage\t{cov:.4f}\n")
            fh.write(f"id_transform\t{transform}\n")
            fh.write(f"recommended_transform\t{recommend}\n")
            fh.write(f"recommended_strip_id_prefix\t{best_pref}\n")
            for pref, (n, hit) in tier_cov.items():
                fh.write(f"tier_{pref or 'none'}_proteins\t{n}\n")
                fh.write(f"tier_{pref or 'none'}_annotated\t{hit}\n")

    if cov < min_cov:
        die(f"only {100*cov:.1f}% of proteins matched (emapper_min_coverage "
            f"is {100*min_cov:.0f}%). Fix the id mismatch rather than lowering "
            "that threshold — every unmatched protein would be misreported as "
            "4_dark. Set emapper_id_transform or emapper_strip_id_prefix; the "
            "WARN lines above rank the candidates. The partial table is left "
            f"at {tmp} for inspection.")
    os.replace(tmp, out)


def stage_emapper(cfg, p):
    pre = cfg.get("emapper_precomputed") or ""
    pre = [pre] if isinstance(pre, str) and pre else list(pre or [])
    tf = cfg["emapper_id_transform"]
    if tf not in ID_TRANSFORMS:
        die(f"unknown emapper_id_transform '{tf}'; "
            f"choose from {sorted(ID_TRANSFORMS)}")
    pref = cfg.get("emapper_strip_id_prefix") or ""
    pref = [pref] if isinstance(pref, str) and pref else list(pref or [])
    pref = [str(x) for x in pref if str(x)]
    if pre:
        prepare_emapper(pre, cfg["proteins_faa"], p.emapper, p.emapper_report,
                        tf, cfg["emapper_min_coverage"],
                        cfg["emapper_warn_coverage"],
                        cfg["emapper_diagnose_rows"], pref)
    else:
        if pref:
            log("emapper_strip_id_prefix only applies when reusing a "
                "precomputed table; ignored for a live emapper.py run", "WARN")
        if not have("emapper.py"):
            die("emapper.py not found and emapper_precomputed is empty. "
                "Either install eggnog-mapper or point emapper_precomputed "
                "at an existing .emapper.annotations file.")
        ram = cfg.get("ram_gb") or 0
        cmd = ["emapper.py", "-i", cfg["proteins_faa"], "--itype", "proteins",
               # -o is a prefix and --output_dir is the directory; putting a
               # path in -o and "." in --output_dir depended on emapper
               # joining them naively.
               "-o", "emapper", "--output_dir", f"{p.R}/eggnog",
               "--temp_dir", f"{p.R}/eggnog",
               "--data_dir", cfg["db"]["eggnog_data"], "-m", "diamond",
               "--sensmode", "more-sensitive",
               "--cpu", cfg["threads"], "--override"]
        need = int(cfg.get("emapper_dbmem_min_gb", 64))
        if ram and ram >= need:
            cmd.append("--dbmem")
            log(f"--dbmem enabled ({ram} GB budget >= {need} GB)")
        elif ram:
            log(f"--dbmem left off: {ram} GB budget is below "
                f"emapper_dbmem_min_gb ({need} GB), and loading the annotation "
                "database into less memory than it needs only causes swapping")
        run_cmd(cmd + tool_args(cfg, "emapper"))


# ======================================================================
# stages: homology and topology
# ======================================================================
def stage_pfam(cfg, p):
    if not have("hmmsearch"):
        die("hmmsearch not found (conda install -c bioconda hmmer)")
    # hmmsearch has no memory flag; its footprint tracks --cpu and the model
    # set, so the RAM budget is honoured by the CPU split rather than a flag.
    with atomic_out(p.pfam) as tmp:
        run_cmd(["hmmsearch", "--cut_ga", "--noali", "--cpu", cfg["threads"],
                 "--tblout", tmp, "-o", os.devnull]
                + tool_args(cfg, "hmmsearch")
                + [cfg["db"]["pfam_hmm"], cfg["proteins_faa"]])


def stage_dbcan(cfg, p):
    if not have("hmmsearch"):
        die("hmmsearch not found (conda install -c bioconda hmmer)")
    with atomic_out(p.dbcan) as tmp:
        run_cmd(["hmmsearch", "--domE", cfg["thresholds"]["dbcan_evalue"],
                 "--noali", "--cpu", cfg["threads"],
                 "--domtblout", tmp, "-o", os.devnull]
                + tool_args(cfg, "hmmsearch")
                + [cfg["db"]["dbcan_hmm"], cfg["proteins_faa"]])


# DIAMOND prints its own scoring constants on every makedb run: "Scoring
# parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)". They are
# what turns a raw alignment score into a bit score, and therefore what decides
# whether a hit of a given length can reach a given e-value at all.
_DMND_LAMBDA, _DMND_K = 0.267, 0.041
# Mean of the BLOSUM62 diagonal (116/20). A perfect self-match of an
# average-composition peptide scores this per residue and nothing scores more
# over a whole sequence, so the estimate below is optimistic on purpose: it is
# the BEST an alignment of that length could ever do, not a typical one.
_BLOSUM62_MEAN_DIAGONAL = 5.8
# Query length for the estimate. Long on purpose - the e-value scales with it,
# so a generous query makes the check err towards saying nothing.
_EVALUE_QUERY_LEN = 1000
# `diamond makedb` on a single 11-residue sequence writes 143 bytes, so a file
# below this is smaller than DIAMOND's own header and cannot be a database at
# all. Observed: a failed makedb left a ZERO-BYTE .dmnd on disk, and searching
# it reports no hits - indistinguishable in the output from a real absence of
# virulence factors.
_DMND_MIN_BYTES = 128


def diamond_evalue_for(cfg, tag):
    """The e-value this DIAMOND database is searched at.

    Per-database because one threshold cannot fit both a virulence-factor
    database of 350-residue proteins and a bacteriocin database whose median
    sequence is 15 residues; see diamond_evalues in the config.
    """
    over = (cfg.get("diamond_evalues") or {})
    raw = over[tag] if tag in over else cfg["thresholds"]["diamond_evalue"]
    try:
        return float(raw)
    except (TypeError, ValueError):
        where = (f"diamond_evalues.{tag}" if tag in over
                 else "thresholds.diamond_evalue")
        die(f"{where} must be a number, not {raw!r}")


def diamond_min_pident_for(cfg, tag):
    """The identity floor this DIAMOND database is filtered at.

    Per-database for the same reason the e-value is: one number cannot fit a
    reference database whose hits are read as "this protein IS that one"
    (CARD, VFDB, BAGEL, floored at 50) and one read as "this protein is in
    that family" (where 30 is the useful setting).
    """
    over = (cfg.get("diamond_min_pidents") or {})
    raw = over[tag] if tag in over else cfg["thresholds"]["diamond_min_pident"]
    try:
        return float(raw)
    except (TypeError, ValueError):
        where = (f"diamond_min_pidents.{tag}" if tag in over
                 else "thresholds.diamond_min_pident")
        die(f"{where} must be a number, not {raw!r}")


def best_possible_evalue(letters, typical_len):
    """The smallest e-value a hit against this database could ever reach.

    E = m*n*2**-S', with S' the bit score. Feeding it a perfect self-match of
    the database's typical sequence length answers a question the pipeline was
    never asking before the run: is this search capable of a hit at all?
    """
    if not letters or not typical_len or typical_len <= 0:
        return None
    bits = (_DMND_LAMBDA * _BLOSUM62_MEAN_DIAGONAL * float(typical_len)
            - math.log(_DMND_K)) / math.log(2.0)
    try:
        return _EVALUE_QUERY_LEN * float(letters) * 2.0 ** -bits
    except OverflowError:
        return 0.0            # long sequences: any e-value is reachable


def _fasta_lengths(path, cap=200000):
    """Sequence lengths in a FASTA, for the median. Capped: the point is the
    shape of the database, and reading a 100 GB UniRef to learn it is not."""
    out, cur = [], 0
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        out.append(cur)
                    if len(out) >= cap:
                        return out
                    cur = 0
                else:
                    cur += len(line.strip())
    except OSError:
        return out
    if cur:
        out.append(cur)
    return out


# BAGEL4 ships two different things and they look alike on disk: the
# bacteriocin SEQUENCE files, and the motif seed set its HMM/regex step is
# built from - headers like LE-nisin, MA-lantibiotic, ggmotif, lasso, with a
# median length of 15 residues. The seed set is what got built into a DIAMOND
# database on a real run, and it returned exactly 0 hits against 38,204
# proteins. The e-value advice below was the wrong answer to that: blastp
# against 15-residue seeds is a search for those fifteen residues, not for the
# molecules they mark, and no threshold makes it into a bacteriocin search.
#
# 25 residues, not 40: mature nisin is 34, so a real bacteriocin database can
# be genuinely short and must not be accused of being a seed set.
MOTIF_SEED_MAX_LEN = 25
MOTIF_SEED_MARKERS = ("ggmotif", "lasso", "motif", "seed")


def diamond_source_fasta(cfg, tag, path):
    """A FASTA on disk this .dmnd was built from, or "".

    .fas is the name `doctor --fix` stages beside the database; a
    sources.diamond entry that is a local path rather than a URL is the other
    way a config records where the sequences came from.
    """
    stem = os.path.splitext(path)[0]
    named = ((cfg.get("sources") or {}).get("diamond") or {}).get(tag, "")
    cands = ([named] if named and os.path.isfile(named) else []) + \
        [stem + ext for ext in (".fas", ".faa", ".fasta")]
    for cand in cands:
        if os.path.isfile(cand):
            return cand
    return ""


def motif_seed_evidence(cfg, tag, path, typical):
    """Why this database looks like a motif seed set rather than sequences.

    Returns a list of human-readable reasons, empty when it looks like a
    normal protein database. Two signals, either sufficient: a typical
    sequence too short to be a protein at all, and headers carrying the words
    a seed set uses. The headers are only readable when a source FASTA is
    beside the .dmnd, since `diamond dbinfo` reports counts and nothing else.
    """
    why = []
    if typical is not None and typical < MOTIF_SEED_MAX_LEN:
        why.append(f"its typical sequence is {typical:.0f} residues, shorter "
                   "than any whole protein")
    fasta = diamond_source_fasta(cfg, tag, path)
    if fasta:
        found, prefixed, n = set(), False, 0
        with opener(fasta) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                n += 1
                tok = (line[1:].split() or [""])[0].lower()
                # LE- and MA- are checked as a PREFIX of the accession, not as
                # a substring: "gamma-haemolysin" contains "ma-".
                prefixed = prefixed or tok.startswith(("le-", "ma-"))
                low = line.lower()
                found |= {m for m in MOTIF_SEED_MARKERS if m in low}
                if n >= 200:
                    break
        if prefixed:
            found.add("LE-/MA- accession prefixes")
        if found:
            why.append(f"{os.path.basename(fasta)} carries "
                       + ", ".join(sorted(found)) + " in its headers")
    return why


def diamond_db_profile(cfg, tag, path):
    """(typical_len, letters, provenance) for a DIAMOND database.

    `diamond dbinfo` is the only thing that can read a .dmnd, and it reports
    Sequences and Letters, so the typical length it yields is the mean. When
    diamond is not installed or the file predates dbinfo, fall back to a source
    FASTA the config names or that `doctor --fix` staged beside the database,
    where the median is available. When neither is, return (None, None, why)
    and say so rather than guessing a length the check would then act on.
    """
    why = ""
    if have("diamond"):
        r = None
        try:
            r = subprocess.run([resolve_tool("diamond"), "dbinfo", "-d", path],
                               capture_output=True, text=True, timeout=120,
                               encoding="utf-8", errors="replace")
        except (OSError, subprocess.SubprocessError) as e:
            why = f"diamond dbinfo failed to run ({e})"
        if r is not None:
            got = {}
            for line in (r.stdout or "").splitlines():
                k, _, v = line.strip().partition("  ")
                v = v.strip()
                if k in ("Sequences", "Letters") and v.isdigit():
                    got[k] = int(v)
            if "Sequences" in got and "Letters" in got:
                n, letters = got["Sequences"], got["Letters"]
                return ((letters / n if n else 0), letters,
                        f"mean of {n} sequences / {letters} letters, from "
                        "`diamond dbinfo`")
            why = "`diamond dbinfo` reported no Sequences/Letters for this file"
    else:
        why = "diamond is not installed, so the .dmnd cannot be read"

    cand = diamond_source_fasta(cfg, tag, path)
    if cand:
        lens = sorted(_fasta_lengths(cand))
        if lens:
            return (lens[len(lens) // 2], sum(lens),
                    f"median of {len(lens)} sequences in {cand}")
    return (None, None,
            why + ", and no source FASTA is on disk beside it or named in "
            f"sources.diamond.{tag}")


def diamond_db_check(cfg, tag, path):
    """(refusal, warning) for one configured DIAMOND database; either may be
    None. A missing file is the caller's business, not this function's.

    All three of these were real. A `diamond makedb` that had failed left a
    zero-byte .dmnd, which the stage would have searched, reporting no hits -
    the same output as a real absence of virulence factors. A database whose
    sequences are too short for the configured e-value cannot produce a hit
    however good the alignment. And the BAGEL database that prompted both
    checks turned out to be neither: it was built from BAGEL4's motif SEED
    set, 262 entries of median length 15, so its 0 hits against 38,204
    proteins were not a threshold problem at all and lowering --evalue would
    have produced meaningless hits instead of meaningless silence.
    """
    try:
        size = os.path.getsize(path)
    except OSError:
        return None, None
    stem = os.path.splitext(path)[0]
    if size < _DMND_MIN_BYTES:
        return (f"{tag}: {path} is {size} bytes, smaller than a DIAMOND "
                "header - a failed `diamond makedb` leaves a file like this, "
                "and searching it reports 0 hits, which reads in the output "
                "exactly like a real absence. Rebuild it: diamond makedb "
                f"--in <fasta> -d {stem}"), None

    typical, letters, prov = diamond_db_profile(cfg, tag, path)
    if letters is not None and (letters == 0 or not typical):
        return (f"{tag}: {path} holds no sequences ({prov}). Searching it "
                "would report 0 hits, which is indistinguishable from a real "
                f"absence. Rebuild it: diamond makedb --in <fasta> -d {stem}"), None
    ev = diamond_evalue_for(cfg, tag)
    weight = (cfg.get("diamond_weights") or {}).get(tag)
    # A database that cannot answer is also holding a scoring weight that says
    # it can contribute to the export ranking.
    note = (f" It also carries diamond_weights {tag}: {weight}, which claims "
            "it can contribute to the score.") if weight else ""
    seed = motif_seed_evidence(cfg, tag, path, typical)
    if seed:
        # Deliberately INSTEAD of the e-value advice below, not alongside it.
        # Telling someone to lower --evalue here sends them to tune a
        # threshold on a database that is the wrong kind of thing, and the
        # tuned search still answers a question nobody asked.
        return None, (
            f"{tag}: this looks like a motif or seed set rather than a "
            f"protein sequence database - {'; and '.join(seed)}. A blastp "
            "against it searches for those residues, not for the molecules "
            "they mark, so its hits and its 0 hits both say nothing about the "
            "biology, and NO e-value makes that a real search. BAGEL in "
            "particular ships both: build the database from its bacteriocin "
            "sequence files, not from the seed set its HMM step uses. If this "
            f"really is a database of very short peptides, set "
            f"sources.diamond.{tag} to the FASTA so this check can see what "
            f"it is.{note}")
    if typical is None:
        return None, (f"{tag}: {prov}, so nothing here can tell whether "
                      f"--evalue {ev:g} is reachable for this database. If it "
                      "returns 0 hits, that may be the threshold rather than "
                      "the biology.")
    best = best_possible_evalue(letters, typical)
    if best is None or best <= ev:
        return None, None
    return None, (
        f"{tag}: sequences here are about {typical:.0f} residues ({prov}), and "
        f"the best e-value even a perfect alignment that long could reach is "
        f"~{best:.0e} - above the configured --evalue {ev:g}. This search is "
        "incapable of a hit before it starts, so 0 hits will say nothing about "
        f"the biology. Set a per-database e-value (diamond_evalues: "
        f"{{{tag}: 1e-3}}), or record that this database needs different "
        f"settings than the rest.{note}")


def stage_diamond(cfg, p):
    if not have("diamond"):
        die("diamond not found (conda install -c bioconda diamond)")
    dbs = cfg["db"].get("diamond") or {}
    if not dbs:
        log("no diamond databases configured, nothing to do", "WARN")
        open(p.diamond_done, "w", encoding="utf-8").close()
        return
    jobs = [(tag, path) for tag, path in dbs.items() if os.path.exists(path)]
    for tag, path in dbs.items():
        if not os.path.exists(path):
            log(f"diamond database missing, skipping: {tag} -> {path}", "WARN")
    if not jobs:
        open(p.diamond_done, "w", encoding="utf-8").close()
        return

    # A database that cannot answer is worse than one that is missing: the
    # missing one is reported above, the broken one returns 0 hits, and 0 hits
    # is what a real absence of virulence factors looks like too. Both halves
    # of this were observed on one run - a zero-byte .dmnd left by a failed
    # makedb, and a BAGEL database built from a motif seed set rather than
    # from bacteriocin sequences.
    refusals = []
    for tag, path in jobs:
        bad, warn = diamond_db_check(cfg, tag, path)
        if bad:
            refusals.append(bad)
        elif warn:
            log(warn, "WARN")
    if refusals:
        die("refusing to search a DIAMOND database that cannot answer:\n  "
            + "\n  ".join(refusals))

    # Several small databases in parallel beat one after another: DIAMOND's
    # thread scaling is sublinear, so 4 jobs at N/4 threads finish sooner than
    # 4 jobs at N threads in sequence.
    base_ev = float(cfg["thresholds"]["diamond_evalue"])
    base_id = float(cfg["thresholds"]["diamond_min_pident"])
    for tag, _path in jobs:
        ev = diamond_evalue_for(cfg, tag)
        if ev != base_ev:
            log(f"{tag}: searching at --evalue {ev:g} from diamond_evalues, "
                f"not the thresholds.diamond_evalue {base_ev:g} the other "
                "databases use")
        pid = diamond_min_pident_for(cfg, tag)
        if pid != base_id:
            log(f"{tag}: searching at --id {pid:g} from diamond_min_pidents, "
                f"not the thresholds.diamond_min_pident {base_id:g} the other "
                "databases use, so nothing below that identity is reported "
                "at all")

    workers = min(len(jobs), max(1, int(cfg.get("diamond_workers", 4))))
    per = max(1, int(cfg["threads"]) // workers)
    log(f"{len(jobs)} database(s), {workers} at a time x {per} threads")

    # DIAMOND's peak memory is roughly 6 GB per unit of block size (-b), so
    # the budget is divided by the concurrent jobs and converted here.
    ram = cfg.get("ram_gb") or 0
    mem = ["-b", f"{max(0.4, min(12.0, (ram / workers) / 6.0)):.2f}",
           "-c", "1"] if ram else []
    if mem:
        log(f"block size {mem[1]} per job from a {ram} GB budget")

    def one(job):
        tag, dbpath = job
        # The .done marker is atomic by construction, but the per-database
        # tables it stands for were not: a killed DIAMOND left a half-written
        # <tag>.tsv that integrate read as that database's complete answer.
        with atomic_out(f"{p.diamond_dir}/{tag}.tsv") as tmp:
            # --id as well as the reader's filter. DIAMOND applies it during
            # the search, so a floored database writes thousands of rows
            # instead of hundreds of thousands and the parse is not the place
            # the promise first takes effect.
            run_cmd(["diamond", "blastp", "-q", cfg["proteins_faa"],
                     "-d", dbpath, "-o", tmp, "--very-sensitive",
                     "-e", f"{diamond_evalue_for(cfg, tag):g}",
                     "--id", f"{diamond_min_pident_for(cfg, tag):g}",
                     "--max-target-seqs", 5, "--threads", per, "--quiet"]
                    + mem + tool_args(cfg, "diamond")
                    + ["--outfmt", "6", "qseqid", "sseqid", "pident",
                       "length", "evalue", "bitscore", "qcovhsp", "scovhsp",
                       "stitle"])

    errs = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(one, j): j[0] for j in jobs}
        for f in concurrent.futures.as_completed(futs):
            try:
                f.result()
            except BaseException as e:                 # noqa: BLE001
                errs.append(f"{futs[f]}: {e}")
    if errs:
        die("diamond failed for " + "; ".join(errs))
    open(p.diamond_done, "w", encoding="utf-8").close()


def stage_signalp(cfg, p):
    if not have("signalp6"):
        die("signalp6 not found. SignalP 6.0 needs an academic licence; "
            "install it, or set run.topology: false (you will lose the single "
            "most informative filter for effector candidates).")
    # --write_procs parallelises output writing, not inference, so
    # --torch_num_threads is what actually uses the CPU allocation.
    dest = os.path.dirname(p.signalp)
    with atomic_out(p.signalp) as tmp:
        # signalp6 writes a directory, not a file, so it is pointed at a
        # scratch directory beside the real one. Everything it produced is
        # moved across afterwards and prediction_results.txt — the file the
        # resume test looks at — lands last, by rename.
        work = f"{tmp}.d"
        _atomic_rm(work)
        os.makedirs(work, exist_ok=True)
        try:
            cmd = ["signalp6", "--fastafile", cfg["proteins_faa"],
                   "--organism", "other", "--output_dir", work,
                   "--format", "none",
                   "--mode", cfg.get("signalp_mode", "fast"),
                   "--write_procs", cfg["threads"],
                   "--torch_num_threads", cfg["threads"]]
            bsize = int(cfg.get("signalp_batch_size", 0) or 0)
            if bsize:
                cmd += ["--bsize", bsize]
            run_cmd(cmd + tool_args(cfg, "signalp6"))
            got = os.path.join(work, "prediction_results.txt")
            if not os.path.exists(got):
                die(f"signalp6 exited 0 but wrote no prediction_results.txt "
                    f"in {work}; nothing to adopt as {p.signalp}")
            # The side files (gff3, processed fasta) are not read here but
            # people do look at them, so they are kept rather than binned.
            for name in sorted(os.listdir(work)):
                if name == "prediction_results.txt":
                    continue
                _atomic_rm(os.path.join(dest, name))
                shutil.move(os.path.join(work, name),
                            os.path.join(dest, name))
            os.replace(got, tmp)
        finally:
            _atomic_rm(work)


def cuda_probe():
    """-> (usable, detail): what the GPU stages will actually find.

    Two questions that are easy to conflate. `nvidia-smi` says a card is
    physically present and its driver is loaded. `torch.cuda.is_available()` is
    what esmfold and tmbed actually test, and a CPU-only torch wheel answers
    "no" on a machine with a perfectly good card - which is the single most
    confusing way for a GPU stage to fail, because the hardware is right there.
    Report both, and say plainly when they disagree.

    Deliberately does not import torch when it is absent: `doctor` must stay
    fast and must run on a machine that has no torch at all.
    """
    smi = shutil.which("nvidia-smi") is not None
    card = ""
    if smi:
        try:
            r = subprocess.run([resolve_tool("nvidia-smi"), "-L"],
                               capture_output=True, text=True, timeout=20)
            card = (r.stdout or "").strip().splitlines()[0] if r.returncode == 0                 and r.stdout.strip() else ""
        except (OSError, subprocess.SubprocessError):
            card = ""
    if "torch" not in sys.modules and not importlib.util.find_spec("torch"):
        return (False, "torch is not installed, so no stage can use a GPU"
                + (f" (a card is present: {card})" if card else ""))
    try:
        import torch
        if torch.cuda.is_available():
            try:
                name = torch.cuda.get_device_name(0)
            except Exception:                               # noqa: BLE001
                name = card or "device 0"
            return True, f"torch reports CUDA available: {name}"
        if card:
            return (False, f"a card is present ({card}) but torch reports CUDA "
                    "unavailable - usually a CPU-only torch build; reinstall "
                    "torch with the CUDA index for your driver")
        return False, "no CUDA device and no nvidia-smi"
    except Exception as e:                                  # noqa: BLE001
        return False, f"torch could not be queried ({type(e).__name__}: {e})"


# TMbed writes nothing until it finishes. A single invocation over a whole
# proteome is therefore an all-or-nothing bet measured in days: the run on
# 455,571 proteins (214.9M residues) was still going after two days with an
# empty output file, and the two before it died at 2 h 36 min with nothing
# recoverable. Splitting the
# input into chunks turns that into a series of checkpoints. The price is one
# ProtT5 load per chunk, which is why the chunk count is bounded from both
# ends: tmbed_chunk_residues sets the floor, TMBED_MAX_PARTS the ceiling.
TMBED_MAX_PARTS = 256


def tmbed_chunk_plan(lengths, budget, max_parts=TMBED_MAX_PARTS):
    """Group (id, length) pairs into length-sorted, residue-budgeted chunks.

    Longest first, for two reasons. ProtT5 pads every sequence in a batch out
    to the longest one in it, so a chunk of similar lengths wastes less work
    than one mixing 30-residue peptides with 3000-residue proteins. And
    whatever is going to exhaust the device is in the first chunk, where it
    costs one chunk to discover instead of the whole stage.

    Returns (chunks, budget). The budget returned can be larger than the one
    asked for: a sequence longer than the budget still needs a chunk, and a
    budget small enough to ask for more than max_parts chunks would spend more
    time loading weights than predicting.
    """
    if not lengths:
        return [], 0
    total = sum(n for _, n in lengths)
    budget = int(budget or 0)
    if budget <= 0:                      # 0 = one invocation, as before
        return [[pid for pid, _ in lengths]], total
    # Greedy packing closes a chunk when the NEXT sequence would overflow it,
    # so every chunk but the last can fall short of the budget by as much as
    # the longest sequence. A floor of total/max_parts alone therefore still
    # overshoots the ceiling: 5000 x 100 residues against 256 parts planned
    # 264 of them. Adding the longest sequence to the floor guarantees every
    # non-final chunk holds more than total/max_parts, and so bounds the
    # count. It is also what gives an over-long sequence a chunk of its own
    # rather than dropping it.
    budget = max(budget, -(-total // max_parts) + max(n for _, n in lengths))
    chunks, cur, cur_n = [], [], 0
    for pid, n in sorted(lengths, key=lambda x: (-x[1], x[0])):
        if cur and cur_n + n > budget:
            chunks.append(cur)
            cur, cur_n = [], 0
        cur.append(pid)
        cur_n += n
    if cur:
        chunks.append(cur)
    return chunks, budget


def iter_tmbed_records(path):
    """(header, sequence, labels) for every COMPLETE record of a 3-line file.

    One definition of "complete", used by both the parser and the chunk
    bookkeeping. The format has no trailer, so a TMbed killed mid-write leaves
    a header and a sequence with no label line; that is not a prediction and
    must never be counted, copied or adopted as one.
    """
    if not os.path.exists(path):
        return
    buf = []
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                buf = [line]
            elif len(buf) in (1, 2):
                buf.append(line)
                if len(buf) == 3:
                    yield tuple(buf)
                    buf = []


def stage_tmbed(cfg, p):
    if not have("tmbed"):
        die("tmbed not found (pip install tmbed && tmbed download)")
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(cfg["gpu_device"]))
    # --use-gpu on its own is fatal on a CPU-only host: TMbed only tolerates a
    # missing/failing GPU when --cpu-fallback is given, so "auto" asks for the
    # GPU and lets it fall back rather than losing the topology evidence.
    want = str(cfg.get("tmbed_use_gpu", "auto")).strip().lower()
    if want in ("auto", "true", "false"):
        gpu = {"auto": ["--use-gpu", "--cpu-fallback"],
               "true": ["--use-gpu", "--no-cpu-fallback"],
               "false": ["--no-use-gpu"]}[want]
    else:
        die(f"unknown tmbed_use_gpu '{want}'; choose auto, true or false")
    if want == "auto":
        log("tmbed: GPU preferred, CPU fallback allowed "
            "(set tmbed_use_gpu: true to make a missing GPU fatal)")
    if want != "false":
        # "Fell back to CPU" is not a detail at this scale. TMbed embeds with
        # ProtT5, so a CPU fallback on a large proteome is hours of work per
        # chunk. Say the number out loud before it starts, not after.
        usable, why = cuda_probe()
        if not usable:
            n = sum(1 for _ in read_fasta(cfg["proteins_faa"]))
            log(f"tmbed: {why}. This stage will run on CPU over {n:,} "
                "protein(s), where it is one to two orders of magnitude "
                "slower than on a GPU. Set tmbed_use_gpu: false to accept "
                "that deliberately, true to make it fatal, or "
                "run.topology: false to skip both topology stages", "WARN")

    faa = cfg["proteins_faa"]
    out_dir = os.path.dirname(p.tmbed) or "."
    os.makedirs(out_dir, exist_ok=True)

    # One over-long protein kills the whole stage, and no device setting saves
    # it. ProtT5's attention score matrix is length-squared x heads: titin, at
    # 34,350 residues, asks for 141 GB in a single allocation. Observed twice
    # on real data - 8.79 GiB refused on a 16 GB card, then 151 GB refused on
    # a 94 GB host under --cpu-fallback. Cap the input instead, and record
    # what was cut rather than letting the exclusion pass unnoticed.
    #
    # Lengths only: a 455,571-protein FASTA does not fit in a list of
    # sequences, so each chunk re-reads the file rather than holding it.
    cap = int(cfg.get("tmbed_max_len", 0) or 0)
    work, excluded = [], []
    for pid, seq in read_fasta(faa):
        (excluded if cap and len(seq) > cap else work).append((pid, len(seq)))
    if excluded:
        excl = f"{out_dir}/tmbed_excluded.tsv"
        with open(excl, "w", encoding="utf-8") as fh:
            fh.write("protein_id\tlength\n")
            for pid, n in sorted(excluded, key=lambda x: -x[1]):
                fh.write(f"{pid}\t{n}\n")
        log(f"tmbed: {len(excluded)} protein(s) longer than "
            f"tmbed_max_len={cap} are excluded; ProtT5 attention is "
            f"length-squared and the longest here "
            f"({max(n for _, n in excluded)} aa) would need more memory than "
            f"any device present. They get no topology evidence and are "
            f"listed in {excl}", "WARN")

    if not work:
        # Nothing to predict is not a failure: a proteome of nothing but
        # over-cap sequences, or an empty FASTA, should leave an empty
        # prediction file rather than handing TMbed an empty input and
        # reporting whatever it does with it. The stage is empty_ok, so the
        # file is adoptable on a rerun.
        with atomic_out(p.tmbed) as tmp:
            open(tmp, "w", encoding="utf-8").close()
        log(f"tmbed: no sequence is left to predict "
            f"({len(excluded)} excluded by tmbed_max_len={cap}), so "
            f"{p.tmbed} is empty and nothing is run", "WARN")
        return

    chunks, budget = tmbed_chunk_plan(work, cfg.get("tmbed_chunk_residues"))
    total_res = sum(n for _, n in work)
    asked = int(cfg.get("tmbed_chunk_residues") or 0)
    log(f"tmbed: {len(work):,} protein(s), {total_res / 1e6:.1f}M residues -> "
        f"{len(chunks)} chunk(s) of at most {budget / 1e6:.1f}M residues, "
        f"longest sequence {max(n for _, n in work):,} aa. TMbed writes "
        "nothing until it finishes, so each chunk is a checkpoint: an "
        "interrupted run resumes from the last one instead of starting over."
        + (f" The budget asked for ({asked / 1e6:.1f}M) was raised to keep "
           f"the plan under {TMBED_MAX_PARTS} chunks, because ProtT5 is "
           "loaded once per chunk." if asked and budget > asked else "")
        + (" Set tmbed_chunk_residues: 0 for a single invocation."
           if len(chunks) > 1 else ""))

    # Every work protein, mapped to its chunk. Entries are removed as
    # predictions are written, so whatever is left at the end is exactly the
    # set that got none - taken from the files, not from the loop's own
    # bookkeeping, because a chunk that exits 0 can still be short.
    where = {}
    for i, ids in enumerate(chunks):
        for pid in ids:
            where[pid] = i
    lengths = dict(work)
    parts = f"{out_dir}/tmbed_parts"
    width = len(str(len(chunks)))

    def part_paths(i):
        tag = f"{i:0{width}d}"
        return (f"{parts}/{tag}.faa", f"{parts}/{tag}.pred",
                f"{parts}/{tag}.pred.part")

    def already(i):
        """A chunk counts as done when its committed output holds a record
        for every sequence that went into it."""
        pred = part_paths(i)[1]
        if not os.path.exists(pred):
            return False
        return sum(1 for _ in iter_tmbed_records(pred)) >= len(chunks[i])

    single = len(chunks) == 1 and not excluded
    os.makedirs(parts, exist_ok=True)
    if single:
        # Nothing was filtered and nothing is split, so TMbed reads the
        # original FASTA. This is the whole of the old behaviour, and it keeps
        # a second copy of the proteome off the disk for ordinary runs.
        inputs = {0: faa}
    else:
        pending = [i for i in range(len(chunks)) if not already(i)]
        inputs = {i: part_paths(i)[0] for i in range(len(chunks))}
        if pending:
            # Rewritten every run rather than reused: a chunk FASTA is written
            # before anything reads it, so a killed writer leaves a short one,
            # and a short input would quietly shrink the chunk.
            fhs = {i: open(inputs[i], "w", encoding="utf-8") for i in pending}
            try:
                for pid, seq in read_fasta(faa):
                    fh = fhs.get(where.get(pid, -1))
                    if fh is not None:
                        fh.write(f">{pid}\n{seq}\n")
            finally:
                for fh in fhs.values():
                    fh.close()

    errors = {}                       # chunk index -> why it failed
    consecutive = 0
    max_consecutive = int(cfg.get("tmbed_max_consecutive_failures", 2) or 0)
    t0, res_done = time.time(), 0
    for i, ids in enumerate(chunks):
        pred, part = part_paths(i)[1], part_paths(i)[2]
        res = sum(lengths[pid] for pid in ids)
        if already(i):
            log(f"tmbed: chunk {i + 1}/{len(chunks)} is already predicted; "
                "skipping")
            res_done += res
            continue
        if len(chunks) > 1:
            eta = ""
            if res_done:
                rate = (time.time() - t0) / res_done
                eta = (f", ~{(total_res - res_done) * rate / 3600:.1f}h left")
            log(f"tmbed: chunk {i + 1}/{len(chunks)}, {len(ids):,} "
                f"sequence(s), {res / 1e6:.1f}M residues{eta}")
        _atomic_rm(part)
        cmd = ["tmbed", "predict", "-f", inputs[i], "-p", part,
               "--out-format", "0"] + gpu
        bs = int(cfg.get("tmbed_batch_size", 0) or 0)
        if bs:
            cmd += ["--batch-size", str(bs)]
        cmd += tool_args(cfg, "tmbed")
        # run_cmd, not subprocess.run. This stage is the one the progress
        # heartbeat was written for - tmbed ran 2 h 36 min and then died,
        # twice, with nothing in the log between the command and the failure.
        # run_cmd takes the env, so CUDA_VISIBLE_DEVICES still pins the
        # device, and it resolves the binary through resolve_tool so tmbed is
        # launched by the path `have` found rather than by a bare name.
        try:
            run_cmd(cmd, env=env)
        except StageError:
            # die() raises StageError, which IS a RuntimeError. A bare
            # `except RuntimeError` below would swallow "tmbed not found" or a
            # record with an empty identifier and report it as a chunk that
            # failed, which is neither true nor retryable.
            raise
        except RuntimeError as e:
            # Whatever TMbed managed to write before it died stays where it
            # is: .pred.part is read by the concatenation below, which copies
            # complete records only. The chunk is NOT committed, so a rerun
            # retries it.
            kept = sum(1 for _ in iter_tmbed_records(part))
            msg = str(e).strip().splitlines()
            errors[i] = msg[-1] if msg else "no message"
            consecutive += 1
            log(f"tmbed: chunk {i + 1}/{len(chunks)} failed after writing "
                f"{kept}/{len(ids)} prediction(s), which are kept: "
                f"{errors[i]}", "WARN")
            if max_consecutive and consecutive >= max_consecutive:
                log(f"tmbed: {consecutive} chunk(s) in a row failed, so the "
                    f"remaining {len(chunks) - i - 1} are not attempted - a "
                    "device that has stopped responding fails all of them the "
                    "same way, slowly", "WARN")
                break
            continue
        consecutive = 0
        os.replace(part, pred)
        res_done += res

    # Concatenated through the record iterator rather than copied byte for
    # byte, so a truncated tail in a salvaged .pred.part cannot reach the
    # committed file. The record ORDER here is by length, not the order of
    # proteins_faa; nothing reads it positionally (parse_tmbed builds a dict).
    with atomic_out(p.tmbed) as tmp:
        written = 0
        with open(tmp, "w", encoding="utf-8") as out:
            for i in range(len(chunks)):
                pred, part = part_paths(i)[1], part_paths(i)[2]
                src = pred if os.path.exists(pred) else part
                for hdr, seq, lab in iter_tmbed_records(src):
                    out.write(f"{hdr}\n{seq}\n{lab}\n")
                    where.pop(hdr[1:].split()[0], None)
                    written += 1
        if where:
            # Name the casualties in a file rather than only in the log, so
            # the shortfall survives into the results directory and can be
            # read back by whoever asks why a protein has no topology.
            miss = f"{out_dir}/tmbed_failed.tsv"
            with open(miss, "w", encoding="utf-8") as fh:
                fh.write("protein_id\tlength\tchunk\terror\n")
                for pid in sorted(where, key=lambda q: -lengths.get(q, 0)):
                    i = where[pid]
                    fh.write(f"{pid}\t{lengths.get(pid, '')}\t{i}\t"
                             f"{errors.get(i, 'not attempted')}\n")
            log(f"tmbed: {len(where):,} of {len(work):,} protein(s) got no "
                f"prediction; listed in {miss}", "WARN")
            if not cfg.get("tmbed_allow_partial"):
                die(f"tmbed finished {written:,} of {len(work):,} "
                    f"prediction(s) across {len(chunks)} chunk(s).\n"
                    f"  Every completed chunk is kept in {parts}, so "
                    "rerunning resumes from there rather than starting "
                    "over.\n"
                    "  A card that has stopped responding usually needs the "
                    "machine or the WSL session restarted, not another "
                    "attempt.\n"
                    f"  The proteins that got nothing are listed in {miss}. "
                    "To go on without them, set tmbed_allow_partial: true.")
            log("tmbed: continuing with a partial topology set because "
                "tmbed_allow_partial is on; an absent helix or strand count "
                "here means not attempted, not absent", "WARN")
        log(f"tmbed: {written:,} prediction(s) -> {p.tmbed}")

    # Only once the committed file exists. Until then the parts ARE the
    # result, and a run that dies between the two must be able to resume.
    shutil.rmtree(parts, ignore_errors=True)


def stage_cluster(cfg, p):
    if not have("mmseqs"):
        die("mmseqs not found (conda install -c bioconda mmseqs2)")
    ram = cfg.get("ram_gb") or 0
    with atomic_out(p.cluster) as tmp:
        pref = f"{tmp}.mm"
        run_cmd(["mmseqs", "easy-cluster", cfg["proteins_faa"],
                 pref, f"{p.R}/cluster/tmp",
                 "--min-seq-id", cfg["thresholds"]["cluster_min_seq_id"],
                 "-c", cfg["thresholds"]["cluster_coverage"],
                 "--cov-mode", 0, "--cluster-mode", 0,
                 "--threads", cfg["threads"], "-v", 1]
                + (["--split-memory-limit", f"{ram}G"] if ram else [])
                + tool_args(cfg, "mmseqs"))
        got = f"{pref}_cluster.tsv"
        if not os.path.exists(got):
            die(f"mmseqs exited 0 but wrote no {got}; nothing to adopt as "
                f"{p.cluster}")
        # The representative/all-sequence fastas keep the names the docs use.
        for suf in ("_rep_seq.fasta", "_all_seqs.fasta"):
            if os.path.exists(pref + suf):
                os.replace(pref + suf, f"{p.R}/cluster/fam{suf}")
        os.replace(got, tmp)


def stage_ncbifam(cfg, p):
    if not have("hmmsearch"):
        die("hmmsearch not found (conda install -c bioconda hmmer)")
    with atomic_out(p.ncbifam) as tmp:
        run_cmd(["hmmsearch",
                 cfg["thresholds"].get("ncbifam_cutoff", "--cut_tc"),
                 "--noali", "--cpu", cfg["threads"], "--tblout", tmp,
                 "-o", os.devnull]
                + tool_args(cfg, "hmmsearch")
                + [cfg["db"]["ncbifam_hmm"], cfg["proteins_faa"]])


def stage_kofam(cfg, p):
    exe = "exec_annotation" if have("exec_annotation") else (
        "kofamscan" if have("kofamscan") else None)
    if not exe:
        die("KOfamScan not found (conda install -c bioconda kofamscan), or set "
            "run.kofam: false")
    os.makedirs(f"{p.R}/kofam", exist_ok=True)
    with atomic_out(p.kofam) as tmp:
        run_cmd([exe, "-o", tmp, "-f", "detail-tsv",
                 "--no-report-unannotated",
                 "-p", cfg["db"]["kofam_profiles"],
                 "-k", cfg["db"]["kofam_ko_list"],
                 "--cpu", cfg["threads"], "--tmp-dir", f"{p.R}/kofam/tmp"]
                + tool_args(cfg, "kofamscan") + [cfg["proteins_faa"]])


# InterProScan only accepts IUPAC residues. Everything else here is mapped to
# X rather than dropped, so residue offsets in the TSV still line up with the
# sequences the rest of the pipeline sees.
_NON_IUPAC = str.maketrans({c: "X" for c in "JUOBZ*"})


def sanitise_faa(src, dst):
    """Copy a protein FASTA, replacing residues InterProScan refuses.

    A single '*' (Prodigal/Macrel keep the terminal stop) aborts InterProScan
    with a Java exception hours into the longest stage of the run, so the
    query gets cleaned instead of being passed through verbatim.
    Returns (records, records_changed).
    """
    n = n_fixed = 0
    with open(dst, "w", encoding="utf-8") as fh:
        for pid, seq in read_fasta(src):
            n += 1
            clean = seq.upper().rstrip("*").translate(_NON_IUPAC)
            if clean != seq:
                n_fixed += 1
            fh.write(f">{pid}\n{clean}\n")
    return n, n_fixed


def stage_interpro(cfg, p):
    exe = cfg["db"].get("interproscan_sh") or "interproscan.sh"
    if not (os.path.exists(exe) or have(exe)):
        die(f"InterProScan not found at '{exe}'; set db.interproscan_sh or "
            "run.interpro: false")
    apps = cfg.get("interpro_applications", "")
    os.makedirs(f"{p.R}/interpro/tmp", exist_ok=True)
    query = f"{p.R}/interpro/query.faa"
    n, n_fixed = sanitise_faa(cfg["proteins_faa"], query)
    if n_fixed:
        log(f"interpro: {n_fixed}/{n} sequences carried '*' or a non-IUPAC "
            f"letter; searching a cleaned copy at {query}", "WARN")
    # -dp disables the precalculated lookup: it only holds UniParc matches, so
    # for novel metagenome ORFs it is a network round trip that finds nothing.
    with atomic_out(p.interpro) as out_tmp:
        cmd = [exe, "-i", query, "-f", "TSV", "-o", out_tmp,
               "-cpu", cfg["threads"], "-iprlookup", "-goterms", "-dp",
               "-T", f"{p.R}/interpro/tmp"]
        if apps:
            cmd += ["-appl", apps]
        ram = cfg.get("ram_gb") or 0
        env = dict(os.environ)
        if ram:
            # _JAVA_OPTIONS is read by the JVM itself, so it survives
            # whichever launcher script this InterProScan build uses.
            heap = f"-Xmx{max(2, ram)}g"
            env["_JAVA_OPTIONS"] = (
                env.get("_JAVA_OPTIONS", "") + " " + heap).strip()
            env["JAVA_OPTS"] = (env.get("JAVA_OPTS", "") + " " + heap).strip()
            log(f"JVM heap {heap} from a {ram} GB budget")
        run_cmd(cmd + tool_args(cfg, "interproscan"), env=env)


def _count_fasta(path):
    """Records in a FASTA, without holding any of it in memory."""
    n = 0
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def profile_query_faa(cfg, p, what):
    """Query FASTA for the profile searches, or None when there is no work.

    dark_all.faa, not dark.faa. dark.faa is the ESMFold work-list, so it
    carries max_dark_structures and max_len_structure — a GPU-hours budget
    and a folding memory limit. Querying it made profile search a function of
    that budget, and turned the 3p_profile_only bin (and any "rescued by
    profile search" fraction) into an artefact of an effector_score ranking
    over an arbitrary 2,000 proteins.

    An empty file means "every protein already carries some annotation", not
    "no work-list": the old fallback read that as a licence to search the
    whole proteome, i.e. one hhblits run per identified protein.
    """
    src = p.dark_all
    if not os.path.exists(src):
        if os.path.exists(p.dark):
            # An older results directory, or the two-machine hand-off where
            # only dark.faa was rsynced across. Searching the capped work-list
            # is wrong but bounded; searching the whole proteome is neither.
            log(f"{what}: {p.dark_all} is missing, so this falls back to the "
                f"CAPPED structure work-list {p.dark} "
                f"(max_dark_structures={cfg.get('max_dark_structures')}, "
                f"max_len_structure={cfg.get('max_len_structure')}). Any "
                "'rescued by profile search' fraction is then computed on "
                "that subset rather than on every unannotated protein — "
                "rerun the integrate stage, or copy dark_all.faa across "
                "as well.", "WARN")
            src = p.dark
        else:
            log(f"{what}: neither {p.dark_all} nor {p.dark} exists, so this "
                "falls back to the WHOLE proteome — one search per identified "
                "protein, which is almost certainly not what you want. Run "
                "the integrate stage first.", "WARN")
            return cfg["proteins_faa"]
    if not nonempty(src):
        log(f"{what}: no unannotated proteins to search, skipping")
        return None
    log(f"{what}: querying {_count_fasta(src)} unannotated proteins "
        f"from {src}")
    return src


def stage_hhblits(cfg, p):
    if not have("hhblits"):
        die("hhblits not found (conda install -c bioconda hhsuite), or set "
            "run.hhblits: false")
    os.makedirs(p.hhr_dir, exist_ok=True)
    src = profile_query_faa(cfg, p, "hhblits")
    if src is None:
        open(p.hhr_done, "w", encoding="utf-8").close()
        return
    todo = [(pid, seq) for pid, seq in read_fasta(src)
            if not os.path.exists(f"{p.hhr_dir}/{pid}.hhr")]
    if not todo:
        log("hhblits: nothing new to search")
        open(p.hhr_done, "w", encoding="utf-8").close()
        return

    workers = max(1, int(cfg.get("hhblits_workers", 4)))
    per = max(1, int(cfg["threads"]) // workers)
    ram = cfg.get("ram_gb") or 0
    # hhblits -maxmem is GB per process and defaults to 3, which silently
    # truncates alignments on large databases.
    hh_mem = ["-maxmem", f"{max(1.0, (ram / workers)):.1f}"] if ram else []
    log(f"{len(todo)} queries, {workers} workers x {per} cpu"
        + (f" x {hh_mem[1]} GB" if hh_mem else ""))

    def one(job):
        pid, seq = job
        # Per-worker query file: a single shared temp path was a race waiting
        # to happen, and hhblits takes one query at a time regardless.
        q = f"{p.hhr_dir}/_q_{os.getpid()}_{threading.get_ident()}.fasta"
        with open(q, "w", encoding="utf-8") as fh:
            fh.write(f">{pid}\n{seq}\n")
        # Write via .part and rename: the resume test is "does <pid>.hhr
        # exist", so a query killed mid-write would otherwise be adopted as
        # finished and that protein's hits lost silently.
        part = f"{p.hhr_dir}/{pid}.hhr.part"
        run_cmd(["hhblits", "-i", q, "-d", cfg["db"]["hhblits_db"],
                 "-o", part,
                 "-n", cfg.get("hhblits_iterations", 2),
                 "-cpu", per, "-v", 0]
                + hh_mem + tool_args(cfg, "hhblits"))
        os.replace(part, f"{p.hhr_dir}/{pid}.hhr")

    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        for _ in ex.map(one, todo):
            done += 1
            if done % 100 == 0:
                log(f"hhblits: {done}/{len(todo)}")
    for f in glob.glob(f"{p.hhr_dir}/_q_*.fasta") + glob.glob(
            f"{p.hhr_dir}/*.hhr.part"):
        os.remove(f)
    log(f"hhblits: {done} new profiles searched into {p.hhr_dir}")
    open(p.hhr_done, "w", encoding="utf-8").close()


def stage_jackhmmer(cfg, p):
    if not have("jackhmmer"):
        die("jackhmmer not found (part of HMMER), or set run.jackhmmer: false")
    src = profile_query_faa(cfg, p, "jackhmmer")
    if src is None or not nonempty(src):
        log("nothing to search, skipping")
        with atomic_out(p.jackhmmer) as tmp:
            open(tmp, "w", encoding="utf-8").close()
        return
    with atomic_out(p.jackhmmer) as tmp:
        run_cmd(["jackhmmer", "-N", cfg.get("jackhmmer_iterations", 3),
                 "--noali", "--cpu", cfg["threads"],
                 "-E", cfg["thresholds"].get("jackhmmer_evalue", 1e-5),
                 "--tblout", tmp, "-o", os.devnull]
                + tool_args(cfg, "jackhmmer")
                + [src, cfg["db"]["jackhmmer_db"]])


def stage_smorf(cfg, p):
    """Small-ORF calling. Runs on CONTIGS, not proteins: Prodigal in meta mode
    discards most ORFs below ~90 nt, so bacteriocins, TA toxins and RiPPs were
    never candidates. The output has to be added to the search database and the
    MS data re-searched — no annotation stage can recover them otherwise."""
    fna = cfg.get("contigs_fna", "")
    if not fna or not os.path.exists(fna):
        die("run.smorf needs contigs_fna pointing at the assembly")
    os.makedirs(f"{p.R}/smorf", exist_ok=True)
    made = []
    # SMORFinder's pip package is called smorfinder but its console script is
    # `smorf`; asking for the package name meant this half never ran. `meta` is
    # the metagenome-assembly mode (and the only one that takes -t); `single`
    # is for an isolate genome.
    exe = "smorf" if have("smorf") else (
        "smorfinder" if have("smorfinder") else None)
    mode = str(cfg.get("smorf_mode", "meta")).strip().lower()
    if mode not in ("meta", "single"):
        die(f"unknown smorf_mode '{mode}'; choose meta (assembly) or "
            "single (isolate genome)")
    if exe:
        cmd = [exe, mode, "-o", f"{p.R}/smorf/smorfinder"]
        if mode == "meta":
            cmd += ["-t", cfg["threads"]]
        run_cmd(cmd + tool_args(cfg, "smorfinder") + [fna])
        made += glob.glob(f"{p.R}/smorf/smorfinder/*.faa")
    else:
        log("smorf/smorfinder not found; skipping that half", "WARN")
    if have("macrel"):
        run_cmd(["macrel", "contigs", "--fasta", fna,
                 "--output", f"{p.R}/smorf/macrel", "-t", cfg["threads"]]
                + tool_args(cfg, "macrel"))
        # Only the smorfs file: macrel also writes macrel.out.all_orfs.faa,
        # the entire Prodigal ORF catalogue, whose ids are identical, so the
        # old `*.faa*` glob quietly made this stage emit every ORF in the
        # assembly under the name "small ORFs".
        made += glob.glob(f"{p.R}/smorf/macrel/*smorfs.faa*")
    else:
        log("macrel not found; skipping that half", "WARN")
    if not made:
        # Warn and skip rather than die: every other stage is independent of
        # this one, and losing a whole run over a missing optional tool is a
        # worse outcome than an empty candidate list.
        log("neither smorf(inder) nor macrel is installed, so no small ORFs "
            "were called; install macrel (pip install macrel) if you want "
            "this bin", "WARN")
        with atomic_out(p.smorf_faa) as tmp:
            open(tmp, "w", encoding="utf-8").close()
        return
    maxlen = int(cfg["thresholds"].get("smorf_max_len", 100))
    seen, over, starred = set(), 0, 0
    with atomic_out(p.smorf_faa) as tmp, open(tmp, "w", encoding="utf-8") as out:
        for f in made:
            kept = 0
            for pid, seq in read_fasta(f):
                # Prodigal-style callers keep the terminal stop as '*'. It is
                # not a residue: it adds 1 to the length the smorf_max_len
                # test sees and InterProScan rejects the whole file over it.
                clean = seq.replace("*", "")
                if clean != seq:
                    starred += 1
                if pid in seen:
                    continue
                if len(clean) > maxlen:
                    over += 1
                    continue
                seen.add(pid)
                out.write(f">{pid}\n{clean}\n")
                kept += 1
            log(f"smorf: {kept} kept from {f}")
    log(f"smorf: {len(seen)} small ORFs (<= {maxlen} aa) -> {p.smorf_faa}; "
        f"{over} dropped as too long, {starred} had a '*' stripped")
    log("ACTION REQUIRED: append these to your MS search database and "
        "re-search. They are absent from the current results by construction, "
        "so their absence is not evidence of absence.", "WARN")


# ======================================================================
# stage: genomic context
# ======================================================================
CONTEXT_PATTERNS = {
    "t6ss": re.compile(
        r"\bVgrG\b|\bHcp\b|\bTss[A-M]\b|IcmF|ImpA|EvpB|PAAR|\bRhs\b|"
        r"type VI secretion", re.I),
    "t3ss": re.compile(
        r"\bSct[A-Z]\b|\bYop[A-Z]\b|\bHrp[A-Z]\b|\bEsc[A-Z]\b|\bSpa[0-9]|"
        r"\bInv[A-J]\b|type III secretion|needle", re.I),
    "secretion_other": re.compile(
        r"type (I|II|IV|V|VII|IX) secretion|\bTad[A-Z]\b|\bVir[BD][0-9]|"
        r"autotransporter|\bTps[AB]\b", re.I),
    "susc_susd": re.compile(
        r"\bSusC\b|\bSusD\b|TonB-dependent|TonB_dep_Rec|SusD-like|"
        r"RagB|Plug_translocon", re.I),
}

# prophage and bgc are not single substrings. Bare `phage` matched inside
# `macrophage`, the 3-letter `Xre` matched every Xre-family regulator, and
# `acyl carrier` / `ketosyn` / `AMP-binding` / `radical SAM` match core
# fatty-acid, CoA-ligase and biotin biosynthesis - present in every genome.
# Since both flags fold into context_mge (+2 on the effector score), a
# near-constant +2 destroyed the score's ability to rank. Each flag now needs
# one family-diagnostic term plus a second, DIFFERENT family in the window.
CONTEXT_FAMILIES = {
    "prophage": {
        "terminase": r"\bterminase\b",
        "head": r"\bcapsid\b|\bportal protein\b|\bprohead\b|\bhead protein\b",
        "tail": r"tail fib|\bbaseplate\b|major tail|tail sheath|tail tape",
        "phage": r"\bphages?\b|\bprophages?\b|\bbacteriophage\b",
        "lysis": r"\bholin\b|\bendolysin\b",
        "integration": r"\bintegrase\b|\bexcisionase\b|\bantirepressor\b",
    },
    "bgc": {
        "condensation": r"\bcondensation domain\b|PF00668|\bNRPS\b|"
                        r"non-?ribosomal peptide synthetase",
        "polyketide": r"\bpolyketide synthase\b|\bPKS\b|\btrans-AT\b",
        # two RiPP families, not one, so that a lanthipeptide synthetase
        # beside a YcaO cyclodehydratase counts as two independent hits
        "lanthipeptide": r"\bLan[BCM]\b|lanthipeptide|\bnisin\b",
        "ripp_other": r"\bYcaO\b|lasso peptide|\bthiopeptide\b|"
                      r"\bsactipeptide\b|\bPqqD\b",
        "ketosynthase": r"beta-ketoacyl|\bketosyn|PF00109|PF02801",
        "carrier": r"\bacyl carrier\b|phosphopantetheine|PF00550",
        "adenylation": r"\bAMP-binding\b|PF00501",
        "thioesterase": r"\bthioesterase\b|PF00975|PF07859",
    },
}
# At least one matched family must be on this list, so that FabB beside AcpP
# (the fatty-acid operon of every bacterium) is not read as a BGC and a lone
# integrase beside a holin is not read as a prophage.
CONTEXT_DIAGNOSTIC = {
    "prophage": ("terminase", "head", "tail", "phage"),
    "bgc": ("condensation", "polyketide", "lanthipeptide", "ripp_other"),
}
CONTEXT_FAMILIES_RE = {
    flag: {fam: re.compile(rx, re.I) for fam, rx in fams.items()}
    for flag, fams in CONTEXT_FAMILIES.items()
}

# Toxins that never see the Sec machinery - T6SS effectors, CDI/CdiA,
# colicin-like bacteriocins - dominate gut anaerobes, so the immunity-pair
# test cannot be gated on a signal peptide alone.
TOXIN_HINT_RE = re.compile(
    r"\btoxin\b|\bbacteriocin\b|\bcolicin\b|\bpyocin\b|\bRhs\b|\bCdi[AB]\b|"
    r"contact-dependent (growth )?inhibition|\bLXG\b|\bMaf[AB]\b|"
    r"\bTs[ei][0-9]\b|\bVgrG\b|\bPAAR\b|ADP-ribosyltransferase", re.I)

GFF_TEXT_ATTRS = ("gene", "product", "Note", "eC_number", "gene_functions")


def _family_flag(flag, txt):
    """True when two different core families of `flag` are present in `txt`
    and at least one of them is diagnostic for that flag."""
    hits = {fam for fam, rx in CONTEXT_FAMILIES_RE[flag].items() if rx.search(txt)}
    return len(hits) >= 2 and bool(hits & set(CONTEXT_DIAGNOSTIC[flag]))


def parse_gff(path):
    """CDS coordinates plus whatever the annotator wrote about each ORF.

    The `product=`/`gene=` text is the point: in a metaproteome only a few
    percent of the ORFs on the assembly are identified by MS, so a
    neighbourhood built solely from the MS-derived evidence tables is empty
    for almost every window and every flag collapses to False while the GFF
    on disk names the neighbours explicitly.
    """
    per_type = {"CDS": defaultdict(list), "gene": defaultdict(list)}
    n_bad = 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in per_type:
                continue
            attrs = {}
            for kv in f[8].rstrip(";").split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    # GFF3 escapes the separators inside attribute values.
                    for esc, ch in (("%2C", ","), ("%3B", ";"), ("%3D", "="),
                                    ("%20", " "), ("%25", "%")):
                        v = v.replace(esc, ch).replace(esc.lower(), ch)
                    attrs[k.strip()] = v.strip()
            # locus_tag before ID: NCBI and Bakta write ID=cds-WP_... /
            # ID=gene-b1, neither of which joins to the protein FASTA, while
            # locus_tag does. Prokka sets ID == locus_tag, so it is unaffected.
            pid = (attrs.get("locus_tag") or attrs.get("ID")
                   or attrs.get("Name") or attrs.get("protein_id"))
            if not pid:
                continue
            pid = pid.split()[0]
            for pre in ("cds-", "gene-", "rna-"):
                if pid.startswith(pre):
                    pid = pid[len(pre):]
                    break
            desc = " ".join(attrs[k] for k in GFF_TEXT_ATTRS if attrs.get(k))
            try:
                per_type[f[2]][f[0]].append(
                    (int(f[3]), int(f[4]), f[6], pid, desc))
            except ValueError:
                n_bad += 1
    if n_bad:
        log(f"gff: {n_bad} feature lines had unparsable coordinates", "WARN")

    # NCBI/Bakta emit a `gene` row AND a `CDS` row per locus. Keeping both
    # doubled every ORF, which silently halved the physical span of the
    # window. CDS is authoritative; `gene` is only a fallback for GFFs that
    # carry nothing else.
    contigs = per_type["CDS"]
    if not contigs and per_type["gene"]:
        contigs = per_type["gene"]
        log("gff: no CDS features found, falling back to `gene` features", "WARN")
    n_dup = 0
    for c in list(contigs):
        seen, keep = set(), []
        for e in sorted(contigs[c]):
            key = (e[0], e[1], e[3])
            if key in seen:
                n_dup += 1
                continue
            seen.add(key)
            keep.append(e)
        contigs[c] = keep
    if n_dup:
        log(f"gff: dropped {n_dup} duplicate features (same id and interval)",
            "WARN")
    return contigs


def stage_context(cfg, p):
    gff = cfg.get("gff") or ""
    if not gff:
        log("run.context is on but `gff` is empty; skipping", "WARN")
        with atomic_out(p.context) as tmp:
            pd.DataFrame(columns=["protein_id"]).to_csv(
                tmp, sep="\t", index=False)
        return
    if not os.path.exists(gff):
        die(f"gff not found: {gff}")

    contigs = parse_gff(gff)
    em = parse_emapper(p.emapper) if os.path.exists(p.emapper) else None
    pf = parse_hmm_tblout(p.pfam) if os.path.exists(p.pfam) else {}
    sp = parse_signalp6(p.signalp) if os.path.exists(p.signalp) else {}
    dc = (parse_hmm_domtblout(p.dbcan, cfg["thresholds"]["dbcan_min_cov"],
                              cfg["thresholds"]["dbcan_evalue"])
          if os.path.exists(p.dbcan) else {})

    # The GFF annotation is the base layer: it covers every ORF on the
    # assembly, not just the MS-identified subset. Evidence from the search
    # stages is overlaid on top of it.
    text = {}
    n_gff_text = 0
    for orfs in contigs.values():
        for o in orfs:
            if o[4]:
                text[o[3]] = o[4]
                n_gff_text += 1
    if em is not None:
        for c in ("Description", "Preferred_name", "PFAMs"):
            if c in em.columns:
                for pid, v in em[c].items():
                    text[pid] = text.get(pid, "") + " " + str(v)
    for pid, hits in pf.items():
        text[pid] = text.get(pid, "") + " " + " ".join(h[0] for h in hits)
    for pid, hits in dc.items():
        text[pid] = text.get(pid, "") + " " + " ".join(h[0] for h in hits)

    win = cfg["context_window"]
    max_len = cfg["immunity_max_len"]
    max_gap = cfg["immunity_max_gap"]
    rows = []
    n = sum(len(v) for v in contigs.values())
    log(f"context: {n} ORFs on {len(contigs)} contigs")

    cand = {}          # partner id -> id of the toxin it sits next to
    n_neigh = n_neigh_text = 0
    for orfs in contigs.values():
        for i, (start, end, strand, pid, _desc) in enumerate(orfs):
            lo, hi = max(0, i - win), min(len(orfs), i + win + 1)
            others = [o for o in orfs[lo:hi] if o[3] != pid]
            n_neigh += len(others)
            n_neigh_text += sum(1 for o in others if text.get(o[3], "").strip())
            neigh = " ".join(text.get(o[3], "") for o in others)
            flags = {k: bool(rx.search(neigh)) for k, rx in CONTEXT_PATTERNS.items()}
            flags["prophage"] = _family_flag("prophage", neigh)
            flags["bgc"] = _family_flag("bgc", neigh)
            immunity, partner = False, ""
            if (sp.get(pid, ("", ""))[0] in ("SP", "LIPO", "TAT", "TATLIPO")
                    or TOXIN_HINT_RE.search(text.get(pid, ""))):
                # The partner is the gene downstream in TRANSCRIPTION order.
                # orfs is sorted by coordinate regardless of strand, so on the
                # minus strand that is orfs[i-1]; only ever looking at
                # orfs[i+1] made every minus-strand cassette invisible.
                j = i + 1 if strand != "-" else i - 1
                if 0 <= j < len(orfs):
                    ns, ne, nstr, npid, _ = orfs[j]
                    gap = (ns - end) if strand != "-" else (start - ne)
                    # -1: the CDS span includes the stop codon.
                    if (nstr == strand and gap <= max_gap
                            and (ne - ns + 1) // 3 - 1 <= max_len
                            and sp.get(npid, ("", ""))[0] in ("", "OTHER")):
                        immunity, partner = True, npid
                        cand[npid] = pid
            # Polysaccharide utilisation locus: several CAZymes plus a
            # SusC/SusD-like importer in the same neighbourhood. This is what
            # dbCAN-PUL looks for, computed from evidence already in hand.
            # The query itself is excluded, exactly as it is for the pattern
            # flags above; counting it let a lone SusC call its own PUL.
            window_ids = [o[3] for o in others]
            n_caz = sum(1 for o in window_ids if dc.get(o))
            has_sus = any(CONTEXT_PATTERNS["susc_susd"].search(text.get(o, ""))
                          for o in window_ids)
            flags["pul"] = n_caz >= cfg.get("pul_min_cazymes", 2) and has_sus
            flags["n_cazymes_in_window"] = n_caz
            flags["immunity_pair"] = immunity
            flags["immunity_partner"] = partner
            rows.append({"protein_id": pid, **flags})

    # The interesting protein in a toxin/immunity cassette is the small,
    # usually dark partner, not the toxin - so it gets a flag of its own
    # instead of appearing only as a string in the toxin's row.
    for r in rows:
        toxin = cand.get(r["protein_id"], "")
        r["immunity_candidate"] = bool(toxin)
        r["immunity_toxin"] = toxin

    with atomic_out(p.context) as tmp:
        pd.DataFrame(rows).to_csv(tmp, sep="\t", index=False)

    # An all-False context.tsv is indistinguishable from "no interesting
    # neighbourhoods found", so say what was actually loaded and what fired.
    loaded = {"gff_product": n_gff_text > 0, "emapper": em is not None,
              "pfam": bool(pf), "signalp": bool(sp), "dbcan": bool(dc)}
    log("context: evidence " + ", ".join(
        f"{k}={'yes' if v else 'absent'}" for k, v in loaded.items()))
    if not loaded["signalp"]:
        log("context: no signalp table -> immunity_pair can only fire on a "
            "toxin-like description", "WARN")
    if not loaded["dbcan"]:
        log("context: no dbcan table -> pul cannot fire", "WARN")
    if n_neigh:
        frac = 100.0 * n_neigh_text / n_neigh
        log(f"context: {frac:.1f}% of window neighbours carry annotation text",
            "WARN" if frac < 25 else "INFO")
    keys = list(CONTEXT_PATTERNS) + ["prophage", "bgc", "pul", "immunity_pair",
                                     "immunity_candidate"]
    pos = {k: sum(1 for r in rows if r.get(k) is True) for k in keys}
    log("context: " + ", ".join(f"{k}={v}" for k, v in pos.items()))
    if rows and not any(pos.values()):
        log(f"context: every flag is False across {len(rows)} ORFs; the "
            "context columns will contribute nothing to the score", "WARN")
    log(f"context: wrote {p.context}")


# ======================================================================
# stage: integrate
# ======================================================================
def ncbifam_family_desc(cfg, p):
    """NCBIfam/TIGRFAM family name -> DESC, cached beside the tblout.

    Returns {} when the library cannot be read, and the caller then counts
    every NCBIfam hit as informative — the behaviour that existed before this
    test, never a silent demotion on missing data.
    """
    cache, lib = p.ncbifam_desc, (cfg.get("db") or {}).get("ncbifam_hmm") or ""
    have_lib = bool(lib) and os.path.exists(lib)
    # Rebuild when the library is newer than the cache: pointing ncbifam_hmm
    # at a different library and silently reusing the old descriptions would
    # mis-call families rather than fail.
    fresh = (have_lib and os.path.exists(cache)
             and os.path.getmtime(cache) >= os.path.getmtime(lib))
    if nonempty(cache) and (not have_lib or fresh):
        out = {}
        with opener(cache) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[0]:
                    out[f[0]] = f[1]
        if out:
            return out
        log(f"{cache}: no family descriptions parsed; rebuilding", "WARN")
    if not have_lib:
        log("db.ncbifam_hmm is not readable, so NCBIfam family descriptions "
            "are unknown and every NCBIfam hit is counted as informative "
            f"({lib or 'unset'})", "WARN")
        return {}
    # Gigabytes of HMM library for a few MB of DESC lines, so this is done
    # once per library and cached.
    log(f"reading NCBIfam family descriptions from {lib} (cached in {cache})")
    d = parse_hmm_lib_desc(lib)
    if not d:
        log(f"{lib}: no NAME/DESC pairs found. That is not an HMMER library "
            "in the expected format; every NCBIfam hit will be counted as "
            "informative.", "WARN")
        return {}
    try:
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        with open(cache, "w", encoding="utf-8") as fh:
            for k in sorted(d):
                fh.write(f"{k}\t{d[k]}\n")
    except OSError as e:
        log(f"could not write {cache} ({e}); the library will be re-read next "
            "run", "WARN")
    log(f"ncbifam: {len(d)} family descriptions available")
    return d


def structure_shortfall_message(cfg, p, requested_lengths, have_pdb):
    """(text, level) for "fewer models than we asked for".

    On the very FIRST structure run this said

        WARN 0/1913 requested structures exist in .../results/structures;
             the rest were skipped (OOM) or never folded

    which is a correct sentence in the wrong situation: the directory was
    empty because nothing had been folded YET, and a reader reasonably
    concluded the run had already lost 1,913 models to the OOM killer. The
    three causes are not the same event and must not share a sentence:

      never folded   no .done marker, no plddt.tsv, no models - the esmfold
                     stage has not run, so nothing has been lost and this is
                     not a warning at all;
      still folding  models exist but no .done marker - an interrupted or
                     not-yet-finished pass, so the rest are pending;
      lost           esmfold finished and models are still missing. Even here
                     OOM is only named for the ones length does not already
                     explain, because max_len_structure excludes a protein
                     from folding before ESMFold sees it.
    """
    n_req = len(requested_lengths)
    missing = [pid for pid in requested_lengths if pid not in have_pdb]
    n_got = n_req - len(missing)
    cap = int(cfg.get("max_len_structure") or 0)
    too_long = sum(1 for pid in missing if cap and requested_lengths[pid] > cap)
    finished = os.path.exists(p.struct_done)
    # Any of the three is proof that folding has happened at least once;
    # .done alone is not, because it is written at the END of the stage.
    started = finished or n_got > 0 or nonempty(f"{p.structures}/plddt.tsv")

    if not started:
        return (f"none of the {n_req} proteins in {p.dark} have been folded "
                f"yet: {p.structures} is empty and the esmfold stage has not "
                f"run (no {p.struct_done}, no plddt.tsv). Nothing has been "
                "lost - this pass simply has no structural evidence, so "
                "structure_attempted is False for every protein.", "INFO")
    if not finished:
        return (f"{n_got}/{n_req} requested structures exist in "
                f"{p.structures}, and esmfold has not finished (no "
                f"{p.struct_done}), so the {len(missing)} missing are still "
                "pending rather than lost; rerun the esmfold stage.", "WARN")

    parts = [f"{n_got}/{n_req} requested structures exist in {p.structures} "
             "although esmfold has finished"]
    if too_long:
        parts.append(f"{too_long} are longer than max_len_structure={cap} and "
                     "were never submitted")
    rest = len(missing) - too_long
    if rest:
        # esmfold_failed.tsv is the record of what it actually gave up on, so
        # the two causes no longer have to share one hedged sentence. Missing
        # proteins that are NOT in it were added to dark.faa after the last
        # fold and were never attempted, which is a different thing entirely.
        tried = set()
        fail_tbl = f"{p.structures}/esmfold_failed.tsv"
        if os.path.exists(fail_tbl):
            with contextlib.suppress(OSError):
                with opener(fail_tbl) as fh:
                    next(fh, None)
                    tried = {l.split("	")[0] for l in fh if l.strip()}
        gave_up = sum(1 for pid in missing if pid in tried)
        if gave_up:
            parts.append(f"{gave_up} were attempted and failed twice, listed "
                         f"with the error in {fail_tbl}")
        never = rest - gave_up
        if never and tried:
            parts.append(f"{never} are absent from that list, so they were "
                         "added to dark.faa after the last fold and never "
                         "attempted")
        elif never:
            parts.append(f"{never} were skipped or never folded - there is no "
                         f"{os.path.basename(fail_tbl)}, so this run predates "
                         "the failure record and the two cannot be separated")
    return ("; ".join(parts), "WARN")


def build_annotation(cfg, p, emit_dark=None, emit_dark_all=None):
    th, w = cfg["thresholds"], cfg["weights"]

    # DUF_RE names Pfam DUF families and nothing else. The descriptions that
    # come back from InterProScan, HHblits and Foldseek need a wider net:
    # Pfam's UPF0xxx families, PANTHER's "FAMILY NOT NAMED" and TrEMBL's
    # "hypothetical protein" are exactly as uninformative as a DUF, and
    # counting them as annotation is what drains 3d_duf_only and 4_dark.
    uninf_re = re.compile(
        r"^DUF\d|^UPF\d|unknown function|uncharacteri|hypothetical|"
        r"family not named|predicted protein|putative protein",
        re.IGNORECASE)

    # One pass over the FASTA computing length and the C-terminal motif.
    # Retaining every sequence cost ~1 GB on a full ORF catalogue and bought
    # nothing: the dark-bin FASTA is written by re-reading the file.
    lo = th["lpxtg_min_offset"]
    hi = th["lpxtg_max_offset"]

    def _lpxtg(sq):
        end = len(sq) - lo
        if end <= 0:
            return False
        win_start = max(0, len(sq) - hi)
        m = LPXTG_RE.search(sq[win_start:end])
        if not m:
            return False
        # A bare LP.TG 4-mer occurs by chance in ~0.03% of proteins, and on its
        # own it used to set surface_or_secreted, which is the shortlist gate.
        # A real sortase substrate follows the motif with a hydrophobic
        # membrane anchor and a short positively charged tail; requiring both
        # is what separates the signal from a coincidence.
        tail = sq[win_start + m.end():]
        if len(tail) < 10:
            return False
        hydro = sum(c in "AVILMFWCGP" for c in tail)
        basic = sum(c in "KR" for c in tail[-7:])
        return hydro >= 0.5 * len(tail) and basic >= 1

    ids, lens, lpx, seen, dups = [], [], [], set(), 0
    for pid, sq in read_fasta(cfg["proteins_faa"]):
        if pid in seen:
            dups += 1
            continue
        seen.add(pid)
        ids.append(pid)
        lens.append(len(sq))
        lpx.append(_lpxtg(sq))
    if not ids:
        die(f"no sequences in {cfg['proteins_faa']}")
    if dups:
        log(f"{dups} duplicate FASTA ids, first kept — duplicates break the "
            "join to the quant table; deduplicate the database first", "WARN")
    idx = pd.Index(ids, name="protein_id")
    df = pd.DataFrame({"length": lens, "lpxtg": lpx}, index=idx)
    log(f"integrate: {len(df)} proteins")

    # ---- orthology
    em = None
    if os.path.exists(p.emapper):
        em = parse_emapper(p.emapper)
        ndup = int(em.index.duplicated().sum())
        if ndup:
            log(f"{ndup} duplicate query ids in emapper output, first kept", "WARN")
            em = em[~em.index.duplicated(keep="first")]
        nhit = len(em.index.intersection(idx))
        log(f"emapper covers {nhit}/{len(idx)} "
            f"({100*nhit/len(idx):.1f}%) of the input proteins")
        if nhit == 0:
            log("zero id overlap between the fasta and the emapper table — "
                "almost always an id formatting mismatch, not missing "
                "annotation", "SEVER")
    else:
        log("no eggnog annotations found", "WARN")

    df["ko"] = col(em, "KEGG_ko", idx)
    df["cog_cat"] = col(em, "COG_category", idx)
    df["description"] = col(em, "Description", idx)
    df["preferred_name"] = col(em, "Preferred_name", idx)
    df["ec"] = col(em, "EC", idx)
    df["pfams_emapper"] = col(em, "PFAMs", idx)
    df["cazy"] = col(em, "CAZy", idx)
    df["og"] = col(em, "eggNOG_OGs", idx).str.split(",").str[0]
    df["seed_taxid"] = col(em, "seed_ortholog", idx).str.split(".").str[0]
    df["n_pathway_specific"] = col(em, "KEGG_Pathway", idx).map(
        lambda s: len(set(MAP_RE.findall(s)) - GLOBAL_MAPS))
    df["has_ko"] = df["ko"].str.contains(KO_RE, na=False)
    df["in_specific_pathway"] = df["n_pathway_specific"] > 0

    # KEGG_ko x KEGG_Pathway pairs across the whole emapper table give a
    # run-local KO -> map lookup for free. Without one, a KO that only
    # KOfamScan found had n_pathway_specific = 0 by construction and could
    # never leave 2_ko_orphan, so the advertised KOfam rescue and the
    # "invisible to KEGG enrichment" number contradicted each other.
    ko2maps = {}
    if em is not None and {"KEGG_ko", "KEGG_Pathway"} <= set(em.columns):
        for kos, paths in zip(em["KEGG_ko"].fillna(""),
                              em["KEGG_Pathway"].fillna("")):
            maps = set(MAP_RE.findall(str(paths))) - GLOBAL_MAPS
            if not maps:
                continue
            for k in KO_RE.findall(str(kos)):
                ko2maps.setdefault(k, set()).update(maps)

    # ---- Pfam
    df["pfam_hits"], df["pfam_accs"], df["duf_only"] = "", "", False
    if os.path.exists(p.pfam):
        pf = parse_hmm_tblout(p.pfam)
        df["pfam_hits"] = from_dict(
            {k: ";".join(sorted({h[0] for h in v})) for k, v in pf.items()}, idx)
        df["pfam_accs"] = from_dict(
            {k: ";".join(sorted({h[1].split(".")[0] for h in v}))
             for k, v in pf.items()}, idx)
        # uninf_re, not DUF_RE: Pfam's "uncharacterised protein family" models
        # are named UPF0xxx and never matched the DUF pattern, so a UPF-only
        # protein was reported as annotated.
        df["duf_only"] = from_dict(
            {k: bool(v) and all(uninf_re.search(h[0]) for h in v)
             for k, v in pf.items()},
            idx, False).astype(bool)
    anchor = set(cfg.get("anchor_pfams") or [])
    df["anchor_domain"] = df["pfam_accs"].map(
        lambda s: bool(anchor & set(s.split(";"))) if s else False)

    # ---- dbCAN
    df["dbcan_hits"] = ""
    if os.path.exists(p.dbcan):
        dc = parse_hmm_domtblout(p.dbcan, th["dbcan_min_cov"], th["dbcan_evalue"])
        df["dbcan_hits"] = from_dict(
            {k: ";".join(sorted({h[0] for h in v})) for k, v in dc.items()}, idx)

    # ---- NCBIfam / TIGRFAM
    # NCBIfam is not uniformly informative: it carries families whose whole
    # DESC is "hypothetical protein" or "DUF1801 domain-containing protein".
    # Counting those as annotation promotes a protein out of 4_dark on a
    # family name that says nothing, which is exactly what the Pfam DUF test
    # exists to prevent. hmmsearch --tblout reports the description of the
    # TARGET (our protein), never of the query HMM, so the family DESC has to
    # come from the library itself.
    df["ncbifam_hits"] = ""
    df["ncbifam_accs"] = ""
    df["ncbifam_uninformative"] = False
    if os.path.exists(p.ncbifam):
        nf = parse_hmm_tblout(p.ncbifam)
        df["ncbifam_hits"] = from_dict(
            {k: ";".join(sorted({h[0] for h in v})) for k, v in nf.items()}, idx)
        # The accession, kept for the same reason pfam_accs is. Without it the
        # only record of an NCBIfam call is its family NAME, and a name cannot
        # be matched against InterProScan's NCBIfam member database, which
        # reports accessions - so two searches of the SAME library looked 97%
        # discordant when they in fact agreed. The tblout already carried it;
        # this column just stops throwing it away. Version suffix stripped, as
        # for Pfam, so NF033709.1 and NF033709 are one family.
        df["ncbifam_accs"] = from_dict(
            {k: ";".join(sorted({h[1].split(".")[0] for h in v
                                 if h[1] and h[1] != "-"}))
             for k, v in nf.items()}, idx)
        if nf and not cfg.get("ncbifam_uninformative_test", True):
            log("ncbifam_uninformative_test is off: every NCBIfam hit counts "
                "as sequence annotation, including 'hypothetical protein' "
                "families", "WARN")
        elif nf:
            fam = ncbifam_family_desc(cfg, p)
            unknown = set()

            def _nf_uninformative(hits):
                seen, all_known = [], True
                for h in hits:
                    d = fam.get(h[0]) or fam.get(h[1]) or ""
                    if not d:
                        # Fail open: a family we cannot describe keeps its
                        # evidence, so a missing library can only leave the
                        # bins as they were, never demote on no information.
                        unknown.add(h[0])
                        all_known = False
                    else:
                        seen.append(d)
                return (all_known and bool(seen)
                        and all(uninf_re.search(d) for d in seen))

            if fam:
                df["ncbifam_uninformative"] = from_dict(
                    {k: _nf_uninformative(v) for k, v in nf.items()},
                    idx, False).astype(bool)
                n_u = int(df["ncbifam_uninformative"].sum())
                log(f"ncbifam: {n_u} protein(s) whose only NCBIfam families "
                    "are uninformative (DUF/hypothetical/unknown function); "
                    "they no longer count as sequence annotation"
                    + (f". {len(unknown)} family name(s) in the tblout have no "
                       "DESC in the library and are counted as informative"
                       if unknown else ""))

    # ---- KOfamScan: a control on eggNOG's KO calls, not just more coverage
    df["kofam_ko"], df["kofam_desc"] = "", ""
    if os.path.exists(p.kofam):
        kf = parse_kofam(p.kofam)
        df["kofam_ko"] = from_dict(
            {k: ";".join(sorted({h[0] for h in v})) for k, v in kf.items()}, idx)
        df["kofam_desc"] = from_dict(
            {k: v[0][3] for k, v in kf.items() if v}, idx)

    egg_ko = df["ko"].str.contains(KO_RE, na=False)
    kof_ko = df["kofam_ko"].str.contains(KO_RE, na=False)
    df["ko_source"] = [
        "both" if a and b else "eggnog" if a else "kofam" if b else ""
        for a, b in zip(egg_ko, kof_ko)]
    # ko_source "both" said nothing about agreement, so a protein eggNOG called
    # K00001 and KOfam called K99999 looked exactly like a confirmed one. The
    # disagreement was only ever a log count; keep it per protein.
    df["ko_conflict"] = [
        bool(e and k and not (set(KO_RE.findall(es or ""))
                              & set(KO_RE.findall(ks or ""))))
        for e, k, es, ks in zip(egg_ko, kof_ko, df["ko"], df["kofam_ko"])]
    if kof_ko.any():
        rescued = int((kof_ko & ~egg_ko).sum())
        disagree = int(pd.Series(df["ko_conflict"]).sum())
        log(f"KOfamScan: {int(kof_ko.sum())} proteins with a KO; {rescued} that "
            f"eggNOG missed; {disagree} where the two disagree entirely")
        if rescued:
            log(f"{rescued} proteins leave the KO-less bins on KOfam evidence "
                "alone. Report this: it is a direct test of whether those bins "
                "are biology or an artefact of eggNOG's DIAMOND search", "WARN")
        if ko2maps:
            def _kofam_maps(s):
                out = set()
                for k in KO_RE.findall(s or ""):
                    out |= ko2maps.get(k, set())
                return len(out)
            extra = df["kofam_ko"].fillna("").map(_kofam_maps)
            gained = int(((df["n_pathway_specific"] == 0) & (extra > 0)).sum())
            df["n_pathway_specific"] = np.maximum(
                df["n_pathway_specific"], extra)
            df["in_specific_pathway"] = df["n_pathway_specific"] > 0
            if gained:
                log(f"{gained} KOfam-only proteins reach a specific KEGG map "
                    "through the run's own KO->map pairs; they now bin as "
                    "1_ko_pathway instead of 2_ko_orphan")

    # ---- InterProScan
    df["interpro_sigs"], df["interpro_ipr"], df["interpro_go"] = "", "", ""
    df["interpro_informative"] = False
    if os.path.exists(p.interpro):
        skip = {x.strip() for x in
                str(th.get("interpro_ignore_analyses", "")).split(",") if x.strip()}
        ip = parse_interproscan(p.interpro)
        # Prefix, not equality: InterProScan reports SignalP as SignalP_EUK /
        # SignalP_GRAM_POSITIVE / SignalP_GRAM_NEGATIVE, none of which equals
        # the configured "SignalP", so signal-peptide-only proteins were being
        # promoted out of the dark bin by the very analysis meant to be ignored.
        def _skipped(a):
            return any(a == s or a.startswith(s) for s in skip)
        keep = {k: [r for r in v if not _skipped(r["analysis"])]
                for k, v in ip.items()}
        df["interpro_sigs"] = from_dict({k: ";".join(sorted({
            f"{r['analysis']}:{r['signature']}" for r in v}))
            for k, v in keep.items()}, idx)
        df["interpro_ipr"] = from_dict({k: ";".join(sorted({
            r["ipr"] for r in v if r["ipr"] and r["ipr"] != "-"}))
            for k, v in keep.items()}, idx)
        df["interpro_go"] = from_dict({k: ";".join(sorted({
            g for r in v for g in str(r["go"]).split("|") if g and g != "-"}))
            for k, v in keep.items()}, idx)
        # About half of a real InterProScan TSV (Gene3D, SUPERFAMILY, SMART,
        # CDD and most PANTHER rows) carries an empty signature description.
        # An empty string matches no DUF pattern, so every such row used to
        # count as informative on its own; fall back to the InterPro
        # description and require a name that actually says something.
        def _ip_informative(r):
            # InterProScan writes "-" for a field it has nothing for, so an
            # empty description and a "-" description mean the same thing.
            for f in ("desc", "ipr_desc"):
                txt = str(r.get(f) or "").strip()
                if txt and txt != "-":
                    return not bool(uninf_re.search(txt))
            return False
        df["interpro_informative"] = from_dict(
            {k: bool(v) and any(_ip_informative(r) for r in v)
             for k, v in keep.items()}, idx, False).astype(bool)

    # ---- profile-profile and iterative profile search
    df["hh_hit"], df["hh_prob"], df["hh_desc"] = "", float("nan"), ""
    if os.path.isdir(p.hhr_dir):
        hh = parse_hhr_dir(p.hhr_dir, th.get("hhblits_min_prob", 90.0))
        if hh:
            df["hh_hit"] = from_dict({k: v[0] for k, v in hh.items()}, idx)
            df["hh_prob"] = from_dict(
                {k: v[1] for k, v in hh.items()}, idx, float("nan"))
            df["hh_desc"] = from_dict({k: v[2] for k, v in hh.items()}, idx)
            log(f"hhblits: {len(hh)} proteins with a profile hit at prob >= "
                f"{th.get('hhblits_min_prob', 90.0)}")
    df["jackhmmer_hit"] = ""
    if os.path.exists(p.jackhmmer):
        # parse_hmm_tblout keys on column 1, which is right for hmmsearch (the
        # HMM is the query and the protein the target) and inverted for
        # jackhmmer, where we search our dark proteins against UniRef50: there
        # column 1 is the UniRef target and column 3 the query. Keyed the
        # hmmsearch way the dict was full of UniRef ids, reindexing onto the
        # protein index produced "" for every protein, and warn_if_no_records
        # stayed quiet because records did exist. Re-key on the query here
        # rather than in the shared parser, which hmmsearch also uses.
        jh_raw = parse_hmm_tblout(p.jackhmmer)
        jh = defaultdict(list)
        for tgt, recs in jh_raw.items():
            for qname, _qacc, ev, sc, _tdesc in recs:
                if qname == tgt:
                    continue      # self-hit, if the DB contains the query
                jh[qname].append((tgt, _qacc, ev, sc, _tdesc))
        df["jackhmmer_hit"] = from_dict(
            {k: sorted(v, key=lambda h: h[2])[0][0] for k, v in jh.items() if v},
            idx)
        n_jh = int(df["jackhmmer_hit"].ne("").sum())
        log(f"jackhmmer: {n_jh} proteins with a homolog "
            f"(from {len(jh)} queries in the tblout)")
        if jh and n_jh == 0:
            log("jackhmmer produced hits but none of the query ids match the "
                "FASTA — check that dark_all.faa and the tblout come from "
                "the same run", "WARN")

    # ---- fold groups from self-clustering the unannotated structures
    df["fold_cluster"], df["fold_cluster_size"] = "", 0
    if os.path.exists(p.fold_clusters):
        fc = pd.read_csv(p.fold_clusters, sep="\t", encoding="utf-8", encoding_errors="replace")
        m = dict(zip(fc["member"].astype(str), fc["rep"].astype(str)))
        sz = fc.groupby("rep").size()
        df["fold_cluster"] = from_dict(m, idx)
        df["fold_cluster_size"] = df["fold_cluster"].map(sz).fillna(0).astype(int)

    # ---- targeted databases (globbed, so config changes need no code change)
    dia_weights = cfg.get("diamond_weights") or {}
    dia_tags = []
    for path in sorted(glob.glob(f"{p.diamond_dir}/*.tsv")):
        tag = os.path.splitext(os.path.basename(path))[0]
        dia_tags.append(tag)
        if tag not in dia_weights:
            log(f"no diamond_weights entry for '{tag}'; hits count as annotation "
                "but contribute 0 to the score", "WARN")
        # The same per-database e-value the search used: filtering the table
        # back down to thresholds.diamond_evalue here would quietly undo a
        # diamond_evalues entry and leave the user's setting doing nothing.
        # The same per-database identity floor the search used, for the same
        # reason as the e-value above: reading the table back at the global
        # thresholds.diamond_min_pident would quietly undo a
        # diamond_min_pidents entry. It also has to be applied HERE and not
        # only on the command line, because a <tag>.tsv adopted from another
        # machine, or written before the floor existed, never saw --id.
        min_pid = diamond_min_pident_for(cfg, tag)
        hits = parse_diamond(path, diamond_evalue_for(cfg, tag),
                             th["diamond_min_qcov"], min_pid)
        df[f"{tag}_hit"] = from_dict({k: v[0] for k, v in hits.items()}, idx)
        df[f"{tag}_pident"] = from_dict(
            {k: v[1] for k, v in hits.items()}, idx, float("nan"))
        df[f"{tag}_desc"] = from_dict({k: v[3] for k, v in hits.items()}, idx)

    # ---- topology
    df["sp_class"], df["sp_cs"] = "", ""
    if os.path.exists(p.signalp):
        sp = parse_signalp6(p.signalp)
        df["sp_class"] = from_dict({k: v[0] for k, v in sp.items()}, idx)
        df["sp_cs"] = from_dict({k: v[1] for k, v in sp.items()}, idx)
    df["n_tmh"], df["n_tmb"] = 0, 0
    if os.path.exists(p.tmbed):
        tb = parse_tmbed(p.tmbed)
        df["n_tmh"] = from_dict({k: v[0] for k, v in tb.items()}, idx, 0).astype(int)
        df["n_tmb"] = from_dict({k: v[1] for k, v in tb.items()}, idx, 0).astype(int)

    df["small_protein"] = df["length"] <= th["smorf_max_len"]

    # ---- de novo families
    df["family_id"] = list(idx)
    if os.path.exists(p.cluster):
        cl = parse_cluster(p.cluster)
        df["family_id"] = from_dict(cl, idx).where(lambda x: x != "",
                                                   pd.Series(idx, index=idx))

    # ---- structure
    df["foldseek_target"], df["foldseek_desc"] = "", ""
    df["foldseek_db"] = ""
    df["foldseek_prob"], df["foldseek_tm"] = float("nan"), float("nan")
    # Which target database wins when a protein hits several. Empty config
    # means the order the user listed their targets in: foldseek_target is
    # the one they called primary, extras follow in their own order. Matched
    # against the basenames stage_foldseek writes into hits.tsv.
    fs_prio = list(cfg.get("foldseek_target_priority") or [])
    if not fs_prio:
        _tgts = [(cfg.get("db") or {}).get("foldseek_target") or ""] + list(
            (cfg.get("db") or {}).get("foldseek_extra_targets") or [])
        fs_prio = [os.path.basename(str(t).rstrip("/\\")) for t in _tgts if t]
    fs = parse_foldseek(p.foldseek, th["foldseek_evalue"],
                        th["foldseek_min_prob"], th["foldseek_min_tmscore"],
                        target_priority=fs_prio)
    if fs:
        df["foldseek_target"] = from_dict({k: v[0] for k, v in fs.items()}, idx)
        df["foldseek_prob"] = from_dict(
            {k: v[1] for k, v in fs.items()}, idx, float("nan"))
        df["foldseek_tm"] = from_dict(
            {k: v[2] for k, v in fs.items()}, idx, float("nan"))
        df["foldseek_desc"] = from_dict({k: v[3] for k, v in fs.items()}, idx)
        # Provenance, so a reader can see whether "3s_structure_only" rests on
        # a described Swiss-Prot entry or an unannotated AFDB50 model.
        df["foldseek_db"] = from_dict({k: v[4] for k, v in fs.items()}, idx)

    pats = cfg.get("toxin_fold_patterns") or []
    # Word boundaries: unanchored, "RTX" matched inside "MRTXase" and "Rhs"
    # inside longer names, each adding the largest single effector weight.
    tox_re = re.compile(r"\b(?:" + "|".join(pats) + r")\b", re.I) if pats \
        else None
    # Read hh_desc as well as foldseek_desc. It is free - both columns are
    # already built - and it lifts the ceiling above the handful of proteins
    # that got a Foldseek hit at all.
    def _tox(col):
        if col not in df.columns:
            return pd.Series(False, index=df.index)
        return df[col].fillna("").astype(str).map(
            lambda s: bool(s) and bool(tox_re.search(s)))
    df["toxin_fold"] = (_tox("foldseek_desc") | _tox("hh_desc")) if tox_re         else False
    # AFDB50 entry headers are AlphaFold accessions (AF-<UniProt>-F1-model_v4)
    # with no protein name, so on the default target the toxin term can never
    # fire and a reader takes "no toxin folds" for biology. Say so once.
    if tox_re and df["foldseek_target"].ne("").any():
        # A header carrying a protein name has whitespace after the accession;
        # a bare AFDB accession has none.
        _named = df["foldseek_desc"].fillna("").str.strip().str.contains(
            r"\s", regex=True, na=False)
        if not bool(_named.any()):
            log("Foldseek target headers carry no description (AFDB50-style "
                "accessions), so toxin_fold cannot fire; add a PDB or "
                "Swiss-Prot target if you want that term", "WARN")

    # Membership of dark.faa is a request, not an outcome: ESMFold skips a
    # protein when both OOM retries fail, and an interrupted overnight run
    # leaves most of dark.faa unfolded. Reporting the request as "attempted"
    # made "we folded it and found nothing" indistinguishable from "we never
    # folded it" — the one distinction the bins exist to keep honest.
    df["structure_requested"] = False
    req_len = {}
    if os.path.exists(p.dark):
        # Lengths, not just ids: max_len_structure excludes a long protein
        # from folding before ESMFold ever sees it, so it is one of the
        # reasons a model can be missing without anything having gone wrong.
        req_len = {pid: len(seq) for pid, seq in read_fasta(p.dark)}
        df["structure_requested"] = idx.isin(set(req_len))
    have_pdb = {os.path.basename(f)[:-4]
                for f in glob.glob(f"{p.structures}/*.pdb")}
    df["structure_attempted"] = idx.isin(have_pdb)
    # Requested AND in this table: p.structures can hold models for proteins
    # that are no longer in dark.faa, and counting those made the shortfall
    # look smaller than it was (len(have_pdb) could even exceed the number
    # requested, silencing the message entirely).
    req_len = {pid: n for pid, n in req_len.items()
               if pid in set(idx[df["structure_requested"]])}
    if req_len and not set(req_len) <= have_pdb:
        log(*structure_shortfall_message(cfg, p, req_len, have_pdb))

    # ---- genomic context
    ctx = parse_context(p.context)

    # `v == 1` was also true for the integer column n_cazymes_in_window, so
    # the literal string "n_cazymes_in_window" appeared in context_flags for
    # every protein with exactly one CAZyme neighbour (but not 0, 2 or 3) and
    # was printed in the report's shortlist as if it were a detected feature.
    # Only genuine booleans belong in a flag string.
    def _is_flag(v):
        if isinstance(v, str):
            return v.strip().lower() == "true"
        return isinstance(v, (bool, np.bool_)) and bool(v)

    df["context_flags"] = from_dict(
        {pid: ";".join(k for k, v in d.items() if _is_flag(v))
         for pid, d in ctx.items()}, idx)
    df["context_mge"] = df["context_flags"].str.contains(
        "prophage|bgc|t3ss|t6ss|secretion", case=False, na=False)
    df["context_pul"] = df["context_flags"].str.contains("pul", case=False,
                                                         na=False)
    df["context_immunity"] = df["context_flags"].str.contains(
        "immunity_pair", case=False, na=False)

    # ---- bins
    # KOfam KOs count towards has_ko, but ko_source records the provenance so
    # the rescue is quantifiable rather than silently absorbed.
    if "kofam_ko" in df.columns:
        df["has_ko"] = df["has_ko"] | df["kofam_ko"].str.contains(KO_RE, na=False)

    # eggNOG's own PFAMs column is domain evidence and has to count. Binning
    # on hmmsearch's pfam_hits alone meant that with the pfam stage off, a
    # protein eggNOG had already assigned a domain to was reported as "no
    # evidence": on the UC metaproteome 6,271 of 17,377 dark proteins (36%)
    # carried an eggNOG Pfam, so the dark bin was inflated by a missing join
    # rather than by biology.
    emapper_pfam = df["pfams_emapper"].fillna("") if "pfams_emapper" in df.columns \
        else pd.Series("", index=idx)
    emapper_pfam_informative = emapper_pfam.ne("") & ~emapper_pfam.str.contains(
        uninf_re, na=False)
    # ncbifam_uninformative is the NCBIfam analogue of duf_only: the family
    # is named, its function is not, so the hit is not annotation.
    informative = ((df["pfam_hits"].ne("") & ~df["duf_only"])
                   | (df["ncbifam_hits"].ne("")
                      & ~df["ncbifam_uninformative"].fillna(False).astype(bool))
                   | emapper_pfam_informative
                   | df["interpro_informative"].astype(bool))
    targeted = pd.Series(False, index=idx)
    for t in dia_tags:
        targeted = targeted | df[f"{t}_hit"].ne("")
    df["has_seq_annotation"] = (informative | df["cazy"].ne("")
                                | df["dbcan_hits"].ne("") | targeted)
    # A protein whose only domain evidence is an eggNOG DUF belongs in the
    # DUF bin, not the dark one, for the same reason as an hmmsearch DUF. An
    # uninformative-only NCBIfam hit is the same case: the family is known,
    # the function is not, so it lands in 3d_duf_only rather than 4_dark.
    df["duf_only"] = (df["duf_only"].fillna(False).astype(bool)
                      | (~df["has_seq_annotation"].fillna(False).astype(bool)
                         & (emapper_pfam.ne("")
                            | df["ncbifam_uninformative"].fillna(False)
                            .astype(bool))))
    # "A homolog exists" is not annotation. The default hhblits DB is Pfam, so
    # a dark protein's best profile hit is often a DUF, and AFDB50 is mostly
    # TrEMBL "Uncharacterized protein" entries — both used to leave 4_dark and
    # be counted as a structure/profile rescue. Demote only when the hit
    # actually carries an uninformative description: an empty description
    # (an AFDB accession, or jackhmmer, which records none) keeps its evidence,
    # so this can only shrink the rescue claim, never inflate it.
    def _uninformative(s):
        s = (s or "").strip()
        return bool(s) and bool(uninf_re.search(s))

    fs_uninf = df["foldseek_desc"].fillna("").map(_uninformative)
    hh_uninf = df["hh_desc"].fillna("").map(_uninformative)
    df["structure_evidence"] = df["foldseek_target"].fillna("").ne("") & ~fs_uninf
    df["has_profile_hit"] = ((df["hh_hit"].ne("") & ~hh_uninf)
                             | df["jackhmmer_hit"].ne(""))
    n_demoted = int((fs_uninf & df["foldseek_target"].ne("")).sum()
                    + (hh_uninf & df["hh_hit"].ne("")).sum())
    if n_demoted:
        log(f"{n_demoted} structure/profile hits are to uncharacterised "
            "targets (DUF/hypothetical) and do not count as a rescue", "WARN")

    # np.select, not apply(axis=1): the row-wise apply built a Series per
    # protein and was the single largest cost in this function.
    hk = df["has_ko"].fillna(False).astype(bool)
    df["bin"] = np.select(
        [hk & df["in_specific_pathway"].fillna(False).astype(bool),
         hk,
         df["has_seq_annotation"].fillna(False).astype(bool),
         df["duf_only"].fillna(False).astype(bool),
         df["structure_evidence"].fillna(False).astype(bool),
         df["has_profile_hit"].fillna(False).astype(bool)],
        # The condition list above is in BIN_ORDER, so the bin names come
        # straight from it; 4_dark is the default, i.e. the last entry.
        list(BIN_ORDER[:-1]),
        default=BIN_ORDER[-1])
    df["kegg_enrichment_visible"] = df["bin"] == "1_ko_pathway"

    # duf_only is tested before the structural and profile conditions, so no
    # amount of fold or profile evidence can move a DUF-only protein out of
    # 3d_duf_only — deliberate (a DUF still names the family), but it means
    # those GPU slots buy a column, not a bin change. Count them so the
    # bin_summary rescue numbers are read with that in mind.
    _duf_rescued = int((df["bin"].eq("3d_duf_only")
                        & (df["structure_evidence"].fillna(False)
                           | df["has_profile_hit"].fillna(False))).sum())
    if _duf_rescued:
        log(f"{_duf_rescued} 3d_duf_only proteins have structure or profile "
            "evidence; the DUF keeps them in that bin, so they are not "
            "counted in the 3s/3p rescue")

    # The bin is the FIRST evidence that applies, so it hides every later one:
    # a 3d_duf_only or 2_ko_orphan protein with a good fold looks, in the bin
    # column alone, exactly like one with none. rescued_by records the
    # structure/profile evidence a protein actually has, for every bin, so the
    # rescue can be counted without adding a bin and rippling that new level
    # through the report, BIN_LEVELS, BIN_COLS and dark.faa's bin list.
    # For 3s_structure_only / 3p_profile_only it restates the bin; that is
    # deliberate, so "has structure evidence" is one test, not two.
    _struct = df["structure_evidence"].fillna(False).astype(bool)
    _prof = df["has_profile_hit"].fillna(False).astype(bool)
    df["rescued_by"] = np.select(
        [_struct & _prof, _struct, _prof],
        ["structure;profile", "structure", "profile"], default="")
    _hidden = int(((_struct | _prof)
                   & ~df["bin"].isin(["3s_structure_only",
                                      "3p_profile_only"])).sum())
    if _hidden:
        log(f"{_hidden} protein(s) carry structure or profile evidence that "
            "their bin does not name; the rescued_by column of the "
            "annotation table lists it")

    # ---- effector priority score
    def B(col_, default=False):
        if col_ not in df.columns:
            return pd.Series(default, index=idx)
        return df[col_].fillna(False).astype(bool)

    sp_c = df["sp_class"].fillna("")
    score = pd.Series(0, index=idx, dtype="int64")
    score += (sp_c == "SP") * w["signal_sec_spi"]
    score += (sp_c == "LIPO") * w["signal_lipo_spii"]
    score += sp_c.isin(["TAT", "TATLIPO"]) * w["signal_tat"]
    # SignalP 6 has six classes; PILIN (Sec/SPIII) was parsed and then ignored
    # everywhere, so type IV pilins and competence proteins — bona fide surface
    # proteins, often smORF-sized — scored below an OTHER-class protein and
    # were invisible to the shortlist gate.
    score += (sp_c == "PILIN") * w.get("signal_pilin_spiii", 2)
    score += B("lpxtg") * w["lpxtg_anchor"]
    score += B("anchor_domain") * w["slh_or_anchor_domain"]
    score += (df["n_tmb"] > 0) * w["tm_beta_barrel"]
    score += ((df["n_tmh"] >= 2) & sp_c.isin(["", "OTHER"])) * w["multi_tm_helix_penalty"]
    # A 31%-identity VFDB hit is not the same evidence as a 90% one, and both
    # used to score the full weight. Below diamond_strong_pident the weight is
    # halved rather than dropped, so a weak hit still ranks above no hit.
    strong_pid = th.get("diamond_strong_pident", 50)
    # A database floored at or above diamond_strong_pident cannot produce a
    # weak hit: everything that survived the filter scores full weight and the
    # halving below is dead for it. Said once, because the alternative is a
    # reader concluding from the config that VFDB hits are being graded by
    # identity when in fact they cannot be.
    floored = [t for t in dia_tags
               if diamond_min_pident_for(cfg, t) >= strong_pid]
    if floored:
        log(f"{', '.join(sorted(floored))}: filtered at or above "
            f"diamond_strong_pident={strong_pid:g}, so every surviving hit "
            "scores the full diamond_weights value and the half-weight rule "
            "for weak hits never applies to them")
    # VFDB is not one kind of evidence. Its own VFC category code says which,
    # and a flat weight throws that away: on a real gut metaproteome, of 3,308
    # VFDB hits the two largest categories were "Immune modulation" (965) and
    # "Nutritional/Metabolic factor" (903, i.e. GroEL, ClpP, GuaA, LPS
    # biosynthesis) - each collecting the largest DIAMOND weight in the config
    # - while the categories that actually name an exported effector,
    # VFC0086 (effector delivery, 232) and VFC0235 (exotoxin, 132), were 11% of
    # the signal. Keyed on the NUMERIC code, not the prose: VFDB can reword a
    # category name, it will not renumber it.
    cat_w = cfg.get("vfdb_category_weights") or {}
    for tag in dia_tags:
        hit_ = df[f"{tag}_hit"].fillna("").ne("")
        pid_ = pd.to_numeric(df[f"{tag}_pident"], errors="coerce").fillna(0)
        wt = dia_weights.get(tag, 0)
        by_category = tag == "vfdb" and cat_w and f"{tag}_desc" in df.columns
        if by_category and not wt:
            # The category map SPLITS diamond_weights.vfdb; it is not a second
            # way in. A database whose weight is 0 or missing contributes
            # nothing - that is what the "no diamond_weights entry" WARN,
            # doctor and the README all promise - and applying the map anyway
            # handed 1-4 points per hit back to a user who had deliberately
            # zeroed VFDB, so the only way to stop VFDB scoring was to empty
            # the category map as well. Said out loud, because a map that is
            # ignored is exactly the kind of thing that otherwise looks like
            # the weighting simply not working.
            log("vfdb: diamond_weights.vfdb is "
                + ("0" if "vfdb" in dia_weights else "absent")
                + ", so VFDB hits score 0 and vfdb_category_weights is not "
                "applied; give vfdb a non-zero weight to weight its "
                "categories", "WARN")
            by_category = False
        if by_category:
            code = df[f"{tag}_desc"].fillna("").astype(str).str.extract(
                r"\((VFC\d+)\)", expand=False)
            n_code = int(code.notna().sum())
            n_hit = int(hit_.sum())
            if n_hit:
                # Said once, so a VFDB format change surfaces as a line rather
                # than as every hit silently taking the fallback weight.
                log(f"vfdb: {n_code}/{n_hit} hit(s) carry a VFC category code; "
                    "those are weighted by category, the rest by "
                    f"diamond_weights.vfdb={wt}")
            wt_s = code.map(lambda c: cat_w.get(c, wt)).fillna(wt).astype(int)
        else:
            wt_s = wt
        score += (hit_ & (pid_ >= strong_pid)) * wt_s
        score += (hit_ & (pid_ < strong_pid)) * (wt_s // 2)
    score += (df["cazy"].fillna("").ne("") | df["dbcan_hits"].fillna("").ne("")) * w["cazy_hit"]
    score += B("small_protein") * w["small_protein"]
    score += B("toxin_fold") * w["foldseek_toxin_fold"]
    score += B("context_mge") * w["context_mge_or_secretion"]
    score += B("context_pul") * w.get("context_pul", 1)
    score += B("context_immunity") * w["context_immunity_pair"]
    score += (~hk) * w["no_ko"]
    # export_score, not effector_score. Every term above asks whether a protein
    # LEAVES THE CELL or sits on its surface - signal peptide class, beta
    # barrel, LPXTG or SLH anchor, CAZy, small size, a toxin-like fold, a
    # mobile-element or secretion neighbourhood, no KO. None of them asks
    # whether it is an effector of a secretion system, and after the effectors
    # stage was removed nothing in the tool does. Naming it effector_score
    # promised a claim the evidence never supported, which on a gut commensal
    # metaproteome is exactly the claim a reviewer would reject.
    df["export_score"] = score
    # The old name, carried one release so existing scripts, notebooks and R
    # code keep working. Same numbers, not a second opinion. Drop it after the
    # next release; the column comment above says which one to prefer.
    df["effector_score"] = score
    df["surface_or_secreted"] = (
        df["sp_class"].isin(["SP", "LIPO", "TAT", "TATLIPO", "PILIN"])
        | df["lpxtg"] | df["anchor_domain"] | (df["n_tmb"] > 0))

    df = df.sort_values(["export_score", "length"], ascending=[False, True])

    if emit_dark:
        # Dark AND DUF-only: a DUF names a family, not a function, so those
        # gain as much from a fold as the fully dark ones.
        sel = df[df["bin"].isin(["4_dark", "3d_duf_only", "3p_profile_only"])]
        # Entrapment, decoy, contaminant and host sequences are in the search
        # database on purpose. Nothing annotates them, so they land in 4_dark
        # with a high effector_score (no_ko, small_protein) and float to the
        # top of this work list — the false-discovery control would eat the
        # structure budget and reappear in the report as a novel dark protein.
        prefixes = tuple(cfg.get("exclude_id_prefixes") or ())
        if prefixes:
            drop = np.asarray(
                [str(s).startswith(prefixes) for s in sel.index], dtype=bool)
            if drop.any():
                hits = {}
                for s in sel.index[drop]:
                    for pre in prefixes:
                        if str(s).startswith(pre):
                            hits[pre] = hits.get(pre, 0) + 1
                            break
                log(f"{int(drop.sum())} unannotated proteins excluded from "
                    f"{emit_dark} by exclude_id_prefixes {hits}. Entrapment or "
                    "contaminant sequences among the identifications are a "
                    "finding in their own right, not fold candidates", "WARN")
                sel = sel[~drop]
        # The uncapped set first: this is what "unannotated" means, and what
        # the profile searches query. Neither cap below belongs to it —
        # max_len_structure is an ESMFold memory limit and
        # max_dark_structures is a GPU-hours budget.
        n_all = len(sel)
        want_all = set(sel.index) if emit_dark_all else set()
        cap = sel[sel["length"] <= cfg["max_len_structure"]]
        n_elig = len(cap)
        cap = cap.head(cfg["max_dark_structures"])
        want = set(cap.index)
        # One pass over the protein FASTA feeding both files; on a full ORF
        # catalogue a second pass is minutes of pointless I/O.
        with contextlib.ExitStack() as st:
            t_dark = st.enter_context(atomic_out(emit_dark))
            fh = st.enter_context(open(t_dark, "w", encoding="utf-8"))
            fh_all = None
            if emit_dark_all:
                t_all = st.enter_context(atomic_out(emit_dark_all))
                fh_all = st.enter_context(open(t_all, "w", encoding="utf-8"))
            for pid, sq in read_fasta(cfg["proteins_faa"]):
                if pid in want:
                    fh.write(f">{pid}\n{sq}\n")
                    want.discard(pid)
                if fh_all is not None and pid in want_all:
                    fh_all.write(f">{pid}\n{sq}\n")
                    want_all.discard(pid)
        if emit_dark_all:
            log(f"{n_all} unannotated proteins -> {emit_dark_all} "
                "(uncapped; this is the query set for hhblits and jackhmmer)")
        log(f"{len(cap)}/{n_elig} unannotated proteins -> {emit_dark} "
            "(structure work-list: highest effector_score first, ties "
            "shortest first, length <= "
            f"max_len_structure={cfg['max_len_structure']})")
        if n_elig > len(cap):
            log(f"capped at max_dark_structures={cfg['max_dark_structures']}: "
                f"{n_elig - len(cap)} eligible proteins are not folded. The "
                "cap bounds STRUCTURES only — hhblits and jackhmmer query "
                f"{os.path.basename(emit_dark_all) if emit_dark_all else 'dark_all.faa'}"
                f", all {n_all} of them", "WARN")

    return df


# Evidence a protein can carry, as (label, column). A column absent from the
# frame is skipped rather than reported as 0%: "this stage did not run" and
# "this stage found nothing" are different answers and must not share a cell.
TIER_EVIDENCE = [
    ("eggnog", "og"), ("ko", "ko"), ("pfam", "pfam_hits"),
    ("ncbifam", "ncbifam_hits"), ("kofam", "kofam_ko"),
    ("interpro", "interpro_sigs"), ("dbcan", "dbcan_hits"),
    ("cazy", "cazy"),
]


def write_tier_coverage(df, path, cfg=None):
    """Coverage split by identifier prefix, for a merged search database.

    bin_summary.tsv answers "what did this proteome look like"; this answers
    "and did its parts look alike", which for a database merged from several
    catalogues plus this study's own assembly is the question the headline
    number hides. Writes nothing and says why when the ids are not tiered, so
    an absent file is never ambiguous.

    exclude_id_prefixes is applied FIRST. A run whose proteins_faa is the
    whole search database rather than the identified subset carries the
    decoys, the entrapment set and the contaminants, and on the database this
    was written for those are the two LARGEST namespaces in the file --
    18,318,713 rev_ and 2,280,823 ent_ against 11,379,230 uhgpL_. A tier table
    whose top row is the decoy set is not a description of the biology, and it
    would also spend the twelve-tier budget on namespaces that are there to be
    ignored. They are dropped and counted out loud, never silently.
    """
    if cfg is not None:
        prefixes = tuple(cfg.get("exclude_id_prefixes") or ())
        if prefixes:
            drop = np.asarray([str(s).startswith(prefixes) for s in df.index],
                              dtype=bool)
            if drop.any():
                hits = {}
                for s in df.index[drop]:
                    for pre in prefixes:
                        if str(s).startswith(pre):
                            hits[pre] = hits.get(pre, 0) + 1
                            break
                log(f"tier coverage: {int(drop.sum()):,} protein(s) excluded "
                    f"by exclude_id_prefixes {hits} before the split. Decoy, "
                    "entrapment and contaminant namespaces are in the search "
                    "database to be ignored, and reporting their coverage "
                    "beside a real catalogue's would invite reading them as "
                    "one")
                df = df[~drop]
    if not len(df):
        log("tier coverage: nothing is left after exclude_id_prefixes, so no "
            f"per-tier table is written; {os.path.basename(path)} is absent "
            "for that reason, not because a stage failed", "WARN")
        return False
    tiers = id_tiers(df.index)
    if not tiers:
        log("protein ids are not split by a source prefix (or carry more "
            f"than {MAX_ID_TIERS} distinct ones), so no per-tier coverage "
            f"table is written; {os.path.basename(path)} is absent for that "
            "reason, not because a stage failed")
        return False
    m = pd.Series(list(df.index), index=df.index).str.extract(
        r"^([^_|:.]*[_|:.])", expand=False).fillna("")
    rows = []
    for pref, _n in sorted(tiers.items(), key=lambda kv: -kv[1]):
        sub = df[m == pref]
        row = {"tier": pref or "(no prefix)", "n": len(sub),
               "pct_of_proteome": round(100.0 * len(sub) / len(df), 1)}
        for label, col in TIER_EVIDENCE:
            if col in sub.columns:
                row[f"pct_{label}"] = round(
                    100.0 * sub[col].fillna("").astype(str).ne("").mean(), 1)
        row["pct_dark"] = round(100.0 * sub["bin"].eq("4_dark").mean(), 1)
        row["median_export_score"] = sub["export_score"].median()
        # The KEY under the tag, not the tag. Reported per tier because two
        # tiers with the same key are one namespace under two labels.
        shapes = pd.Series([id_key_shape(q) for q in sub.index]).value_counts()
        row["key_shape"] = shapes.index[0] if len(shapes) else ""
        row["key_shape_pct"] = round(
            100.0 * shapes.iloc[0] / len(sub), 1) if len(shapes) else 0.0
        rows.append(row)
    t = pd.DataFrame(rows).set_index("tier")
    with atomic_out(path) as tmp:
        t.to_csv(tmp, sep="\t")
    log(f"{len(t)} identifier tier(s) in the protein set; coverage per tier "
        f"-> {path}")
    shared = {}
    for tier, shape in t["key_shape"].items():
        shared.setdefault(shape, []).append(tier)
    for shape, tiers in shared.items():
        if len(tiers) > 1:
            log(f"identifier key {shape} is shared by {', '.join(tiers)}: "
                "these are one namespace under several tags, so the same "
                "protein can appear once per tag, one row of a precomputed "
                "annotation table annotates all of them, and EVERY one of "
                "those tags has to be in emapper_strip_id_prefix or its tier "
                "loses that table entirely", "WARN")
    for line in t.to_string().splitlines():
        log(line)
    return True


def write_summary(df, path):
    g = df.groupby("bin")
    s = pd.DataFrame({
        "n": g.size(),
        "pct": (g.size() / len(df) * 100).round(1),
        "n_secreted_or_surface": g["surface_or_secreted"].sum(),
        "n_small": g["small_protein"].sum(),
        "median_export_score": g["export_score"].median(),
    })
    # Every bin gets a row, in the order the classifier assigns them.
    # groupby drops a bin with no proteins, so an empty 3s_structure_only was
    # indistinguishable from a Foldseek stage that never ran, and the rows came
    # out alphabetically (3p before 3s), which is not the order they are
    # decided in.
    absent = [b for b in BIN_ORDER if b not in s.index]
    unknown = [b for b in s.index if b not in BIN_ORDER]
    if unknown:
        # A bin name that BIN_ORDER does not know would be dropped by the
        # reindex, and the report would never show those proteins at all.
        die(f"bin(s) {unknown} are not in BIN_ORDER {list(BIN_ORDER)}; the "
            "classifier and the bin vocabulary have drifted apart.")
    s = s.reindex(list(BIN_ORDER))
    s.index.name = "bin"
    for c in ("n", "pct", "n_secreted_or_surface", "n_small"):
        s[c] = s[c].fillna(0)
    if absent:
        log("no proteins in " + ", ".join(absent)
            + "; reported as zero rows so an empty bin is not mistaken for a "
              "stage that never ran")
    s.loc["TOTAL"] = [len(df), 100.0, int(df["surface_or_secreted"].sum()),
                      int(df["small_protein"].sum()), df["export_score"].median()]
    # Assigning the TOTAL row upcast the count columns to float, so the table
    # the tutorial tells people to read reported "92.0 proteins". pct and the
    # median stay float on purpose.
    for c in ("n", "n_secreted_or_surface", "n_small"):
        s[c] = s[c].astype("int64")
    with atomic_out(path) as tmp:
        s.to_csv(tmp, sep="\t")
    # log(), not a bare stderr write: this is the Phase 4c check the tutorial
    # says to watch for with `tail -f metaannot.log`, and it never got there.
    for line in s.to_string().splitlines():
        log(line)
    log(f"{(~df['kegg_enrichment_visible']).mean()*100:.1f}% of proteins are "
        "invisible to KEGG pathway enrichment")


def stage_integrate_pass1(cfg, p):
    # dark_all.faa is a new output, so a results directory written before the
    # work-list was split has dark.faa and nothing else, exists_all() fails
    # and this stage reruns. Say why, or the rerun looks like cache damage.
    if os.path.exists(p.dark) and not os.path.exists(p.dark_all):
        log(f"{p.dark_all} is absent: this results directory predates the "
            "split of the structure work-list from the profile-search query "
            "set, so integrate runs once more to write it. No threshold, "
            "weight or bin changes.", "WARN")
    df = build_annotation(cfg, p, emit_dark=p.dark, emit_dark_all=p.dark_all)
    with atomic_out(p.pass1) as tmp:
        df.to_csv(tmp, sep="\t")
    log(f"wrote {p.pass1}")


# Outputs that only exist after the first integrate pass. If none of them are
# present, the second pass would recompute a table identical to the first.
def post_integrate_evidence(p):
    """Evidence that can only appear after the first integrate pass.

    Tests the evidence itself, never the stage sentinels: hhr_done and the
    other .done markers are zero-byte files by design, so a nonempty() check
    on them reports "no evidence" exactly when the stage has just produced
    some, and the second pass would silently drop it.
    """
    return (nonempty(p.foldseek)
            or nonempty(p.fold_clusters)
            or nonempty(p.jackhmmer)
            or bool(glob.glob(f"{p.hhr_dir}/*.hhr")))


# Every column build_annotation writes on every run, in the order it writes
# them. The per-DIAMOND-database columns (<tag>_hit, <tag>_pident, <tag>_desc)
# are deliberately absent: which of those exist depends on which databases the
# config names, so their absence is a difference of configuration rather than
# of version. ANN_CORE_COLS is what stage_integrate_final checks an
# annotation_pass1.tsv against before adopting one it did not write itself;
# keep it in step with build_annotation, which
# test_ann_core_cols_lists_every_column_build_annotation_writes enforces.
ANN_CORE_COLS = (
    "protein_id", "length", "lpxtg", "ko", "cog_cat", "description",
    "preferred_name", "ec", "pfams_emapper", "cazy", "og", "seed_taxid",
    "n_pathway_specific", "has_ko", "in_specific_pathway", "pfam_hits",
    "pfam_accs", "duf_only", "anchor_domain", "dbcan_hits", "ncbifam_hits",
    "ncbifam_accs", "ncbifam_uninformative", "kofam_ko", "kofam_desc",
    "ko_source", "ko_conflict", "interpro_sigs", "interpro_ipr",
    "interpro_go", "interpro_informative", "hh_hit", "hh_prob", "hh_desc",
    "jackhmmer_hit", "fold_cluster", "fold_cluster_size", "sp_class", "sp_cs",
    "n_tmh", "n_tmb", "small_protein", "family_id", "foldseek_target",
    "foldseek_desc", "foldseek_db", "foldseek_prob", "foldseek_tm",
    "toxin_fold", "structure_requested", "structure_attempted",
    "context_flags", "context_mge", "context_pul", "context_immunity",
    "has_seq_annotation", "structure_evidence", "has_profile_hit", "bin",
    "kegg_enrichment_visible", "rescued_by", "export_score", "effector_score",
    "surface_or_secreted",
)


def stage_integrate_final(cfg, p):
    have_extra = post_integrate_evidence(p)
    if not have_extra and os.path.exists(p.pass1):
        # Nothing has changed since pass one, so re-deriving every evidence
        # column would burn the same seconds for the same answer. Only
        # structure_attempted can differ, because pass one writes dark.faa
        # after it reads it, and that column feeds nothing else.
        _head = pd.read_csv(p.pass1, sep="\t", nrows=0, encoding="utf-8", encoding_errors="replace")
        # The pass1 on disk was not necessarily written by THIS metaannot.
        # `--only finalise` over a v0.2 results directory landed here with a
        # pass1 carrying effector_score and no export_score, wrote it straight
        # back out as annotation_final.tsv, and only then did write_summary
        # ask for the median export_score and raise KeyError('export_score').
        # What survived was a half-applied upgrade: a v0.2-shaped
        # annotation_final.tsv with no export_score at all, a stale
        # bin_summary.tsv, no source_agreement.tsv, and an R report whose
        # col_or_na filled export_score with NA and ordered the shortlist by
        # it. The check is over the whole column set rather than export_score
        # alone because v0.2 also lacks ncbifam_accs, whose absence merely
        # drops a comparison from source_agreement.tsv and says nothing.
        # Refuse before a single byte is written.
        _missing = [c for c in ANN_CORE_COLS if c not in _head.columns]
        if _missing:
            die(f"{p.pass1} is missing {_missing[:8]}"
                + (f" and {len(_missing) - 8} more column(s)"
                   if len(_missing) > 8 else "")
                + ". It was written by an older metaannot whose "
                "build_annotation produced a different set of columns, so "
                "finalise cannot reuse it and nothing has been written. "
                "Rebuild it with this version first: `--from integrate` does "
                "integrate and finalise in one go, or `--force --only "
                "integrate` then finalise.")
        log("no structure or profile evidence, reusing the first pass")
        # Every identifier-like column must be read as str. Pinning only
        # protein_id let a numeric seed_taxid round-trip through float64 and
        # be written back as "821.0", which matches nothing downstream.
        _str_cols = {c: str for c in _head.columns
                     if c in set(ANN_STR_COLS) or c == "protein_id"}
        df = pd.read_csv(p.pass1, sep="\t", low_memory=False,
                         dtype=_str_cols, encoding="utf-8", encoding_errors="replace").set_index("protein_id")
        if os.path.exists(p.dark):
            tried = {pid for pid, _ in read_fasta(p.dark)}
            df["structure_requested"] = df.index.isin(tried)
        # structure_attempted means a PDB exists, not that we asked for one.
        have_pdb = {os.path.basename(f)[:-4]
                    for f in glob.glob(f"{p.structures}/*.pdb")}
        df["structure_attempted"] = df.index.isin(have_pdb)
    else:
        df = build_annotation(cfg, p)
    with atomic_out(p.final) as tmp:
        df.to_csv(tmp, sep="\t")
    write_summary(df, p.summary)
    write_tier_coverage(df, p.tier_coverage, cfg)
    write_source_agreement(df, p.agreement)
    log(f"wrote {p.final}")


# ======================================================================
# stage: structure
# ======================================================================
# Pairs of columns that live in the SAME identifier namespace, and can
# therefore be compared directly rather than merely counted together. Each
# entry is (label, how to read side A, how to read side B, what it tests).
# A tuple ("interpro_sigs", "Pfam") means "the Pfam: entries inside
# interpro_sigs"; a bare string is a whole column.
AGREEMENT_PAIRS = [
    ("pfam: hmmsearch vs interproscan",
     "pfam_accs", ("interpro_sigs", "Pfam"),
     "the same Pfam-A library searched by two implementations"),
    ("ncbifam: hmmsearch vs interproscan",
     "ncbifam_accs", ("interpro_sigs", "NCBIfam"),
     "the same NCBIfam library searched by two implementations"),
    ("ko: eggnog vs kofamscan",
     "ko", "kofam_ko",
     "orthology by DIAMOND search against orthology by per-family HMM"),
    ("pfam names: hmmsearch vs eggnog",
     "pfam_hits", "pfams_emapper",
     "domains found directly against domains carried by the eggNOG ortholog"),
]


def _agreement_sets(df, spec):
    """One column, or one member database inside interpro_sigs, as id sets.

    Version suffixes are dropped and the `ko:` prefix eggNOG writes is
    stripped, because an identifier that differs only in its decoration is the
    same identifier and counting it as a disagreement would be an artefact of
    formatting rather than a finding.
    """
    if isinstance(spec, tuple):
        col, prefix = spec
        if col not in df.columns:
            return None
        pref = prefix + ":"

        def read(v):
            if not isinstance(v, str) or not v.strip():
                return frozenset()
            return frozenset(x.strip()[len(pref):].split(".")[0]
                             for x in v.split(";")
                             if x.strip().startswith(pref))
        return df[col].map(read)
    if spec not in df.columns:
        return None

    def read(v):
        if not isinstance(v, str) or not v.strip():
            return frozenset()
        out = set()
        for x in re.split(r"[;,]", v):
            x = x.strip()
            if x.lower().startswith("ko:"):
                x = x[3:]
            x = x.split(".")[0]
            if x and x != "-":
                out.add(x)
        return frozenset(out)
    return df[spec].map(read)


def source_agreement(df):
    """Do the sources that reach the same protein AGREE about it?

    Coverage overlap and concordance are different questions, and only the
    second one tells you whether the evidence is corroborated. Two searches of
    the same library reaching 34,000 proteins in common says nothing until you
    ask whether they name the same families. This computes that, per pair.

    The distinction that matters most in the output is `disjoint`: both sides
    called something and they share nothing. A handful of those is ordinary
    (a paralogue boundary, a threshold near a family edge). A rate near 100%
    is almost never real disagreement - it means the two columns are in
    different namespaces and nothing is being compared at all. That is not
    hypothetical: NCBIfam family NAMES were being compared against
    InterProScan NCBIfam ACCESSIONS, which read as 97% conflict between two
    searches of one library that in fact agreed.
    """
    rows = []
    for label, aspec, bspec, about in AGREEMENT_PAIRS:
        a, b = _agreement_sets(df, aspec), _agreement_sets(df, bspec)
        if a is None or b is None:
            continue
        ha, hb = a.map(bool), b.map(bool)
        both = ha & hb
        n = int(both.sum())
        if not n:
            continue
        ident = part = disj = a_sup = b_sup = mutual = 0
        for x, y in zip(a[both], b[both]):
            if x == y:
                ident += 1
            elif x & y:
                part += 1
                if x > y:
                    a_sup += 1
                elif x < y:
                    b_sup += 1
                else:
                    mutual += 1
            else:
                disj += 1
        rows.append(dict(
            comparison=label, tests=about,
            a=aspec if isinstance(aspec, str) else ":".join(aspec),
            b=bspec if isinstance(bspec, str) else ":".join(bspec),
            a_only=int((ha & ~hb).sum()), b_only=int((hb & ~ha).sum()),
            both=n, identical=ident, overlapping=part, disjoint=disj,
            a_superset=a_sup, b_superset=b_sup, mutually_exclusive=mutual,
            pct_agree=round(100.0 * (ident + part) / n, 1),
            pct_disjoint=round(100.0 * disj / n, 1)))
    return pd.DataFrame(rows)


def write_source_agreement(df, path):
    """Write the concordance table and say what it found."""
    tab = source_agreement(df)
    cols = ["comparison", "tests", "a", "b", "a_only", "b_only", "both",
            "identical", "overlapping", "disjoint", "a_superset",
            "b_superset", "mutually_exclusive", "pct_agree", "pct_disjoint"]
    if tab.empty:
        # Written anyway, with its header. "No two sources here share a
        # namespace" is a real answer, and a declared output that is only
        # sometimes created makes the stage look unfinished and rerun forever.
        tab = pd.DataFrame(columns=cols)
    with atomic_out(path) as tmp:
        tab.to_csv(tmp, sep="	", index=False)
    if tab.empty:
        log("agreement | no two sources share an identifier namespace in this "
            f"run, so there is nothing to cross-check; wrote {path} empty")
        return
    for _, r in tab.iterrows():
        line = (f"agreement | {r['comparison']}: {r['both']} protein(s) called "
                f"by both, {r['pct_agree']}% agree "
                f"({r['identical']} identical, {r['overlapping']} overlapping), "
                f"{r['disjoint']} disjoint")
        # Near-total disjointness between two views of one library is a
        # namespace mismatch far more often than it is a real disagreement,
        # so say which it looks like rather than reporting a number that
        # invites the wrong conclusion.
        if r["both"] >= 50 and r["pct_disjoint"] >= 90:
            log(f"{line}. {r['pct_disjoint']}% disjoint between two views of "
                f"{r['tests']} is not a credible rate of real disagreement - "
                f"check that {r['a']} and {r['b']} carry the same kind of "
                "identifier (names against accessions compare as total "
                "conflict)", "WARN")
        elif r["pct_disjoint"] >= 20:
            log(line + " — high enough to be worth reading before trusting "
                "either source alone", "WARN")
        else:
            log(line)
    log(f"wrote {path}")


def mean_plddt(pdb_text):
    """Mean CA B-factor of an ESMFold model, which is where it stores pLDDT.
    None when the model carries no CA atoms."""
    b = [float(l[60:66]) for l in pdb_text.splitlines()
         if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    return sum(b) / len(b) if b else None


def read_plddt(structures):
    """protein_id -> mean pLDDT (None for a recorded "NA") from plddt.tsv.
    Empty when the table is missing, e.g. structures rsynced in by hand."""
    out = {}
    path = f"{structures}/plddt.tsv"
    if not nonempty(path):
        return out
    with open(path, encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                out[f[0]] = float(f[2])
            except ValueError:
                out[f[0]] = None
    return out


def vram_fit_length(cfg, torch):
    """The longest sequence this card can fold without oversubscribing, or None.

    Call this AFTER the weights are resident, so what it measures is the
    headroom that is actually left rather than the size of the card.

    Why this exists. ESMFold's cost does not rise smoothly with length; it
    rises smoothly until the working set stops fitting in VRAM and then falls
    off a cliff. Measured over 1,819 folds on one 16 GB card: 22.1 s median at
    470-478 aa, then 140 s at 481 aa and up to 2,053 s by 491 aa. A 0.6%
    increase in length cost 6x, and shortly after that 90x. Nothing reports an
    OOM, because the driver quietly pages device memory to host RAM instead of
    failing - so the run does not stop, it just stops being finishable, and the
    only outward sign is that a stage which was going to take an hour is now
    going to take a week.

    It gets worse than slow. On the machine this was measured on, folding above
    that cliff also produced repeated `CUDA driver error: device not ready`
    faults and then took the whole host down twice with a hypervisor bugcheck
    inside eleven minutes - the GPU there is reached through a virtualisation
    layer, and sustained paging across it is what broke. That is a defect in
    somebody else's code and not something this tool can fix, but staying under
    the cliff avoids it entirely, which is the point of this function.

    The estimate. Peak footprint above the resident weights is dominated by
    terms quadratic in length (the pair representation and the triangular
    attention that runs over it), so

        max_len = sqrt(headroom / bytes_per_residue_pair)

    `esmfold_bytes_per_residue_pair` is EMPIRICAL, calibrated so that the
    measured 478 aa / 4.8 GB observation comes out exactly at 478, and it is a
    config key precisely because one card is not a law of nature. Raise it to
    be more conservative, lower it if your card demonstrably folds longer
    sequences at a smooth rate, or set `esmfold_vram_cap: false` to switch the
    whole check off and go back to `max_len_structure` alone.
    """
    if not cfg.get("esmfold_vram_cap", True):
        return None
    try:
        free, total = torch.cuda.mem_get_info()
    except Exception as e:                                  # noqa: BLE001
        # Old torch, or a build where mem_get_info is absent. Not a reason to
        # fail the stage - just say the guard is off rather than pretend it ran.
        log(f"esmfold: cannot read free VRAM ({type(e).__name__}), so the "
            "length cap falls back to max_len_structure alone", "WARN")
        return None
    reserve = float(cfg.get("esmfold_vram_reserve_gb", 0.5)) * 1024 ** 3
    per_pair = float(cfg.get("esmfold_bytes_per_residue_pair", 21000))
    headroom = free - reserve
    if headroom <= 0 or per_pair <= 0:
        log(f"esmfold: only {free / 1024 ** 3:.1f} GB of VRAM is free after "
            f"loading the weights, which is under the {reserve / 1024 ** 3:.1f} "
            "GB reserve; folding anything at all will oversubscribe this card",
            "WARN")
        return 0
    n = int((headroom / per_pair) ** 0.5)
    log(f"esmfold: {free / 1024 ** 3:.1f} of {total / 1024 ** 3:.1f} GB VRAM "
        f"free with the weights resident, so sequences up to {n} aa fit "
        "without paging to host memory. Past that, folds slow by one to two "
        "orders of magnitude with no OOM raised (esmfold_vram_cap: false "
        "turns this off; esmfold_bytes_per_residue_pair recalibrates it)")
    return n


def folded_already(path):
    """A committed structure, as opposed to what a killed writer left behind.

    The resume rule for this stage is "the .pdb is there, so it is folded",
    and open() truncates the moment it is called, so a run interrupted between
    the open and the flush leaves a 0-byte or header-only file that
    os.path.exists cannot tell from a finished one. Every later run then
    counted it, and the protein was permanently absent from Foldseek with
    nothing anywhere to say why. Structures are renamed into place now, so
    this cannot happen again, but the files already on disk from before it
    still can.

    Stops at the first coordinate line, so checking a whole directory costs
    one short read per file rather than a full pass over gigabytes.
    """
    try:
        if os.path.getsize(path) == 0:
            return False
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    return True
    except OSError:
        return False
    return False


def stage_esmfold(cfg, p):
    os.makedirs(p.structures, exist_ok=True)
    if not nonempty(p.dark):
        log("no unannotated proteins to fold, skipping")
        open(p.struct_done, "w", encoding="utf-8").close()
        return

    # Decide whether there is ANY work before touching the GPU. Loading the
    # weights puts ~11 GB on the card, and on a machine where the GPU is
    # reached through a virtualisation layer that upload is itself a risk: a
    # resumed run whose whole remaining work-list is already folded, or is
    # excluded by max_len_structure, used to load the model in full and then
    # discover it had nothing to do. The VRAM cap below still needs the
    # weights resident to measure headroom, so only the STATIC limit can be
    # applied this early - which is exactly the one that answers "is this
    # stage already finished?".
    static_cap = int(cfg["max_len_structure"])
    pending = [q for q, t in read_fasta(p.dark)
               if len(t) <= static_cap
               and not folded_already(f"{p.structures}/{q}.pdb")]
    if not pending:
        n_have = len(glob.glob(f"{p.structures}/*.pdb"))
        log(f"esmfold: every sequence at or under max_len_structure="
            f"{static_cap} is already folded ({n_have} model(s) on disk), so "
            "the GPU is not touched at all")
        open(p.struct_done, "w", encoding="utf-8").close()
        return
    # CUDA_VISIBLE_DEVICES is read once, when the CUDA runtime initialises.
    # Setting it after `import torch` (or after the first torch.cuda call)
    # is a no-op, and gpu_device silently folds on physical GPU 0.
    os.environ["CUDA_VISIBLE_DEVICES"] = str(cfg["gpu_device"])
    if "torch" in sys.modules:
        log("torch was already imported before gpu_device could be applied; "
            "folding may not land on the requested GPU", "WARN")
    try:
        import torch
    except ImportError as e:
        die(f"ESMFold needs torch ({e}). Install it, or set "
            "run.structure: false.")
    if not torch.cuda.is_available():
        die(f"no CUDA device visible (gpu_device={cfg['gpu_device']}); ESMFold "
            "on CPU is impractical at this scale. Set run.structure: false, or "
            "run this stage on the GPU box and rsync results/structures back.")

    # Two ways to get the same weights. fair-esm is the reference
    # implementation, but its esmfold extra needs openfold built from a pinned
    # 2022 commit whose CUDA kernels do not compile against a current toolkit
    # at all — on an sm_120 card there is no version of it that works.
    # transformers ships the same facebook/esmfold_v1 weights with openfold's
    # needed parts vendored, exposes the same infer_pdb(seq) -> pdb string,
    # and builds on modern torch. Prefer fair-esm when it is importable so
    # nothing changes for an existing install; fall back rather than fail.
    model, backend = None, ""
    try:
        import esm
        model = esm.pretrained.esmfold_v1().eval().cuda()
        backend = "fair-esm"
    except Exception as e_esm:
        try:
            from transformers import EsmForProteinFolding
            model = EsmForProteinFolding.from_pretrained(
                "facebook/esmfold_v1", low_cpu_mem_usage=True).eval().cuda()
            backend = "transformers"
        except Exception as e_hf:
            die("ESMFold needs either fair-esm[esmfold] or transformers.\n"
                f"  fair-esm:     {type(e_esm).__name__}: {e_esm}\n"
                f"  transformers: {type(e_hf).__name__}: {e_hf}\n"
                "Install one, or set run.structure: false.")
    log(f"esmfold: {backend} backend on "
        f"{torch.cuda.get_device_name(0)}")

    base_chunk = cfg["esmfold_chunk_size"]

    def set_chunk(n):
        """Chunked attention is what keeps a long sequence inside VRAM.

        fair-esm puts set_chunk_size on the model; transformers puts it on the
        folding trunk. If neither exists the OOM retry below still works, it
        just cannot make the second attempt cheaper, so say so once.
        """
        for obj in (model, getattr(model, "trunk", None),
                    getattr(model, "esm_folding_trunk", None)):
            f = getattr(obj, "set_chunk_size", None) if obj is not None else None
            if callable(f):
                f(n)
                return True
        return False

    if not set_chunk(base_chunk):
        log("this ESMFold build exposes no set_chunk_size, so esmfold_chunk_size "
            "has no effect and an OOM retry cannot lower it", "WARN")

    def normalise_plddt(pdb):
        """Put pLDDT on the 0-100 scale every other part of this tool assumes.

        fair-esm writes pLDDT into the B-factor column as 0-100; transformers
        writes the same quantity as 0-1. thresholds.esmfold_min_plddt is 70,
        so under the transformers backend every model would score below the
        gate and the entire structure stage would be silently discarded — the
        worst kind of failure, because the PDBs exist and look fine. Detect
        the scale from the data rather than from the backend name, so a build
        that changes convention is still handled.
        """
        vals = [l[60:66] for l in pdb.splitlines() if l.startswith(("ATOM", "HETATM"))]
        try:
            hi = max(float(v) for v in vals if v.strip())
        except ValueError:
            return pdb
        if hi > 1.5:
            return pdb
        out = []
        for l in pdb.splitlines(True):
            if l.startswith(("ATOM", "HETATM")) and len(l) >= 66:
                try:
                    l = f"{l[:60]}{float(l[60:66]) * 100:6.2f}{l[66:]}"
                except ValueError:
                    pass
            out.append(l)
        return "".join(out)

    # Same trap the interpro stage already handles, in a stage that did not.
    # A Prodigal/Prokka ORF can carry a trailing '*' for the stop codon, and
    # ESMFold rejects the sequence outright: "Invalid character in the
    # sequence: *". Observed on real data at protein 350 of 1912, three hours
    # into a run, killing the stage over 7 sequences out of 1912. A stop is
    # not a residue, so it is dropped; anything else outside the standard 20
    # becomes X, which ESMFold accepts and which says "unknown residue"
    # rather than inventing one.
    _STD = set("ACDEFGHIKLMNPQRSTVWY")

    def foldable(s):
        body = s[:-1] if s.endswith("*") else s
        return "".join(c if c in _STD else "X" for c in body.upper())

    cap, cap_why = int(cfg["max_len_structure"]), "max_len_structure"
    vram_cap = vram_fit_length(cfg, torch)
    if vram_cap is not None and vram_cap < cap:
        cap, cap_why = vram_cap, "the card's free VRAM"

    raw = [(q, s) for q, s in read_fasta(p.dark)
           if len(s) <= cap]
    seqs, n_fixed = [], 0
    for q, s in raw:
        f = foldable(s)
        if f != s:
            n_fixed += 1
        if f:
            seqs.append((q, f))
    if n_fixed:
        log(f"esmfold: {n_fixed}/{len(raw)} sequence(s) carried a stop codon or "
            "a non-standard residue; folding a cleaned copy (stop dropped, "
            "anything else outside the standard 20 replaced with X). ESMFold "
            "refuses the raw sequence and the stage would die on it", "WARN")
    seqs.sort(key=lambda x: len(x[1]))

    over = [(q, len(t)) for q, t in read_fasta(p.dark) if len(t) > cap]
    if over:
        longest = max(l for _, l in over)
        log(f"esmfold: {len(over)} sequence(s) are longer than {cap} aa "
            f"(up to {longest}) and will NOT be folded; the limit came from "
            f"{cap_why}. They are reported as never attempted, not as "
            "failures - fold them on a card with more memory, in the cloud, "
            "or on CPU, and drop the models into "
            f"{p.structures} before rerunning foldseek", "WARN")
    log(f"esmfold: folding {len(seqs)} sequences, shortest first, at "
        f"chunk_size={base_chunk}")
    t_fold = time.time()

    plddt_path = f"{p.structures}/plddt.tsv"
    # Resuming is the normal mode for this stage, so the table is appended to,
    # never truncated: a rerun that folds the last 50 of 2000 must not throw
    # away the confidence recorded for the other 1950. Anything already folded
    # but missing a row gets its pLDDT recomputed from the PDB, so the table
    # always covers results/structures/.
    seen = set(read_plddt(p.structures))
    fresh = not nonempty(plddt_path)
    done = skipped = consecutive = 0
    failed, wedged = [], False
    max_consecutive = cfg["esmfold_max_consecutive_failures"]
    with open(plddt_path, "a", encoding="utf-8") as ph:
        if fresh:
            ph.write("protein_id\tlength\tmean_plddt\n")
        for pid, seq in seqs:
            out_pdb = f"{p.structures}/{pid}.pdb"
            if os.path.exists(out_pdb):
                if not folded_already(out_pdb):
                    log(f"esmfold: {out_pdb} holds no atom, which is what an "
                        "interrupted writer leaves behind; refolding it "
                        "rather than counting it", "WARN")
                else:
                    skipped += 1
                    if pid not in seen:
                        with open(out_pdb, encoding="utf-8") as fh:
                            mp = mean_plddt(fh.read())
                        ph.write(f"{pid}\t{len(seq)}\t"
                                 + (f"{mp:.1f}\n" if mp is not None
                                    else "NA\n"))
                        ph.flush()
                        seen.add(pid)
                    continue
            pdb, chunk = None, base_chunk
            why = ""
            for attempt in (0, 1):
                try:
                    with torch.no_grad():
                        pdb = normalise_plddt(model.infer_pdb(seq))
                    break
                except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
                    # Two different faults land here and both answer to a
                    # smaller chunk. OOM is the obvious one: torch before
                    # 1.13, and OOM raised inside cuDNN/cuBLAS, report it as a
                    # plain RuntimeError. The other is the driver watchdog. A
                    # long sequence at a large chunk runs one attention kernel
                    # for long enough that the display driver resets the card
                    # mid-fold, which surfaces as "CUDA driver error: device
                    # not ready", not as an OOM. Halving the chunk shortens
                    # each kernel and clears it, so both faults get the same
                    # treatment and neither is re-raised: one hard sequence
                    # must not cost the whole stage.
                    oom = isinstance(e, torch.cuda.OutOfMemoryError)                         or "out of memory" in str(e).lower()
                    kind = "OOM" if oom else "CUDA fault"
                    why = str(e).splitlines()[0][:120]
                    try:
                        torch.cuda.empty_cache()
                    except Exception:
                        # A wedged device refuses even this. The consecutive
                        # counter below is what stops the run, not this.
                        pass
                    if attempt == 0:
                        chunk = max(8, chunk // 2)
                        set_chunk(chunk)
                        log(f"{kind} on {pid} (len {len(seq)}), retry at "
                            f"chunk_size={chunk}", "WARN")
                    else:
                        log(f"skipping {pid} (len {len(seq)}): {why}", "WARN")
            set_chunk(base_chunk)
            # Hand the cached blocks back between sequences. Folding is a
            # variable-shape workload - the work-list runs from tens of
            # residues to hundreds - so the caching allocator ends up holding
            # blocks shaped for the last sequence rather than the next one.
            # Releasing them costs one synchronisation per protein, which is
            # nothing beside a fold that takes minutes, and it keeps the
            # reserve from sitting at the card's capacity where the next OOM
            # retry has no room to drop the chunk size into.
            with contextlib.suppress(Exception):
                torch.cuda.empty_cache()
            if pdb is None:
                failed.append((pid, len(seq), why))
                consecutive += 1
                if consecutive >= max_consecutive:
                    # Every sequence failing in a row is a card that has
                    # stopped working, not a run of hard proteins. Walking the
                    # rest of the list would take hours to produce nothing.
                    log(f"esmfold: {consecutive} sequences failed in a row, "
                        "which reads as a wedged GPU rather than hard "
                        "sequences; stopping here", "ERROR")
                    wedged = True
                    break
                continue
            consecutive = 0
            # Renamed into place, never written in place. This loop runs for
            # hours and is interrupted often — a wedged card, a bugcheck, a
            # Ctrl-C — and the resume rule is "the file is there, so it is
            # folded". atomic_out's temp is dot-prefixed and keeps the
            # extension, so glob("*.pdb") never sees it and a leftover cannot
            # come back as a structure named ".P0001.9134.7.part".
            with atomic_out(out_pdb) as tmp_pdb:
                with open(tmp_pdb, "w", encoding="utf-8") as fh:
                    fh.write(pdb)
            mp = mean_plddt(pdb)
            ph.write(f"{pid}\t{len(seq)}\t"
                     + (f"{mp:.1f}\n" if mp is not None else "NA\n"))
            ph.flush()
            seen.add(pid)
            done += 1
            # Rate and ETA, not just a count. Folding time rises steeply with
            # length and the work-list is sorted shortest first, so a count
            # alone says nothing about how long the rest will take - and when
            # a run slowed by a factor of ten partway through, the log gave no
            # sign of it and the slowdown had to be read off file timestamps.
            # Every 10, because at the long end of a real dark set one protein
            # can take minutes and 25 of them is hours between lines.
            if done % 10 == 0:
                rate = (time.time() - t_fold) / done
                left = len(seqs) - (done + skipped)
                log(f"esmfold: {done + skipped}/{len(seqs)}, "
                    f"{len(seq)} aa at chunk_size={base_chunk}, "
                    f"{rate:.0f}s each, ~{left * rate / 3600:.1f}h for the "
                    f"remaining {left}")
    # Both counters, because a resumed run folds nothing and "0 structures"
    # reads like a failure.
    log(f"esmfold: {done} new, {skipped} already present, "
        f"{done + skipped} structures in {p.structures}")

    if failed:
        # Name the casualties in a file rather than only in the log, so the
        # shortfall survives into the results directory and can be read back
        # by whoever asks why a protein has no structure evidence.
        miss = f"{p.structures}/esmfold_failed.tsv"
        with open(miss, "w", encoding="utf-8") as fh:
            fh.write("protein_id\tlength\terror\n")
            for pid, ln, why in failed:
                fh.write(f"{pid}\t{ln}\t{why}\n")
        log(f"esmfold: {len(failed)} sequence(s) could not be folded; "
            f"listed in {miss}", "WARN")

    if wedged and not cfg["esmfold_allow_partial"]:
        die(f"esmfold stopped after {max_consecutive} consecutive failures "
            f"with {done + skipped} of {len(seqs)} structures written.\n"
            "  Every PDB already written is kept, so rerunning resumes from "
            "there rather than starting over.\n"
            "  A card that has stopped responding usually needs the machine "
            "or the WSL session restarted, not another attempt.\n"
            "  To go on with the structures already folded instead, set "
            "esmfold_allow_partial: true.")
    if wedged:
        log("esmfold: continuing with a partial structure set because "
            "esmfold_allow_partial is on; Foldseek will search only what was "
            "folded, so absent structure evidence here means not attempted, "
            "not absent", "WARN")

    open(p.struct_done, "w", encoding="utf-8").close()


def stage_foldseek(cfg, p):
    os.makedirs(os.path.dirname(p.foldseek), exist_ok=True)
    pdbs = sorted(glob.glob(f"{p.structures}/*.pdb"))
    if not pdbs:
        log("no structures, writing empty hit table")
        with atomic_out(p.foldseek) as tmp:
            open(tmp, "w", encoding="utf-8").close()
        return
    if not have("foldseek"):
        die("foldseek not found (conda install -c bioconda foldseek)")
    target = cfg["db"]["foldseek_target"]
    if not glob.glob(target + "*"):
        die(f"foldseek target database not found: {target}")
    # Deliberately serial: a Foldseek target index is tens to hundreds of GB
    # and two searches at once will thrash or OOM. The parallelism worth having
    # here is inside Foldseek's own --threads.
    ram = cfg.get("ram_gb") or 0
    fs_mem = ["--split-memory-limit", f"{ram}G"] if ram else []
    th = cfg["thresholds"]

    # The pLDDT gate. A ~45-pLDDT model of a short dark ORF matching a fold at
    # TM 0.5 is noise, and a hit here is what promotes a protein out of
    # 4_dark into 3s_structure_only, so low-confidence models are never
    # searched at all rather than filtered afterwards.
    query, queried = p.structures, pdbs
    min_plddt = float(th.get("esmfold_min_plddt") or 0)
    if min_plddt > 0:
        pl = read_plddt(p.structures)
        if not pl:
            log(f"no plddt.tsv in {p.structures}, so the pLDDT >= "
                f"{min_plddt:g} gate cannot be applied and every model is "
                "searched", "WARN")
        else:
            keep, ungated = [], 0
            for f in pdbs:
                v = pl.get(os.path.basename(f)[:-len(".pdb")])
                if v is None:
                    ungated += 1
                    keep.append(f)
                elif v >= min_plddt:
                    keep.append(f)
            if not keep:
                log(f"all {len(pdbs)} models are below pLDDT {min_plddt:g}; "
                    "writing an empty hit table", "WARN")
                with atomic_out(p.foldseek) as tmp:
                    open(tmp, "w", encoding="utf-8").close()
                return
            if len(keep) < len(pdbs):
                # Foldseek takes a directory, so the survivors are linked into
                # one rather than passed as thousands of argv entries.
                query = f"{p.R}/foldseek/query_hq"
                shutil.rmtree(query, ignore_errors=True)
                os.makedirs(query, exist_ok=True)
                for f in keep:
                    dst = f"{query}/{os.path.basename(f)}"
                    try:
                        os.link(f, dst)
                    except OSError:
                        shutil.copyfile(f, dst)
                queried = keep
                log(f"pLDDT gate: searching {len(keep)}/{len(pdbs)} models at "
                    f"or above {min_plddt:g}"
                    + (f", {ungated} with no pLDDT recorded" if ungated else ""),
                    "WARN")

    targets = [target] + list(cfg["db"].get("foldseek_extra_targets") or [])
    # The per-target result files are concatenated into one table, and Foldseek
    # writes nothing that says which database a row came from — the accession
    # cannot do it either, since AFDB50 and AFDB-SwissProt both look like
    # AF-<acc>-F1-model_v4. Without that column the merged table is ranked by
    # raw bitscore alone and an unannotated AFDB50 model beats a described
    # Swiss-Prot or PDB hit. Label every row with its target database instead.
    # The field list is named once so the header written below always describes
    # the columns actually requested.
    # Ask for the columns the analysis actually needs, not the legacy ten.
    # FOLDSEEK_COLS says as much in its own comment and the stage was ignoring
    # it: without qtmscore the TM gate silently falls back to alntmscore,
    # which is normalised by the ALIGNMENT, so a 40-residue local match inside
    # a 300-residue query can score 0.6 while the two proteins share almost no
    # fold; and without qlen the coverage backstop cannot be applied at all.
    # finalise even printed "re-run the foldseek stage to get qtmscore and
    # qlen", which was advice the stage could not take, because re-running
    # asked for the same ten columns again.
    fs_fields = ",".join(FOLDSEEK_COLS)
    parts, seen = [], {}
    for i, tgt in enumerate(targets):
        if not glob.glob(str(tgt) + "*"):
            log(f"foldseek target missing, skipping: {tgt}", "WARN")
            continue
        # The basename of the target path is the name a user recognises and is
        # what foldseek_target_priority is matched against.
        lab = os.path.basename(str(tgt).rstrip("/\\")) or str(tgt)
        if lab in seen:
            # Two targets with the same basename would be indistinguishable in
            # the output, which is the defect this column exists to fix.
            lab = f"{lab}#{i}"
            log(f"two foldseek targets share the basename "
                f"'{os.path.basename(str(tgt))}'; the second is labelled "
                f"'{lab}' in {p.foldseek} and will not match "
                "foldseek_target_priority by its bare name", "WARN")
        seen[lab] = tgt
        outp = f"{p.R}/foldseek/hits_{i}.tsv"
        tmpd = f"{p.R}/foldseek/tmp{i}"
        def search(fields):
            run_cmd(["foldseek", "easy-search", query, tgt, outp,
                     tmpd, "--format-output", fields,
                     "-e", cfg["thresholds"]["foldseek_evalue"],
                     "--max-seqs", 300, "--threads", cfg["threads"], "-v", 1]
                    + fs_mem + tool_args(cfg, "foldseek"))

        # try/finally so the scratch tree goes on EVERY exit, not only the
        # successful one. Against AFDB50 it is tens to hundreds of GB, and
        # CLAUDE.md already records `results/foldseek/tmp*` as a thing that is
        # never cleaned up; each re-raise below used to add another way to
        # leave one behind, and the one that matters is a genuine failure - a
        # full disk or an OOM kill - where the tree is largest and a leak
        # hurts most.
        try:
            try:
                search(fs_fields)
            except RuntimeError as e:
                # RuntimeError, not StageError: run_cmd raises a plain
                # RuntimeError on a non-zero exit, and StageError is a
                # SUBCLASS of it, so `except StageError` could never catch the
                # one failure this fallback exists for. It caught nothing and
                # the branch was unreachable; the test that covered it
                # injected a StageError no tool ever raises. StageError is
                # still caught here, being a subclass.
                #
                # A rejected format code costs SECONDS, not the search.
                # EasyStructureSearch.cpp validates --format-output in
                # getOutputFormat (line 42) while the temporary directory is
                # not created until line 59, so foldseek exits before it
                # prefilters anything and leaves no tree behind. The retry is
                # therefore cheap, which is why the test below is permissive:
                # when the message says a format code was rejected but the
                # field cannot be named, fall back anyway rather than lose the
                # structural evidence to a wording change.
                # Both halves of the phrase must be on ONE line.
                # "<path> does not exist" is stock MMseqs2 wording for a
                # missing database or temp dir, so tested across the whole
                # multi-line stderr tail this would splice an unrelated line
                # onto the words "format code" and read a path as a rejected
                # column - falling back on a failure that has nothing to do
                # with the format list.
                hit = next((ln for ln in str(e).splitlines()
                            if "format code" in ln.lower()
                            and "does not exist" in ln.lower()), "")
                if fs_fields == ",".join(FOLDSEEK_COLS_LEGACY) or not hit:
                    raise
                m = re.search(r"format code\W*(\S+?)\s+does not exist",
                              hit, re.I)
                # Defensive, not a claim about upstream: both emitters print
                # the name bare (foldseek LocalParameters.cpp, MMseqs2
                # Parameters.cpp), so this only covers wording drift.
                bad = m.group(1).strip("\"'`") if m else ""
                # Matched case-insensitively above, so compare that way too -
                # being half case-insensitive would classify QTMSCORE as
                # undroppable and say something false about it.
                droppable = set(FOLDSEEK_COLS) - set(FOLDSEEK_COLS_LEGACY)
                lower = {c.lower() for c in droppable}
                if bad and bad.lower() not in lower:
                    # Falling back cannot remove this one: FOLDSEEK_COLS_LEGACY
                    # is a strict SUBSET of FOLDSEEK_COLS, so the legacy list
                    # asks for it too and the retry would fail identically.
                    log(f"foldseek rejected the format code '{bad}', which the "
                        "legacy column list asks for as well, so falling back "
                        "cannot help: the only columns it drops are "
                        + ", ".join(sorted(droppable)) + ". Not retrying; the "
                        "error below is the real one.", "WARN")
                    raise
                named = f"the format code '{bad}'" if bad else \
                    "a format code it could not name"
                log(f"foldseek rejected {named}, so this build predates "
                    "qtmscore/ttmscore (Foldseek 5 and earlier); falling back "
                    "to the legacy columns, which drop "
                    + ", ".join(sorted(droppable)) + ". The TM gate will use "
                    "alntmscore, which is normalised by the alignment rather "
                    "than the query, and no coverage filter can be applied - "
                    "a short local match can pass it. Upgrade Foldseek to "
                    "restore the stricter gate.", "WARN")
                fs_fields = ",".join(FOLDSEEK_COLS_LEGACY)
                search(fs_fields)
        finally:
            shutil.rmtree(tmpd, ignore_errors=True)
        parts.append((lab, outp))
    n_rows = 0
    with atomic_out(p.foldseek) as tmp, open(tmp, "w", encoding="utf-8") as out:
        # A named header, because with target_db appended a legacy 10-column
        # row and a new 11-column one cannot be told apart by field count.
        out.write("#" + "\t".join(fs_fields.split(",") + [FOLDSEEK_DB_COL])
                  + "\n")
        for lab, f in parts:
            with open(f, encoding="utf-8") as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line:
                        out.write(line + "\t" + lab + "\n")
                        n_rows += 1
    if not n_rows:
        # A header on its own is a non-empty file with no records, and
        # warn_if_no_records would read that as a format failure. No hits is a
        # real outcome here, so say it as an empty table instead.
        open(p.foldseek, "w", encoding="utf-8").close()
        log("foldseek: no rows in any target result file; writing an empty "
            "hit table", "WARN")
    log("foldseek: searched " + str(len(parts)) + " target database(s)"
        + (": " + ", ".join(lab for lab, _ in parts) if parts else ""))

    # Cluster the unannotated structures against each other. Fifty dark
    # proteins sharing a fold is a far stronger signal than fifty singletons,
    # and it needs no reference database at all.
    if cfg.get("foldseek_self_cluster", True):
        try:
            run_cmd(["foldseek", "easy-cluster", query,
                     f"{p.R}/foldseek/selfclu", f"{p.R}/foldseek/tmpc",
                     "-e", th.get("foldseek_cluster_evalue", 1e-2),
                     # Without an explicit TM-score and coverage requirement a
                     # "fold group" can be built from loose local 3Di matches,
                     # over-merging unrelated dark proteins into the headline
                     # "recurrent unknown folds" number.
                     "--tmscore-threshold",
                     th.get("foldseek_cluster_tmscore", 0.50),
                     "-c", th.get("foldseek_cluster_coverage", 0.80),
                     "--cov-mode", 0,
                     "--threads", cfg["threads"], "-v", 1]
                    + fs_mem + tool_args(cfg, "foldseek"))
            shutil.rmtree(f"{p.R}/foldseek/tmpc", ignore_errors=True)
            cl = f"{p.R}/foldseek/selfclu_cluster.tsv"
            if nonempty(cl):
                d = pd.read_csv(cl, sep="\t", header=None,
                                names=["rep", "member"], encoding="utf-8", encoding_errors="replace")
                # Same normalisation as parse_foldseek: Foldseek may append a
                # chain name (`X.pdb_A`), and stripping only the extension
                # there leaves ids that match nothing in the protein index.
                _chain = r"\.(pdb|cif)(_[A-Za-z0-9]+)?$"
                d["rep"] = d["rep"].str.replace(_chain, "", regex=True)
                d["member"] = d["member"].str.replace(_chain, "", regex=True)
                with atomic_out(p.fold_clusters) as _fc:
                    d.to_csv(_fc, sep="\t", index=False)
                sizes = d.groupby("rep").size()
                log(f"foldseek self-clustering: {len(sizes)} fold groups, "
                    f"largest {int(sizes.max())} members, "
                    f"{int((sizes > 1).sum())} non-singleton")
                # An id format the normalisation missed shows up downstream as
                # "self-clustering was not run", so say it here instead.
                on_disk = {os.path.basename(f)[:-len(".pdb")] for f in queried}
                hit = len(set(d["member"]) & on_disk)
                if hit < len(on_disk):
                    log(f"fold_clusters: {hit}/{len(on_disk)} ids matched a "
                        "folded protein — check the id format", "WARN")
        except Exception as e:
            # Not just RuntimeError: a zero-byte or truncated selfclu table
            # raises pandas' EmptyDataError/ParserError, and this cosmetic
            # sub-step must never fail the stage that just spent hours
            # searching AFDB50.
            log(f"self-clustering failed (non-fatal): {e}", "WARN")


# ======================================================================
# stage: join to the quant matrix
# ======================================================================
def group_members(row, member_cols):
    """Ordered, de-duplicated group members. FragPipe keeps the extra members
    in `Indistinguishable Proteins`, not in the id column, so both are read."""
    ids, seen = [], set()
    for c in member_cols:
        v = row.get(c, "")
        if not isinstance(v, str) or not v:
            continue
        for pid in re.split(r"[;,]", v):
            pid = pid.strip()
            if pid and pid not in seen:
                seen.add(pid)
                ids.append(pid)
    return ids


# ======================================================================
# feature-level input: FragPipe peptide/ion tables and MSstats exports
# ======================================================================
# ======================================================================
# FragPipe manifest: the experimental design, straight from the run table
# ======================================================================
# A .fp-manifest is tab-separated with no header:
#     path <TAB> experiment <TAB> bioreplicate <TAB> data_type
# It already carries condition (experiment) and blocking (bioreplicate) per
# raw file, so when one is supplied nothing else needs to be hand-written.

def manifest_basename(path):
    """Basename of a raw-file path, whatever separator wrote it.

    FragPipe on Windows writes 'D:\\runs\\EX1.mzML'; os.path.basename on Linux
    does not split backslashes, so the whole path would become the sample name
    when the pipeline is run on the analysis machine.
    """
    return os.path.splitext(re.split(r"[\\/]", str(path))[-1])[0]


def read_manifest(path):
    rows = []
    with opener(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.split("\t")
            if len(f) == 1:
                # A hand-edited manifest may use commas or aligned columns.
                # Single spaces are NOT a separator: raw-file paths routinely
                # contain them, and guessing there silently produced a design
                # with one nameless group.
                f = re.split(r"\s{2,}|,", line)
            if len(f) == 1:
                die(f"{path}:{lineno}: only one field in "
                    f"'{line.strip()[:80]}'. A .fp-manifest is tab separated: "
                    "path <TAB> experiment <TAB> bioreplicate <TAB> data_type. "
                    "Single spaces are not accepted because file paths "
                    "contain them.")
            f = [x.strip() for x in f] + [""] * (4 - len(f))
            rows.append({"file": f[0], "experiment": f[1],
                         "bioreplicate": f[2], "data_type": f[3]})
    if not rows:
        die(f"{path}: no rows; expected a FragPipe .fp-manifest "
            "(path/experiment/bioreplicate/data_type, tab separated)")
    m = pd.DataFrame(rows)
    m["basename"] = [manifest_basename(x) for x in m["file"]]
    # FragPipe names output columns after the experiment, appending the
    # bioreplicate when one is set, and falls back to the file basename.
    m["sample"] = [
        (f"{e}_{b}" if e and b else e if e else bn)
        for e, b, bn in zip(m["experiment"], m["bioreplicate"], m["basename"])]
    # Repeated sample names are how FragPipe denotes FRACTIONS: every fraction
    # file of one sample carries the same experiment and bioreplicate, and
    # FragPipe writes a single quant column for the group. Rejecting them made
    # the tool unusable on any fractionated run, which is most offline-
    # fractionated proteomics. Collapse them to one row per sample instead,
    # keeping the fraction count, and only refuse when the rows genuinely
    # disagree about the design.
    m["n_fractions"] = m.groupby("sample")["file"].transform("size")
    if (m["n_fractions"] > 1).any():
        # Fractions of one sample must agree about the design. If they do not,
        # the rows are not fractions at all and collapsing them would invent a
        # sample that never existed, so refuse loudly instead.
        bad = [s for s, g in m.groupby("sample")
               if g["data_type"].str.strip().str.upper().nunique() > 1]
        if bad:
            die(f"{path}: sample(s) {sorted(bad)[:5]} have rows with different "
                "data_type. Fractions of one sample must share "
                "experiment, bioreplicate and data_type.")
        frac = m.loc[m["n_fractions"] > 1, "sample"].nunique()
        log(f"manifest: {frac} sample(s) are split across multiple fraction "
            "files; each group is one quant column in FragPipe's output and is "
            "treated here as one sample")
        m = (m.sort_values("basename")
               .drop_duplicates(subset="sample", keep="first")
               .reset_index(drop=True))
    log(f"manifest: {len(m)} sample(s), "
        f"{m['experiment'].nunique()} experiment(s), "
        f"{m['bioreplicate'].replace('', pd.NA).nunique()} bioreplicate(s)")
    return m


def manifest_candidates(row):
    """Names FragPipe might have used for this run, best guess first."""
    e, b, bn = row["experiment"], row["bioreplicate"], row["basename"]
    out = []
    if e and b:
        out += [f"{e}_{b}", f"{e}-{b}"]
    if e:
        out.append(e)
    out += [bn, row["file"]]
    return [x for x in dict.fromkeys(out) if x]


def map_manifest_to_columns(m, columns, suffix=""):
    """-> ({column: sample}, unmatched_manifest_rows, unmatched_columns).

    `columns` are the detected sample columns; `suffix` is what FragPipe
    appends to the sample name in a wide table (e.g. " Intensity").
    """
    stripped = {}
    for c in columns:
        base = c[:-len(suffix)].strip() if suffix and c.endswith(suffix) else c
        stripped[c] = base
    by_base = {}
    for c, base in stripped.items():
        by_base.setdefault(base, c)
        by_base.setdefault(manifest_basename(base), c)

    mapping, unmatched, claimed = {}, [], defaultdict(list)
    for _, row in m.iterrows():
        hit = next((by_base[cand] for cand in manifest_candidates(row)
                    if cand in by_base), None)
        if hit is None:
            unmatched.append(row["sample"])
        else:
            claimed[hit].append(row["sample"])
            mapping[hit] = row["sample"]
    # Two manifest samples resolving to one quant column used to be silent
    # (last one wins), which produced a design with fewer samples than the
    # manifest and the wrong replicate labels. Fractions are already collapsed
    # in read_manifest, so a collision here is a genuine ambiguity.
    dup = {c: s for c, s in claimed.items() if len(s) > 1}
    if dup:
        det = "; ".join(f"{c!r} <- {s}" for c, s in sorted(dup.items())[:5])
        die("manifest: several samples map to the same quant column: " + det +
            ". Rename the experiment/bioreplicate pairs so each run resolves "
            "to its own column, or drop the extra manifest rows.")
    return mapping, unmatched, [c for c in columns if c not in mapping]


def die_manifest_unmatched(miss_rows, path, columns):
    """One wording for "the manifest and the quant table do not line up".

    Both the peptide-level and the protein-level reader hit this, and they
    used to say it differently: one listed the columns present, the other did
    not, so the same mistake was harder to diagnose from one input than from
    the other.
    """
    die(f"{len(miss_rows)} manifest run(s) match no column in {path}: "
        f"{list(miss_rows)[:5]}\n  columns present: {list(columns)[:8]}\n"
        "  metaannot tried experiment_bioreplicate, experiment and the file "
        "basename. Fix the run names in the manifest so they match the quant "
        "columns, or unset `manifest:` and let the column names stand as the "
        "sample names.")


# ======================================================================
# taxonomy: Unipept peptide LCA, and comparison with eggNOG seed taxonomy
# ======================================================================
# eggNOG's seed_ortholog taxid is the taxon of the best-matching *reference*
# protein, not of the organism in the sample. Unipept's peptide LCA is an
# independent estimate with different biases (UniProt coverage, tryptic-only
# index, deliberately conservative LCA). Where they agree, the taxon-based
# steps downstream are on firm ground; where they disagree, they are not.

# NCBI renamed the top cellular rank from 'superkingdom' to 'domain' in March
# 2025, and Unipept's API/CLI followed (domain_id/domain_name). Old taxdumps
# and old pept2lca exports still say superkingdom, so both spellings are
# accepted everywhere and normalised to 'domain'.
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_ALIASES = {"superkingdom": "domain"}
TRYPTIC_RE = re.compile(r"(?<=[KR])(?!P)")


def strip_modifications(seq):
    """Peptide sequence as Unipept wants it: bare residues, upper case."""
    seq = str(seq).strip()
    # FragPipe's MSstats export marks the peptide TERMINI with a lowercase n/c
    # in front of the modification mass ("n[230]GEQ...", "PEPTIDEc[-0.98]").
    # Deleting only the bracket group left the marker behind, and the upper()
    # below then turned it into a real N or C residue, so the peptide sent to
    # Unipept was never found.
    seq = re.sub(r"^[nc](?=[\[({])", "", seq)
    seq = re.sub(r"c(?=[\[({][^\])}]*[\])}]$)", "", seq)
    seq = re.sub(r"\[[^\]]*\]|\([^)]*\)|\{[^}]*\}", "", seq)
    # "K.PEPTIDER.S" notation: the flanking residues belong to the protein,
    # not to the peptide.
    seq = re.sub(r"^[A-Za-z-]\.|\.[A-Za-z-]$", "", seq)
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def tryptic_split(pep):
    """Unipept indexes fully tryptic peptides, so a peptide with a missed
    cleavage is often absent. Splitting after K/R (not before P) recovers it,
    at the cost of shorter, less specific peptides."""
    parts = [x for x in TRYPTIC_RE.split(pep) if len(x) >= 5]
    return parts or ([pep] if len(pep) >= 5 else [])


def unipept_http(peptides, cfg):
    """POST a batch to the Unipept API. Convenience path only — the offline
    route (unipept_result in the config) is the one that does not depend on a
    reachable server or on the API version matching."""
    import urllib.request
    u = cfg["unipept"]
    body = json.dumps({"input": list(peptides),
                       "equate_il": bool(u.get("equate_il", True)),
                       "extra": True, "names": True}).encode()
    req = urllib.request.Request(
        u.get("api_url", "https://api.unipept.ugent.be/api/v2/pept2lca.json"),
        data=body, headers={"Content-Type": "application/json",
                            "Accept": "application/json",
                            "User-Agent": f"metaannot/{__version__}"})
    last = None
    for attempt in range(int(u.get("retries", 3))):
        try:
            with urllib.request.urlopen(req, timeout=int(u.get("timeout", 180))) as r:
                return json.loads(r.read().decode())
        except Exception as e:                      # noqa: BLE001
            last = e
            time.sleep(float(u.get("sleep", 1.0)) * (attempt + 1))
    raise RuntimeError(
        f"Unipept request failed after retries: {last}. Run it yourself "
        "instead (the CLI is an npm package now, the old Ruby gem is "
        "unmaintained):\n  npm install -g unipept-cli\n  unipept pept2lca "
        "--equate --all -i peptides.txt -o pept2lca.csv\nthen set "
        "unipept.result to that file. An export from the unipept.ugent.be web "
        "interface is NOT accepted: it has no taxon-id columns.")


def read_unipept_result(path):
    """Ingest pept2lca output from the CLI or a previous run. Tolerant of
    csv/tsv and of the v1/v2 column spellings."""
    df = pd.read_csv(path, sep=None, engine="python", encoding="utf-8", encoding_errors="replace")
    ren = {"peptide": "peptide", "Peptide": "peptide",
           "taxon_id": "taxon_id", "taxon_name": "taxon_name",
           "taxon_rank": "taxon_rank"}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    if "peptide" not in df.columns:
        die(f"{path}: no 'peptide' column; is this a pept2lca result?")
    for r in RANKS:
        # *_id only. A bare `genus` column holds the NAME, and accepting it
        # would put species names into the taxid comparison. The alias covers
        # the pre-2025 spelling of the top rank (superkingdom_id).
        olds = [k for k, v in RANK_ALIASES.items() if v == r]
        for base in [r] + olds:
            for cand in (f"{base}_id", base.capitalize() + "_id",
                         f"{base}_taxon_id"):
                if cand in df.columns and f"{r}_id" not in df.columns:
                    df[f"{r}_id"] = df[cand]
    keep = ["peptide", "taxon_id", "taxon_name", "taxon_rank"] + \
           [f"{r}_id" for r in RANKS]
    # A table with only a peptide column (a web-interface export, or a CLI run
    # without --all) used to pass silently and give every protein an empty
    # lineage, which the verdict logic then reported as disagreement. Refuse.
    if "taxon_id" not in df.columns:
        die(f"{path}: no 'taxon_id' column. This is not a pept2lca result; "
            "run 'unipept pept2lca --equate --all' (npm install -g "
            "unipept-cli), not the web interface.")
    if not any(f"{r}_id" in df.columns for r in RANKS):
        die(f"{path}: no rank taxid columns ({', '.join(r + '_id' for r in RANKS)}). "
            "The lineage columns come from the --all flag: "
            "unipept pept2lca --equate --all -i peptides.txt -o pept2lca.csv")
    return df[[c for c in keep if c in df.columns]]


class NCBITaxonomy:
    """Minimal nodes.dmp / names.dmp reader, so eggNOG seed taxids can be
    given a lineage offline rather than through another web service."""

    def __init__(self, d):
        self.parent, self.rank, self.name = {}, {}, {}
        self.merged, self.deleted = {}, set()
        self.unresolved = set()
        with open(os.path.join(d, "nodes.dmp")) as fh:
            for line in fh:
                f = [x.strip() for x in line.split("|")]
                self.parent[f[0]], self.rank[f[0]] = f[1], f[2]
        # eggNOG 5 seed taxids come from a 2018 taxonomy; many strain ids have
        # since been merged into another node or deleted. Without merged.dmp
        # those proteins looked like taxonomic conflicts, and a *newer*
        # taxdump made it worse, not better.
        mp = os.path.join(d, "merged.dmp")
        if os.path.exists(mp):
            with open(mp, encoding="utf-8") as fh:
                for line in fh:
                    f = [x.strip() for x in line.split("|")]
                    if len(f) > 1 and f[0] and f[1]:
                        self.merged[f[0]] = f[1]
        dp = os.path.join(d, "delnodes.dmp")
        if os.path.exists(dp):
            with open(dp, encoding="utf-8") as fh:
                for line in fh:
                    t = line.split("|")[0].strip()
                    if t:
                        self.deleted.add(t)
        with open(os.path.join(d, "names.dmp")) as fh:
            for line in fh:
                # Cheap substring test first: names.dmp has several million
                # lines and only a fraction are scientific names.
                if "scientific name" not in line:
                    continue
                f = [x.strip() for x in line.split("|")]
                if len(f) > 3 and f[3] == "scientific name":
                    self.name[f[0]] = f[1]
        log(f"NCBI taxonomy: {len(self.parent)} nodes, "
            f"{len(self.merged)} merged, {len(self.deleted)} deleted")

    def current(self, taxid):
        """The live taxid for an id that may have been merged away, or ''
        when it no longer exists at all."""
        t = str(taxid).strip()
        seen = set()
        while t in self.merged and t not in seen:
            seen.add(t)
            t = self.merged[t]
        return t if t in self.parent else ""

    def sci_name(self, taxid):
        return self.name.get(self.current(taxid), "")

    def lineage(self, taxid):
        """-> {rank: taxid} for the ranks in RANKS. Empty if the taxid does not
        resolve, i.e. it was deleted rather than merged."""
        t = self.current(taxid)
        if not t:
            self.unresolved.add(str(taxid))
            return {}
        out, seen = {}, set()
        while t and t in self.parent and t not in seen and t != "1":
            seen.add(t)
            # A taxdump older than March 2025 still calls the top rank
            # 'superkingdom'; RANK_ALIASES folds it onto 'domain'.
            r = RANK_ALIASES.get(self.rank.get(t), self.rank.get(t))
            if r in RANKS and r not in out:
                out[r] = t
            t = self.parent[t]
        return out


def stage_unipept(cfg, p):
    u = cfg.get("unipept") or {}
    os.makedirs(f"{p.R}/unipept", exist_ok=True)

    if u.get("result"):
        if not os.path.exists(u["result"]):
            die(f"unipept.result not found: {u['result']}")
        df = read_unipept_result(u["result"])
        with atomic_out(p.unipept_lca) as tmp:
            df.to_csv(tmp, sep="\t", index=False)
        log(f"unipept: ingested {len(df)} peptide assignments from "
            f"{u['result']}")
        return

    fmt = cfg["quant_format"]
    if fmt not in FEATURE_FORMATS:
        die(f"Unipept needs peptide sequences, but quant_format is '{fmt}'. "
            "Use fragpipe_peptide / fragpipe_ion / msstats_csv / "
            "msstats_feature, or set run.unipept: false.")
    feats = peptide_features(cfg, "unipept")
    raw = {strip_modifications(x) for x in feats["peptide"]}
    raw = {x for x in raw if len(x) >= 5}

    peps = set()
    if u.get("split_missed_cleavages", True):
        for x in raw:
            peps.update(tryptic_split(x))
        log(f"unipept: {len(raw)} peptides -> {len(peps)} fully tryptic "
            "fragments (split_missed_cleavages)")
    else:
        peps = raw
        log(f"unipept: {len(peps)} peptides")

    # Peptide-keyed cache: a rerun only queries what is new.
    cache, cached = {}, 0
    if os.path.exists(p.unipept_cache) and os.path.getsize(p.unipept_cache) > 0:
        try:
            c = pd.read_csv(p.unipept_cache, sep="\t", encoding="utf-8", encoding_errors="replace")
        except pd.errors.EmptyDataError:
            c = None
        if c is None or "peptide" not in c.columns:
            # A previous run that assigned nothing used to leave a zero-byte
            # cache here, and every later run then died in read_csv with no
            # explanation. Start over instead.
            log(f"unipept: cache {p.unipept_cache} is empty or has no "
                "'peptide' column; ignoring it", "WARN")
        else:
            cache = {r["peptide"]: r for _, r in c.iterrows()}
        cached = len(cache)
    with atomic_out(p.unipept_peptides) as tmp, open(tmp, "w", encoding="utf-8") as fh:
        fh.write("\n".join(sorted(peps)) + "\n")
    log(f"unipept: peptide list -> {p.unipept_peptides}")

    todo = sorted(peps - set(cache))
    log(f"unipept: {cached} cached, {len(todo)} to query")

    if todo:
        if not u.get("allow_http", False):
            die("unipept.allow_http is false and the cache is incomplete.\n"
                "Run Unipept yourself and point unipept.result at the output "
                "(npm install -g unipept-cli):\n"
                f"  unipept pept2lca --equate --all -i {p.unipept_peptides} "
                "-o pept2lca.csv\n"
                "then set unipept.result to pept2lca.csv")
        bs = int(u.get("batch_size", 100))
        rows = []
        for i in range(0, len(todo), bs):
            batch = todo[i:i + bs]
            rows.extend(unipept_http(batch, cfg))
            time.sleep(float(u.get("sleep", 0.3)))
            if (i // bs) % 20 == 0:
                log(f"unipept: {min(i + bs, len(todo))}/{len(todo)}")
        # The API simply omits peptides it cannot place, so without an explicit
        # "asked, not found" row every rerun re-sent them — which in a gut
        # metaproteome is most of the list. Cache the misses too.
        found = {r.get("peptide") for r in rows if isinstance(r, dict)}
        if not found and not cache:
            die(f"unipept: the API returned no LCA for any of {len(todo)} "
                "peptides. Not writing a cache (an empty one poisons every "
                "later run). Check unipept.api_url and your network, or run "
                "the CLI yourself and set unipept.result.")
        missing = [x for x in todo if x not in found]
        if missing:
            log(f"unipept: {len(missing)}/{len(todo)} queried peptides have no "
                "LCA; cached as misses so reruns do not re-query them")
            rows += [{"peptide": x, "taxon_id": ""} for x in missing]
        new = pd.DataFrame(rows)
        allrows = pd.concat([pd.DataFrame(list(cache.values())), new],
                            ignore_index=True) if cache else new
        with atomic_out(p.unipept_cache) as tmp:
            allrows.to_csv(tmp, sep="\t", index=False)
    df = read_unipept_result(p.unipept_cache)
    with atomic_out(p.unipept_lca) as tmp:
        df.to_csv(tmp, sep="\t", index=False)
    log(f"unipept: {len(df)} peptide assignments -> {p.unipept_lca}")


def consensus_taxon(lineages, min_frac, min_peptides):
    """Deepest rank at which a majority of a protein's peptides agree.

    Not the LCA of the LCAs: one spurious peptide would drag that to the root.
    """
    # A peptide whose LCA is root carries no rank at all, and one that stops at
    # phylum says nothing about genus. Counting those in the denominator made
    # every conserved (well-quantified) protein lose its taxon, so the vote at
    # each rank is taken only over the peptides that are resolved that deep.
    lineages = [l for l in lineages if l]
    n = len(lineages)
    if n < min_peptides:
        return None, None, 0.0, n
    for i in range(len(RANKS) - 1, -1, -1):
        r = RANKS[i]
        deep = RANKS[i:]
        informative = [l for l in lineages if any(l.get(rk) for rk in deep)]
        vals = [l.get(r) for l in informative if l.get(r)]
        if not vals or len(informative) < min_peptides:
            continue
        # Counter is O(n); the previous list.count inside a loop over set(vals)
        # was O(n^2) for proteins with many peptides. Sorted for a
        # deterministic tie-break — set order varies with PYTHONHASHSEED.
        from collections import Counter
        counts = sorted(((c, v) for v, c in Counter(vals).items()),
                        key=lambda x: (-x[0], x[1]))
        cnt, top = counts[0]
        # A 2-vs-2 split used to pass min_fraction 0.5 and was resolved by
        # string order. Require a strict winner and fall back to the shared
        # parent rank instead.
        if len(counts) > 1 and counts[1][0] == cnt:
            continue
        d = len(informative)
        if cnt / d >= min_frac:
            return top, r, cnt / d, n
    return None, None, 0.0, n


def stage_taxonomy(cfg, p):
    if not os.path.exists(p.unipept_lca):
        die(f"no Unipept result at {p.unipept_lca}; run the unipept stage first")
    lca = read_unipept_result(p.unipept_lca)
    rank_cols = [rk for rk in RANKS if f"{rk}_id" in lca.columns]
    peps = lca["peptide"].astype(str).tolist()
    cols = {rk: lca[f"{rk}_id"].astype("object").where(
        lca[f"{rk}_id"].notna()).tolist() for rk in rank_cols}
    lin_of = {}
    for i, pep in enumerate(peps):
        d = {}
        for rk in rank_cols:
            v = cols[rk][i]
            if v is not None and v == v:
                d[rk] = str(v).split(".")[0]
        lin_of[pep] = d

    fmt = cfg["quant_format"]
    if fmt not in FEATURE_FORMATS:
        die("the taxonomy comparison needs peptide-level input")
    ann = pd.read_csv(p.final, sep="\t", dtype=str, low_memory=False, encoding="utf-8", encoding_errors="replace")
    ann = ann.rename(columns={ann.columns[0]: "protein_id"})
    # dtype=str still leaves float NaN in empty cells, and str(nan) is the
    # truthy string 'nan'. Without this every protein with no eggNOG hit — most
    # of a Prokka/smORF database — was scored 'conflict' rather than
    # 'eggnog_missing', and 'nan' turned up in the unresolved-taxid warning.
    seed = ann["seed_taxid"] if "seed_taxid" in ann.columns else pd.Series(dtype=str)
    taxon_of = dict(zip(ann["protein_id"], seed.fillna("").astype(str)))

    feats = peptide_features(cfg, "taxonomy")
    u = cfg.get("unipept") or {}
    split = u.get("split_missed_cleavages", True)

    # One vote per source peptide, not per tryptic fragment and not per charge
    # state: an ion table repeats a peptide once per charge/modification, and
    # a single missed-cleavage peptide yields two fragments, which was enough
    # to satisfy consensus_min_peptides=2 on the strength of one peptide.
    per_protein = defaultdict(dict)
    n_pep, n_hit = 0, 0
    for pep_raw, cands in zip(feats["peptide"], feats["candidates"]):
        pep = strip_modifications(pep_raw)
        n_pep += 1
        frags = tryptic_split(pep) if split else [pep]
        lins = [lin_of[f] for f in frags if f in lin_of]
        lins = [l for l in lins if l]
        if not lins:
            continue
        n_hit += 1
        # Collapse the fragments of one peptide into the deepest lineage they
        # all agree on; a fragment that resolves less deeply simply stops.
        merged = {}
        for rk in RANKS:
            vals = {l[rk] for l in lins if l.get(rk)}
            if len(vals) != 1:
                break
            merged[rk] = vals.pop()
        if not merged:
            continue
        for c in dict.fromkeys(cands):
            if c:
                per_protein[c][pep] = merged
    frac_hit = n_hit / n_pep if n_pep else 0.0
    msg = (f"taxonomy: {n_hit}/{n_pep} feature rows ({frac_hit:.1%}) have a "
           f"Unipept LCA, over {len(per_protein)} proteins")
    log(msg, "WARN" if frac_hit < 0.1 else "INFO")
    if frac_hit < 0.1:
        log("taxonomy: fewer than 10% of peptides matched the LCA table. Check "
            "that unipept.split_missed_cleavages matches how the peptide list "
            "was produced, and that the pept2lca file is for THIS dataset.",
            "WARN")

    minf = float(u.get("consensus_min_fraction", 0.5))
    minp = int(u.get("consensus_min_peptides", 2))
    rows = []
    for prot, by_pep in per_protein.items():
        tid, rank, frac, n = consensus_taxon(list(by_pep.values()), minf, minp)
        # '' rather than None: pandas >= 3 turns a None in a string column into
        # NaN, and the `uni is None` test below then never fired, so every
        # protein without a consensus was reported as a taxonomic conflict.
        rows.append({"protein_id": prot, "unipept_taxid": tid or "",
                     "unipept_rank": rank or "",
                     "unipept_agreement": round(frac, 3),
                     "n_peptides_with_lca": n})
    prot_tax = pd.DataFrame(rows)
    os.makedirs(f"{p.R}/unipept", exist_ok=True)
    with atomic_out(p.protein_taxonomy) as tmp:
        prot_tax.to_csv(tmp, sep="\t", index=False)
    log(f"taxonomy: consensus for {len(prot_tax)} proteins -> {p.protein_taxonomy}")

    # ---- compare with eggNOG ------------------------------------------
    tdir = (cfg.get("db") or {}).get("ncbi_taxonomy", "")
    tax = None
    if tdir and os.path.exists(os.path.join(tdir, "nodes.dmp")):
        tax = NCBITaxonomy(tdir)
    else:
        log("db.ncbi_taxonomy not set (needs nodes.dmp/names.dmp), so eggNOG "
            "seed taxids cannot be resolved to a lineage; the comparison is "
            "limited to exact taxid identity", "WARN")

    def _taxid(v):
        """A taxid, or '' for every flavour of missing. Both sides arrive as
        str, None or float NaN depending on pandas version and column dtype."""
        if v is None or (isinstance(v, float) and v != v):
            return ""
        s = str(v).strip()
        # '742765.0' happens when a taxid has been through a float column.
        if s.endswith(".0") and s[:-2].isdigit():
            s = s[:-2]
        return "" if s.lower() in ("", "nan", "none", "na") else s

    comp = []
    for _, r in prot_tax.iterrows():
        pid = r["protein_id"]
        egg = _taxid(taxon_of.get(pid, ""))
        uni = _taxid(r["unipept_taxid"])
        row = {"protein_id": pid, "eggnog_taxid": egg, "unipept_taxid": uni,
               "unipept_rank": r["unipept_rank"],
               "unipept_agreement": r["unipept_agreement"],
               "n_peptides_with_lca": r["n_peptides_with_lca"]}
        # Missing is not disagreement, and two missings are certainly not
        # "identical" — which is what str(nan) == str(nan) used to produce.
        if not egg or not uni:
            row["verdict"] = ("eggnog_missing" if not egg else
                              "unipept_missing")
        elif tax is None:
            row["verdict"] = "identical" if egg == uni else "differ_no_lineage"
        elif not tax.current(egg):
            # A seed taxid that is neither current nor merged (eggNOG 5 ships a
            # 2018 taxonomy) is an unusable reference, not a contradiction.
            tax.unresolved.add(egg)
            row["deepest_agreement"] = "none"
            row["verdict"] = "eggnog_unresolved"
        else:
            el, ul = tax.lineage(egg), tax.lineage(uni)
            deepest, differs = None, False
            for rk in RANKS:
                if el.get(rk) and ul.get(rk):
                    row[f"agree_{rk}"] = el[rk] == ul[rk]
                    if el[rk] == ul[rk]:
                        deepest = rk
                    else:
                        differs = True
            row["deepest_agreement"] = deepest or "none"
            row["eggnog_name"] = tax.sci_name(egg)
            row["unipept_name"] = tax.sci_name(uni)
            # 'conflict' now means a rank defined on BOTH sides actually
            # disagrees. Unipept's LCA stopping at family is a limit of the
            # evidence, not a contradiction, and used to be counted as one.
            row["verdict"] = ("identical" if egg == uni else
                              "conflict" if differs else
                              "concordant" if deepest in ("genus", "species") else
                              "concordant_above_genus" if deepest else
                              "no_common_rank")
        comp.append(row)
    comp = pd.DataFrame(comp)
    if tax is not None and tax.unresolved:
        log(f"{len(tax.unresolved)} taxid(s) resolve to nothing even after "
            f"merged.dmp, e.g. {sorted(tax.unresolved)[:5]}; they were deleted "
            "from NCBI Taxonomy (eggNOG 5 seeds are from a 2018 taxdump), so "
            "those proteins get verdict 'eggnog_unresolved' rather than a "
            "lineage comparison. Adding merged.dmp/delnodes.dmp to "
            "db.ncbi_taxonomy is what recovers them, not a newer taxdump.",
            "WARN")
    with atomic_out(p.taxonomy_comparison) as tmp:
        comp.to_csv(tmp, sep="\t", index=False)
    log("taxonomy: eggNOG vs Unipept verdicts:\n" +
        comp["verdict"].value_counts().to_string())
    log(f"taxonomy: wrote {p.taxonomy_comparison}")


FEATURE_FORMATS = {"fragpipe_peptide", "fragpipe_ion", "fragpipe_tmt",
                   "msstats_csv",
                   "msstats_feature"}
PROTEIN_FORMATS = {"diann", "fragpipe", "msstats_protein"}
ALL_FORMATS = FEATURE_FORMATS | PROTEIN_FORMATS


def join_cols(df, cols, sep="_"):
    """Row-wise string join that survives missing values.

    pandas >= 3 gives astype(str) a StringDtype whose NA is a real float nan,
    so "_".join(row) raises TypeError on any column with a missing value —
    FragmentIon and ProductCharge in an MSstats export are usually all-NA.
    """
    parts = [df[c].astype(str).fillna("").tolist() for c in cols]
    return pd.Series(["_".join(t) for t in zip(*parts)], index=df.index)


def first_token(v):
    """The identifier a fasta header would yield for this value.

    read_fasta keys every protein on the first whitespace-delimited token, so
    an id carrying a description ("CDPNAMPK_339076 hypothetical protein",
    which is what FragPipe writes into 'Protein ID' for a non-UniProt
    database) has to be reduced the same way or nothing joins.
    """
    s = str(v).strip()
    return s.split()[0] if s else ""


def split_ids(v):
    if not isinstance(v, str) or not v or v.lower() == "nan":
        return []
    return [first_token(x) for x in re.split(r"[;,]", v) if x.strip()]


def read_delim_table(path, **kw):
    """pd.read_csv with the delimiter taken from the header line only.

    sep=None forces pandas onto the pure-Python parser. On a 400k-row
    combined_peptide.tsv that is ~18x slower and needs 4-6x the memory of the
    C parser for no benefit whatever: the first line already says which
    delimiter this is.
    """
    with open(path, "r", newline="", errors="replace", encoding="utf-8") as fh:
        head = fh.readline()
    sep = max(("\t", ",", ";", "|"), key=head.count)
    if head.count(sep) == 0:
        sep = "," if path.lower().endswith(".csv") else "\t"
    # setdefault, not a literal keyword: **kw is the caller's, and passing
    # encoding twice is a TypeError rather than a preference.
    kw.setdefault("encoding", "utf-8")
    kw.setdefault("encoding_errors", "replace")
    return pd.read_csv(path, sep=sep, low_memory=False, **kw)


def excluded_prefixes(cfg):
    """Decoy/contaminant id prefixes, as a tuple ready for str.startswith."""
    return tuple(cfg.get("exclude_id_prefixes")
                 or ["rev_", "decoy_", "contam_", "Cont_", "CON__"])


# ======================================================================
# FragPipe TMT (isobaric): the per-plex tables
# ======================================================================
# Read TMTn/{ion,peptide}.tsv, and deliberately NOT tmt-report/. The
# tmt-report matrices are already log2 and median-centred, so the log2 step
# downstream would take the log of a log; they are protein level, which
# deletes the shared-peptide rule, peptide_assignment, peptide_evidence.tsv
# and the peptide assay of the R object — the layer this tool exists for
# against a strain-redundant metagenome database; and they carry
# TMT-Integrator's own protein inference, which is exactly the inference such
# a database makes least trustworthy. Per-plex reporter intensities are
# LINEAR, so the existing log2 path, the roll-up and the median-of-ratios
# size factor apply to them unchanged.

TMT_LEVEL_FILES = {"ion": "ion.tsv", "peptide": "peptide.tsv"}


def natural_key(s):
    """Sort key that puts TMT10 after TMT9 rather than after TMT1."""
    return [(1, int(t)) if t.isdigit() else (0, t.lower())
            for t in re.split(r"(\d+)", str(s)) if t != ""]


def header_columns(path):
    """Column names from the header line alone.

    "What is this file?" has to be answerable without parsing its body: the
    TMT flavour of msstats.csv carries unquoted commas in Protein.Description
    and dies in the C parser with a tokenising error that names neither TMT
    nor the file, so the recogniser has to run before pandas does.
    """
    with open(path, "r", newline="", errors="replace", encoding="utf-8") as fh:
        head = fh.readline()
    sep = max(("\t", ",", ";", "|"), key=head.count)
    if head.count(sep) == 0:
        sep = "," if path.lower().endswith(".csv") else "\t"
    return [c.strip().strip('"') for c in head.rstrip("\r\n").split(sep)]


def refuse_isobaric_matrix(path, cols):
    """Refuse the two TMT files that no reader here can honestly read.

    Both are quantitative, both look plausible, and both are wrong in a way
    that produces numbers instead of an error: the tmt-report matrices are
    already log2 and median-centred (this tool would log them a second time)
    and carry TMT-Integrator's protein inference, and the TMT msstats.csv
    holds one row per PSM with the channels in 'Channel <mass>' columns, which
    no format here maps to samples. Naming what was found is the point: the
    user reached for the file that looked most like a matrix.
    """
    cols = [str(c) for c in cols]
    if "ReferenceIntensity" in cols:
        die(f"{path}: this is a TMT-Integrator tmt-report matrix (it has a "
            "'ReferenceIntensity' column). Its values are already log2 and "
            "median-centred, so quantifying them here would log-transform "
            "them a second time; it is protein level, so the shared-peptide "
            "rule, peptide_assignment and peptide_evidence.tsv have nothing "
            "to work on; and its protein inference is TMT-Integrator's, which "
            "is the inference a strain-redundant metagenome database makes "
            "least trustworthy. Set quant_format: fragpipe_tmt and point "
            "quant_table at the run directory holding the per-plex TMTn/ "
            "folders, whose reporter intensities are linear.")
    chan = [c for c in cols if re.fullmatch(r"Channel[ _.]\S+", c)]
    if chan:
        die(f"{path}: this is the TMT flavour of msstats.csv — one row per "
            f"PSM, with {len(chan)} reporter channel(s) in columns named "
            f"like {chan[:3]}, which carry the label mass and not a sample "
            "name. quant_format 'msstats_csv' expects the label-free export "
            "(a single 'Intensity' column per run) and would either fail to "
            "parse this file or quantify the wrong column. Set "
            "quant_format: fragpipe_tmt and point quant_table at the run "
            "directory holding the per-plex TMTn/ folders.")


def refuse_per_plex_reporter_table(path, cols, token="Intensity"):
    """Refuse a PROTEIN-level table whose intensities are reporter channels.

    Only for the protein-level formats. A per-plex TMTn/protein.tsv is the one
    isobaric file nothing else catches: it has no ReferenceIntensity and no
    'Channel <mass>' column, and the protein-level column detector takes every
    numeric column that is not declared metadata — so it quantifies ONE plex
    as if it were the experiment, with 'Length', 'Protein Qvalue' and 'Razor
    Intensity' sitting in the matrix beside the channels, and says nothing.

    The signature is the PREFIX form: FragPipe names a reporter column
    '<token> <sample>' and a label-free column '<sample> <token>', so a
    combined_protein.tsv can never match this and label-free input is
    untouched.
    """
    rep = [str(c) for c in cols
           if re.fullmatch(re.escape(token) + r" \S.*", str(c))]
    if not rep:
        return
    die(f"{path}: this is a per-plex FragPipe TMT table, not a label-free "
        f"protein table. {len(rep)} reporter-ion column(s) named "
        f"'{token} <sample>' are present, e.g. {rep[:3]}, which is how "
        "FragPipe names an isobaric channel. Quantifying it at protein level "
        "would report a SINGLE plex as the whole experiment, sweep the "
        f"numeric metadata beside it (Length, Protein Qvalue, Razor {token}) "
        "into the matrix as if those were samples, and lose the peptide layer "
        "the shared-peptide rule needs. Set quant_format: fragpipe_tmt and "
        "point quant_table at the run directory holding the per-plex TMTn/ "
        "folders; it reads every plex and joins them at feature level.")


def read_tmt_annotation(path, plex):
    """-> [(channel, sample), ...] in file order.

    FragPipe writes '<channel> <sample>' per line, e.g. '131C Pool01'. The
    sample name is what the reporter COLUMNS are named after, so it, not the
    channel, is the identifier the rest of this reader keys on.
    """
    rows = []
    with opener(path) as fh:
        for lineno, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split(None, 1)
            if len(parts) < 2:
                die(f"{path}:{lineno}: '{s[:60]}' names a channel but no "
                    "sample. A FragPipe TMT annotation is '<channel> "
                    "<sample>' per line, e.g. '131C Pool01'.")
            rows.append((parts[0], parts[1].strip()))
    if not rows:
        die(f"{path}: empty; expected one '<channel> <sample>' line per "
            f"channel of plex {plex}")
    for what, vals in (("channel", [c for c, _ in rows]),
                       ("sample name", [s for _, s in rows])):
        dup = sorted({v for v in vals if vals.count(v) > 1})
        if dup:
            die(f"{path}: {what} {dup} appears more than once, so the "
                f"channels of plex {plex} cannot be mapped to samples")
    return rows


def tmt_plex_dirs(root, cfg):
    """-> [(plex, directory), ...], naturally sorted."""
    t = cfg.get("tmt") or {}
    pat = str(t.get("plex_glob") or "TMT*")
    hits = glob.glob(pat if os.path.isabs(pat) else os.path.join(root, pat))
    # fnmatchcase on top of glob: glob follows the FILESYSTEM's case rules, so
    # on Windows (and on a case-insensitive mac volume) 'TMT*' also matched
    # the sibling 'tmt-report' directory and the run died on the plex that
    # never existed. The pattern has to mean the same thing on every platform.
    base = lambda h: os.path.basename(h.rstrip("/\\"))
    dirs = sorted((h for h in hits
                   if os.path.isdir(h)
                   and fnmatch.fnmatchcase(base(h), os.path.basename(pat))),
                  key=lambda h: natural_key(base(h)))
    if not dirs:
        here = sorted(os.listdir(root))[:12] if os.path.isdir(root) else []
        die(f"no plex directory matches tmt.plex_glob '{pat}' under {root}. "
            f"{root} contains: {here}. quant_format 'fragpipe_tmt' reads the "
            "per-plex FragPipe output directories (TMT1/, TMT2/, ...), so "
            "quant_table must be the run directory that holds them.")
    return [(base(d), d) for d in dirs]


def tmt_annotation_path(plex, pdir, cfg):
    """The annotation file of one plex.

    FragPipe names it <PLEX>_annotation.txt (TMT1/TMT1_annotation.txt) and
    never a plain annotation.txt, which is what a first guess reaches for.
    """
    t = cfg.get("tmt") or {}
    ann = t.get("annotation") or "{plex}_annotation.txt"
    if isinstance(ann, dict):
        # An explicit map is an explicit statement: a plex missing from it is
        # a mistake to report, not a reason to guess the pattern back.
        if plex not in ann:
            die(f"tmt.annotation is a map and has no entry for plex '{plex}' "
                f"(it lists {sorted(ann)}). Add one, or use the "
                "'{plex}_annotation.txt' pattern form.")
        cand = str(ann[plex])
    else:
        cand = str(ann).replace("{plex}", plex)
    # Relative to the PLEX DIRECTORY, which is what the README documents, and
    # not to the process's working directory. Resolving against the cwd made
    # the file that was read depend on where the run was started from: a
    # 'ann/{plex}.txt' pattern, or any map entry, picked up a same-named file
    # beside the shell rather than the plex's own annotation, and mapping the
    # reporter columns through the wrong plex's annotation is silent — the
    # channels line up, the sample names do not. An absolute path is still
    # taken exactly as written.
    rel = cand
    if not os.path.isabs(cand):
        cand = os.path.join(pdir, cand)
    if not os.path.exists(cand):
        near = sorted(f for f in os.listdir(pdir) if "annotation" in f.lower())
        die(f"plex {plex} has no annotation file: {cand} does not exist. "
            f"{pdir} contains {near or 'no annotation-like file'}. FragPipe "
            "writes '<PLEX>_annotation.txt', not 'annotation.txt'; set "
            "tmt.annotation to the right pattern, or to a {plex: path} map. "
            + (f"'{rel}' does exist relative to the working directory "
               f"{os.getcwd()}, but tmt.annotation is resolved against the "
               "plex directory; give an absolute path if that is the file "
               "you mean. " if rel != cand and os.path.exists(rel) else "")
            + "Without it the reporter columns cannot be mapped to samples.")
    return cand


def quant_inputs(cfg):
    """The quant file(s) a stage signature must digest: [quant_table] for
    every format but fragpipe_tmt, whose input is a directory of them.

    The signature digests a directory by its size and mtime, and neither
    changes when FragPipe rewrites a table inside it — so listing
    quant_table would leave the join stage cached across a re-search. Listing
    the tables and annotations themselves makes the cache honest. Tolerant by
    design: this runs before any stage, on a config that may not point
    anywhere yet.
    """
    root = cfg.get("quant_table") or ""
    if cfg.get("quant_format") != "fragpipe_tmt":
        return [root]
    t = cfg.get("tmt") or {}
    fname = TMT_LEVEL_FILES.get(str(t.get("level") or "ion").lower(), "ion.tsv")
    try:
        # Inside the try, and ValueError caught below: this float() ran before
        # read_fragpipe_tmt could validate anything, so a min_purity written
        # as '90%' came out of the signature as a bare traceback instead of
        # the reader's own message naming the key and its range. The signature
        # is not the place to judge the value; falling back to the directory
        # leaves the reader to die properly.
        #
        # psm.tsv is an input only when min_purity actually reads it; listing
        # it unconditionally would invalidate every cached join the moment
        # FragPipe rewrote a file the run never opened.
        want_psm = float(t.get("min_purity") or 0) > 0
        out = []
        for plex, pdir in tmt_plex_dirs(root, cfg):
            out.append(os.path.join(pdir, fname))
            if want_psm:
                out.append(os.path.join(pdir, "psm.tsv"))
            with contextlib.suppress(StageError):
                out.append(tmt_annotation_path(plex, pdir, cfg))
        return sorted(out)
    except (StageError, OSError, TypeError, ValueError):
        return [root]


def _tmt_map_reporter_columns(path, df, ann, token):
    """-> (columns, sample names, 'sample name'|'channel').

    FragPipe names the reporter columns after the ANNOTATED SAMPLE
    ('Intensity Pool01'), which is the prefix form the label-free suffix rule
    cannot see. A run annotated after the fact can still carry the channel
    ('Intensity 131C'), so both are accepted — but which one was used is
    logged, because the sample names in the output come from it.
    """
    prefix = token + " "
    present = [c for c in df.columns if str(c).startswith(prefix)]
    for how, keys in (("sample name", [s for _, s in ann]),
                      ("channel", [c for c, _ in ann])):
        want = [prefix + k for k in keys]
        if all(w in df.columns for w in want):
            extra = [c for c in present if c not in want]
            if extra:
                die(f"{path}: {len(extra)} reporter column(s) {extra[:4]} are "
                    f"not in the annotation, which lists {len(keys)} "
                    f"{how}(s) {keys[:4]}. The annotation does not describe "
                    "this file, so its channels cannot be mapped to samples.")
            return want, [s for _, s in ann], how
    die(f"{path}: the reporter columns cannot be mapped to samples. The "
        f"annotation lists {len(ann)} channels "
        f"{[f'{c}={s}' for c, s in ann][:4]}, and the columns starting "
        f"'{prefix}' are {present[:6] or 'none at all'}. FragPipe names them "
        f"'{token} <sample>' after the annotated sample; check that the "
        "annotation belongs to this plex.")
    return None, None, None                       # unreachable; die() raises


# psm.tsv names the same three things as ion.tsv, with two different column
# names. Kept as a map rather than hard-coded so the key the purity join uses
# is literally the key the feature ids were built from.
_TMT_PSM_COLUMNS = {"Peptide Sequence": "Peptide",
                    "Modified Sequence": "Modified Peptide"}


def tmt_psm_purity(pdir, plex, key):
    """-> Series feature_id -> MEDIAN precursor purity of its PSMs.

    Purity exists only in psm.tsv (ion.tsv, peptide.tsv and protein.tsv have
    no such column), so a purity filter at feature level is a join, and the
    median is the aggregate that fits what the join produces: a feature's
    reporter intensities are a SUM over its PSMs, so no single PSM's purity
    describes it, and the median says whether the typical contributing
    spectrum was clean. It is an approximation of the per-PSM filter
    TMT-Integrator would apply before summarising — FragPipe has already
    summed by the time this reader sees the file — and it is reported as one.
    """
    path = os.path.join(pdir, "psm.tsv")
    if not os.path.exists(path):
        die(f"tmt.min_purity is set, but plex {plex} has no psm.tsv: {path} "
            "does not exist. Purity is written only into psm.tsv, so the "
            "filter cannot be applied without it. Remove tmt.min_purity, or "
            "point quant_table at a run directory that still has its PSM "
            "tables.")
    df = read_delim_table(path)
    if "Purity" not in df.columns:
        die(f"{path}: no 'Purity' column, so tmt.min_purity cannot be "
            f"applied. It has: {list(df.columns)[:12]}...")
    cols = []
    for k in key:
        c = _TMT_PSM_COLUMNS.get(k, k)
        if c not in df.columns:
            die(f"{path}: no '{c}' column, so its purities cannot be keyed "
                f"onto the {k} of the feature table. psm.tsv has: "
                f"{list(df.columns)[:12]}...")
        cols.append(c)
    pur = pd.to_numeric(df["Purity"], errors="coerce")
    return pur.groupby(join_cols(df, cols)).median()


def read_fragpipe_tmt(root, cfg):
    """-> (features, int_cols, design), from the per-plex FragPipe TMT output.

    Exactly the shape read_feature_table returns for a label-free table —
    feature_id, peptide, razor_protein, candidates, one column per sample —
    so the roll-up, the shared-peptide rule and everything downstream are
    untouched. design carries a plex column as well.
    """
    t = cfg.get("tmt") or {}
    level = str(t.get("level") or "ion").lower()
    if level not in TMT_LEVEL_FILES:
        die(f"tmt.level must be 'ion' or 'peptide', got '{level}'")
    fname = TMT_LEVEL_FILES[level]
    try:
        min_purity = float(t.get("min_purity") or 0)
    except (TypeError, ValueError):
        die(f"tmt.min_purity must be a number between 0 and 1, got "
            f"{t.get('min_purity')!r}")
    if not 0 <= min_purity <= 1:
        die(f"tmt.min_purity must be between 0 and 1 (FragPipe's Purity is a "
            f"fraction), got {min_purity}")
    raw_norm = t.get("within_plex_normalise", "median")
    norm = str("median" if raw_norm is None else raw_norm).lower()
    if norm not in ("median", "none"):
        die(f"tmt.within_plex_normalise must be 'median' or 'none', got "
            f"{raw_norm!r}")
    if os.path.isfile(root):
        refuse_isobaric_matrix(root, header_columns(root))
        die(f"quant_format 'fragpipe_tmt' reads the per-plex FragPipe "
            f"directories, so quant_table must be the run directory that "
            f"holds them (the one with {t.get('plex_glob') or 'TMT*'} "
            f"subdirectories), not the single file {root}.")
    if not os.path.isdir(root):
        die(f"quant_table not found: {root}. quant_format 'fragpipe_tmt' "
            "expects the FragPipe run directory holding the per-plex TMTn/ "
            "folders.")
    if cfg.get("manifest"):
        # The TMT manifest lists LC-MS RUNS, and its experiment column is the
        # PLEX. Renaming channels from it is impossible and taking its
        # experiment as the condition would be inferring the condition from
        # the plex, which is exactly the wrong answer.
        log("tmt: `manifest` is ignored for quant_format 'fragpipe_tmt'. A "
            "FragPipe TMT manifest names LC-MS runs, not reporter channels, "
            "and its experiment column is the plex — using it as a condition "
            "would infer the condition from the batch. Sample names come from "
            "each plex's annotation file.", "WARN")

    plexes = tmt_plex_dirs(root, cfg)
    pref = excluded_prefixes(cfg)
    # The same word FragPipe uses in the label-free tables, used here as a
    # PREFIX ("Intensity Pool01") instead of a suffix ("Pool01 Intensity").
    token = cfg.get("feature_intensity_suffix", "Intensity")
    ref_name = str(t.get("reference_name") or "")
    ref_chan = str(t.get("reference_channel") or "")
    if ref_name and ref_chan:
        die(f"tmt.reference_name ('{ref_name}') and tmt.reference_channel "
            f"('{ref_chan}') are both set, and they can disagree per plex. "
            "Set one: the sample name is the stable signal when the reference "
            "moves between channels.")
    use_ratios = bool(t.get("use_reference_ratios", False))
    if use_ratios and not (ref_name or ref_chan):
        die("tmt.use_reference_ratios is on but no reference is named. Set "
            "tmt.reference_name (a glob on the annotated sample name, e.g. "
            "'Pool*') or tmt.reference_channel (e.g. '131C').")
    drop_empty = bool(t.get("drop_empty_channels", True))
    zero_missing = bool(cfg.get("zero_intensity_is_missing", True))

    meta = {}          # feature_id -> dict(peptide, razor, candidates, plexes)
    frames, design_rows, owner = [], [], {}
    int_cols, n_zero, n_rows = [], 0, 0
    refs = []                       # (plex, channel, sample) of the reference
    ratio_had = ratio_lost = 0      # values before / values lost to the divide
    pur_drop = pur_unknown = pur_seen = 0    # tmt.min_purity bookkeeping
    norm_lo, norm_hi = [], []       # log2 of the within-plex channel factors
    for plex, pdir in plexes:
        path = os.path.join(pdir, fname)
        if not os.path.exists(path):
            # No fallback to tmt-report/: a missing per-plex table means this
            # run is not the one being described, and quantifying a different
            # file to fill the hole would be undetectable in the output.
            die(f"plex {plex} has no {fname}: {path} does not exist "
                f"(tmt.level is '{level}'). {pdir} contains "
                f"{sorted(f for f in os.listdir(pdir) if f.endswith('.tsv'))}.")
        refuse_isobaric_matrix(path, header_columns(path))
        ann = read_tmt_annotation(tmt_annotation_path(plex, pdir, cfg), plex)
        df = read_delim_table(path)
        cols, samples, how = _tmt_map_reporter_columns(path, df, ann, token)
        nonnum = [c for c in cols if not pd.api.types.is_numeric_dtype(df[c])]
        if nonnum:
            die(f"{path}: reporter column(s) {nonnum[:4]} are not numeric, so "
                "they are not reporter intensities")

        # ---- identity: same rule as the label-free reader -------------
        prot = "Protein" if "Protein" in df.columns else "Protein ID"
        if prot not in df.columns:
            die(f"{path}: no 'Protein' or 'Protein ID' column; is this a "
                "FragPipe per-plex ion.tsv/peptide.tsv?")
        bad = pd.Series(False, index=df.index)
        for c in ("Is Decoy", "Is Contaminant"):
            if c in df.columns:
                bad |= df[c].astype(str).str.lower().isin(["true", "1", "yes"])
        if pref:
            bad |= df[prot].astype(str).map(first_token).str.startswith(pref)
        if bool(bad.any()):
            log(f"tmt {plex}: {int(bad.sum())} decoy/contaminant row(s) "
                f"dropped from {path}", "WARN")
            df = df.loc[~bad].reset_index(drop=True)

        if level == "ion":
            # All three, not "modified sequence or sequence": the feature id
            # has to be comparable ACROSS plexes, and only ~45% of ion keys
            # are shared between two plexes, so the outer join is only as
            # good as this key.
            key = [c for c in ("Peptide Sequence", "Modified Sequence",
                               "Charge") if c in df.columns]
        else:
            key = [c for c in ("Peptide Sequence", "Peptide")
                   if c in df.columns][:1]
        if not key:
            die(f"{path}: no peptide sequence column")
        fid = join_cols(df, key)
        pep_col = [c for c in ("Peptide Sequence", "Peptide")
                   if c in df.columns][:1]
        pep = (df[pep_col[0]].astype(str).fillna("") if pep_col
               else pd.Series([""] * len(df), index=df.index))
        # FragPipe leaves 'Peptide Sequence' empty on some ion rows while
        # 'Modified Sequence' is filled. Those rows still identify a peptide,
        # and the taxonomy stages read this column, so recover it rather than
        # sending an empty string to Unipept.
        if "Modified Sequence" in df.columns:
            empty = pep.eq("") | pep.eq("nan")
            if bool(empty.any()):
                ms = df["Modified Sequence"].astype(str).fillna("")
                pep = pep.mask(empty, ms.map(strip_modifications))
                log(f"tmt {plex}: {int(empty.sum())} row(s) have no "
                    "'Peptide Sequence'; the sequence was recovered from "
                    "'Modified Sequence'", "WARN")
        # A row with neither a sequence nor a modified sequence identifies
        # nothing, and every such row seen so far is all-zero filler. Left in,
        # they all collapse onto one feature id and merge into each other.
        blank = fid.str.replace("_", "", regex=False).str.strip().eq("")
        if bool(blank.any()):
            log(f"tmt {plex}: {int(blank.sum())} row(s) carry no peptide "
                "identity at all (no sequence and no modified sequence) and "
                "were dropped", "WARN")
            keep = ~blank
            df, fid, pep = (df.loc[keep].reset_index(drop=True),
                            fid[keep].reset_index(drop=True),
                            pep[keep].reset_index(drop=True))

        # ---- co-isolation: the only place purity exists is psm.tsv -----
        if min_purity > 0:
            med = fid.map(tmt_psm_purity(pdir, plex, key))
            unknown = med.isna()
            low = med.lt(min_purity).fillna(False)
            pur_seen += len(fid)
            pur_drop += int(low.sum())
            pur_unknown += int(unknown.sum())
            if bool(unknown.any()):
                # Kept, not dropped: a feature no PSM row matches is a join
                # failure (FragPipe leaves 'Modified Peptide' empty on rows
                # ion.tsv writes a modified sequence for), and deleting IDs
                # for that would look exactly like a purity filter working.
                log(f"tmt {plex}: {int(unknown.sum())} of {len(fid)} feature(s) "
                    "match no row of psm.tsv, so their purity is unknown; "
                    "they are KEPT — an unmatched key is a join failure, not "
                    "a co-isolated precursor", "WARN")
            log(f"tmt {plex}: min_purity={min_purity} drops "
                f"{int(low.sum())} of {len(fid)} feature(s) whose MEDIAN PSM "
                "purity is below it")
            keep = ~low
            df, fid, pep = (df.loc[keep].reset_index(drop=True),
                            fid[keep].reset_index(drop=True),
                            pep[keep].reset_index(drop=True))

        razor = df[prot].astype(str).map(first_token)
        mapped = (df["Mapped Proteins"] if "Mapped Proteins" in df.columns
                  else pd.Series([""] * len(df), index=df.index))
        if "Mapped Proteins" not in df.columns:
            log(f"tmt {plex}: no 'Mapped Proteins' column, so every feature "
                "is treated as unique to its razor protein and shared-peptide "
                "filtering is inactive", "WARN")

        vals = df[cols].copy()
        vals.columns = samples
        if zero_missing:
            nz = int((vals == 0).sum().sum())
            n_zero += nz
            n_rows += vals.size
            vals = vals.where(vals != 0)

        # ---- one row per feature within the plex ----------------------
        vals.index = pd.Index(fid, name="feature_id")
        dup = int(fid.duplicated().sum())
        if dup:
            log(f"tmt {plex}: {dup} row(s) repeat a feature id and were "
                "summed; the reindex the outer join needs cannot carry a "
                "duplicated key", "WARN")
            vals = vals.groupby(level=0, sort=False).sum(min_count=1)
        for f, r, m, pp in zip(fid, razor, mapped, pep):
            rec = meta.get(f)
            cand = [r] + [x for x in split_ids(m) if not x.startswith(pref)]
            if rec is None:
                meta[f] = {"peptide": pp, "razor": r,
                           "cand": dict.fromkeys(cand),
                           "razors": {r: 1}, "plexes": {plex: 1}}
            else:
                rec["cand"].update(dict.fromkeys(cand))
                rec["razors"][r] = 1
                rec["plexes"][plex] = 1

        # ---- channels that are not samples ----------------------------
        # FragPipe names an unassigned channel <PLEX>_<CHANNEL> in the
        # annotation. Its signal is isotope carry-over, not a sample.
        empty_ch = [(c, s) for c, s in ann if s == f"{plex}_{c}"]
        keep_samples = [s for _, s in ann]
        if empty_ch and drop_empty:
            log(f"tmt {plex}: {len(empty_ch)} channel(s) "
                f"{[c for c, _ in empty_ch]} carry the placeholder name "
                f"{[s for _, s in empty_ch]}, which is how FragPipe writes an "
                "unassigned channel; dropped (set tmt.drop_empty_channels "
                "false to keep them)", "WARN")
            keep_samples = [s for s in keep_samples
                            if s not in {s2 for _, s2 in empty_ch}]
        elif empty_ch:
            log(f"tmt {plex}: {len(empty_ch)} unassigned channel(s) "
                f"{[s for _, s in empty_ch]} kept as samples "
                "(tmt.drop_empty_channels is false); their signal is isotope "
                "carry-over, not a sample", "WARN")

        chan_of = {s: c for c, s in ann}

        # ---- within-plex normalisation, before anything is joined -----
        # The channels of one plex are the same LC-MS run, so what differs
        # between them is how much peptide was loaded and how completely it
        # was labelled: a per-channel constant with no biology in it, which
        # the roll-up would otherwise sum straight into the protein.
        # Centring on the plex's own median channel rather than on 1 keeps the
        # values linear and leaves the BETWEEN-plex difference untouched —
        # that one is the batch the plex term (or the report's own median
        # normalisation) is there to absorb, and removing it here would hide
        # it from both. median() commutes with log2, so this is exactly the
        # per-channel median centring the report would do, applied one plex at
        # a time and before the roll-up rather than after it.
        if norm == "median" and keep_samples:
            med = vals[keep_samples].median(axis=0, skipna=True)
            usable = med[med.gt(0) & med.notna()]
            dead = [s for s in keep_samples if s not in usable.index]
            if dead:
                log(f"tmt {plex}: channel(s) {dead} have no positive median "
                    "and were left unscaled by the within-plex normalisation",
                    "WARN")
            if len(usable):
                fac = float(usable.median()) / usable
                vals[usable.index] = vals[usable.index].mul(fac, axis=1)
                lg = np.log2(fac.astype(float))
                norm_lo.append(float(lg.min()))
                norm_hi.append(float(lg.max()))
                log(f"tmt {plex}: within-plex median centring applied to "
                    f"{len(usable)} channel(s); scale factors log2 "
                    f"{lg.min():+.2f}..{lg.max():+.2f} (set "
                    "tmt.within_plex_normalise: none to keep FragPipe's "
                    "numbers)")
                # A channel that is both far off scale AND much emptier than
                # its neighbours is the case this step handles WORST: the
                # median is taken over OBSERVED values only, so a channel
                # whose low end went missing has a median sitting above its
                # true centre and is scaled up too little. An unequal but
                # complete load is exactly what median centring is for and is
                # not worth a warning; this combination is.
                gaps = vals[usable.index].isna().mean()
                thin = [s for s in usable.index
                        if abs(float(lg[s])) > 1.0
                        and float(gaps[s]) > float(gaps.median()) + 0.10]
                if thin:
                    log(f"tmt {plex}: channel(s) " + ", ".join(
                        f"{s} ({chan_of.get(s, '?')}) scaled "
                        f"{2 ** float(lg[s]):.2f}x with "
                        f"{100 * float(gaps[s]):.0f}% missing"
                        for s in thin) + f" are far off the plex scale AND "
                        f"much emptier than the rest (typical "
                        f"{100 * float(gaps.median()):.0f}%). The median is "
                        "taken over observed values, so these are "
                        "UNDER-corrected; check the loading before trusting "
                        "them, or drop them from the annotation", "WARN")
        elif norm == "none" and len(plexes) > 1:
            log(f"tmt {plex}: tmt.within_plex_normalise is 'none', so the "
                "channels of this plex keep whatever loading difference they "
                "were labelled with; the roll-up sums across them", "WARN")

        # ---- the reference channel, resolved per plex -----------------
        ref = ""
        if ref_name or ref_chan:
            if ref_name:
                hit = [s for s in keep_samples
                       if fnmatch.fnmatchcase(s, ref_name)]
                what = f"tmt.reference_name '{ref_name}'"
            else:
                hit = [s for s in keep_samples
                       if fnmatch.fnmatchcase(chan_of[s], ref_chan)]
                what = f"tmt.reference_channel '{ref_chan}'"
            if len(hit) != 1:
                die(f"plex {plex}: {what} matches {len(hit)} of its channels "
                    f"{hit or ''}, not exactly one. Its channels are "
                    f"{[f'{c}={s}' for c, s in ann]}. The reference is not at "
                    "a fixed position in every plex, so name it by the sample "
                    "name (tmt.reference_name, e.g. 'Pool*') when the channel "
                    "moves.")
            ref = hit[0]
            refs.append((plex, chan_of[ref], ref))
        else:
            # Nothing configured. Say so where it can be acted on rather than
            # quietly quantifying a bridge channel as if it were a sample.
            pool = [s for s in keep_samples if s.lower().startswith("pool")]
            if len(pool) == 1:
                log(f"tmt {plex}: channel {chan_of[pool[0]]} is named "
                    f"'{pool[0]}', which looks like a reference/bridge "
                    "channel. It is being quantified as an ordinary sample; "
                    "set tmt.reference_name: 'Pool*' to mark it, and "
                    "tmt.use_reference_ratios: true to divide by it", "WARN")

        if use_ratios:
            # A reference of 0 is not a reference. FragPipe writes 0 for "not
            # quantified", and with zero_intensity_is_missing false it reaches
            # here as a number: x/0 wrote +inf into every other channel of the
            # plex and 0/0 wrote NaN, while notna() counts an inf as an
            # observed value, so both the per-plex line below and the run's
            # "use_reference_ratios cost ..." summary reported nothing lost
            # over a poisoned matrix. The inf then survived sum(min_count=1)
            # into the protein matrix, stayed inf through log2, and turned the
            # taxon size factors of those samples into NaN. A feature whose
            # reference is 0 (or non-finite) has no reference in this plex,
            # which is the case already handled and counted.
            denom = vals[ref]
            denom = denom.where(denom.ne(0) & np.isfinite(denom))
            n_bad = int(denom.isna().sum())
            others = [s for s in keep_samples if s != ref]
            # Counted before and after the divide, not from the feature count:
            # the features with no reference are mostly sparse ones, so the
            # share of FEATURES lost and the share of VALUES lost differ by an
            # order of magnitude, and only the second says whether this
            # treatment quietly emptied the matrix.
            was = int(vals[others].notna().sum().sum())
            vals = vals.div(denom, axis=0)
            now = int(vals[others].notna().sum().sum())
            ratio_had += was
            ratio_lost += was - now
            keep_samples = others
            log(f"tmt {plex}: every channel divided by the reference "
                f"'{ref}' ({chan_of[ref]}); the reference column itself is "
                f"dropped, and {n_bad} feature(s) with no usable reference "
                "value (missing, 0 or non-finite) became missing in this "
                f"plex, costing {was - now} of {was} value(s)")
        elif ref:
            # The covariate treatment. A pooled bridge is not a biological
            # sample: left in the sample columns it acquires a condition in
            # the design, joins a group's mean, and shifts the size factors
            # towards a pool that is in every plex by construction.
            keep_samples = [s for s in keep_samples if s != ref]
            log(f"tmt {plex}: reference channel is '{ref}' ({chan_of[ref]}); "
                "dropped from the sample columns because it is a pooled "
                "bridge, not a biological sample (the covariate treatment: "
                "the plex stays in the model). Set tmt.use_reference_ratios "
                "true to divide every channel by it instead")

        for s in keep_samples:
            if s in owner:
                die(f"plexes {owner[s]} and {plex} both claim the sample name "
                    f"'{s}'. Sample names are the columns of the joined "
                    "matrix, so two plexes cannot share one; fix the "
                    "annotation files.")
            owner[s] = plex
            # No reference row: under either treatment the reference is not a
            # sample column, and the report matches every design row to a
            # column of the quant matrix and stops when one is missing.
            design_rows.append({"sample": s, "plex": plex,
                                "channel": chan_of[s]})
        int_cols += keep_samples
        frames.append(vals[keep_samples])
        log(f"tmt {plex}: {len(vals)} {level} feature(s), {len(ann)} channels "
            f"mapped by {how} -> {len(keep_samples)} sample column(s)")

    # ---- the outer join across plexes --------------------------------
    # reindex, not merge: a feature not identified in a plex must be NA for
    # every sample of that plex, never 0. Zero is a measurement here (and
    # FragPipe writes plenty of them), so filling one in for "not identified"
    # would turn structured, plex-shaped missingness into fold change.
    idx = pd.Index(list(meta), name="feature_id")
    mat = pd.concat([f.reindex(idx) for f in frames], axis=1)
    feats = pd.DataFrame({
        "feature_id": list(meta),
        "peptide": [m["peptide"] for m in meta.values()],
        "razor_protein": [m["razor"] for m in meta.values()]})
    feats["candidates"] = [list(m["cand"]) for m in meta.values()]
    seen = pd.Series([len(m["plexes"]) for m in meta.values()])
    hist = " ".join(f"{n}:{c}" for n, c in sorted(seen.value_counts().items()))
    log(f"tmt: {len(feats)} feature(s) over {len(plexes)} plexes and "
        f"{len(int_cols)} samples; features identified in N plexes -> {hist}")
    if zero_missing and n_zero:
        log(f"tmt: {n_zero} of {n_rows} reporter cell(s) "
            f"({100.0 * n_zero / max(n_rows, 1):.1f}%) are 0, which FragPipe "
            "writes for 'not quantified'; treated as missing (set "
            "zero_intensity_is_missing false to keep them)", "WARN")
    conflict = sum(1 for m in meta.values() if len(m["razors"]) > 1)
    if conflict:
        log(f"tmt: {conflict} feature(s) have a different razor protein in "
            "different plexes; the first plex that saw the feature wins and "
            "the candidate lists are unioned, so the shared-peptide rule sees "
            "every protein any plex mapped the feature to", "WARN")

    if min_purity > 0:
        log(f"tmt: min_purity={min_purity} dropped {pur_drop} of {pur_seen} "
            f"per-plex feature row(s) ({100.0 * pur_drop / max(pur_seen, 1):.1f}%) "
            f"on the median purity of their PSMs; {pur_unknown} matched no "
            "PSM row and were kept. This limits co-isolation; it does not "
            "correct the ratio compression co-isolation causes")

    minp = int(t.get("min_plexes") or 1)
    if minp > 1:
        keep = (seen >= minp).to_numpy()
        log(f"tmt: min_plexes={minp} drops {int((~keep).sum())} of "
            f"{len(feats)} feature(s) identified in fewer plexes")
        feats, mat = feats.loc[keep].reset_index(drop=True), mat.loc[keep]
    feats = pd.concat([feats, mat.reset_index(drop=True)], axis=1)

    if use_ratios:
        pct = 100.0 * ratio_lost / max(ratio_had, 1)
        log(f"tmt: use_reference_ratios cost {ratio_lost} of {ratio_had} "
            f"non-reference value(s) ({pct:.2f}%), which had no reference in "
            "their own plex and so became missing for that whole plex",
            "WARN" if pct >= 5 else "INFO")
        if pct >= 25:
            # Not fatal — the user asked for ratios — but at this rate the
            # ratio matrix is a different, much sparser experiment than the
            # intensity matrix, and that has to be said before the roll-up
            # rather than inferred from a thin result.
            log(f"tmt: {pct:.1f}% of the measured values are gone, so the "
                "ratio matrix is substantially sparser than the intensities. "
                "The covariate treatment (tmt.use_reference_ratios: false, "
                "plex in design_formula) keeps them and models the plex "
                "instead", "WARN")

    design = pd.DataFrame(design_rows)
    # The condition is NOT in these files and is never taken from the plex,
    # which is a batch: a plex-versus-plex contrast is a batch effect
    # presented as a hypothesis. It is either derivable from the sample names
    # or the user's to write down, and which of the two happened is recorded.
    design, cond_note = _tmt_add_condition(design, cfg)
    notes = ["input:            FragPipe TMT, "
             f"{len(plexes)} plex(es), {len(int_cols)} sample column(s)",
             f"condition source: {cond_note}"]
    # The normalisation belongs beside the numbers, not only in a log that
    # scrolls away: a matrix that has been median-centred per channel and one
    # that has not are different data, and nothing downstream can tell them
    # apart by looking.
    notes.append("within-plex norm: " + (
        f"median centring per channel, log2 factors "
        f"{min(norm_lo):+.2f}..{max(norm_hi):+.2f} over {len(plexes)} plex(es)"
        if norm == "median" and norm_lo else
        "median centring per channel (no channel could be scaled)"
        if norm == "median" else
        "none (tmt.within_plex_normalise: none) — channels carry their "
        "loading differences into the roll-up"))
    if min_purity > 0:
        notes.append(
            f"purity filter:    median PSM purity >= {min_purity} "
            f"(from psm.tsv); {pur_drop} of {pur_seen} feature(s) dropped, "
            f"{pur_unknown} unjudged and kept")
    notes.append(f"feature min_plexes: {int(t.get('min_plexes') or 1)} "
                 "(feature level, before the roll-up)")
    if refs:
        notes.append("reference:        " + (
            "ratios (every channel divided by its plex reference; "
            f"{ratio_lost}/{ratio_had} value(s) lost to a missing reference)"
            if use_ratios else
            "covariate (dropped from the design; plex stays in the model)"))
        notes.append("reference channel: " + ", ".join(
            f"{p}={c}/{s}" for p, c, s in refs))
    else:
        notes.append("reference:        none named (tmt.reference_name / "
                     "tmt.reference_channel are unset)")
    # attrs, not a return value: every caller of read_feature_table unpacks a
    # 3-tuple, and widening that signature for one format would touch every
    # label-free path.
    design.attrs["design_notes"] = notes
    return feats, int_cols, design


# Separators a sample name might carry its condition in front of. "." is
# included because FragPipe run names often use it, but a bare digit after the
# split is a replicate index, never a condition.
_TMT_NAME_SEPS = ("_", "-", ".")


def _tmt_split_condition(samples):
    """-> (mapping, how) for an UNAMBIGUOUS name split, else (None, why).

    Only a split that partitions every sample into at least two levels of at
    least two samples each is accepted, and only when no other separator gives
    a DIFFERENT partition. Anything looser invents a hypothesis: "MF0030" and
    "MF0071" would become one condition per sample, and "resp-1_a" would mean
    two different things depending on which separator was tried first.
    """
    found = {}
    for sep in _TMT_NAME_SEPS:
        if not all(sep in s for s in samples):
            continue
        m = {s: s.split(sep)[0] for s in samples}
        lv = sorted(set(m.values()))
        if len(lv) < 2 or any(not x or x.isdigit() for x in lv):
            continue
        if any(sum(1 for v in m.values() if v == x) < 2 for x in lv):
            continue
        found[sep] = m
    if not found:
        return None, ("no separator (" + ", ".join(_TMT_NAME_SEPS) + ") splits "
                      "every sample name into two or more conditions of two "
                      "or more samples")
    # Compared as PARTITIONS, not as label maps: "a-x_1" and "a-x_2" fall
    # together whichever separator is used, and only a split that groups the
    # samples differently is a real ambiguity about what the condition is.
    parts = {frozenset(frozenset(s for s in m if m[s] == lv)
                       for lv in set(m.values())) for m in found.values()}
    if len(parts) > 1:
        return None, ("the sample names group differently on " +
                      ", ".join(f"'{s}'" for s in found) +
                      ", so which part of the name is the condition is "
                      "ambiguous")
    sep = next(s for s in _TMT_NAME_SEPS if s in found)
    return found[sep], f"the sample name before the first '{sep}'"


def _tmt_add_condition(design, cfg):
    """Add a `group` column to the TMT design when it can be had honestly.

    -> (design, note). The note goes into design_record.txt, because a
    condition that was DERIVED and one that was WRITTEN DOWN are not the same
    kind of claim and the table on disk cannot tell them apart afterwards.
    """
    spec = str((cfg.get("tmt") or {}).get("condition_from_name", "auto"))
    samples = [str(s) for s in design["sample"]]

    def _metadata_note():
        """What to tell the reader about supplying the condition themselves.

        The advice used to be 'supply it in analysis.metadata' whether or not
        analysis.metadata was already set and already covered every sample -
        which reads as a defect when it is only a division of labour: the
        design recovered from the input carries what the input knows, and the
        metadata is merged later, at report time. Saying so is the difference
        between a warning and a false alarm.
        """
        a = cfg.get("analysis") or {}
        path = a.get("metadata") or ""
        if not path:
            return ("supply it in analysis.metadata, keyed on sample")
        if not os.path.exists(path):
            return (f"analysis.metadata is set to {path}, which does not "
                    "exist; the report will have no condition either")
        col = a.get("sample_col") or "sample"
        try:
            md = read_delim_table(path)
        except Exception:                                   # noqa: BLE001
            return (f"analysis.metadata ({path}) could not be read here, so "
                    "whether it supplies the condition is unknown")
        if col not in md.columns:
            return (f"analysis.metadata ({path}) has no '{col}' column "
                    f"(analysis.sample_col), only {list(md.columns)[:6]}")
        have = set(md[col].astype(str))
        miss = [s for s in samples if s not in have]
        if miss:
            return (f"analysis.metadata ({path}) is keyed on '{col}' but "
                    f"does not name {len(miss)} of these samples "
                    f"({miss[:4]}), so the report will drop them")
        return (f"analysis.metadata ({path}) does name every sample, so the "
                "report supplies the condition; this affects only "
                "design_from_input.tsv, which records what the INPUT knew")
    if not spec:
        log("tmt: tmt.condition_from_name is empty, so no condition is "
            f"derived from the sample names; {_metadata_note()}", "WARN")
        return design, "not derived (tmt.condition_from_name is empty)"
    if spec == "auto":
        mapping, how = _tmt_split_condition(samples)
        if mapping is None:
            log(f"tmt: the condition could not be derived from the sample "
                f"names ({how}), so the design has no group column. It is NOT "
                f"taken from the plex, which is a batch: {_metadata_note()}",
                "WARN")
            return design, f"not derived ({how})"
    else:
        try:
            rx = re.compile(spec)
        except re.error as e:
            die(f"tmt.condition_from_name is neither 'auto', empty, nor a "
                f"valid regular expression: {e}")
        if rx.groups != 1:
            die(f"tmt.condition_from_name '{spec}' has {rx.groups} capture "
                "groups; it needs exactly one, and that group is the "
                "condition.")
        hit = {s: rx.search(s) for s in samples}
        miss = [s for s, m in hit.items() if not m or not m.group(1)]
        if miss:
            die(f"tmt.condition_from_name '{spec}' captures nothing in "
                f"{len(miss)} of {len(samples)} sample name(s), e.g. "
                f"{miss[:5]}. Every sample needs a condition, so fix the "
                "pattern or write analysis.metadata by hand.")
        mapping = {s: hit[s].group(1) for s in samples}
        how = f"tmt.condition_from_name '{spec}'"
    design = design.copy()
    design["group"] = [mapping[s] for s in samples]
    sizes = design.groupby("group").size().to_dict()
    log(f"tmt: condition derived from {how}: {sizes}. This is a GUESS from "
        "the annotation's sample names — check it, or set analysis.metadata "
        "to state the condition explicitly", "WARN")
    tab = design.groupby(["plex", "group"]).size().unstack(fill_value=0)
    if len(tab) > 1 and (tab > 0).sum(axis=1).max() == 1:
        die("the condition derived from the sample names is perfectly "
            "confounded with the plex: each of the "
            f"{len(tab)} plexes contains exactly one condition "
            f"({dict(zip(tab.index, tab.idxmax(axis=1)))}). The plex is a TMT "
            "batch, so no model can tell the batch from the biology and any "
            "fold change would be both. Either the derivation is wrong (set "
            "tmt.condition_from_name or analysis.metadata), or the experiment "
            "cannot answer this question — a TMT design needs each condition "
            "spread over several plexes.")
    return design, f"derived from {how}"


def read_fragpipe_tmt_peptides(root, cfg):
    """The identification half of read_fragpipe_tmt: peptides + candidates.

    Same contract as read_feature_peptides — no intensity is touched — so the
    taxonomy stages still run on a TMT project whose channels this reader
    would refuse (an annotation that does not match, a plex without one).
    """
    t = cfg.get("tmt") or {}
    level = str(t.get("level") or "ion").lower()
    if level not in TMT_LEVEL_FILES:
        die(f"tmt.level must be 'ion' or 'peptide', got '{level}'")
    if not os.path.isdir(root):
        die(f"quant_table not found: {root}. quant_format 'fragpipe_tmt' "
            "expects the FragPipe run directory holding the per-plex TMTn/ "
            "folders.")
    pref = excluded_prefixes(cfg)
    out = []
    for plex, pdir in tmt_plex_dirs(root, cfg):
        path = os.path.join(pdir, TMT_LEVEL_FILES[level])
        if not os.path.exists(path):
            die(f"plex {plex} has no {TMT_LEVEL_FILES[level]}: {path}")
        refuse_isobaric_matrix(path, header_columns(path))
        df = read_delim_table(path)
        prot = "Protein" if "Protein" in df.columns else "Protein ID"
        if prot not in df.columns:
            die(f"{path}: no 'Protein' or 'Protein ID' column")
        bad = pd.Series(False, index=df.index)
        for c in ("Is Decoy", "Is Contaminant"):
            if c in df.columns:
                bad |= df[c].astype(str).str.lower().isin(["true", "1", "yes"])
        if pref:
            bad |= df[prot].astype(str).map(first_token).str.startswith(pref)
        df = df.loc[~bad].reset_index(drop=True)
        pep_col = [c for c in ("Peptide Sequence", "Peptide")
                   if c in df.columns][:1]
        if not pep_col:
            die(f"{path}: no peptide sequence column")
        pep = df[pep_col[0]].astype(str).fillna("")
        if "Modified Sequence" in df.columns:
            empty = pep.eq("") | pep.eq("nan")
            if bool(empty.any()):
                ms = df["Modified Sequence"].astype(str).fillna("")
                pep = pep.mask(empty, ms.map(strip_modifications))
        razor = df[prot].astype(str).map(first_token)
        mapped = (df["Mapped Proteins"] if "Mapped Proteins" in df.columns
                  else pd.Series([""] * len(df), index=df.index))
        one = pd.DataFrame({"peptide": pep, "razor_protein": razor})
        one["candidates"] = [[r] + [x for x in split_ids(m)
                                    if not x.startswith(pref)]
                             for r, m in zip(razor, mapped)]
        out.append(one[one["peptide"].ne("")])
    res = pd.concat(out, ignore_index=True)
    # One row per (peptide, protein): the same peptide is identified in
    # several plexes, and the callers count rows.
    res["_k"] = res["peptide"] + "\t" + res["razor_protein"]
    res = res[~res["_k"].duplicated()].drop(columns="_k").reset_index(drop=True)
    log(f"tmt: {len(res)} distinct (peptide, razor protein) pair(s) across "
        f"{len(out)} plexes")
    return res


def read_feature_table(path, fmt, cfg):
    """-> (features, int_cols, design)

    features: feature_id, razor_protein, candidates (list of str), + int cols
    design:   sample/condition/replicate table recovered from the input, or None
    """
    if fmt == "fragpipe_tmt":
        # Not a single table: one per plex, joined here rather than by a
        # search engine. Same return shape all the same.
        return read_fragpipe_tmt(path, cfg)
    # Recognised from the header alone, before pandas parses the body: the TMT
    # msstats.csv dies in the C parser on unquoted commas in
    # Protein.Description, and a tmt-report matrix parses perfectly and is
    # already log2. Both have to be named, not guessed at downstream.
    refuse_isobaric_matrix(path, header_columns(path))
    df = read_delim_table(path)
    design = None

    if fmt in ("fragpipe_peptide", "fragpipe_ion"):
        # "Protein" first, not "Protein ID". For a UniProt database the two
        # differ only in decoration, but for a metagenome database FragPipe
        # fills "Protein ID" with "<id> <description>" — a value that can
        # never match the fasta, whose ids are the first whitespace token.
        # Preferring "Protein ID" silently left every protein unannotated.
        prot = "Protein" if "Protein" in df.columns else "Protein ID"
        if prot not in df.columns:
            die(f"{path}: no 'Protein ID' or 'Protein' column; is this a "
                "FragPipe combined_peptide/combined_ion table?")
        if df[prot].astype(str).str.contains(" ").any():
            n = int(df[prot].astype(str).str.contains(" ").sum())
            log(f"{n} value(s) in '{prot}' contain a space, so they carry a "
                "description rather than a bare identifier; using the first "
                "token, which is what the fasta and the eggNOG table are "
                "keyed on", "WARN")
        # Decoys and contaminants are rows of the same table. Left in, bovine
        # trypsin and keratin become quantified, "unannotated" proteins, and
        # in a smORF study their short peptides are exactly what gets called
        # a small protein.
        bad = pd.Series(False, index=df.index)
        for c in ("Is Decoy", "Is Contaminant"):
            if c in df.columns:
                bad |= df[c].astype(str).str.lower().isin(["true", "1", "yes"])
        pref = excluded_prefixes(cfg)
        if pref:
            bad |= df[prot].astype(str).map(first_token).str.startswith(pref)
        if bool(bad.any()):
            log(f"{int(bad.sum())} decoy/contaminant row(s) dropped from {path} "
                f"(prefixes {list(pref)} or an 'Is Decoy'/'Is Contaminant' "
                "flag); set exclude_id_prefixes to change the list", "WARN")
            df = df.loc[~bad].reset_index(drop=True)
        suffix = cfg.get("feature_intensity_suffix", "Intensity")
        drop = tuple(cfg.get("feature_exclude_suffixes",
                             ["MaxLFQ Intensity", "Spectral Count"]))
        int_cols = [c for c in df.columns
                    if c.endswith(suffix) and not c.endswith(drop)
                    and pd.api.types.is_numeric_dtype(df[c])]
        if not int_cols:
            die(f"{path}: no columns ending in '{suffix}'. Columns look like: "
                f"{list(df.columns)[:12]}")
        # Isobaric output is not supported, and the way it fails is the
        # dangerous part: FragPipe names TMT channels "Intensity <sample>",
        # which this suffix rule cannot see, so the only column it matches is
        # the bare MS1 "Intensity" — one precursor value pooled over every
        # channel. Quantifying that as the sole sample is silently wrong, so
        # refuse instead of returning it.
        reporter = [c for c in df.columns
                    if c.startswith(suffix + " ") and c not in int_cols
                    and pd.api.types.is_numeric_dtype(df[c])]
        if reporter and int_cols == [suffix]:
            die(f"{path}: this looks like isobaric (TMT/iTRAQ) output. "
                f"{len(reporter)} reporter-ion column(s) named "
                f"'{suffix} <sample>' are present, e.g. {reporter[:3]}, and "
                f"the only column matching '<sample> {suffix}' is the bare "
                f"MS1 '{suffix}', which is one precursor value pooled over "
                f"all channels. quant_format '{fmt}' does not read "
                "reporter-ion channels: quantifying this file would report a "
                "single sample and be silently wrong. Set quant_format: "
                "fragpipe_tmt and point quant_table at the run directory "
                "holding the per-plex TMTn/ folders, which reads the channels "
                "of every plex and joins them; or use a label-free or DIA-NN "
                "input.")
        rename = {}
        if cfg.get("manifest"):
            m = read_manifest(cfg["manifest"])
            rename, miss_rows, miss_cols = map_manifest_to_columns(
                m, int_cols, " " + suffix)
            if miss_rows:
                die_manifest_unmatched(miss_rows, path, int_cols)
            if miss_cols:
                log(f"{len(miss_cols)} quant column(s) are not in the manifest "
                    f"and will be dropped: {miss_cols[:5]}", "WARN")
            int_cols = [c for c in int_cols if c in rename]
            log(f"manifest: mapped {len(int_cols)} quant columns to sample names")
        if fmt == "fragpipe_ion":
            key = [c for c in ("Modified Sequence", "Peptide Sequence")
                   if c in df.columns][:1] + \
                  [c for c in ("Charge",) if c in df.columns]
        else:
            key = [c for c in ("Peptide Sequence", "Peptide") if c in df.columns][:1]
        if not key:
            die(f"{path}: no peptide sequence column")
        fid = join_cols(df, key)
        mapped = df["Mapped Proteins"] if "Mapped Proteins" in df.columns else ""
        if "Mapped Proteins" not in df.columns:
            log("no 'Mapped Proteins' column: every feature is treated as unique "
                "to its razor protein, so shared-peptide filtering is inactive",
                "WARN")
        razor = df[prot].astype(str).map(first_token)
        # A decoy/contaminant/entrapment candidate can never carry eggNOG
        # taxonomy, so leaving it in the candidate list makes every peptide it
        # touches shared_unknown_taxon — i.e. it vetoes real quantification.
        cand = [[r] + [x for x in split_ids(m) if not x.startswith(pref)]
                for r, m in
                zip(razor,
                    mapped if isinstance(mapped, pd.Series) else [""] * len(df))]
        pep_col = [c for c in ("Peptide Sequence", "Peptide") if c in df.columns]
        feats = pd.DataFrame({
            "feature_id": fid,
            "peptide": (df[pep_col[0]].astype(str) if pep_col
                        else fid.map(lambda x: x.split("_")[0])),
            "razor_protein": razor})
        feats["candidates"] = cand
        vals = df[int_cols]
        if cfg.get("zero_intensity_is_missing", True):
            nz = int((vals == 0).sum().sum())
            if nz:
                log(f"{nz} intensity cell(s) are 0, which FragPipe writes for "
                    "'not quantified'; treated as missing, because summing "
                    "them as real zeros turns missingness into fold change "
                    "(set zero_intensity_is_missing false to keep them)",
                    "WARN")
            vals = vals.where(vals != 0)
        if rename:
            vals = vals.rename(columns=rename)
            int_cols = [rename[c] for c in int_cols]
        feats = pd.concat([feats, vals], axis=1)
        design = None
        if cfg.get("manifest"):
            mm = read_manifest(cfg["manifest"])
            design = mm[["sample", "experiment", "bioreplicate"]].rename(
                columns={"experiment": "group", "bioreplicate": "replicate"})
            design = design[design["sample"].isin(int_cols)]
        return feats, int_cols, design

    # ---- long formats ------------------------------------------------
    if fmt == "msstats_csv":
        need_cols = ["ProteinName", "PeptideSequence", "Run", "Intensity"]
        miss = [c for c in need_cols if c not in df.columns]
        if miss:
            die(f"{path}: MSstats input missing {miss}")
        fkey = ["PeptideSequence"] + [
            c for c in ("PrecursorCharge", "FragmentIon", "ProductCharge")
            if c in df.columns and df[c].notna().any()]
        df["_feature"] = join_cols(df, fkey)
        prot_col, run_col, val_col = "ProteinName", "Run", "Intensity"
        cond = "Condition" if "Condition" in df.columns else None
        rep = "BioReplicate" if "BioReplicate" in df.columns else None
    else:  # msstats_feature (FeatureLevelData)
        ren = {"PROTEIN": "ProteinName", "FEATURE": "_feature",
               "PEPTIDE": "PeptideSequence", "RUN": "Run"}
        df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
        if "_feature" not in df.columns and "PeptideSequence" in df.columns:
            df["_feature"] = df["PeptideSequence"].astype(str)
        val_col = "INTENSITY" if "INTENSITY" in df.columns else "ABUNDANCE"
        if val_col not in df.columns:
            die(f"{path}: no INTENSITY or ABUNDANCE column in FeatureLevelData")
        if val_col == "ABUNDANCE":
            log("FeatureLevelData ABUNDANCE is log2 and already normalised by "
                "MSstats; de-logging it for roll-up", "WARN")
            df[val_col] = 2 ** df[val_col]
        prot_col, run_col = "ProteinName", "Run"
        cond = "GROUP_ORIGINAL" if "GROUP_ORIGINAL" in df.columns else (
            "GROUP" if "GROUP" in df.columns else None)
        rep = "SUBJECT_ORIGINAL" if "SUBJECT_ORIGINAL" in df.columns else (
            "SUBJECT" if "SUBJECT" in df.columns else None)

    if cond:
        cols = [run_col, cond] + ([rep] if rep else [])
        design = (df[cols].drop_duplicates()
                  .rename(columns={run_col: "sample", cond: "group",
                                   **({rep: "replicate"} if rep else {})}))
    log("MSstats-format input carries only the razor protein per feature, so "
        "shared-peptide and taxon-unique filtering are unavailable; use "
        "combined_peptide.tsv or combined_ion.tsv if you need them", "WARN")

    dup = int(df.duplicated(subset=["_feature", prot_col, run_col]).sum())
    if dup:
        log(f"{dup} duplicated (feature, run) row(s) in {path} were summed — "
            "fractions, or charge states the feature key does not separate; "
            "the scale then differs between runs with different fraction "
            "counts", "WARN")
    # groupby().sum(min_count=1), not pivot_table(aggfunc="sum"): pandas sums
    # an all-NaN cell to 0.0, so a feature never seen in a run was written as
    # a real zero and added as one by the roll-up.
    wide = (df.groupby(["_feature", prot_col, run_col], sort=False)[val_col]
              .sum(min_count=1).unstack(run_col))
    int_cols = [str(c) for c in wide.columns]
    wide.columns = int_cols
    wide = wide.reset_index().rename(columns={"_feature": "feature_id",
                                              prot_col: "razor_protein"})
    pep_map = dict(zip(df["_feature"], df["PeptideSequence"].astype(str))) \
        if "PeptideSequence" in df.columns else {}
    wide["peptide"] = [pep_map.get(f, str(f).split("_")[0])
                       for f in wide["feature_id"]]
    wide["razor_protein"] = wide["razor_protein"].astype(str)
    wide["candidates"] = [[x] for x in wide["razor_protein"]]
    return wide, int_cols, design


def read_feature_peptides(path, fmt, cfg):
    """-> DataFrame(peptide, razor_protein, candidates). No intensities.

    The identification half of read_feature_table. stage_unipept and
    stage_taxonomy ask what peptides were seen and which proteins they could
    have come from; neither reads a single intensity. Routing them through the
    full reader meant a table whose QUANTIFICATION this tool cannot use —
    isobaric channels, no column matching the intensity suffix, a manifest run
    that maps to nothing — aborted the taxonomy work too, although its peptide
    column was perfectly readable. Decoy/contaminant filtering is kept: those
    rows must not become taxon votes.
    """
    if fmt == "fragpipe_tmt":
        return read_fragpipe_tmt_peptides(path, cfg)
    refuse_isobaric_matrix(path, header_columns(path))
    df = read_delim_table(path)
    pref = excluded_prefixes(cfg)

    if fmt in ("fragpipe_peptide", "fragpipe_ion"):
        prot = "Protein" if "Protein" in df.columns else "Protein ID"
        if prot not in df.columns:
            die(f"{path}: no 'Protein ID' or 'Protein' column; is this a "
                "FragPipe combined_peptide/combined_ion table?")
        bad = pd.Series(False, index=df.index)
        for c in ("Is Decoy", "Is Contaminant"):
            if c in df.columns:
                bad |= df[c].astype(str).str.lower().isin(["true", "1", "yes"])
        if pref:
            bad |= df[prot].astype(str).map(first_token).str.startswith(pref)
        if bool(bad.any()):
            log(f"peptide-only reader: {int(bad.sum())} decoy/contaminant "
                f"row(s) dropped from {path}", "WARN")
            df = df.loc[~bad].reset_index(drop=True)
        pep_col = [c for c in ("Peptide Sequence", "Peptide")
                   if c in df.columns][:1]
        if not pep_col:
            die(f"{path}: no peptide sequence column")
        razor = df[prot].astype(str).map(first_token)
        mapped = (df["Mapped Proteins"] if "Mapped Proteins" in df.columns
                  else pd.Series([""] * len(df), index=df.index))
        cand = [[r] + [x for x in split_ids(m) if not x.startswith(pref)]
                for r, m in zip(razor, mapped)]
        out = pd.DataFrame({"peptide": df[pep_col[0]].astype(str),
                            "razor_protein": razor})
        out["candidates"] = cand
        return out

    # ---- MSstats long formats: one razor protein per row, no shared list ---
    ren = {"PROTEIN": "ProteinName", "PEPTIDE": "PeptideSequence"}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    if "PeptideSequence" not in df.columns:
        die(f"{path}: no PeptideSequence column; MSstats input cannot supply "
            "peptides for the taxonomy stages")
    if "ProteinName" not in df.columns:
        die(f"{path}: no ProteinName column")
    razor = df["ProteinName"].astype(str).map(first_token)
    keep = ~razor.str.startswith(pref) if pref else pd.Series(True, index=df.index)
    if not bool(keep.all()):
        # The full reader does not screen decoys out of MSstats input, but a
        # decoy peptide can never carry a real LCA, so letting it vote here
        # would only add noise. Logged because the two readers differ on it.
        log(f"peptide-only reader: {int((~keep).sum())} row(s) with a "
            f"decoy/contaminant protein dropped from {path}", "WARN")
    df, razor = df.loc[keep], razor.loc[keep]
    out = pd.DataFrame({"peptide": df["PeptideSequence"].astype(str),
                        "razor_protein": razor.values})
    # One row per (feature, run) in a long table; the callers count rows, so
    # collapse to one row per (peptide, protein) to keep the tallies honest.
    out = out.drop_duplicates().reset_index(drop=True)
    out["candidates"] = [[x] for x in out["razor_protein"]]
    return out


def peptide_features(cfg, stage):
    """Peptides + candidate proteins for a taxonomy stage, either way round.

    Governed by peptide_only_reader (auto | always | never). 'auto' keeps the
    full reader's behaviour whenever it works and only falls back when it
    refuses, so nothing about an already-working project changes. The path
    taken is always logged: the two readers agree on peptides and candidates
    but not on row counts for MSstats input, and the reader used has to be
    visible in the log when the numbers are read back.
    """
    fmt = cfg["quant_format"]
    mode = str(cfg.get("peptide_only_reader", "auto")).lower()
    if mode not in ("auto", "always", "never"):
        die(f"unknown peptide_only_reader '{mode}'; choose auto, always or never")
    if mode == "always":
        log(f"{stage}: reading peptides only (peptide_only_reader: always); "
            "quantification columns are not touched")
        return read_feature_peptides(cfg["quant_table"], fmt, cfg)
    try:
        feats, _, _ = read_feature_table(cfg["quant_table"], fmt, cfg)
        log(f"{stage}: peptides taken from the full quant reader")
        return feats
    except StageError as e:
        if mode == "never":
            raise
        log(f"{stage}: the full quant reader refused this table ({e}); "
            "falling back to the peptide-only reader, because this stage "
            "needs peptides and not intensities. Quantification is still "
            "unavailable — only the taxonomy work continues.", "WARN")
        return read_feature_peptides(cfg["quant_table"], fmt, cfg)


ASSIGNMENT_MODES = ("protein_unique", "taxon_unique",
                    "taxon_or_family_unique", "razor")
ROLLUP_METHODS = ("sum", "median_polish")


def _nanmed(a, axis):
    """Median over the observed entries; 0 where a whole row/column is NaN.

    np.nanmedian emits a RuntimeWarning per all-NaN slice, and an all-NaN row
    (a peptide missing from every sample) is the normal case here, so the
    empty slices are excluded rather than warned about. Returning 0 for them
    is what the polish needs: no observation, no effect to remove.
    """
    ok = ~np.isnan(a)
    if axis == 1:
        out = np.zeros(a.shape[0])
        idx = ok.any(axis=1)
        if idx.any():
            out[idx] = np.nanmedian(a[idx, :], axis=1)
    else:
        out = np.zeros(a.shape[1])
        idx = ok.any(axis=0)
        if idx.any():
            out[idx] = np.nanmedian(a[:, idx], axis=0)
    return out


def median_polish(mat, max_iter=10, tol=1e-4):
    """Tukey median polish of a log2 (feature x sample) matrix with NaNs.

    Returns the fitted sample profile (overall + column effects), i.e. the
    protein log2 abundance with each feature's own response level removed.
    That is the whole point: a peptide that ionises ten times better than its
    neighbour, or that is observed in only half the samples, shifts its own
    row and not the between-sample comparison.

    Samples where no feature was observed come back NaN - nothing was measured
    there, and the polish must not invent a value.
    """
    r = np.array(mat, dtype="float64")
    nrow, ncol = r.shape
    row_eff = np.zeros(nrow)
    col_eff = np.zeros(ncol)
    overall = 0.0
    seen = ~np.isnan(r)
    dead_col = ~seen.any(axis=0)          # sample with no observed feature
    prev = None
    for _ in range(max_iter):
        rm = _nanmed(r, 1)
        r = r - rm[:, None]
        row_eff += rm
        m = float(np.median(row_eff)) if nrow else 0.0
        row_eff -= m
        overall += m
        cm = _nanmed(r, 0)
        r = r - cm[None, :]
        col_eff += cm
        m = float(np.median(col_eff)) if ncol else 0.0
        col_eff -= m
        overall += m
        fit = overall + col_eff
        if prev is not None and float(np.max(np.abs(fit - prev))) < tol:
            break
        prev = fit
    fit = overall + col_eff
    fit[dead_col] = np.nan
    return fit


def rollup_median_polish(use, int_cols):
    """Protein-level matrix by median polish, rescaled to preserve the sum.

    The polish fixes ratios, not magnitudes: its output is defined only up to
    an additive constant per protein. Each protein is therefore shifted so its
    total over the samples where it was seen equals the plain sum's total,
    exactly as MaxLFQ does. Absolute numbers stay on the scale users already
    read; only the way they are apportioned between samples changes.
    """
    mat = use[int_cols].to_numpy(dtype="float64", na_value=np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log2(np.where(mat > 0, mat, np.nan))
    keys = use["_assigned"].to_numpy()
    order = np.argsort(keys, kind="stable")
    keys_s, lg_s, mat_s = keys[order], lg[order], mat[order]
    if not len(keys_s):
        return pd.DataFrame(columns=int_cols,
                            index=pd.Index([], name="_assigned"))
    bounds = np.flatnonzero(np.r_[True, keys_s[1:] != keys_s[:-1], True])
    ids, rows, n_single, n_polished = [], [], 0, 0
    for a, b in zip(bounds[:-1], bounds[1:]):
        block, raw = lg_s[a:b], mat_s[a:b]
        if b - a == 1:
            # One feature: the polish of a single row is that row. Skipping it
            # keeps the common case fast and identical to the sum.
            prof = block[0].copy()
            n_single += 1
        else:
            prof = median_polish(block)
            n_polished += 1
        lin = np.where(np.isnan(prof), np.nan, np.exp2(prof))
        tot_sum = np.nansum(raw)
        tot_pol = np.nansum(lin)
        if tot_pol > 0 and tot_sum > 0:
            lin = lin * (tot_sum / tot_pol)
        ids.append(keys_s[a])
        rows.append(lin)
    out = pd.DataFrame(rows, columns=int_cols,
                       index=pd.Index(ids, name="_assigned"))
    log(f"roll-up: median polish over {n_polished} multi-feature protein(s); "
        f"{n_single} single-feature protein(s) pass through unchanged")
    return out


def rollup_features(feats, int_cols, taxon_of, mode, min_features,
                    family_of=None, rollup_method="sum"):
    """Peptide/ion -> protein, with an explicit rule for shared features.

    protein_unique   keep only features matching exactly one protein
    taxon_unique     keep features whose candidate proteins all share one
                     taxon, and assign them to the razor protein
    taxon_or_family_unique
                     as taxon_unique, plus features whose candidates cannot be
                     compared by taxon (at least one has none) but do all share
                     one MMseqs family_id. A family is a sequence cluster, not
                     an organism - opt-in only, and recorded as its own class.
    razor            keep everything, assign to the razor protein (FragPipe's
                     own behaviour; arbitrary in a strain-redundant metagenome)

    rollup_method decides how the kept features become a number: "sum" (the
    default, and what every previously published metaannot number was computed
    with) or "median_polish" (log-space, ratio-preserving).
    """
    if mode not in ASSIGNMENT_MODES:
        die(f"peptide_assignment must be one of {list(ASSIGNMENT_MODES)}, "
            f"got '{mode}'")
    if rollup_method not in ROLLUP_METHODS:
        die(f"rollup_method must be one of {list(ROLLUP_METHODS)}, "
            f"got '{rollup_method}'")
    fget = (family_of or {}).get
    assign, klass = [], []
    tget = taxon_of.get
    for razor, cands in zip(feats["razor_protein"], feats["candidates"]):
        if len(cands) == 1:              # the common case, no set work needed
            klass.append("unique")
            assign.append(razor)
            continue
        cands = [c for c in dict.fromkeys(cands) if c]
        taxa = {str(tget(c, "") or "") for c in cands}
        unknown = "" in taxa or "nan" in taxa
        taxa.discard(""); taxa.discard("nan")
        if len(cands) <= 1:
            k = "unique"
        elif unknown:
            # Cannot be shown to be taxon-unique: at least one candidate has no
            # taxonomy. Treating it as taxon-unique would readmit exactly the
            # unannotated proteins this pipeline exists to scrutinise.
            # Under taxon_or_family_unique the user has said that for their
            # question one sequence family is close enough to one unit. The
            # class is then named family_unique and never taxon_unique, so the
            # choice stays visible in peptide_evidence.tsv and feature_quant.
            k = "shared_unknown_taxon"
            if family_of is not None:
                fams = {str(fget(c, "") or "") for c in cands}
                if len(fams) == 1 and not (fams & {"", "nan", "None"}):
                    k = "family_unique"
        elif len(taxa) == 1:
            k = "taxon_unique"
        else:
            k = "shared"
        klass.append(k)
        if mode == "protein_unique":
            assign.append(razor if k == "unique" else None)
        elif mode == "taxon_unique":
            assign.append(razor if k in ("unique", "taxon_unique") else None)
        elif mode == "taxon_or_family_unique":
            assign.append(razor if k in ("unique", "taxon_unique",
                                         "family_unique") else None)
        elif mode == "razor":
            assign.append(razor)
        else:
            die(f"peptide_assignment must be one of {list(ASSIGNMENT_MODES)}, "
                f"got '{mode}'")

    feats = feats.assign(_assigned=assign, _class=klass)
    n = len(feats)
    kept = feats["_assigned"].notna()
    n_fam = int((feats["_class"] == "family_unique").sum())
    log(f"features: {n} total, {int(kept.sum())} assigned under '{mode}' "
        f"(unique {int((feats['_class']=='unique').sum())}, "
        f"taxon-unique {int((feats['_class']=='taxon_unique').sum())}, "
        f"family-unique {n_fam}, "
        f"shared {int((feats['_class']=='shared').sum())}, "
        f"shared-unknown-taxon {int((feats['_class']=='shared_unknown_taxon').sum())})")
    if mode == "taxon_or_family_unique":
        log(f"peptide_assignment=taxon_or_family_unique: {n_fam} feature(s) "
            "kept because their candidates share an MMseqs family_id rather "
            "than a taxon. A family is a sequence cluster at "
            "cluster_min_seq_id, NOT an organism, so that intensity is "
            "attributed to one representative of the cluster",
            "WARN" if n_fam else "INFO")

    use = feats[kept]
    if rollup_method == "median_polish":
        log("roll-up: rollup_method=median_polish - protein values are a "
            "Tukey median polish in log2 space rescaled to the summed total, "
            "not the plain sum. They are not comparable with numbers from a "
            "rollup_method=sum run", "WARN")
        quant = rollup_median_polish(use, int_cols)
    else:
        log("roll-up: rollup_method=sum - protein values are the plain sum "
            "over observed features")
        quant = use.groupby("_assigned", sort=False)[int_cols].sum(min_count=1)

    # Booleans computed once and summed, rather than a Python lambda per
    # group: with ~10^5 proteins the lambdas dominated this stage.
    use = use.assign(one=1,
                     _is_unique=(use["_class"] == "unique").astype("int64"),
                     _is_taxu=(use["_class"] == "taxon_unique").astype("int64"),
                     _is_famu=(use["_class"] == "family_unique").astype("int64"))
    ev = use.groupby("_assigned", sort=False).agg(
        n_features_used=("one", "sum"),
        n_unique=("_is_unique", "sum"),
        n_taxon_unique=("_is_taxu", "sum"),
        n_family_unique=("_is_famu", "sum"))
    dropped = (feats.loc[~kept].assign(one=1)
               .groupby("razor_protein", sort=False)["one"].sum()
               .rename("n_features_dropped"))
    ev = ev.join(dropped, how="outer").fillna({"n_features_dropped": 0})
    for c in ("n_features_used", "n_unique", "n_taxon_unique",
              "n_family_unique", "n_features_dropped"):
        ev[c] = ev[c].fillna(0).astype(int)
    # Named, not renamed afterwards: the outer join above takes its index name
    # from whichever side is non-empty, so with NO feature assigned at all
    # (every peptide shared across taxa — the case a strain-redundant database
    # produces) the index arrived as 'razor_protein' and the rename raised
    # KeyError('protein_id') instead of writing an empty evidence table.
    ev.index.name = "protein_id"
    ev = ev.reset_index()
    ev["protein_id"] = ev["protein_id"].astype(str)
    # A protein quantified mostly from features it shares with same-taxon (or,
    # under taxon_or_family_unique, same-family) neighbours is a weaker
    # measurement than one carried by its own peptides: the intensity is real,
    # but which member of the group it belongs to is an assumption. Flagged
    # rather than dropped, because in a strain-redundant metagenome that is
    # most of the small unannotated proteins this tool exists to find.
    ev["taxon_unique_dominated"] = ((ev["n_taxon_unique"]
                                     + ev["n_family_unique"])
                                    > ev["n_unique"])
    n_dom = int(ev["taxon_unique_dominated"].sum())
    if n_dom:
        log(f"{n_dom}/{len(ev)} protein(s) rest more on shared-but-taxon- "
            "(or family-) unique features than on their own unique ones "
            "(taxon_unique_dominated in peptide_evidence.tsv and "
            "annotated_quant.tsv)", "WARN")

    # Said out loud whatever the threshold is: most proteins <= 100 aa are
    # single-peptide by nature, so this number is the size of the population
    # min_features_per_protein would remove.
    single = int((ev["n_features_used"] == 1).sum())
    if single:
        log(f"{single} protein(s) rest on a single assigned feature "
            f"(min_features_per_protein = {min_features})")
    if min_features > 1:
        enough = ev.loc[ev["n_features_used"].fillna(0) >= min_features, "protein_id"]
        before = len(quant)
        quant = quant.loc[quant.index.isin(set(enough))]
        log(f"{len(quant)}/{before} proteins retained with >= {min_features} "
            f"assigned features; the {before - len(quant)} dropped are absent "
            "from annotated_quant.tsv and from every report table", "WARN")
    return (quant.reset_index().rename(columns={"_assigned": "group_id"}),
            ev, feats)


FRAGPIPE_META = [
    "Protein", "Protein ID", "Entry Name", "Gene", "Gene Names", "Description",
    "Organism", "Protein Length", "Coverage", "Protein Existence",
    "Indistinguishable Proteins", "Protein Probability",
    "Top Peptide Probability", "Combined Total Peptides",
    "Combined Spectral Count", "Combined Unique Spectral Count",
    "Combined Total Spectral Count", "Total Peptides", "Unique Peptides",
    "Razor Peptides", "Total Spectral Count", "Unique Spectral Count",
    "Razor Spectral Count",
    # Booleans, and pandas calls a bool column numeric: without these two,
    # "Is Decoy"/"Is Contaminant" were auto-detected as sample intensities.
    "Is Decoy", "Is Contaminant",
]

# Columns that must stay text when the annotation table is read back.
# Left to pandas, a numeric-looking taxid becomes a float: 821 -> "821.0",
# which silently breaks every downstream taxonomy join.
ANN_STR_COLS = ["protein_id", "seed_taxid", "ko", "family_id", "pfam_hits",
                "bin", "og", "ec", "cazy", "dbcan_hits", "cog_cat",
                "rescued_by"]


def collapse_taxon_rank(cfg, mapping):
    """Collapse strain-level taxids to cfg['taxon_rank'] through the taxdump.

    A seed_ortholog taxid names a reference *strain*, so two ORFs of the same
    gut organism routinely carry different ones: their shared peptides are
    then classed 'shared' and dropped, and the median-of-ratios reference is
    computed over singleton pseudo-taxa. Collapsing first makes 'taxon' mean
    an organism. Off by default ("") because it needs a taxdump.
    """
    rank = str(cfg.get("taxon_rank", "") or "").strip()
    if not rank:
        # Said out loud on every run. The default is not "no taxonomy": it is
        # one taxon per eggNOG REFERENCE GENOME, which is what "taxon-unique"
        # and the taxon size factors are then computed over. Nobody reading
        # "1,842 distinct taxa" should have to open the config to learn that
        # several of those can be one gut organism.
        n = len({str(v) for v in mapping.values()
                 if v and str(v) not in ("nan", "None")})
        log(f"taxon_rank='' (the default): {n} distinct taxa, kept as raw "
            "eggNOG seed_ortholog taxids. Those name reference GENOMES, not "
            "organisms, so two ORFs of one gut species can carry different "
            "ones — their shared peptides then class as 'shared' and drop, "
            "and each pseudo-taxon gets its own size factor. Set taxon_rank "
            "to species/genus/family (needs db.ncbi_taxonomy) to collapse "
            "them first.")
        return mapping
    if rank not in RANKS:
        die(f"taxon_rank must be one of {RANKS} or empty, got '{rank}'")
    tdir = (cfg.get("db") or {}).get("ncbi_taxonomy", "")
    if not (tdir and os.path.exists(os.path.join(tdir, "nodes.dmp"))):
        die(f"taxon_rank='{rank}' needs db.ncbi_taxonomy (nodes.dmp/names.dmp). "
            "Set it, or clear taxon_rank to use the raw seed taxids.")
    tax = NCBITaxonomy(tdir)
    out, lost, cache = {}, 0, {}
    for pid, tid in mapping.items():
        t = str(tid or "")
        if not t or t in ("nan", "None"):
            out[pid] = ""
            continue
        if t not in cache:
            cache[t] = tax.lineage(t).get(rank, "")
        out[pid] = cache[t]
        if not cache[t]:
            lost += 1
    before = len({v for v in mapping.values() if v})
    after = len({v for v in out.values() if v})
    log(f"taxon_rank={rank}: {before} seed taxids collapsed to {after} {rank} "
        f"taxa; {lost} protein(s) have no {rank} in the taxdump and are "
        "excluded from taxon-based steps", "WARN" if lost else "INFO")
    return out


def resolve_taxonomy(cfg, p, ann):
    """-> {protein_id: taxid} under the configured taxonomy_source."""
    src = cfg.get("taxonomy_source", "eggnog")
    egg = dict(zip(ann["protein_id"],
                   (ann["seed_taxid"] if "seed_taxid" in ann.columns
                    else pd.Series([""] * len(ann))).fillna("").astype(str)))
    if src == "eggnog":
        return collapse_taxon_rank(cfg, egg)
    if not os.path.exists(p.taxonomy_comparison):
        log(f"taxonomy_source is '{src}' but {p.taxonomy_comparison} is absent; "
            "falling back to eggnog", "WARN")
        return collapse_taxon_rank(cfg, egg)
    comp = pd.read_csv(p.taxonomy_comparison, sep="\t", dtype=str, encoding="utf-8", encoding_errors="replace")
    if src == "unipept":
        out = dict(zip(comp["protein_id"], comp["unipept_taxid"].fillna("")))
        n = sum(1 for v in out.values() if v and v != "nan")
        log(f"taxonomy_source=unipept: {n} proteins carry a peptide-LCA consensus")
        return collapse_taxon_rank(cfg, {
            k: ("" if v in ("", "nan", "None") else v) for k, v in out.items()})
    if src == "concordant":
        ok = comp["verdict"].isin(["identical", "concordant"])
        keep = set(comp.loc[ok, "protein_id"])
        # Spelled out, because "concordant" sounds like it only removes
        # disagreements: everything the comparison could not decide is blanked
        # too — unipept_missing, eggnog_missing, a merged/obsolete taxid, and
        # every protein with no row at all because it never reached
        # consensus_min_peptides, which is most single-peptide small proteins.
        dropped = {k: int(v) for k, v in
                   comp["verdict"].value_counts().to_dict().items()
                   if k not in ("identical", "concordant")}
        absent = len([k for k in egg if k not in set(comp["protein_id"])])
        log(f"taxonomy_source=concordant: {len(keep)}/{len(comp)} proteins agree "
            "at genus or below and keep their eggNOG taxid; excluded from "
            f"taxon-based steps: {dropped} plus {absent} protein(s) with no "
            "row in the comparison at all",
            "WARN" if (dropped or absent) else "INFO")
        return collapse_taxon_rank(
            cfg, {k: (v if k in keep else "") for k, v in egg.items()})
    die(f"taxonomy_source must be eggnog, unipept or concordant, got '{src}'")


def taxon_size_factors(df, tax_col, int_cols, min_proteins):
    """Per-taxon, per-sample reference by median of ratios.

    A summed reference is not robust: one strongly changing protein inflates
    the sum that every other member of its taxon is then divided by, giving
    the passengers a bias in the opposite direction. Leave-one-out protects
    the changing protein itself but not its neighbours. The median over a
    taxon's proteins is unmoved by a minority that changes, which is exactly
    the DESeq2 size-factor argument applied within a taxon.

    Taxa with fewer than `min_proteins` fall back to the plain sum, since a
    median over two or three proteins is not robust either.

    The reference is DESeq2's: computed over proteins observed in *every*
    sample, so it does not move with the sample set. A per-protein mean over
    whichever samples happened to be observed centres a protein seen in 2 of
    36 samples on those 2, and a protein detected in one group only drags that
    group's median — with 60% missing values that is the common case, not the
    exception. When too few proteins are complete, the poscounts variant is
    used (mean over observed) but a sample's factor is only trusted when its
    median rests on at least `min_proteins` ratios; otherwise the taxon falls
    back to the sum.
    """
    rows = []
    for taxid, g in df.groupby(tax_col, sort=True):
        m = g[int_cols].apply(pd.to_numeric, errors="coerce")
        m = m.where(m > 0)
        logm = np.log(m)
        complete = logm.dropna(axis=0, how="any")
        fac, method = None, None
        if len(g) >= min_proteins and len(complete) >= min_proteins:
            ref = complete.mean(axis=1)
            fac = np.exp((complete.sub(ref, axis=0)).median(axis=0, skipna=True))
            method = "median_of_ratios"
        elif len(g) >= min_proteins and m.notna().sum(axis=1).ge(2).sum() >= min_proteins:
            ref = logm.mean(axis=1)
            ratios = logm.sub(ref, axis=0)
            fac = np.exp(ratios.median(axis=0, skipna=True))
            # A median over one or two ratios is not a size factor; leave the
            # sample NaN rather than reporting a number the report cannot
            # tell apart from a well-supported one.
            thin = ratios.notna().sum(axis=0) < min_proteins
            fac = fac.where(~thin)
            method = "median_of_ratios_poscounts"
            if bool(fac.notna().sum() == 0):
                fac, method = None, None
        if fac is None:
            tot = m.sum(axis=0, min_count=1)
            fac = tot / (tot[tot.notna()].mean() if tot.notna().any() else 1.0)
            method = "sum_fallback"
        rec = {tax_col: taxid, "n_proteins": len(g), "method": method}
        rec.update({c: (float(fac[c]) if pd.notna(fac.get(c, np.nan))
                        else float("nan")) for c in int_cols})
        rows.append(rec)
    return pd.DataFrame(rows)


def tmt_size_factor_plex_exposure(tx, tax_col, int_cols, plex_of, min_proteins):
    """-> (n_exposed, n_taxa): taxa whose size factor rests on plex-confined
    proteins.

    The taxon size factor is a median of ratios over all samples, so the
    question worth asking of an isobaric run is not whether it carries the
    plex — a per-sample loading shift is exactly what a size factor is for,
    and it is shared by every taxon — but whether plex-shaped MISSINGNESS
    changes it taxon by taxon. It can: with fewer than `min_proteins` members
    observed in every plex, taxon_size_factors() loses its complete-case
    reference and falls to the poscounts variant, whose median then mixes
    proteins whose reference was computed inside one plex with proteins whose
    reference spans them all. Measured on a fixture (see the tests): the
    taxon-specific part of the factor stays within 0.03 log2 for well-observed
    taxa, and is displaced 1.40 log2 for a taxon with 7 of its 10 proteins
    confined to one plex — in that plex only, since it has no factor at all in
    the other two. Those are the numbers README.md and the v0.3.0 CHANGELOG
    entry quote; this docstring carried a different pair, which left three
    places disagreeing about one fixture.
    """
    by_plex = {}
    for c in int_cols:
        by_plex.setdefault(plex_of.get(c, ""), []).append(c)
    if len(by_plex) < 2:
        return 0, 0
    exposed = total = 0
    for _, g in tx.groupby(tax_col, sort=False):
        if len(g) < min_proteins:
            continue                      # already on the sum fallback
        total += 1
        m = g[int_cols].apply(pd.to_numeric, errors="coerce")
        m = m.where(m > 0)
        complete = np.ones(len(m), dtype=bool)
        for cs in by_plex.values():
            complete &= m[cs].notna().any(axis=1).to_numpy()
        if int(complete.sum()) < min_proteins:
            exposed += 1
    return exposed, total


def stage_join(cfg, p):
    qpath, fmt = cfg["quant_table"], cfg["quant_format"]
    if not os.path.exists(qpath):
        log(f"quant table not found, skipping join: {qpath}", "WARN")
        return
    if fmt not in ALL_FORMATS:
        die(f"quant_format must be one of {sorted(ALL_FORMATS)}, got '{fmt}'")
    # A tmt-report matrix reads perfectly as a wide protein table and is
    # already log2, so nothing downstream would ever notice. Name it here.
    if fmt in PROTEIN_FORMATS and os.path.isfile(qpath):
        refuse_isobaric_matrix(qpath, header_columns(qpath))
        # And the per-plex protein.tsv, which carries neither of that
        # function's two markers and would otherwise be quantified as a
        # label-free protein table — one plex reported as the experiment.
        refuse_per_plex_reporter_table(
            qpath, header_columns(qpath),
            cfg.get("feature_intensity_suffix", "Intensity"))
    # A roll-up method only means something where there is something to roll
    # up. Ignored quietly, a config saying median_polish next to a DIA-NN
    # protein matrix would describe numbers the search engine produced.
    _rm = str(cfg.get("rollup_method", "sum") or "sum")
    # Checked here as well as in rollup_features, so a typo fails in the first
    # second rather than after the quant table has been read.
    if _rm not in ROLLUP_METHODS:
        die(f"rollup_method must be one of {list(ROLLUP_METHODS)}, got '{_rm}'")
    if _rm != "sum" and fmt not in FEATURE_FORMATS:
        die(f"rollup_method='{_rm}' applies only to feature-level input "
            f"({sorted(FEATURE_FORMATS)}). quant_format='{fmt}' is already "
            "rolled up to proteins by the search engine, so nothing here "
            "would use it. Set rollup_method: sum, or point quant_table at a "
            "peptide/ion table.")

    head = pd.read_csv(p.final, sep="\t", nrows=0, encoding="utf-8", encoding_errors="replace")
    dtypes = {c: str for c in ANN_STR_COLS if c in head.columns}
    ann = pd.read_csv(p.final, sep="\t", dtype=dtypes, low_memory=False, encoding="utf-8", encoding_errors="replace")
    ann = ann.rename(columns={ann.columns[0]: "protein_id"})
    required = ["protein_id", "bin"]
    missing_cols = [c for c in required if c not in ann.columns]
    if missing_cols:
        die(f"{p.final} is missing {missing_cols}. It does not look like a "
            "metaannot annotation table — it may be truncated or from another "
            "run. Delete it and rerun the finalise stage.")
    ann["protein_id"] = ann["protein_id"].astype(str)
    ann = ann[~ann["protein_id"].duplicated(keep="first")]

    # One taxonomy, resolved once: computing it twice logged its summary
    # twice and invited the two uses to drift apart.
    eff_taxonomy = resolve_taxonomy(cfg, p, ann)

    # ---- feature-level input: roll up to protein first -----------------
    evidence = None
    if fmt in FEATURE_FORMATS:
        feats, int_cols_f, design = read_feature_table(qpath, fmt, cfg)
        taxon_of = eff_taxonomy
        mode = cfg.get("peptide_assignment", "taxon_unique")
        # Only built when the user asked for it: family_of is what turns a
        # shared_unknown_taxon feature into a family_unique one, and building
        # it unconditionally would make an opt-in rule look active in the log
        # of every run.
        family_of = None
        if mode == "taxon_or_family_unique":
            if "family_id" not in ann.columns:
                die("peptide_assignment='taxon_or_family_unique' needs the "
                    "family_id column, which the cluster stage writes. Enable "
                    "run.cluster and rerun finalise, or use taxon_unique.")
            family_of = {k: ("" if str(v) in ("", "nan", "None") else str(v))
                         for k, v in zip(ann["protein_id"], ann["family_id"])}
            # With the cluster stage off, build_annotation fills family_id with
            # the protein's own id, so every family is a singleton and no two
            # candidates can ever share one. Silently running a fallback that
            # cannot fire would look like "the family rule found nothing".
            sizes = defaultdict(int)
            for _f in family_of.values():
                if _f:
                    sizes[_f] += 1
            n_multi = sum(1 for c in sizes.values() if c > 1)
            if not n_multi:
                die("peptide_assignment='taxon_or_family_unique' but no "
                    f"family_id in {p.final} has more than one member: the "
                    "cluster stage did not run, so family_id is just the "
                    "protein id and the family fallback can never fire. "
                    "Enable run.cluster and rerun finalise, or use "
                    "taxon_unique.")
            log(f"peptide_assignment=taxon_or_family_unique: {len(sizes)} "
                f"MMseqs families, {n_multi} of them with more than one "
                "member, available as the fallback unit for features whose "
                "candidates have no taxonomy")
        needs_taxa = mode in ("taxon_unique", "taxon_or_family_unique")
        if needs_taxa and not any(taxon_of.values()):
            if family_of:
                log("no seed_taxid in the annotation, so every shared feature "
                    "falls to the family rule under "
                    "taxon_or_family_unique", "WARN")
            else:
                log("no seed_taxid in the annotation, so taxon_unique cannot "
                    "be evaluated; falling back to protein_unique", "WARN")
                mode = "protein_unique"
        rollup_method = str(cfg.get("rollup_method", "sum") or "sum")
        q, evidence, feat_class = rollup_features(
            feats, int_cols_f, taxon_of, mode,
            cfg.get("min_features_per_protein", 1),
            family_of=family_of, rollup_method=rollup_method)
        # Recorded per protein, not only in the log: a table on disk must say
        # how its numbers were made, or a median_polish run and a sum run are
        # indistinguishable once the log is gone.
        evidence["rollup_method"] = rollup_method
        evidence["peptide_assignment"] = mode
        os.makedirs(p.quant_dir, exist_ok=True)
        if cfg.get("export_feature_quant", True):
            fq = feat_class[["feature_id", "peptide", "razor_protein",
                             "_class", "_assigned"] + int_cols_f].copy()
            fq["candidates"] = [";".join(c) for c in feat_class["candidates"]]
            fq = fq.rename(columns={"_class": "assignment_class",
                                    "_assigned": "assigned_protein"})
            with atomic_out(p.feature_quant) as tmp:
                fq.to_csv(tmp, sep="\t", index=False)
            log(f"feature-level quantification -> {p.feature_quant} "
                f"({len(fq)} features)")
        with atomic_out(f"{p.quant_dir}/peptide_evidence.tsv") as tmp:
            evidence.to_csv(tmp, sep="\t", index=False)
        log(f"join: wrote {p.quant_dir}/peptide_evidence.tsv")
        if design is not None:
            with atomic_out(f"{p.quant_dir}/design_from_input.tsv") as tmp:
                design.to_csv(tmp, sep="\t", index=False)
            log(f"join: recovered a design from the input -> "
                f"{p.quant_dir}/design_from_input.tsv")
            # How the design was arrived at, next to the design itself. The
            # report copies these lines into design_record.txt, because a
            # table of samples cannot say where its condition came from or
            # what happened to a reference channel that is no longer in it.
            notes = list(design.attrs.get("design_notes") or ())
            if notes:
                with atomic_out(f"{p.quant_dir}/design_notes.txt") as tmp:
                    with open(tmp, "w", encoding="utf-8") as fh:
                        fh.write("\n".join(notes) + "\n")
                log(f"join: how that design was made -> "
                    f"{p.quant_dir}/design_notes.txt")
        id_col, member_cols = "group_id", ["group_id"]
        meta, int_cols = ["group_id"], int_cols_f
    elif fmt == "msstats_protein":
        d = read_delim_table(qpath)
        pc = "Protein" if "Protein" in d.columns else "PROTEIN"
        rc = "originalRUN" if "originalRUN" in d.columns else "RUN"
        vc = "LogIntensities" if "LogIntensities" in d.columns else "ABUNDANCE"
        for c in (pc, rc, vc):
            if c not in d.columns:
                die(f"{qpath}: ProteinLevelData needs Protein, "
                    f"originalRUN/RUN and LogIntensities; missing '{c}'")
        gcol = "GROUP_ORIGINAL" if "GROUP_ORIGINAL" in d.columns else (
            "GROUP" if "GROUP" in d.columns else None)
        if gcol:
            scol = "SUBJECT_ORIGINAL" if "SUBJECT_ORIGINAL" in d.columns else (
                "SUBJECT" if "SUBJECT" in d.columns else None)
            cols = [rc, gcol] + ([scol] if scol else [])
            os.makedirs(p.quant_dir, exist_ok=True)
            with atomic_out(f"{p.quant_dir}/design_from_input.tsv") as tmp:
                (d[cols].drop_duplicates()
                 .rename(columns={rc: "sample", gcol: "group",
                                  **({scol: "replicate"} if scol else {})})
                 .to_csv(tmp, sep="\t", index=False))
            log(f"join: recovered a design from the input -> "
                f"{p.quant_dir}/design_from_input.tsv")
        # De-log so the taxon sums downstream are sums of intensities, not of
        # logs. MSstats' own normalisation is preserved by the transform.
        w = d.pivot_table(index=pc, columns=rc, values=vc, aggfunc="mean")
        w = (2 ** w).reset_index().rename(columns={pc: "group_id"})
        w["group_id"] = w["group_id"].astype(str)
        q = w
        id_col, member_cols = "group_id", ["group_id"]
        meta, int_cols = ["group_id"], [c for c in q.columns if c != "group_id"]
        log("MSstats ProteinLevelData: summarisation and any imputation were "
            "done by dataProcess(); check MBimpute/censoredInt there, because "
            "imputed low-abundance values are exactly where the KO-less "
            "fraction sits", "WARN")
    else:
        q = pd.read_csv(qpath, sep="\t", low_memory=False, encoding="utf-8", encoding_errors="replace")

    if fmt == "diann":
        id_col = "Protein.Group"
        # Protein.Ids too: a group's other members were never considered, so
        # any conflict between them was invisible.
        member_cols = [c for c in ["Protein.Group", "Protein.Ids"]
                       if c in q.columns]
        meta = [c for c in ["Protein.Group", "Protein.Ids", "Protein.Names",
                            "Genes", "First.Protein.Description"] if c in q.columns]
    elif fmt == "fragpipe":
        # "Protein" first, for the same reason as in read_feature_table: on a
        # metagenome database FragPipe writes "<id> <description>" into
        # "Protein ID", which can never match a fasta id.
        id_col = "Protein" if "Protein" in q.columns else "Protein ID"
        if ("Protein" in q.columns and "Protein ID" in q.columns
                and not q["Protein"].astype(str).map(first_token).equals(
                    q["Protein ID"].astype(str).map(first_token))):
            log("'Protein' and 'Protein ID' disagree after taking the first "
                f"token in {qpath}; using 'Protein', which is the fasta "
                "header, so the ids match the annotation table", "WARN")
        member_cols = [c for c in [id_col, "Indistinguishable Proteins"]
                       if c in q.columns]
        # An explicit metadata list, not "everything non-numeric": FragPipe
        # emits numeric metadata (Protein Length, Coverage, spectral counts)
        # that would otherwise be summed as if it were sample intensity.
        meta = [c for c in FRAGPIPE_META if c in q.columns]
    if id_col not in q.columns:
        die(f"{qpath}: expected column '{id_col}' not found")

    if fmt in ("diann", "fragpipe"):
        # Same rule as the feature tables: a decoy or a bovine contaminant is
        # not a protein of this metagenome and must not be quantified,
        # annotated or summed into a taxon.
        pref = excluded_prefixes(cfg)
        drop_rows = pd.Series(False, index=q.index)
        for c in ("Is Decoy", "Is Contaminant"):
            if c in q.columns:
                drop_rows |= q[c].astype(str).str.lower().isin(
                    ["true", "1", "yes"])
        if pref:
            drop_rows |= q[id_col].astype(str).map(first_token).str.startswith(pref)
        if bool(drop_rows.any()):
            log(f"{int(drop_rows.sum())} decoy/contaminant row(s) dropped from "
                f"{qpath} (exclude_id_prefixes {list(pref)} or an "
                "'Is Decoy'/'Is Contaminant' flag)", "WARN")
            q = q.loc[~drop_rows].reset_index(drop=True)

    # Detection is only needed for wide protein tables. The feature-level and
    # ProteinLevelData branches already know exactly which columns are samples,
    # and re-detecting here would overwrite that with a guess.
    if fmt in ("diann", "fragpipe"):
        explicit = cfg.get("intensity_columns") or []
        regex = cfg.get("intensity_regex") or ""
        if explicit:
            int_cols = [c for c in explicit if c in q.columns]
            missing = [c for c in explicit if c not in q.columns]
            if missing:
                log(f"intensity_columns not present in {qpath}: {missing}", "WARN")
        elif regex:
            rx = re.compile(regex)
            int_cols = [c for c in q.columns
                        if rx.search(c) and pd.api.types.is_numeric_dtype(q[c])]
        else:
            # not is_bool_dtype: pandas calls a bool column numeric, so
            # "Is Decoy"/"Is Contaminant" were auto-detected as samples.
            int_cols = [c for c in q.columns
                        if c not in meta and pd.api.types.is_numeric_dtype(q[c])
                        and not pd.api.types.is_bool_dtype(q[c])]
    # A manifest overrides column detection for protein-level tables too.
    if fmt in ("diann", "fragpipe") and cfg.get("manifest"):
        m = read_manifest(cfg["manifest"])
        sfx = " Intensity" if fmt == "fragpipe" else ""
        mapping, miss_rows, _ = map_manifest_to_columns(m, int_cols, sfx)
        if miss_rows:
            die_manifest_unmatched(miss_rows, qpath, int_cols)
        q = q.rename(columns=mapping)
        int_cols = [mapping[c] for c in int_cols if c in mapping]
        os.makedirs(p.quant_dir, exist_ok=True)
        with atomic_out(f"{p.quant_dir}/design_from_input.tsv") as tmp:
            (m[["sample", "experiment", "bioreplicate"]]
             .rename(columns={"experiment": "group",
                              "bioreplicate": "replicate"})
             .to_csv(tmp, sep="\t", index=False))
        log(f"manifest: design -> {p.quant_dir}/design_from_input.tsv")

    if fmt == "fragpipe" and int_cols and cfg.get("zero_intensity_is_missing", True):
        # FragPipe's 0 means "not quantified in this run", not "zero".
        nz = int((q[int_cols] == 0).sum().sum())
        if nz:
            log(f"{nz} intensity cell(s) are 0 in {qpath}; treated as missing "
                "(set zero_intensity_is_missing false to keep them)", "WARN")
            q[int_cols] = q[int_cols].where(q[int_cols] != 0)

    if not int_cols:
        log(f"no numeric intensity columns found in {qpath}", "WARN")
    shown = ", ".join(int_cols[:8]) + (" ..." if len(int_cols) > 8 else "")
    log(f"join: {len(q)} groups x {len(int_cols)} intensity columns: {shown}")
    if fmt in ("diann", "fragpipe") and not (cfg.get("intensity_columns")
                                             or cfg.get("intensity_regex")):
        log("if any of those are not sample intensities, set intensity_columns "
            "or intensity_regex in the config", "WARN")

    q = q.copy()
    # first_token, not the raw cell: a DIA-NN "P1;P2" survives unchanged, but
    # a FragPipe "<id> <description>" is reduced to the id the fasta and the
    # eggNOG table are keyed on.
    q["group_id"] = q[id_col].astype(str).map(first_token)
    q["_row"] = range(len(q))

    # Explode with an explicit rank. A merge reorders rows, so "first member
    # after the join" is not the leading protein of the group.
    # Vectorised explode. iterrows() built a pandas Series per protein group
    # and was over half of this stage's runtime.
    joined = q[member_cols[0]].fillna("").astype(str)
    for c in member_cols[1:]:
        joined = joined.str.cat(q[c].fillna("").astype(str), sep=";")
    long = pd.DataFrame({"_row": q["_row"], "group_id": q["group_id"],
                         "protein_id": joined.str.split(r"[;,]")})
    long = long.explode("protein_id", ignore_index=True)
    # Same normalisation for every member: "Indistinguishable Proteins" and
    # "Protein ID" carry descriptions on a metagenome database, and splitting
    # those on commas produced ids like "3-aminomutase".
    long["protein_id"] = (long["protein_id"].str.strip()
                          .str.split(n=1).str[0].fillna(""))
    long = long[long["protein_id"].ne("")]

    # A join that matches nothing used to end in "N wholly unannotated" and
    # exit 0, which reads like a biology result rather than a broken id space.
    known = set(ann["protein_id"])
    ids = long["protein_id"].drop_duplicates()
    if len(ids):
        absent = ids[~ids.isin(known)]
        hit = len(ids) - len(absent)
        if hit == 0:
            die(f"none of the {len(ids)} protein ids in {qpath} appear in "
                f"{p.final}, e.g. {absent.head(3).tolist()}. The quant table "
                "and the annotation were not built from the same fasta, or "
                "the ids carry a prefix the annotation does not have.")
        if hit / len(ids) < 0.5:
            log(f"only {hit}/{len(ids)} ({100*hit/len(ids):.1f}%) protein ids "
                f"from {qpath} are present in {p.final}, e.g. "
                f"{absent.head(3).tolist()}; the rest join as unannotated. "
                "Check that both were built from the same fasta", "WARN")
    # Order is preserved by explode, so rank is a within-group cumcount, and
    # duplicates inside one group are dropped exactly as before.
    long = long[~long.duplicated(subset=["_row", "protein_id"])]
    long["member_rank"] = long.groupby("_row", sort=False).cumcount() + 1
    long = long.merge(ann, on="protein_id", how="left")

    def joined_unique(s):
        return "|".join(sorted({str(x) for x in s.dropna()
                                if str(x) not in ("", "nan")}))

    sizes = long.groupby("_row", sort=False)["protein_id"].transform("size")
    if bool((sizes == 1).all()):
        # Single-member groups (always the case after a feature roll-up):
        # there is nothing to disagree about, so skip the per-group lambdas.
        conf = long[["_row", "group_id"]].copy()
        conf["n_members"] = 1
        conf["n_annotated"] = long["bin"].notna().astype(int).values
        conf["bins"] = long["bin"].fillna("").astype(str).values
        conf["kos"] = (long["ko"].fillna("").astype(str).values
                       if "ko" in long.columns else "")
        conf["taxa"] = (long["seed_taxid"].fillna("").astype(str).values
                        if "seed_taxid" in long.columns else "")
    else:
        aggs = {
            "n_members": ("protein_id", "size"),
            "n_annotated": ("bin", lambda s: int(s.notna().sum())),
            "bins": ("bin", joined_unique),
        }
        if "ko" in long.columns:
            aggs["kos"] = ("ko", joined_unique)
        if "seed_taxid" in long.columns:
            aggs["taxa"] = ("seed_taxid", joined_unique)
        conf = long.groupby(["_row", "group_id"], sort=False).agg(**aggs).reset_index()
    for c in ("kos", "taxa"):
        if c not in conf.columns:
            conf[c] = ""
    conf["bin_conflict"] = conf["bins"].str.contains(r"\|", na=False)
    conf["ko_conflict"] = conf["kos"].str.contains(r"\|", na=False)
    conf["taxon_conflict"] = conf["taxa"].str.contains(r"\|", na=False)
    conf["unannotated_group"] = conf["n_annotated"] == 0

    os.makedirs(p.quant_dir, exist_ok=True)
    bad = conf[conf["bin_conflict"] | conf["ko_conflict"]
               | conf["taxon_conflict"] | conf["unannotated_group"]]
    with atomic_out(f"{p.quant_dir}/group_conflicts.tsv") as tmp:
        bad.drop(columns=["_row"]).to_csv(tmp, sep="\t", index=False)
    log(f"join: {int((conf['bin_conflict'] | conf['ko_conflict']).sum())}/{len(conf)} "
        f"groups internally inconsistent, {int(conf['unannotated_group'].sum())} "
        "wholly unannotated")

    lead = long[long["member_rank"] == 1]
    # Kept before protein_id is dropped: effective_taxid is keyed by protein,
    # and a DIA-NN group_id ("P1;P2") is not a protein id, so mapping the
    # group_id straight through gave every multi-member group no taxon at all.
    rep_protein = dict(zip(lead["group_id"], lead["protein_id"]))
    rep = lead.drop(columns=["protein_id", "member_rank"])
    # suffixes=("", "_ann") keeps the quant table's own column names intact.
    # Without it a shared name (length, Description, ...) becomes length_x /
    # length_y and every later reference to int_cols raises KeyError.
    out = q.merge(rep, on=["_row", "group_id"], how="left", suffixes=("", "_ann"))
    out = out.merge(conf[["_row", "group_id", "n_members", "bin_conflict",
                          "ko_conflict", "taxon_conflict", "unannotated_group"]],
                    on=["_row", "group_id"], how="left", suffixes=("", "_conf"))
    out = out.sort_values("_row").drop(columns=["_row"])
    if evidence is not None:
        out = out.merge(evidence.rename(columns={"protein_id": "group_id"}),
                        on="group_id", how="left", suffixes=("", "_ev"))
    if os.path.exists(p.taxonomy_comparison):
        cmp_cols = ["protein_id", "unipept_taxid", "unipept_rank",
                    "unipept_agreement", "n_peptides_with_lca",
                    "deepest_agreement", "verdict", "unipept_name", "eggnog_name"]
        # dtype=str, as in resolve_taxonomy: without it unipept_taxid comes
        # back float and lands in annotated_quant.tsv as "821.0" next to a
        # string effective_taxid, so the two never join.
        c = pd.read_csv(p.taxonomy_comparison, sep="\t", dtype=str, encoding="utf-8", encoding_errors="replace")
        c = c[[x for x in cmp_cols if x in c.columns]].rename(
            columns={"protein_id": "group_id", "verdict": "taxonomy_verdict"})
        out = out.merge(c, on="group_id", how="left", suffixes=("", "_tax"))
    # Assigned before the file is written. These were previously set further
    # down, after to_csv, so the column never reached disk and the report
    # silently fell back to the eggNOG taxid whatever taxonomy_source said.
    out["effective_taxid"] = (out["group_id"].map(rep_protein)
                              .fillna(out["group_id"])
                              .map(eff_taxonomy).fillna("")
                              .astype(str).replace({"nan": "", "None": ""}))
    out["taxonomy_source"] = cfg.get("taxonomy_source", "eggnog")
    with atomic_out(f"{p.quant_dir}/annotated_quant.tsv") as tmp:
        out.to_csv(tmp, sep="\t", index=False)

    # Which columns of that table are samples, written down rather than left
    # to be re-derived. The report and build_object.R both read this file
    # first, and without it they fall back to design_from_input.tsv — or, with
    # no manifest, to guessing from column types, which is the one path that
    # can sweep an annotation column into the assay. This stage is the only
    # place that KNOWS, because it is what renamed and selected them.
    # One name per line: readLines() at the other end, trimmed, blanks
    # dropped. An empty file is read as "no list", the same as no file.
    odd = [c for c in int_cols if c != c.strip()]
    if odd:
        log(f"sample column name(s) {odd[:3]} have leading or trailing "
            "whitespace. The report and the R object trim what they read from "
            "quant/sample_columns.txt, so those names will not match the "
            "columns of annotated_quant.tsv and will be reported as absent "
            "there; fix the header or the manifest", "WARN")
    with atomic_out(f"{p.quant_dir}/sample_columns.txt") as tmp:
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("".join(c + "\n" for c in int_cols))
    log(f"join: recorded {len(int_cols)} sample column(s) -> "
        f"{p.quant_dir}/sample_columns.txt")

    fam_col = "family_id" if "family_id" in rep.columns else None
    if fam_col:
        t2g = rep[[fam_col, "group_id"]].dropna()
        t2g = t2g[t2g[fam_col].astype(str).str.len() > 0]
        t2g.columns = ["term", "gene"]
        with atomic_out(f"{p.quant_dir}/term2gene_family.tsv") as tmp:
            t2g.drop_duplicates().to_csv(tmp, sep="\t", index=False)

    if "pfam_hits" in rep.columns:
        rows = []
        for _, r in rep.iterrows():
            v = r["pfam_hits"]
            if isinstance(v, str) and v:
                for t in v.split(";"):
                    if t:
                        rows.append((t, r["group_id"]))
        with atomic_out(f"{p.quant_dir}/term2gene_pfam.tsv") as tmp:
            (pd.DataFrame(rows, columns=["term", "gene"])
             .drop_duplicates().to_csv(tmp, sep="\t", index=False))

    # Taxon-level intensity: a "differentially abundant" KO-less protein is far
    # more likely to be tracking its source organism than a metabolic enzyme is.
    # One taxonomy for everything: the same resolve_taxonomy() result that
    # governed the roll-up also governs the taxon reference.
    tax_col = "effective_taxid"
    if tax_col and int_cols:
        tx = out[out[tax_col].notna() & out[tax_col].astype(str).ne("")
                 & out[tax_col].astype(str).ne("nan")]
        if len(tx):
            with atomic_out(f"{p.quant_dir}/taxon_intensity.tsv") as tmp:
                tx.groupby(tax_col)[int_cols].sum().reset_index().to_csv(
                    tmp, sep="\t", index=False)
            sf = taxon_size_factors(tx, tax_col, int_cols,
                                    int(cfg.get("taxon_min_proteins_for_factor", 4)))
            with atomic_out(f"{p.quant_dir}/taxon_size_factors.tsv") as tmp:
                sf.to_csv(tmp, sep="\t", index=False)
            # Counted from the method column, not len(sf): the old line said
            # "median-of-ratios for N" while most of those N were the sum.
            meth = sf["method"].value_counts().to_dict()
            _rank = str(cfg.get("taxon_rank", "") or "").strip()
            log(f"{tx[tax_col].nunique()} distinct taxa "
                f"(taxon_rank={_rank or 'unset, so reference-genome taxids'}); "
                f"size-factor method: { {k: int(v) for k, v in meth.items()} }")
            thin_rows = sf.loc[sf["method"].eq("median_of_ratios_poscounts"),
                               int_cols] if int_cols else sf.iloc[:0]
            nthin = int(thin_rows.isna().any(axis=1).sum()) if len(thin_rows) else 0
            if nthin:
                log(f"{nthin} taxon(s) have no size factor in at least one "
                    "sample, because that sample's median would have rested "
                    "on fewer than taxon_min_proteins_for_factor ratios; the "
                    "ratio model has no adjustment there", "WARN")
            nfall = int(meth.get("sum_fallback", 0))
            if nfall:
                log(f"{nfall}/{len(sf)} taxa have too few proteins for a "
                    "within-taxon median and fall back to the plain sum, where "
                    "one changing protein sets its own reference; raise "
                    "taxon_rank or analysis.taxon_min_proteins if the ratio "
                    "model matters", "WARN")
            # The size factor is computed across ALL samples, so an isobaric
            # run has to be asked whether the plex got into it. A common
            # loading shift does, and should — that is what a size factor is.
            # What must not pass unremarked is the taxon-by-taxon part, which
            # plex-shaped missingness can move (see
            # tmt_size_factor_plex_exposure).
            dpath = f"{p.quant_dir}/design_from_input.tsv"
            if fmt == "fragpipe_tmt" and os.path.exists(dpath):
                dz = pd.read_csv(dpath, sep="\t", dtype=str, encoding="utf-8",
                                 encoding_errors="replace")
                if "plex" in dz.columns:
                    n_exp, n_tot = tmt_size_factor_plex_exposure(
                        tx, tax_col, int_cols,
                        dict(zip(dz["sample"], dz["plex"])),
                        int(cfg.get("taxon_min_proteins_for_factor", 4)))
                    if n_exp:
                        log(f"{n_exp}/{n_tot} taxon(s) have fewer than "
                            "taxon_min_proteins_for_factor protein(s) "
                            "observed in EVERY plex, so their size factor "
                            "falls back to the poscounts variant and mixes "
                            "plex-confined proteins with cross-plex ones; "
                            "that is where a plex effect can reach the "
                            "taxon-specific part of the factor. Raise "
                            "analysis.min_plexes (protein level) or "
                            "tmt.min_plexes (feature level) if the ratio "
                            "model matters for those taxa", "WARN")
                    else:
                        log(f"all {n_tot} taxon(s) with enough proteins have "
                            "taxon_min_proteins_for_factor of them observed "
                            "in every plex, so the size factors rest on a "
                            "plex-complete reference")

    vis_col = "kegg_enrichment_visible" if "kegg_enrichment_visible" in out.columns \
        else ("kegg_enrichment_visible_ann" if "kegg_enrichment_visible_ann"
              in out.columns else None)
    if vis_col:
        vis = out[vis_col].fillna(False).astype(bool)
        log(f"join: {100*(~vis).mean():.1f}% of quantified protein groups are "
            "invisible to KEGG pathway enrichment")
    log(f"join: wrote {p.quant_dir}/")


# ======================================================================
# stage registry, signatures, state
# ======================================================================
_HASH_CACHE = {}
_HASH_FULL_MAX = 64 * 2**20      # digest whole files up to this size
_HASH_EDGE = 8 * 2**20           # larger ones: their first and last 8 MB


def _content_digest(path, size, mtime, full=False):
    """sha1 of a file's content, memoised on (path, size, mtime, mode).

    Big databases are digested at their two ends plus their size instead of
    end to end: reading 36 GB on every signature would cost more than the
    stage it protects, and a rebuilt database changes both ends anyway. That
    is blind to an in-place edit in the middle of a multi-GB file that keeps
    its size, so `full_content_digest: true` hashes every byte instead. The
    mode is part of the memo key, because the two answers differ.
    """
    key = (path, size, mtime, bool(full))
    hit = _HASH_CACHE.get(key)
    if hit is not None:
        return hit
    h = hashlib.sha1()
    try:
        with open(path, "rb") as fh:
            if full or size <= _HASH_FULL_MAX:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            else:
                h.update(fh.read(_HASH_EDGE))
                fh.seek(-_HASH_EDGE, os.SEEK_END)
                h.update(fh.read(_HASH_EDGE))
    except OSError:
        return None
    _HASH_CACHE[key] = h.hexdigest()
    return _HASH_CACHE[key]


def _stat(path, full=False):
    """Identity of one input: its content, not the text of its path.

    Hashing the path string meant that spelling the same directory as
    'results' on one run and as an absolute path on the next recomputed
    everything; hashing the mtime meant that a cp/scp of an unchanged
    database, or a stage rewriting byte-identical output, did the same.
    Position in the stage's input list already says which file this is.
    """
    try:
        st = os.stat(path)
    except OSError:
        return [None, None]
    if not os.path.isfile(path):      # directories, hhblits/foldseek prefixes
        return [st.st_size, int(st.st_mtime)]
    return [st.st_size,
            _content_digest(path, st.st_size, int(st.st_mtime), full)]


def _dig(cfg, dotted):
    cur = cfg
    for k in dotted.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


# Only the parts of the unipept block that change the answer. Hashing the
# whole block meant that tuning the HTTP back-off re-hit the rate-limited API
# for every peptide.
UNIPEPT_KEYS = ["unipept.result", "unipept.allow_http", "unipept.api_url",
                "unipept.equate_il", "unipept.split_missed_cleavages",
                "unipept.consensus_min_fraction",
                "unipept.consensus_min_peptides"]

STAGES = [
    dict(name="emapper", cost=3, enabled="eggnog",
         out=lambda p: [p.emapper],
         inp=lambda c, p: [c["proteins_faa"]] + (
             [c["emapper_precomputed"]] if isinstance(c["emapper_precomputed"], str)
             and c["emapper_precomputed"] else list(c["emapper_precomputed"] or [])),
         keys=["emapper_precomputed", "emapper_id_transform",
               "emapper_strip_id_prefix",
               "emapper_min_coverage", "db.eggnog_data"],
         deps=[], fn=stage_emapper),
    dict(name="pfam", cost=3, enabled="pfam", out=lambda p: [p.pfam],
         inp=lambda c, p: [c["proteins_faa"], c["db"]["pfam_hmm"]],
         keys=["db.pfam_hmm"], deps=[], fn=stage_pfam),
    dict(name="dbcan", cost=2, enabled="dbcan", out=lambda p: [p.dbcan],
         inp=lambda c, p: [c["proteins_faa"], c["db"]["dbcan_hmm"]],
         keys=["db.dbcan_hmm", "thresholds.dbcan_evalue"], deps=[], fn=stage_dbcan),
    dict(name="diamond", cost=2, empty_ok=True, enabled="diamond",
         out=lambda p: [p.diamond_done],
         inp=lambda c, p: [c["proteins_faa"]] + list((c["db"].get("diamond") or {}).values()),
         keys=["db.diamond", "thresholds.diamond_evalue", "diamond_evalues",
               "thresholds.diamond_min_pident", "diamond_min_pidents"],
         deps=[], fn=stage_diamond),
    dict(name="signalp", cost=3, enabled="topology", out=lambda p: [p.signalp],
         inp=lambda c, p: [c["proteins_faa"]], keys=["signalp_mode"],
         deps=[], fn=stage_signalp),
    # gpu=True: this stage takes an exclusive lease on gpu_device. tmbed held
    # 15.5 GB of a 16 GB card; see gpu_workers.
    # empty_ok: a proteome whose every sequence is over tmbed_max_len leaves
    # an empty prediction file on purpose, and a rerun must be able to adopt
    # it rather than re-deciding that there is nothing to do.
    #
    # tmbed_chunk_residues is deliberately NOT a key. It changes how the work
    # is divided and therefore the ORDER of the records, but not one
    # prediction in them, and listing it would throw away a 30-hour stage
    # because someone tuned a checkpoint size. The two keys that DO change
    # what is in the file - how much of a failure is tolerated - are listed.
    dict(name="tmbed", cost=3, enabled="topology", gpu=True, empty_ok=True,
         out=lambda p: [p.tmbed],
         inp=lambda c, p: [c["proteins_faa"]],
         keys=["gpu_device", "tmbed_use_gpu", "tmbed_max_len",
               "tmbed_batch_size", "tmbed_allow_partial",
               "tmbed_max_consecutive_failures"], deps=[], fn=stage_tmbed),
    dict(name="cluster", cost=1, enabled="cluster", out=lambda p: [p.cluster],
         inp=lambda c, p: [c["proteins_faa"]],
         keys=["thresholds.cluster_min_seq_id", "thresholds.cluster_coverage"],
         deps=[], fn=stage_cluster),
    dict(name="ncbifam", cost=3, enabled="ncbifam", out=lambda p: [p.ncbifam],
         inp=lambda c, p: [c["proteins_faa"], c["db"].get("ncbifam_hmm", "")],
         keys=["db.ncbifam_hmm", "thresholds.ncbifam_cutoff"], deps=[], fn=stage_ncbifam),
    dict(name="kofam", cost=3, enabled="kofam", out=lambda p: [p.kofam],
         inp=lambda c, p: [c["proteins_faa"], c["db"].get("kofam_ko_list", "")],
         keys=["db.kofam_profiles", "db.kofam_ko_list"], deps=[], fn=stage_kofam),
    dict(name="interpro", cost=3, enabled="interpro", out=lambda p: [p.interpro],
         inp=lambda c, p: [c["proteins_faa"]],
         keys=["interpro_applications", "db.interproscan_sh"], deps=[], fn=stage_interpro),
    dict(name="smorf", cost=1, empty_ok=True, enabled="smorf", out=lambda p: [p.smorf_faa],
         inp=lambda c, p: [c.get("contigs_fna", "")],
         keys=["smorf_mode", "thresholds.smorf_max_len"], deps=[], fn=stage_smorf),
    dict(name="context", cost=1, enabled="context", out=lambda p: [p.context],
         inp=lambda c, p: [c.get("gff") or "", p.emapper, p.pfam, p.signalp,
                           p.dbcan],
         keys=["gff", "context_window", "immunity_max_len", "immunity_max_gap",
               "pul_min_cazymes", "thresholds.dbcan_min_cov",
               "thresholds.dbcan_evalue"],
         deps=['emapper', 'pfam', 'signalp', 'dbcan'], fn=stage_context),
    dict(name="integrate", cost=2, enabled=None,
         out=lambda p: [p.pass1, p.dark, p.dark_all],
         inp=lambda c, p: [c["proteins_faa"], p.emapper, p.pfam, p.dbcan,
                           p.signalp, p.tmbed, p.cluster, p.context,
                           p.ncbifam, p.kofam, p.interpro,
                           p.diamond_done] + sorted(
                               glob.glob(f"{p.diamond_dir}/*.tsv")),
         # Every config key build_annotation reads UNCONDITIONALLY has to be
         # here, because this stage is what runs it. vfdb_category_weights and
         # foldseek_target_priority were missing: finalise listed both, but on
         # a run with no structure or profile evidence finalise takes the
         # "reusing the first pass" branch and copies annotation_pass1.tsv
         # verbatim, so the only stage that could act on the change was the
         # one the change did not invalidate. Re-weighting VFDB was a silent
         # no-op on every such re-run. The three emit_dark keys below belong
         # to this stage alone, since finalise never writes dark.faa.
         keys=["thresholds", "weights", "diamond_weights",
               "vfdb_category_weights", "foldseek_target_priority",
               "diamond_evalues",
               "diamond_min_pidents", "anchor_pfams",
               "max_dark_structures", "max_len_structure",
               "exclude_id_prefixes", "toxin_fold_patterns",
               "ncbifam_uninformative_test"],
         deps=['emapper', 'pfam', 'dbcan', 'diamond', 'signalp', 'tmbed', 'cluster', 'ncbifam', 'kofam', 'interpro', 'context'], fn=stage_integrate_pass1),
    # dark_all.faa, not dark.faa: the profile searches query the whole
    # unannotated set, the structure work-list is a GPU budget.
    dict(name="jackhmmer", cost=3, empty_ok=True, enabled="jackhmmer", out=lambda p: [p.jackhmmer],
         inp=lambda c, p: [p.dark_all, c["db"].get("jackhmmer_db", "")],
         keys=["db.jackhmmer_db", "jackhmmer_iterations",
               "thresholds.jackhmmer_evalue"], deps=['integrate'], fn=stage_jackhmmer),
    dict(name="hhblits", cost=3, empty_ok=True, enabled="hhblits", out=lambda p: [p.hhr_done],
         inp=lambda c, p: [p.dark_all, c["db"].get("hhblits_db", "")],
         keys=["db.hhblits_db", "hhblits_iterations"], deps=['integrate'], fn=stage_hhblits),
    # gpu=True: ESMFold peaked at 13.3 GB on a single short sequence, so it
    # cannot share a 16 GB card with tmbed; see gpu_workers.
    dict(name="esmfold", cost=3, empty_ok=True, enabled="structure", gpu=True,
         out=lambda p: [p.struct_done],
         inp=lambda c, p: [p.dark],
         keys=["max_len_structure", "esmfold_chunk_size",
               "esmfold_allow_partial", "esmfold_max_consecutive_failures",
               "esmfold_vram_cap", "esmfold_bytes_per_residue_pair",
               "esmfold_vram_reserve_gb"],
         deps=['integrate'], fn=stage_esmfold),
    dict(name="foldseek", cost=2, empty_ok=True, enabled="structure", out=lambda p: [p.foldseek],
         inp=lambda c, p: [p.struct_done],
         keys=["db.foldseek_target", "db.foldseek_extra_targets",
               "thresholds.foldseek_evalue", "foldseek_self_cluster",
               "thresholds.esmfold_min_plddt",
               "thresholds.foldseek_cluster_evalue",
               "thresholds.foldseek_cluster_tmscore",
               "thresholds.foldseek_cluster_coverage"],
         deps=['esmfold'], fn=stage_foldseek),
    dict(name="finalise", cost=2, enabled=None,
         out=lambda p: [p.final, p.summary, p.agreement],
         inp=lambda c, p: [p.pass1, p.foldseek, p.context, p.fold_clusters,
                           p.ncbifam, p.kofam, p.interpro,
                           p.hhr_done, p.jackhmmer],
         keys=["thresholds", "weights", "diamond_weights",
               "vfdb_category_weights",
               "diamond_evalues",
               "diamond_min_pidents", "anchor_pfams",
               "toxin_fold_patterns", "ncbifam_uninformative_test",
               "foldseek_target_priority"],
         deps=['integrate', 'jackhmmer', 'hhblits', 'foldseek', 'context'], fn=stage_integrate_final),
    dict(name="unipept", cost=3, enabled="unipept", out=lambda p: [p.unipept_lca],
         inp=lambda c, p: quant_inputs(c) + [(c.get("unipept") or {}).get("result", "")],
         keys=UNIPEPT_KEYS + ["quant_table", "quant_format", "tmt",
               "peptide_only_reader", "exclude_id_prefixes"],
         deps=[], fn=stage_unipept),
    dict(name="taxonomy", cost=1, enabled="taxonomy", out=lambda p: [p.taxonomy_comparison],
         inp=lambda c, p: [p.unipept_lca, p.final] + quant_inputs(c),
         keys=UNIPEPT_KEYS + ["db.ncbi_taxonomy", "peptide_only_reader",
               "exclude_id_prefixes", "quant_format", "tmt"],
         deps=['unipept', 'finalise'], fn=stage_taxonomy),
    dict(name="join", cost=2, enabled="join",
         out=lambda p: [f"{p.quant_dir}/annotated_quant.tsv"],
         # taxonomy_comparison is a real input: join merges its columns and
         # resolves effective_taxid from it. Omitting it left annotated_quant
         # stale whenever the taxonomy stage rebuilt its verdicts.
         inp=lambda c, p: [p.final, *quant_inputs(c), p.taxonomy_comparison,
                           c.get("manifest") or ""],
         keys=["quant_table", "quant_format", "tmt", "manifest",
               "taxonomy_source",
               "taxon_rank", "peptide_assignment", "rollup_method",
               "min_features_per_protein",
               "zero_intensity_is_missing", "exclude_id_prefixes",
               "intensity_columns", "intensity_regex",
               "feature_intensity_suffix", "feature_exclude_suffixes",
               "export_feature_quant", "taxon_min_proteins_for_factor"],
         deps=['finalise', 'taxonomy'], fn=stage_join),
]
STAGE_NAMES = [s["name"] for s in STAGES]


# How long a stage runs, coarsely. 3 = hours, 2 = minutes, 1 = seconds, and
# the numbers come off two real runs rather than intuition: on 38k proteins
# interproscan took 2.8 h, signalp 56 min, tmbed 52 min, kofam 29 min, pfam
# 27 min, ncbifam 24 min, dbcan 41 s, cluster 14 s; on 455,571 proteins
# kofam took 14.3 h, pfam 7.2 h, ncbifam 5.6 h, dbcan 10 min, diamond 5 min,
# cluster 109 s, and interproscan was still running after two days. Three
# ranks is all the resolution the scheduler can use: it decides which ready
# stage claims a worker first, not when anything finishes.
STAGE_COSTS = {st["name"]: st["cost"] for st in STAGES}


def stage_priority(name):
    """Sort key for one round's ready stages: the longest one goes first.

    Every stage with no dependencies is ready in the first round, and
    stage_workers is 4, so the first four IN TABLE ORDER started and the rest
    waited. That put cluster (109 s) and dbcan (10 min) on the box while
    interproscan — the longest stage in the pipeline by an order of magnitude
    — sat in the queue behind them. Longest-processing-time-first is the
    standard greedy answer to that, and here it costs one sort of a list that
    is never longer than 21.

    Two things this must NOT do. It must not reach the cache: scheduling
    order cannot change a stage's output, so no signature and no keys list
    mentions cost. And it must not default: indexing STAGE_COSTS raises
    KeyError on an unknown name, where a .get(name, 1) would quietly rank a
    stage added without a cost as trivial and reintroduce the exact problem.
    """
    return STAGE_COSTS[name]


def gpu_lease(ready, running, slots, needs_gpu):
    """Which of one round's ready stages may start, given what is running.

    Independent stages run concurrently and the CPU and RAM budgets are split
    between them, but the GPU was not modelled at all: tmbed held 15.5 GB of a
    16 GB card and ESMFold peaked at 13.3 GB on a single short sequence, so
    with run.topology and run.structure both on they cannot fit together and
    whichever loses dies of a CUDA OOM that names no cause. On the run this
    comes from they only avoided each other by accident, because they happened
    to be in separate invocations.

    This leases the device, not the machine: every CPU-only stage passes
    through untouched, so the pipeline is not serialised to achieve it.
    Returns (dispatch, waiting).
    """
    free = max(0, slots - sum(1 for n in running if needs_gpu(n)))
    dispatch, waiting = [], []
    for name in ready:
        if not needs_gpu(name):
            dispatch.append(name)
        elif free:
            free -= 1
            dispatch.append(name)
        else:
            waiting.append(name)
    return dispatch, waiting


def detect_ram_gb():
    """Physical RAM in GB, or 0 when it cannot be determined.

    Three probes, because no one of them covers the platforms this runs on.
    sysconf works on Linux; macOS DEFINES _SC_PHYS_PAGES but sysconf returns
    EINVAL for it, and Darwin has no /proc, so both of the first two fall
    through there and the budget silently became 0 - meaning every stage ran
    with no memory allocation at all on a Mac. hw.memsize is the Darwin answer.
    """
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        if pages and pages > 0:
            return int(os.sysconf("SC_PAGE_SIZE") * pages / 2**30)
    except (ValueError, OSError, AttributeError):
        pass
    try:
        with open("/proc/meminfo", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    return int(int(line.split()[1]) / 2**20)
    except OSError:
        pass
    if sys.platform == "darwin":
        try:
            out = subprocess.run(["sysctl", "-n", "hw.memsize"],
                                 capture_output=True, text=True, timeout=5)
            if out.returncode == 0 and out.stdout.strip().isdigit():
                return int(int(out.stdout.strip()) / 2**30)
        except (OSError, subprocess.SubprocessError):
            pass
    if os.name == "nt":
        # The pipeline stages need POSIX, but doctor, report and object all run
        # on Windows directly, and they read the same budget.
        try:
            import ctypes

            class _MEMSTAT(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong),
                            ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong),
                            ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong),
                            ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong),
                            ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

            st = _MEMSTAT()
            st.dwLength = ctypes.sizeof(_MEMSTAT)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
                return int(st.ullTotalPhys / 2**30)
        except Exception:                                   # noqa: BLE001
            pass
    return 0


def parse_ram(v):
    """Accept 64, '64', '64G', '64GB', '65536M'."""
    if v in (None, "", 0):
        return 0
    t = str(v).strip().upper().rstrip("B")
    mult = 1.0
    if t.endswith("T"):
        mult, t = 1024.0, t[:-1]
    elif t.endswith("G"):
        mult, t = 1.0, t[:-1]
    elif t.endswith("M"):
        mult, t = 1 / 1024.0, t[:-1]
    try:
        gb = float(t) * mult
    except ValueError:
        die(f"could not read a memory size from '{v}'; use e.g. 64, 64G, 512M")
    if gb < 0:
        die(f"negative memory budget '{v}'; use 0 for auto-detect")
    if gb <= 0:
        return 0
    # Round up, never down to zero: 0 means "auto-detect", so truncating
    # 512M to 0 would silently discard an explicit budget.
    return max(1, int(round(gb)))


def stage_cfg(cfg, threads, ram_gb):
    """Shallow copy with this stage's CPU and memory allocation.
    Shallow is enough: both keys are top level and nothing else is mutated."""
    c = dict(cfg)
    c["threads"] = max(1, int(threads))
    c["ram_gb"] = max(1, int(ram_gb)) if ram_gb else 0
    return c


def tool_args(cfg, tool):
    """User-supplied extra flags for a tool, appended verbatim.

    An escape hatch on purpose: this code is pinned to the CLI of a dozen
    programs whose flags move between versions, and a wrong or missing flag
    should be fixable from the config rather than by editing the tool.
    """
    extra = (cfg.get("tool_args") or {}).get(tool) or []
    if isinstance(extra, str):
        extra = shlex.split(extra)
    elif not isinstance(extra, (list, tuple)):
        die(f"tool_args.{tool} must be a command line, either a string or a "
            f"list of flags, not {type(extra).__name__}")
    return [str(x) for x in extra]


_DIGEST_MODE_LOGGED = [False]


def signature(stage, cfg, p):
    full = bool(cfg.get("full_content_digest", False))
    # Said once per run, not once per stage: which mode is in force decides
    # whether an in-place database edit can invalidate a cached stage at all.
    if not _DIGEST_MODE_LOGGED[0]:
        _DIGEST_MODE_LOGGED[0] = True
        log("input digests: " + (
            "full_content_digest is ON — every byte of every input is "
            "hashed. Large databases make this slow, and signatures differ "
            "from those written with it off, so stages will recompute once."
            if full else
            f"files over {_HASH_FULL_MAX // 2**20} MB are digested at their "
            f"first and last {_HASH_EDGE // 2**20} MB plus their size. An "
            "in-place edit to the middle of such a file is not seen; set "
            "full_content_digest: true if that matters."))
    payload = {
        "signature_version": SIGNATURE_VERSION,
        "inputs": [_stat(x, full) for x in stage["inp"](cfg, p) if x],
        "config": {k: _dig(cfg, k) for k in stage["keys"]},
        # tool_args is documented as THE way to correct a wrong flag and is
        # appended to nearly every command line, so a change to it has to
        # invalidate the cache. Hashing the whole block over-invalidates (a
        # new hmmsearch flag also reruns diamond) but never under-invalidates,
        # and a stage silently keeping the output of the flag you just fixed
        # is the worse failure.
        "tool_args": cfg.get("tool_args") or {},
    }
    return hashlib.sha1(
        json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _windows_pid_alive(pid):
    """Whether a pid names a live process on Windows, without signalling it.

    Unprovable means alive, as everywhere else in the lock: a missing API, a
    refused handle or an ambiguous exit code all answer True, and
    --force-unlock is the escape. Only ERROR_INVALID_PARAMETER -- the answer
    Windows gives for a pid that does not exist at all -- is taken as proof of
    death.
    """
    try:
        import ctypes
        from ctypes import wintypes
    except Exception:                                   # pragma: no cover
        return True
    PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
    ERROR_INVALID_PARAMETER = 87
    STILL_ACTIVE = 259
    try:
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.OpenProcess.restype = wintypes.HANDLE
        k32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL,
                                    wintypes.DWORD]
        h = k32.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False,
                            int(pid))
        if not h:
            return ctypes.get_last_error() != ERROR_INVALID_PARAMETER
        try:
            code = wintypes.DWORD()
            if not k32.GetExitCodeProcess(h, ctypes.byref(code)):
                return True
            # STILL_ACTIVE is ambiguous with a process that exited WITH code
            # 259, which is why it errs towards alive rather than away.
            return code.value == STILL_ACTIVE
        finally:
            k32.CloseHandle(h)
    except Exception:                                   # pragma: no cover
        return True


class ResultsLock:
    """One writer per results directory.

    Two concurrent runs write the same files and the same state, and the
    damage is silent: both report success and the outputs are interleaved.
    """

    def __init__(self, path, force=False, empty_grace=60.0):
        self.path = path
        self.force = force
        self.empty_grace = float(empty_grace)
        self.held = False

    def _holder_is_alive(self, info):
        """Whether the process named in a lock file may still be running.

        Unprovable means alive: a lock we cannot disprove is treated as a live
        run, because trampling one is silent corruption while refusing is a
        message and --force-unlock. That covers a lock written on another node
        of the shared array (we cannot see its process table), a pid owned by
        another user (os.kill raises PermissionError, which is an OSError and
        used to read as 'dead'), and a garbled lock file.
        """
        pid, host = info.get("pid"), info.get("host", "")
        if not isinstance(pid, int) or host != socket.gethostname():
            return True
        if os.name == "nt":
            # os.kill(pid, 0) on Windows calls TerminateProcess: asking
            # whether a process is alive that way would KILL it. OpenProcess
            # asks without touching it. Answering "alive" unconditionally, as
            # this did, meant a lock left by a crashed run on Windows could
            # never be reclaimed and every resume needed --force-unlock -- and
            # a crash is exactly when reclaiming has to work.
            return _windows_pid_alive(pid)
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except OSError:
            return True                   # PermissionError: it exists
        return True

    def _empty_and_settled(self):
        """A zero-byte lock old enough that no live writer could still owe it.

        `lock_empty_grace_s` is the width of the only window in which an empty
        lock is legitimate: between another run's O_EXCL create and its write.
        That is microseconds of work, so a minute is already enormous slack,
        and anything past it is a corpse.
        """
        try:
            st = os.stat(self.path)
        except OSError:
            return False
        return (st.st_size == 0
                and time.time() - st.st_mtime > self.empty_grace)

    def __enter__(self):
        os.makedirs(os.path.dirname(self.path) or ".", exist_ok=True)
        payload = json.dumps({"pid": os.getpid(), "host": socket.gethostname(),
                              "started": time.strftime("%Y-%m-%dT%H:%M:%S")})
        # O_EXCL, not exists()-then-write: two runs launched in the same second
        # both used to see no lock and both proceed.
        for _ in range(3):
            try:
                fd = os.open(self.path,
                             os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            except FileExistsError:
                try:
                    with open(self.path, encoding="utf-8") as fh:
                        info = json.load(fh)
                except (OSError, ValueError):
                    info = {}
                pid, host = info.get("pid"), info.get("host", "")
                if self._empty_and_settled() and not self.force:
                    # A zero-byte lock is not a garbled lock, it is one whose
                    # writer died between the O_EXCL create and the write. A
                    # real writer closes that gap in microseconds, so an empty
                    # file that has sat unchanged for minutes cannot be a live
                    # run - it is what a power loss, an OOM kill or a host
                    # bugcheck leaves behind. Treating it as "unprovable, so
                    # assume alive" made a crashed run permanently
                    # unresumable without --force-unlock, which is exactly
                    # backwards: a hard crash is when reclaiming has to work.
                    # The age check is what keeps the genuine race safe.
                    log(f"removing a zero-byte lock left at "
                        f"{time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(os.path.getmtime(self.path)))}"
                        "; its writer died between creating the file and "
                        "writing to it, which is what a crash or a power loss "
                        "leaves behind", "WARN")
                    try:
                        os.remove(self.path)
                    except OSError:
                        pass
                    continue
                if self._holder_is_alive(info) and not self.force:
                    die(f"another metaannot is already running here "
                        f"(pid {pid} on {host or '?'}, started "
                        f"{info.get('started','?')}).\n"
                        "Two runs sharing a results directory corrupt each "
                        "other. Wait for it, or use --force-unlock if you are "
                        "certain it is gone.")
                log(f"removing a stale lock from pid {pid} on {host or '?'}",
                    "WARN")
                try:
                    os.remove(self.path)
                except OSError:
                    pass
                continue
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(payload)
            self.held = True
            return self
        die(f"could not take the results lock {self.path}: it keeps being "
            "recreated, so another metaannot is starting at the same moment.")

    def __exit__(self, *exc):
        if self.held:
            try:
                os.remove(self.path)
            except OSError:
                pass
        return False


class _State(dict):
    """The stage records, plus whether we actually managed to read them.

    A run that could not read the state file knows nothing about the outputs
    lying in the results directory, so it must not adopt them as its own.
    """
    unreadable = False


def load_state(path):
    if not os.path.exists(path):
        return _State()
    try:
        with open(path, encoding="utf-8") as fh:
            return _State(json.load(fh))
    except (OSError, ValueError) as e:
        log(f"state file unreadable ({e}); every stage will be recomputed, "
            "and no existing output will be adopted, because there is no "
            "record of how it was made", "WARN")
        st = _State()
        st.unreadable = True
        return st


def _mtime(path):
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0.0


def _file_note(path):
    """Path, size and mtime — what adoption is being asked to trust."""
    try:
        st = os.stat(path)
    except OSError:
        return f"{path} (missing)"
    return (f"{path} ({st.st_size} bytes, "
            + time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime))
            + ")")


def _stamp_epoch(s):
    """Seconds since the epoch for a '%Y-%m-%dT%H:%M:%S' state stamp."""
    try:
        return time.mktime(time.strptime(str(s), "%Y-%m-%dT%H:%M:%S"))
    except (ValueError, TypeError):
        return None


_STATELOCK = threading.Lock()


def _rmd_bin_vocabulary(text):
    """(levels, colour keys) as the embedded Rmd spells them."""
    def names(var, pat):
        m = re.search(var + r" <- c\((.*?)\)", text, re.S)
        return tuple(re.findall(pat, m.group(1))) if m else ()
    return (names("BIN_LEVELS", r'"([^"]+)"'),
            names("BIN_COLS", r'"([^"]+)"\s*='))


def _check_bin_vocabulary():
    """The report carries its own copy of the bin list. Nothing at run time
    would notice a drift: a bin missing from BIN_LEVELS becomes NA in every
    figure, and a bin missing from BIN_COLS drops out of the legend. Both are
    silent, so the copy is checked here, once, at import.
    """
    lv, cl = _rmd_bin_vocabulary(RMD_TEMPLATE)
    if lv != BIN_ORDER or cl != BIN_ORDER:
        raise AssertionError(
            "the embedded report's bin vocabulary has drifted from "
            f"BIN_ORDER {BIN_ORDER}: BIN_LEVELS={lv}, BIN_COLS keys={cl}. "
            "Edit both to match BIN_ORDER (and the module docstring).")


def save_state(path, state):
    # Locked and written via a temp file: stages finish concurrently, and a
    # torn state file would silently invalidate every cached stage.
    with _STATELOCK:
        tmp = f"{path}.{os.getpid()}.{threading.get_ident()}.tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(dict(state), fh, indent=1, sort_keys=True)
        os.replace(tmp, path)


# ======================================================================
# the R report, carried inside this file
# ======================================================================
ROBJECT_SCRIPT = r"""#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# build_object.R — assemble metaannot's output into one R object.
#
#   QFeatures if available: a peptide assay and a protein assay, linked by
#   the assignment metaannot made, with the full annotation as rowData and
#   the design as colData.
#   SummarizedExperiment if QFeatures is not installed, or if the input was
#   protein-level and there is no peptide assay to link.
#   A plain list if neither is installed, so the data is never trapped.
#
# Everything is read from the TSVs the pipeline wrote. Nothing here recomputes
# anything: if a number differs from the tables, that is a bug, not a variant.
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(utils); library(stats)
})

args <- commandArgs(trailingOnly = TRUE)
RD  <- if (length(args) >= 1) args[1] else "results"
OUT <- if (length(args) >= 2) args[2] else file.path(RD, "metaannot.rds")
# The report writes its tables to results/<analysis.out_subdir>, so the object
# has to look for them in the same place. Hard-coding "analysis" here made a
# non-default out_subdir drop the DE fold-in and the risk table in silence.
SUB <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "analysis"

msg  <- function(...) cat(sprintf(...), "\n", sep = "")
rd   <- function(...) file.path(RD, ...)
read_if <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  x <- try(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                             ...), silent = TRUE)
  if (inherits(x, "try-error")) { msg("  could not read %s", path); return(NULL) }
  x
}

# ---- inputs ----------------------------------------------------------
aq <- read_if(rd("quant", "annotated_quant.tsv"))
if (is.null(aq))
  stop("no annotated_quant.tsv in ", RD, "; run the join stage first",
       call. = FALSE)
# annotation_final.tsv is deliberately not read: annotated_quant.tsv already
# carries the annotation of every quantified group, and reading the full table
# cost minutes and gigabytes for a value nothing consumed.
design  <- read_if(rd("quant", "design_from_input.tsv"))
feat    <- read_if(rd("quant", "feature_quant.tsv"))
evid    <- read_if(rd("quant", "peptide_evidence.tsv"))
sizef   <- read_if(rd("quant", "taxon_size_factors.tsv"))
taxint  <- read_if(rd("quant", "taxon_intensity.tsv"))
conf    <- read_if(rd("quant", "group_conflicts.tsv"))
taxcmp  <- read_if(rd("unipept", "taxonomy_comparison.tsv"))
binsum  <- read_if(rd("bin_summary.tsv"))
t2gfam  <- read_if(rd("quant", "term2gene_family.tsv"))
t2gpfam <- read_if(rd("quant", "term2gene_pfam.tsv"))
risk    <- read_if(rd(SUB, "taxon_normalisation_risk.tsv"))
state   <- if (file.exists(rd(".metaannot_state.json")))
  paste(readLines(rd(".metaannot_state.json"), warn = FALSE), collapse = "\n") else NULL

# differential abundance, if the report has been run
de_files <- list.files(rd(SUB), pattern = "^differential_abundance_.*\\.tsv$",
                       full.names = TRUE)
if (!length(de_files))
  msg("no differential_abundance_*.tsv in %s; the object will carry no DE columns",
      rd(SUB))
de <- setNames(lapply(de_files, read_if),
               sub("^differential_abundance_", "", tools::file_path_sans_ext(basename(de_files))))

# ---- which columns are samples ---------------------------------------
# Positively identified, never guessed from "which columns look numeric".
# The annotation carries dozens of numeric and logical columns, and one of
# them being mistaken for a sample would silently corrupt every assay.
FEATMETA <- c("feature_id", "peptide", "razor_protein", "assignment_class",
              "assigned_protein", "candidates")
SFMETA   <- c("effective_taxid", "seed_taxid", "n_proteins", "method")
# Numeric columns metaannot or FragPipe add that are NOT sample intensities.
# Kept identical to the report's list: two components disagreeing about what a
# sample is would build an assay and a model on different matrices.
ANNOT_NUM <- c("length", "export_score", "effector_score", "n_tmh", "n_tmb",
               "n_pathway_specific", "foldseek_prob", "foldseek_tm",
               "n_members", "member_rank", "fold_cluster_size",
               "n_features_used", "n_unique", "n_taxon_unique",
               # n_family_unique is written only under
               # peptide_assignment=taxon_or_family_unique, but it must be
               # declared unconditionally: samples are guessed from column
               # types whenever quant/sample_columns.txt is absent, and an
               # undeclared integer counter would be built into the assay.
               "n_family_unique",
               "n_features_dropped", "n_peptides_with_lca",
               "unipept_agreement", "hh_prob", "taxon_factor_log2FC",
               "taxon_frac_deviating",
               # taxids: read.delim/readr re-type "816" as an integer even
               # though the pipeline wrote them as text.
               "taxid", "seed_taxid", "effective_taxid", "unipept_taxid",
               # FragPipe protein-level metadata, numeric but not intensity.
               "Protein Length", "Coverage", "Protein Probability",
               "Top Peptide Probability", "Combined Total Peptides",
               "Combined Spectral Count", "Combined Unique Spectral Count",
               "Combined Total Spectral Count", "Total Peptides",
               "Unique Peptides", "Razor Peptides", "Total Spectral Count",
               "Unique Spectral Count", "Razor Spectral Count",
               "Length", "Protein Qvalue", "Total Intensity",
               "Unique Intensity", "Razor Intensity",
               "NumberPSM", "MaxPepProb", "ReferenceIntensity")

# Written by the join stage when it knows exactly which columns it treated as
# intensities. Nothing else is as authoritative, so it comes first.
sc_path <- rd("quant", "sample_columns.txt")
sample_cols <- if (file.exists(sc_path))
  trimws(readLines(sc_path, warn = FALSE)) else character(0)
sample_cols <- sample_cols[nzchar(sample_cols)]

sample_source <- NULL
samples <- NULL
if (length(sample_cols)) {
  samples <- intersect(sample_cols, names(aq))
  sample_source <- "quant/sample_columns.txt"
  gone <- setdiff(sample_cols, names(aq))
  if (length(gone))
    msg(paste("WARNING: %d column(s) named in sample_columns.txt are absent",
              "from annotated_quant.tsv: %s"), length(gone),
        paste(utils::head(gone, 8), collapse = ", "))
} else if (!is.null(design) && "sample" %in% names(design)) {
  samples <- intersect(as.character(design$sample), names(aq))
  sample_source <- "design_from_input.tsv"
  if (!length(samples))
    stop("design_from_input.tsv names ", nrow(design), " sample(s) and none of ",
         "them is a column of annotated_quant.tsv.\n  design: ",
         paste(utils::head(as.character(design$sample), 8), collapse = ", "),
         "\n  quant:  ", paste(utils::head(names(aq), 8), collapse = ", "),
         "\n  The two came from the same run, so this is a naming mismatch ",
         "(channel names\n  against sample names, or a stale design file), not ",
         "something to guess around.", call. = FALSE)
} else if (!is.null(sizef)) {
  samples <- intersect(setdiff(names(sizef), SFMETA), names(aq))
  sample_source <- "taxon_size_factors.tsv"
} else if (!is.null(feat)) {
  samples <- intersect(setdiff(names(feat), FEATMETA), names(aq))
  sample_source <- "feature_quant.tsv"
}
if (is.null(samples) || !length(samples)) {
  # Last resort. Warn loudly: this is the only path that can misfire.
  num <- names(aq)[vapply(aq, function(x) is.numeric(x) && !is.logical(x),
                          logical(1))]
  samples <- setdiff(num, c(ANNOT_NUM, grep("_pident$", names(aq), value = TRUE)))
  sample_source <- "inferred from column types"
  msg("WARNING: no sample-column list, design, size-factor or feature file to")
  msg("  name the samples. Falling back to a type heuristic, which can")
  msg("  misclassify an annotation column. Check the list below, and prefer")
  msg("  supplying a manifest.")
}
if (!length(samples))
  stop("could not identify any sample columns in annotated_quant.tsv",
       call. = FALSE)
# Silently dropping a channel the user expects is worse than a noisy object.
if (!is.null(design) && "sample" %in% names(design)) {
  unmatched <- setdiff(as.character(design$sample), samples)
  if (length(unmatched))
    msg("WARNING: %d design row(s) match no quant column and are dropped: %s",
        length(unmatched), paste(utils::head(unmatched, 8), collapse = ", "))
}
cand_num <- setdiff(names(aq)[vapply(aq, function(x) is.numeric(x) &&
                                       !is.logical(x), logical(1))],
                    c(ANNOT_NUM, grep("_pident$", names(aq), value = TRUE)))
extra <- setdiff(cand_num, samples)
if (length(extra))
  msg(paste("WARNING: %d numeric column(s) look like intensities but are not",
            "in the sample list, so they stay in rowData: %s"), length(extra),
      paste(utils::head(extra, 8), collapse = ", "))
msg("samples (%s): %s", sample_source,
    paste(utils::head(samples, 8), collapse = ", "))
msg("proteins: %d   samples: %d", nrow(aq), length(samples))

prot_id <- if ("group_id" %in% names(aq)) "group_id" else names(aq)[1]
rn <- make.unique(as.character(aq[[prot_id]]))

prot_mat <- as.matrix(aq[, samples, drop = FALSE])
storage.mode(prot_mat) <- "double"
rownames(prot_mat) <- rn

prot_row <- aq[, setdiff(names(aq), samples), drop = FALSE]
rownames(prot_row) <- rn
# attach differential abundance as extra rowData columns, prefixed per contrast
for (nm in names(de)) {
  d <- de[[nm]]
  if (is.null(d) || !"group_id" %in% names(d)) next
  # Only the statistics. The differential table repeats the annotation it was
  # built from, and copying that back in would give every column a second,
  # prefixed twin with no indication which is authoritative.
  keep <- setdiff(names(d), c("group_id", names(prot_row)))
  if (!length(keep)) next
  add <- d[match(rn, as.character(d$group_id)), keep, drop = FALSE]
  names(add) <- paste0("DE.", nm, ".", names(add))
  prot_row <- cbind(prot_row, add)
  msg("  %s: %d statistic column(s) added", nm, length(keep))
}

col_df <- if (!is.null(design)) {
  design[match(samples, as.character(design$sample)), , drop = FALSE]
} else data.frame(sample = samples, stringsAsFactors = FALSE)
rownames(col_df) <- samples

meta <- list(
  metaannot     = list(results_dir = normalizePath(RD, mustWork = FALSE),
                       built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                       state = state),
  bin_summary   = binsum,
  taxon_size_factors = sizef,
  taxon_intensity    = taxint,
  group_conflicts    = conf,
  taxonomy_comparison = taxcmp,
  peptide_evidence   = evid,
  normalisation_risk = risk,
  term2gene     = list(family = t2gfam, pfam = t2gpfam),
  differential  = de
)

have <- function(pkg) requireNamespace(pkg, quietly = TRUE)

obj <- NULL
if (have("SummarizedExperiment")) {
  SE <- SummarizedExperiment::SummarizedExperiment
  se_prot <- SE(assays = list(intensity = prot_mat),
                rowData = S4Vectors::DataFrame(prot_row, check.names = FALSE),
                colData = S4Vectors::DataFrame(col_df, check.names = FALSE))
  # meta is attached once, to whatever object is finally saved. Attaching it
  # here as well doubled the size of every .rds, since the QFeatures container
  # holds se_prot and carries its own copy of the same list.

  se_pep <- NULL
  if (!is.null(feat)) {
    fsamp <- intersect(samples, names(feat))
    if (length(fsamp)) {
      fm <- as.matrix(feat[, fsamp, drop = FALSE])
      storage.mode(fm) <- "double"
      frn <- make.unique(as.character(feat$feature_id))
      rownames(fm) <- frn
      # pad to the protein assay's samples so the two assays share colData
      miss <- setdiff(samples, fsamp)
      if (length(miss)) {
        pad <- matrix(NA_real_, nrow(fm), length(miss),
                      dimnames = list(frn, miss))
        fm <- cbind(fm, pad)[, samples, drop = FALSE]
      }
      frow <- feat[, setdiff(names(feat), fsamp), drop = FALSE]
      rownames(frow) <- frn
      se_pep <- SE(assays = list(intensity = fm),
                   rowData = S4Vectors::DataFrame(frow, check.names = FALSE),
                   colData = S4Vectors::DataFrame(col_df, check.names = FALSE))
      msg("peptide assay: %d features", nrow(fm))
    }
  }

  if (!is.null(se_pep) && have("QFeatures")) {
    obj <- try({
      q <- QFeatures::QFeatures(list(peptides = se_pep, proteins = se_prot),
                                colData = S4Vectors::DataFrame(col_df,
                                                               check.names = FALSE))
      # Record the aggregation metaannot actually performed, so the link in
      # the object is the one used, not one re-derived here.
      q <- try(QFeatures::addAssayLink(q, from = "peptides", to = "proteins",
                                       varFrom = "assigned_protein",
                                       varTo = prot_id), silent = TRUE)
      if (inherits(q, "try-error")) {
        # The real condition, not a guess: addAssayLink fails for several
        # reasons (no assigned_protein matches a group_id gives "No match
        # found"), and printing "version differs" for all of them sent people
        # after the wrong problem.
        msg("  assay link not added (%s); assays are present but unlinked",
            trimws(conditionMessage(attr(q, "condition"))))
        q <- QFeatures::QFeatures(list(peptides = se_pep, proteins = se_prot),
                                  colData = S4Vectors::DataFrame(col_df,
                                                                 check.names = FALSE))
      }
      S4Vectors::metadata(q) <- meta
      q
    }, silent = TRUE)
    if (inherits(obj, "try-error")) {
      msg("  QFeatures failed (%s), using SummarizedExperiment",
          trimws(conditionMessage(attr(obj, "condition"))))
      obj <- NULL
    }
  }
  if (is.null(obj)) {
    obj <- se_prot
    S4Vectors::metadata(obj) <- meta
    if (!is.null(se_pep)) S4Vectors::metadata(obj)$peptides <- se_pep
    msg("class: SummarizedExperiment%s",
        if (!is.null(se_pep)) " (peptide assay in metadata()$peptides)" else "")
  } else {
    msg("class: QFeatures with assays %s",
        paste(names(obj), collapse = ", "))
  }
} else {
  msg("SummarizedExperiment is not installed; writing a plain list.")
  msg("  install it with BiocManager::install(c('QFeatures','SummarizedExperiment'))")
  obj <- list(protein_intensity = prot_mat, protein_annotation = prot_row,
              design = col_df, features = feat, metadata = meta)
  class(obj) <- c("metaannot_list", "list")
}

saveRDS(obj, OUT)
msg("wrote %s (%.1f MB)", OUT, file.info(OUT)$size / 2^20)
if (!is.null(binsum)) { msg("bins:"); print(binsum) }
"""


RMD_TEMPLATE = r"""---
title: "metaannot: differential abundance in the non-KEGG metaproteome"
date: "`r format(Sys.Date(), '%d %B %Y')`"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 3
    code_folding: hide
    df_print: paged
    fig_width: 7
    fig_height: 4.5
params:
__METAANNOT_PARAMS__
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE,
                      fig.align = "center", dpi = 120)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(ggplot2); library(purrr); library(limma)
})
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

# BIN_ORDER in metaannot.py, verbatim: the order the classifier tests the
# bins in, so a figure legend reads in the same order as bin_summary.tsv.
# Each bin keeps the colour it has always had; only the order changed.
BIN_LEVELS <- c("1_ko_pathway", "2_ko_orphan", "3_annotated_no_ko",
                "3d_duf_only", "3s_structure_only", "3p_profile_only", "4_dark")
BIN_COLS <- c("1_ko_pathway" = "#4C72B0", "2_ko_orphan" = "#55A868",
              "3_annotated_no_ko" = "#C44E52", "3d_duf_only" = "#8172B2",
              "3s_structure_only" = "#CCB974", "3p_profile_only" = "#DD8452",
              "4_dark" = "#64676B")

RD  <- params$results_dir
OUT <- file.path(RD, params$out_subdir)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# opts_chunk suppresses warning() and message(), so QC output uses these.
note <- function(...) cat(paste0("NOTE  ", sprintf(...), "\n"))
gate <- function(...) cat(paste0("GATE  ", sprintf(...), "\n"))

# metaannot suffixes annotation columns with _ann when a quant column of the
# same name exists, so never address a column by a bare literal name.
# Named col_or_na, not pick: dplyr exports pick() as a tidyselect helper and
# intercepts it inside mutate()/transmute()/summarise(), so a helper by that
# name is silently never called.
# Returning NULL would be worse than useless: transmute(x = NULL) silently
# DROPS the column and the next filter() fails with "object not found".
col_or_na <- function(df, name, fill = NA) {
  if (name %in% names(df)) return(df[[name]])
  alt <- paste0(name, "_ann")
  if (alt %in% names(df)) return(df[[alt]])
  note("column absent from metaannot output, filled with NA: %s", name)
  rep(fill, nrow(df))
}
need <- function(path) {
  if (!file.exists(path)) stop("missing metaannot output: ", path, call. = FALSE)
  path
}
# readr guesses a column's type from 1000 rows and types an all-empty sample
# as logical; every later value is then parsed to NA. The rare-evidence
# columns this whole report is about (fold_cluster, hh_hit, jackhmmer_hit,
# foldseek_desc, pred_*) are non-empty in far fewer than 1 row in 1000, so
# they would silently arrive blank — and opts_chunk warning=FALSE hides the
# parsing warning that would have said so. Guess from every row instead, and
# print readr's own problem table.
# chr_cols is typed up front rather than converted afterwards: a taxid guessed
# numeric renders as "1e+06" once as.character() gets to it, and every
# taxonomy join then misses silently.
read_tsv_full <- function(path, chr_cols = character(0)) {
  hdr <- names(readr::read_tsv(path, n_max = 0, show_col_types = FALSE))
  present <- intersect(chr_cols, hdr)
  ct <- if (length(present))
    do.call(readr::cols, setNames(rep(list(readr::col_character()),
                                      length(present)), present))
  else readr::cols()
  x <- readr::read_tsv(path, show_col_types = FALSE, col_types = ct,
                       guess_max = .Machine$integer.max)
  pb <- readr::problems(x)
  if (nrow(pb)) {
    note("%d parsing problem(s) in %s; first rows below", nrow(pb), basename(path))
    print(utils::head(pb, 5))
  }
  x
}
csv_param <- function(s) {
  v <- trimws(unlist(strsplit(s %||% "", ",")))
  v[nzchar(v)]
}
`%||%` <- function(a, b) if (is.null(a)) b else a
```

# Experimental design

Everything about conditions, covariates and blocking comes from **your**
metadata file and the `design_formula`, `factor_cols`, `numeric_cols`,
`block_col` and `contrasts` parameters at the top of this document. The
annotation pipeline knows nothing about your experiment; this is where it
gets told.

## Quant matrix and metadata

```{r read-quant}
ann_path   <- need(file.path(RD, "annotation_final.tsv"))
quant_path <- need(file.path(RD, "quant", "annotated_quant.tsv"))

as_chr <- c("protein_id", "group_id", "seed_taxid", "effective_taxid",
            "ko", "family_id", "pfam_hits", "bin", "og", "ec", "cazy",
            "dbcan_hits", "kofam_ko", "ncbifam_hits", "interpro_ipr",
            "taxonomy_verdict", "unipept_taxid")
ann <- read_tsv_full(ann_path,   as_chr) %>%
  mutate(across(any_of(as_chr), as.character))
aq  <- read_tsv_full(quant_path, as_chr) %>%
  mutate(across(any_of(as_chr), as.character))

# Columns metaannot or FragPipe add that are numeric but are not sample
# intensities. Kept in step with build_object.R's ANNOT_NUM: if the two lists
# disagree, the object and the model are built on different matrices.
ANNOT_NUMERIC <- c("length", "export_score", "effector_score", "n_tmh", "n_tmb",
                   "n_pathway_specific", "foldseek_prob", "foldseek_tm",
                   "n_members", "member_rank", "fold_cluster_size",
                   "n_features_used", "n_unique", "n_taxon_unique",
                   # see build_object.R: written only under
                   # peptide_assignment=taxon_or_family_unique, declared
                   # always, because samples can be guessed from column types.
                   "n_family_unique",
                   "n_features_dropped", "n_peptides_with_lca",
                   "unipept_agreement", "hh_prob", "taxon_factor_log2FC",
                   "taxon_frac_deviating",
                   "taxid", "seed_taxid", "effective_taxid", "unipept_taxid",
                   "Protein Length", "Coverage", "Protein Probability",
                   "Top Peptide Probability", "Combined Total Peptides",
                   "Combined Spectral Count", "Combined Unique Spectral Count",
                   "Combined Total Spectral Count", "Total Peptides",
                   "Unique Peptides", "Razor Peptides", "Total Spectral Count",
                   "Unique Spectral Count", "Razor Spectral Count",
                   "Length", "Protein Qvalue", "Total Intensity",
                   "Unique Intensity", "Razor Intensity",
                   "NumberPSM", "MaxPepProb", "ReferenceIntensity")

# Which columns are samples is answered positively wherever metaannot recorded
# it, and only guessed from column types as a last resort: a metadata skeleton
# that lists ReferenceIntensity or Protein Length as a "sample" invites the
# user to give it a group, and nothing downstream would object.
sample_src <- "column types"
sc_path <- file.path(RD, "quant", "sample_columns.txt")
sf_path0 <- file.path(RD, "quant", "taxon_size_factors.tsv")
feat_path0 <- file.path(RD, "quant", "feature_quant.tsv")
design_path0 <- file.path(RD, "quant", "design_from_input.tsv")
candidate_samples <- character(0)
if (file.exists(sc_path)) {
  candidate_samples <- intersect(
    trimws(readLines(sc_path, warn = FALSE)), names(aq))
  sample_src <- "quant/sample_columns.txt"
} else if (file.exists(design_path0)) {
  candidate_samples <- intersect(
    as.character(read_tsv_full(design_path0, "sample")$sample), names(aq))
  sample_src <- "quant/design_from_input.tsv"
} else if (file.exists(sf_path0)) {
  candidate_samples <- intersect(
    setdiff(names(read_tsv_full(sf_path0, c("effective_taxid", "seed_taxid"))),
            c("effective_taxid", "seed_taxid", "n_proteins", "method")),
    names(aq))
  sample_src <- "quant/taxon_size_factors.tsv"
} else if (file.exists(feat_path0)) {
  candidate_samples <- intersect(
    setdiff(names(read_tsv_full(feat_path0, "assigned_protein")),
            c("feature_id", "peptide", "razor_protein", "assignment_class",
              "assigned_protein", "candidates")),
    names(aq))
  sample_src <- "quant/feature_quant.tsv"
}
if (!length(candidate_samples)) {
  candidate_samples <- setdiff(
    names(aq)[vapply(aq, is.numeric, logical(1))],
    c(ANNOT_NUMERIC, grep("_pident$", names(aq), value = TRUE)))
  sample_src <- "column types (last resort)"
  note(paste("no sample-column list, design, size-factor or feature file:",
             "the sample columns below were guessed from column types.\n",
             "     Check them before filling in any metadata."))
}
note("sample columns identified from %s: %d", sample_src, length(candidate_samples))

# design_notes.txt is written by the isobaric reader and by nothing else, so
# its presence is what tells the report that these columns are TMT reporter
# channels rather than LC-MS runs. That changes three things below: the
# missingness is plex-shaped (see min_plexes), the plex has to be named as a
# batch wherever the design is printed, and the fold changes are compressed.
DESIGN_NOTES <- file.path(RD, "quant", "design_notes.txt")
IS_ISOBARIC  <- file.exists(DESIGN_NOTES)
```

```{r metadata-template}
design_auto <- file.path(RD, "quant", "design_from_input.tsv")
# Path comparison, not string comparison: write_report_rmd absolutises the
# parameter, and on Windows that gives backslashes against file.path()'s
# forward ones.
same_path <- function(a, b)
  identical(normalizePath(a, winslash = "/", mustWork = FALSE),
            normalizePath(b, winslash = "/", mustWork = FALSE))
# A typo in analysis.metadata used to fall through to the manifest-derived
# design and produce a complete, plausible report of the wrong conditions.
# Only that silent substitution is refused; with no design to fall back on the
# skeleton is still written below, which is how a metadata file gets started.
if (!file.exists(params$metadata) && !same_path(params$metadata, design_auto) &&
    file.exists(design_auto))
  stop("analysis.metadata points at a file that does not exist:\n  ",
       params$metadata,
       "\n\nIt was set explicitly, so quietly analysing the design metaannot ",
       "recovered from\nyour input instead would report conditions you did not ",
       "ask for. Either fix the\npath, or start from the recovered design:\n  ",
       "cp ", design_auto, " ", params$metadata, call. = FALSE)
meta_path <- if (file.exists(params$metadata)) params$metadata else
  if (file.exists(design_auto)) design_auto else NA_character_

if (!is.na(meta_path) && identical(meta_path, design_auto)) {
  note(paste("no %s, using the design metaannot recovered from your input: %s.",
             "\n      It has sample/group/replicate only — to add covariates,",
             "copy it to\n      %s, add columns, and set design_formula."),
       params$metadata, design_auto, params$metadata)
}

if (is.na(meta_path)) {
  dir.create(dirname(params$metadata), recursive = TRUE, showWarnings = FALSE)
  tmpl <- tibble(
    sample = candidate_samples,
    group  = NA_character_,      # your condition — required
    batch  = NA_character_,      # optional covariate, delete if unused
    cage   = NA_character_,      # optional block, delete if unused
    sex    = NA_character_,
    age    = NA_real_)
  write_tsv(tmpl, params$metadata)
  stop("no metadata found, so a skeleton was written to:\n  ", params$metadata,
       "\n\nIt lists the ", length(candidate_samples),
       " sample columns detected in annotated_quant.tsv.",
       "\nFill in `group` (and any covariates you want to model), delete the ",
       "columns you\ndo not need, then set design_formula / factor_cols / ",
       "numeric_cols / contrasts\nin the YAML header and knit again.",
       call. = FALSE)
}

meta <- read_tsv_full(meta_path, params$sample_col)
if (!params$sample_col %in% names(meta))
  stop(meta_path, " has no column '", params$sample_col, "'; it has: ",
       paste(names(meta), collapse = ", "), call. = FALSE)
meta <- rename(meta, .sample = all_of(params$sample_col)) %>%
  mutate(.sample = as.character(.sample))

# Typing is explicit on purpose. A numeric batch code left untyped becomes a
# continuous covariate and quietly fits a slope through your batches.
for (cc in csv_param(params$factor_cols)) {
  if (!cc %in% names(meta)) stop("factor_cols names a missing column: ", cc, call. = FALSE)
  meta[[cc]] <- factor(meta[[cc]])
}
for (cc in csv_param(params$numeric_cols)) {
  if (!cc %in% names(meta)) stop("numeric_cols names a missing column: ", cc, call. = FALSE)
  meta[[cc]] <- as.numeric(meta[[cc]])
}
untyped <- setdiff(names(meta),
                   c(".sample", csv_param(params$factor_cols),
                     csv_param(params$numeric_cols)))
untyped <- untyped[untyped %in% all.vars(as.formula(params$design_formula))]
if (length(untyped))
  note("used in the formula but not typed, taken as-is: %s",
       paste(untyped, collapse = ", "))

# Match metadata samples to quant columns. DIA-NN writes full file paths as
# column names, so fall back to basename-without-extension.
idx <- match(meta$.sample, names(aq))
if (anyNA(idx)) {
  key <- tools::file_path_sans_ext(basename(names(aq)))
  idx[is.na(idx)] <- match(meta$.sample[is.na(idx)], key)
}
if (anyNA(idx))
  stop("these metadata samples match no column in annotated_quant.tsv:\n  ",
       paste(meta$.sample[is.na(idx)], collapse = "\n  "),
       "\n\ndetected sample columns:\n  ",
       paste(candidate_samples, collapse = "\n  "), call. = FALSE)
int_cols <- names(aq)[idx]
meta$.col <- int_cols

not_numeric <- int_cols[!vapply(aq[int_cols], is.numeric, logical(1))]
if (length(not_numeric))
  stop("these matched columns are not numeric, so the intensity matrix would ",
       "be coerced to character:\n  ", paste(not_numeric, collapse = "\n  "),
       call. = FALSE)

cat(sprintf("%d protein groups, %d samples\n", nrow(aq), nrow(meta)))
meta %>% select(-.col) %>% head(20)
```

## The isobaric design

```{r isobaric-design}
# The per-sample plex, from ONE source for the whole document: the model, the
# plex filter and this printout must not be able to disagree about which
# batch a sample was in. meta first (the user may have written it themselves),
# then the design metaannot recovered from the input.
PLEX_OF <- NULL
if (IS_ISOBARIC) {
  if ("plex" %in% names(meta)) {
    PLEX_OF <- setNames(as.character(meta$plex), meta$.col)
    plex_src <- meta_path
  } else if (file.exists(design_path0)) {
    d0 <- read_tsv_full(design_path0, c("sample", "plex"))
    if (all(c("sample", "plex") %in% names(d0))) {
      PLEX_OF <- setNames(as.character(d0$plex),
                          as.character(d0$sample))[meta$.sample]
      names(PLEX_OF) <- meta$.col
      plex_src <- design_path0
    }
  }
  if (!is.null(PLEX_OF) && anyNA(PLEX_OF)) {
    note(paste("%d sample(s) have no plex in %s, so the plex is not usable",
               "as a batch here"), sum(is.na(PLEX_OF)), plex_src)
    PLEX_OF <- NULL
  }
}

if (!IS_ISOBARIC) {
  cat("Not an isobaric run: no quant/design_notes.txt, so the sample columns",
      "are\nLC-MS runs and there is no plex.\n")
} else {
  cat("metaannot read isobaric (TMT) input. What it did, verbatim from\n",
      DESIGN_NOTES, ":\n\n", sep = "")
  # Printed, not summarised: the reference treatment, where the condition came
  # from and whether the channels were median-centred are choices the numbers
  # cannot be read without, and the log they were made in is long gone by the
  # time anyone reads the report.
  cat(paste0("    ", readLines(DESIGN_NOTES, warn = FALSE)), sep = "\n")
  cat("\n")
  ref_line <- grep("^reference", readLines(DESIGN_NOTES, warn = FALSE),
                   value = TRUE)
  if (length(ref_line))
    note("reference treatment: %s",
         trimws(sub("^reference[^:]*:", "", ref_line[1])))

  if (is.null(PLEX_OF)) {
    note(paste("no plex column in %s, so this report cannot name the batch.",
               "Copy 'plex' across from\n      %s, keyed on sample."),
         meta_path, design_path0)
  } else {
    cat("\nsamples per plex:\n")
    print(table(plex = PLEX_OF[meta$.col]))
    if ("group" %in% names(meta)) {
      cat("\ncondition x plex (a diagonal-only table means the two cannot be",
          "\nseparated; the design diagnostics below stop on it):\n")
      print(table(group = meta$group, plex = PLEX_OF[meta$.col]))
    }
    note(paste("the plex is a COVARIATE here: a TMT batch, and never a condition.",
               "Contrasts are formed over\n      the condition, and metaannot",
               "does not infer the condition from the plex."))
  }
}
```

## The model, before anything is fitted

```{r design}
form <- as.formula(params$design_formula)
vars <- all.vars(form)
missing_vars <- setdiff(vars, names(meta))
if (length(missing_vars))
  stop("design_formula references columns not in the metadata: ",
       paste(missing_vars, collapse = ", "), call. = FALSE)

md <- meta[, vars, drop = FALSE]
if (any(!complete.cases(md)))
  stop(sum(!complete.cases(md)), " sample(s) have NA in a model variable; ",
       "fill them in or drop those samples from the metadata", call. = FALSE)
md <- droplevels(md)

# An isobaric run models the plex as a batch. With one plex there is no batch
# to model, and model.matrix() would fail with "contrasts can be applied only
# to factors with 2 or more levels", which does not say which term it means.
if ("plex" %in% vars && is.factor(md$plex) && nlevels(md$plex) < 2)
  stop("design_formula models 'plex' but every sample is in plex '",
       levels(md$plex)[1], "'. One plex is not a batch effect: drop '+ plex' ",
       "from design_formula.", call. = FALSE)

design <- model.matrix(form, data = md)
rownames(design) <- meta$.col
colnames(design) <- make.names(colnames(design))

cat("coefficient names available to `contrasts`:\n")
print(colnames(design))
cat("\n")
if ("plex" %in% vars)
  note(paste("'plex' is in the model as a TMT BATCH term, not as a",
             "hypothesis: the condition is\n      estimated within plex, and",
             "the plex coefficients are nuisance parameters.\n      Contrasts",
             "are formed over the condition; a plex-versus-plex contrast",
             "would be\n      a batch effect reported as biology."))
as_tibble(design) %>% head(12)
```

### Design diagnostics

A confounded or rank-deficient design fails deep inside `makeContrasts` with
an unhelpful error, or worse, fits and returns nonsense. Checked here instead.

```{r design-diagnostics}
# Checked BEFORE the rank, because rank deficiency is how this shows up and
# "coefficient plexTMT8 is not estimable" names the symptom, not the mistake.
# A TMT plex is a batch: when every plex holds exactly one level of the
# condition, the batch and the biology are the same vector and no model
# separates them.
plex_confounding <- function(md, batch = "plex") {
  if (!batch %in% names(md) || !is.factor(md[[batch]])) return(invisible(NULL))
  b <- droplevels(md[[batch]])
  if (nlevels(b) < 2) return(invisible(NULL))
  fc <- setdiff(names(md)[vapply(md, is.factor, logical(1))], batch)
  for (f in fc) {
    tt <- table(b, droplevels(md[[f]]), dnn = c(batch, f))
    if (nlevels(droplevels(md[[f]])) > 1 && all(rowSums(tt > 0) == 1))
      stop("'", batch, "' is perfectly confounded with '", f, "': each of the ",
           nlevels(b), " ", batch, "es holds exactly one level of '", f,
           "'.\n", paste(capture.output(print(tt)), collapse = "\n"),
           "\n\nThe ", batch, " is a batch, so the batch effect and the ",
           "condition are the same\nvector and no model can separate them: ",
           "every fold change would be both.\nA TMT design needs each ",
           "condition spread over several plexes. Dropping '+ ", batch,
           "'\nfrom design_formula does not fix it — it reports the batch ",
           "effect as biology.", call. = FALSE)
  }
  invisible(NULL)
}
plex_confounding(md)

r <- qr(design)$rank
gate("design: %d samples, %d coefficients, rank %d%s",
     nrow(design), ncol(design), r,
     if (r < ncol(design)) "  <-- RANK DEFICIENT" else "")
if (r < ncol(design)) {
  ne <- limma::nonEstimable(design)
  stop("design is rank deficient; these coefficients are not estimable: ",
       paste(ne, collapse = ", "),
       "\nUsually a covariate is perfectly confounded with the condition ",
       "(every treated\nsample in batch 1, every control in batch 2). No model ",
       "can separate those.",
       if ("plex" %in% vars)
         paste0("\n'plex' is in this model: in a TMT run the plex is that ",
                "batch, and a condition\nthat does not cross plexes cannot be ",
                "adjusted for it.") else "",
       call. = FALSE)
}
if (nrow(design) - r < 1)
  stop("0 residual degrees of freedom: ", nrow(design), " samples and ", r,
       " coefficients leaves nothing to estimate variance from, so no test is\n",
       "  possible. You need replication, or a simpler model.", call. = FALSE)
if (nrow(design) - r < 2)
  gate("only %d residual degree(s) of freedom — the moderated test has very little to work with",
       nrow(design) - r)

fac <- vars[vapply(md, is.factor, logical(1))]
if (length(fac)) {
  cat("\ngroup sizes:\n")
  for (f in fac) print(table(md[[f]], dnn = f))
  # A single-level factor is not rank deficient, so the check above passes and
  # the failure surfaces much later as "no contrasts given".
  one <- fac[vapply(fac, function(f) nlevels(md[[f]]) < 2, logical(1))]
  if (length(one))
    stop("factor(s) with only one level: ", paste(one, collapse = ", "),
         ".\n  Nothing can be contrasted against anything. Check that the ",
         "manifest\n  or metadata really distinguishes your conditions.",
         call. = FALSE)
  small <- lapply(fac, function(f) names(which(table(md[[f]]) < 2)))
  names(small) <- fac
  for (f in fac) if (length(small[[f]]))
    gate("%s has level(s) with a single sample: %s — those get no within-group variance",
         f, paste(small[[f]], collapse = ", "))
}
if (length(fac) > 1) {
  cat("\ncross-tabulation of factors (a zero cell off the diagonal means the\n",
      "two are partly confounded; a diagonal-only table means fully):\n", sep = "")
  for (i in seq_along(fac)) for (j in seq_along(fac)) if (i < j) {
    tt <- table(md[[fac[i]]], md[[fac[j]]], dnn = c(fac[i], fac[j]))
    print(tt)
    if (any(tt == 0)) gate("%s x %s has empty cells — partial confounding",
                           fac[i], fac[j])
  }
}
num <- vars[vapply(md, is.numeric, logical(1))]
if (length(num) && length(fac)) {
  cat("\ncontinuous covariates vs the first factor (large differences between\n",
      "levels mean the covariate carries condition information):\n", sep = "")
  for (n in num) print(tapply(md[[n]], md[[fac[1]]], function(x)
    round(c(mean = mean(x), sd = sd(x)), 3)))
}
```

### Blocking

```{r block}
blk <- NULL
if (nzchar(params$block_col)) {
  if (!params$block_col %in% names(meta))
    stop("block_col names a missing column: ", params$block_col, call. = FALSE)
  if (any(is.na(meta[[params$block_col]])))
    stop("block_col '", params$block_col, "' has NA values; every sample needs ",
         "a block", call. = FALSE)
  # Ordered to the matrix columns, not to metadata row order.
  blk <- droplevels(factor(meta[[params$block_col]]))[match(int_cols, meta$.col)]
  gate("blocking on '%s': %d levels, sizes %s", params$block_col, nlevels(blk),
       paste(table(blk), collapse = "/"))
  if (nlevels(blk) == nrow(meta))
    stop("every sample has a unique ", params$block_col,
         " — that is not a blocking structure", call. = FALSE)
  note(paste("fitted with duplicateCorrelation, not as a fixed effect, so a",
             "block nested\n      inside the condition (cage within diet) does",
             "not consume the contrast"))
}
```

### Contrasts

```{r contrasts}
parse_contrasts <- function(s) {
  parts <- trimws(unlist(strsplit(s, ";")))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) stop("no contrasts given", call. = FALSE)
  named <- grepl("=", parts, fixed = TRUE)
  ex <- ifelse(named, trimws(sub("^[^=]*=", "", parts)), parts)
  nm <- ifelse(named, trimws(sub("=.*$", "", parts)), make.names(parts))
  setNames(ex, nm)
}
CTR <- parse_contrasts(params$contrasts)

used <- unique(unlist(strsplit(CTR, "[-+*/() ]+")))
used <- used[nzchar(used) & is.na(suppressWarnings(as.numeric(used)))]
bad  <- setdiff(used, colnames(design))
if (length(bad))
  stop("contrast term(s) not among the design coefficients: ",
       paste(bad, collapse = ", "),
       "\n  available: ", paste(colnames(design), collapse = ", "),
       "\n  (model.matrix pastes the variable name onto the level, so factor ",
       "`group`\n   with level `fiber` gives the coefficient `groupfiber`)",
       call. = FALSE)

ctr <- limma::makeContrasts(contrasts = unname(CTR), levels = design)
colnames(ctr) <- names(CTR)
PRIMARY <- if (nzchar(params$primary_contrast)) params$primary_contrast else names(CTR)[1]
if (!PRIMARY %in% names(CTR)) stop("primary_contrast not among the contrasts", call. = FALSE)
# PRIMARY always names a limma contrast — it indexes topTable() and the
# contrast coefficient matrix. PRIMARY_LABEL is only ever printed or pasted
# into a filename, and the MSstats branch below overrides that one alone.
PRIMARY_LABEL <- PRIMARY
gate("%d contrast(s); primary = %s  (%s)", length(CTR), PRIMARY, CTR[[PRIMARY]])
ctr
```

# Annotation QC

## Bin composition

```{r bins}
ann <- mutate(ann, bin = factor(bin, levels = BIN_LEVELS))
aq  <- mutate(aq,  bin = factor(col_or_na(aq, "bin"), levels = BIN_LEVELS))

bin_tab <- aq %>% count(bin, .drop = FALSE) %>%
  mutate(pct = round(100 * n / sum(n), 1))
bin_tab
# Groups quantified but absent from annotation_final.tsv (human host,
# contaminants, entrapment sequences) have no bin. Left in the mean, one NA
# turned the headline number into "NA%", so they are counted separately.
n_unbinned <- sum(is.na(aq$bin))
cat(sprintf("\n%.1f%% of annotated quantified protein groups are invisible to KEGG pathway enrichment.\n",
            100 * mean(aq$bin[!is.na(aq$bin)] != "1_ko_pathway")))
if (n_unbinned)
  note(paste("%d of %d quantified group(s) have no annotation row and no bin",
             "(host, contaminant or entrapment sequences?);\n      they are",
             "excluded from the percentage above and appear as the NA row in",
             "the table"), n_unbinned, nrow(aq))

ggplot(bin_tab, aes(bin, n, fill = bin)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)), vjust = -0.15, size = 3) +
  scale_fill_manual(values = BIN_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "protein groups",
       title = "Quantified proteome by annotation evidence") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
```

## Feature support per protein

With peptide- or ion-level input, each protein's quantification rests on a
countable number of features, and metaannot records how the shared ones were
handled. This matters more for the KO-less bins than for anything else: those
proteins are strain-specific, so they carry fewer peptides, and a protein
quantified from one shared peptide is not evidence of anything.

```{r feature-support, fig.height=3.4}
n_used <- col_or_na(aq, "n_features_used")
if (all(is.na(n_used))) {
  note(paste("no peptide_evidence columns in annotated_quant.tsv:",
             "the input was a protein-level table, so feature support is unknown"))
} else {
  aq$n_features_used   <- n_used
  aq$n_unique          <- col_or_na(aq, "n_unique")
  aq$n_features_dropped<- col_or_na(aq, "n_features_dropped")

  print(aq %>% group_by(bin) %>%
          summarise(n = n(),
                    median_features = median(n_features_used, na.rm = TRUE),
                    pct_single_feature = round(100 * mean(n_features_used <= 1, na.rm = TRUE), 1),
                    median_unique = median(n_unique, na.rm = TRUE),
                    features_dropped = sum(n_features_dropped, na.rm = TRUE),
                    .groups = "drop"))

  print(ggplot(aq, aes(bin, n_features_used, fill = bin)) +
          geom_boxplot(outlier.size = 0.4, alpha = 0.85) + scale_y_log10() +
          scale_fill_manual(values = BIN_COLS, guide = "none") +
          labs(x = NULL, y = "assigned features (log10)") +
          theme(axis.text.x = element_text(angle = 25, hjust = 1)))

  if (params$min_features > 0) {
    drop <- aq$n_features_used < params$min_features
    drop[is.na(drop)] <- FALSE
    gate("dropping %d protein(s) with fewer than %d assigned features",
         sum(drop), params$min_features)
    aq <- aq[!drop, , drop = FALSE]
    # bin_tab was computed before this filter; recompute so the composition
    # table describes the same protein set as everything below it.
    bin_tab <- aq %>% count(bin, .drop = FALSE) %>%
      mutate(pct = round(100 * n / sum(n), 1))
    cat("\nbin composition after the feature-support filter:\n"); print(bin_tab)
  }
  ev_path <- file.path(RD, "quant", "peptide_evidence.tsv")
  if (file.exists(ev_path))
    note(paste("assignment rule and per-protein counts: %s.",
               "A protein whose features were mostly dropped as cross-taxon",
               "shared\n      has no defensible quantification, whatever its",
               "annotation says."), ev_path)
}
```

## Where each KO came from

eggNOG assigns KOs by DIAMOND search against seed orthologs; KOfamScan uses
per-family HMMs with adaptive thresholds. Running both is a control, not just
extra coverage: if KOfam rescues a large slice of the KO-less bins, those bins
were partly an artefact of eggNOG's search rather than a fact about the
proteins.

```{r ko-provenance}
src <- col_or_na(aq, "ko_source")
if (all(is.na(src)) || all(src == "", na.rm = TRUE)) {
  note("no ko_source column: KOfamScan was not run, so eggNOG's KO calls are unchecked")
} else {
  aq$ko_source <- factor(ifelse(is.na(src) | src == "", "none", src),
                         levels = c("both", "eggnog", "kofam", "none"))
  print(count(aq, ko_source, .drop = FALSE) %>%
          mutate(pct = round(100 * n / sum(n), 1)))
  n_resc <- sum(aq$ko_source == "kofam")
  gate("%d protein(s) carry a KO from KOfamScan that eggNOG missed (%.1f%% of all)",
       n_resc, 100 * mean(aq$ko_source == "kofam"))
  if (n_resc > 0)
    note(paste("those proteins would have been counted as KO-less on eggNOG",
               "alone.\n      Quote this number alongside any claim about the",
               "size of the non-KEGG fraction."))
  print(ggplot(aq, aes(bin, fill = ko_source)) + geom_bar(position = "fill") +
          scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
          labs(x = NULL, y = NULL, fill = NULL, title = "KO provenance by bin") +
          theme(axis.text.x = element_text(angle = 25, hjust = 1)))
}
```

## Do the sources agree, or merely overlap?

Two sources reaching the same protein is not the same as the two of them
agreeing about it, and only the second tells you the evidence is corroborated.
Each row below compares a pair that shares an identifier namespace, so the
comparison is between the calls themselves rather than between coverage counts.

Read `disjoint` first: both sides called something and they share nothing. A
few of those are ordinary — a paralogue boundary, a threshold near the edge of
a family. A rate near 100% is almost never real disagreement; it means the two
columns hold different kinds of identifier and nothing is being compared.

```{r source-agreement}
ag_path <- file.path(RD, "source_agreement.tsv")
agree <- if (file.exists(ag_path)) read_tsv_full(ag_path) else NULL
if (is.null(agree) || nrow(agree) == 0) {
  note(paste("no two sources in this run share an identifier namespace,",
             "so none of them can be cross-checked against another"))
} else {
  print(agree[, c("comparison", "both", "identical", "overlapping",
                  "disjoint", "pct_agree")])
  for (i in seq_len(nrow(agree))) {
    r <- agree[i, ]
    if (as.integer(r$both) >= 50 && r$pct_disjoint >= 90) {
      gate(paste("%s: %.1f%% disjoint. That is not a credible rate of real",
                 "disagreement between two views of %s -- check that %s and %s",
                 "carry the same kind of identifier"),
           r$comparison, r$pct_disjoint, r$tests, r$a, r$b)
    } else if (r$pct_disjoint >= 20) {
      gate("%s: %.1f%% of shared calls disagree outright; read those before trusting either source alone",
           r$comparison, r$pct_disjoint)
    } else {
      # note() is itself a sprintf wrapper, so pre-formatting here would
      # format twice and the % in "97.4% agree" would have no argument left.
      note("%s: %.1f%% agree over %d shared protein(s)",
           r$comparison, r$pct_agree, as.integer(r$both))
    }
  }
  one_sided <- agree[agree$overlapping > 0 &
                     (agree$a_superset == agree$overlapping |
                      agree$b_superset == agree$overlapping), , drop = FALSE]
  if (nrow(one_sided))
    note(paste("where the two differ, one side is consistently the more",
               "sensitive rather than the two contradicting each other:",
               paste(one_sided$comparison, collapse = "; ")))
}
```

## Recurrent unknown folds

Fifty dark proteins sharing a fold is a far stronger signal than fifty
singletons, and it needs no reference database. Non-singleton clusters are the
ones worth a structure figure.

```{r fold-clusters}
fcs <- col_or_na(aq, "fold_cluster_size")
if (all(is.na(fcs)) || sum(fcs, na.rm = TRUE) == 0) {
  note("no fold_cluster columns: Foldseek self-clustering was not run")
} else {
  aq$fold_cluster <- col_or_na(aq, "fold_cluster")
  aq$fold_cluster_size <- fcs
  big <- aq %>% filter(!is.na(fold_cluster), fold_cluster != "",
                       fold_cluster_size > 1) %>%
    group_by(fold_cluster) %>%
    summarise(size = dplyr::first(fold_cluster_size),
              bins = paste(sort(unique(as.character(bin))), collapse = ", "),
              members = paste(head(group_id, 5), collapse = ", "),
              .groups = "drop") %>% arrange(desc(size))
  gate("%d structural cluster(s) with more than one member", nrow(big))
  print(head(big, 20))
}
```

## Is the KO-less fraction identified as well?

```{r id-rate}
identified <- unique(unlist(strsplit(aq$group_id, ";")))
if (nrow(ann) > nrow(aq) * 1.2) {
  ann %>%
    mutate(found = protein_id %in% identified) %>%
    group_by(bin) %>%
    summarise(in_database = n(), identified = sum(found),
              pct = round(100 * mean(found), 2), .groups = "drop") %>% print()
  tt <- with(mutate(ann, no_ko = !(col_or_na(ann, "has_ko") %in% TRUE),
                    found = protein_id %in% identified),
             table(no_ko, found))
  if (all(dim(tt) == 2)) {
    ft <- fisher.test(tt)
    gate("identification rate, KO-less vs KO-bearing: OR = %.2f (95%% CI %.2f-%.2f), p = %.3g",
         ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)
    note(paste("an OR below 1 means KO-less proteins are identified less often,",
               "which biases every per-bin count below"))
  }
} else {
  note(paste("annotation_final.tsv has %d rows against %d quantified groups:",
             "the annotated fasta looks like the identified subset already,",
             "so the identification-rate comparison is unavailable"),
       nrow(ann), nrow(aq))
}
```

```{r id-quality, fig.height=3.6}
M <- as.matrix(aq[, int_cols]); rownames(M) <- aq$group_id
n_samples <- length(int_cols)
qc <- tibble(group_id = aq$group_id, bin = aq$bin,
             aa_length = col_or_na(aq, "length"),
             n_valid = rowSums(is.finite(M) & M > 0)) %>%
  mutate(frac_valid = n_valid / n_samples)

p1 <- ggplot(qc, aes(bin, aa_length, fill = bin)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) + scale_y_log10() +
  scale_fill_manual(values = BIN_COLS, guide = "none") +
  labs(x = NULL, y = "length (aa, log10)") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
p2 <- ggplot(qc, aes(bin, frac_valid, fill = bin)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  scale_fill_manual(values = BIN_COLS, guide = "none") +
  labs(x = NULL, y = "fraction of samples with signal") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
if (requireNamespace("patchwork", quietly = TRUE))
  patchwork::wrap_plots(p1, p2, nrow = 1) else { print(p1); print(p2) }
```

## Taxonomy: eggNOG seed ortholog vs Unipept peptide LCA

Two independent estimates with different failure modes. eggNOG's `seed_taxid`
is the taxon of the best-matching *reference* protein, which for a divergent
gut organism can be a different genus entirely. Unipept's LCA is computed from
peptides against UniProt, is deliberately conservative, and misses peptides
its tryptic index does not contain.

This is not a side quest: the `taxon_unique` peptide assignment and the ratio
model below both stand on whichever taxonomy metaannot was told to use.

```{r taxonomy-comparison, fig.height=3.8}
verdict <- col_or_na(aq, "taxonomy_verdict")
if (all(is.na(verdict))) {
  note(paste("no taxonomy comparison in the quant table: run metaannot's",
             "unipept and taxonomy stages to produce it"))
} else {
  # "differ (no lineage)" is what the taxonomy stage emits when no NCBI
  # taxonomy dump was configured, which is the default. Missing from the level
  # list it became NA, and the disagreement gate below then reported 0%.
  aq$taxonomy_verdict <- factor(
    verdict, levels = c("identical", "concordant", "concordant_above_genus",
                        "conflict", "differ (no lineage)",
                        "unipept_missing", "eggnog_missing"))
  print(count(aq, taxonomy_verdict, .drop = FALSE) %>%
          mutate(pct = round(100 * n / sum(n), 1)))

  print(ggplot(aq, aes(bin, fill = taxonomy_verdict)) +
          geom_bar(position = "fill") +
          scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
          labs(x = NULL, y = NULL, fill = NULL,
               title = "eggNOG / Unipept agreement by annotation bin") +
          theme(axis.text.x = element_text(angle = 25, hjust = 1)))

  disc <- aq$taxonomy_verdict %in% c("concordant_above_genus", "conflict",
                                     "differ (no lineage)")
  gate("%d/%d proteins (%.1f%%) disagree below family level",
       sum(disc), nrow(aq), 100 * mean(disc))
  n_nolin <- sum(aq$taxonomy_verdict %in% "differ (no lineage)")
  if (n_nolin)
    note(paste("%d of those are 'differ (no lineage)': db.ncbi_taxonomy was",
               "not set, so the two taxids could only be compared for exact",
               "equality.\n      The rank of the disagreement is unknown —",
               "point db.ncbi_taxonomy at a taxdump to resolve it"), n_nolin)
  note(paste("if the KO-less bins show more disagreement than 1_ko_pathway,",
             "that is expected — eggNOG has no close reference for them —\n",
             "     and it is precisely why the ratio model should be read with",
             "the verdict column beside it"))

  dr <- col_or_na(aq, "deepest_agreement")
  if (!all(is.na(dr))) {
    cat("\ndeepest rank of agreement:\n")
    print(sort(table(dr), decreasing = TRUE))
  }
  agr <- col_or_na(aq, "unipept_agreement")
  if (!all(is.na(agr)))
    note(paste("unipept_agreement is the fraction of a protein's peptides",
               "supporting its consensus taxon;\n      low values usually mean",
               "short fragments from missed-cleavage splitting"))

  if (isTRUE(params$require_taxonomy_concordance)) {
    ok <- aq$taxonomy_verdict %in% c("identical", "concordant")
    gate("require_taxonomy_concordance: blanking the taxon of %d discordant protein(s); they drop out of the ratio model",
         sum(!ok))
    if ("effective_taxid" %in% names(aq)) aq$effective_taxid[!ok] <- NA
  }
}
```

## QC gates

```{r gates}
rp <- file.path(RD, "eggnog", "reuse_report.tsv")
if (file.exists(rp)) {
  rr <- read_tsv(rp, show_col_types = FALSE)
  cov <- suppressWarnings(as.numeric(rr$value[rr$metric == "coverage"]))
  rec <- rr$value[rr$metric == "recommended_transform"]
  if (length(cov) && !is.na(cov)) {
    if (cov < 0.90)
      gate("emapper covers only %.1f%% of proteins%s", 100 * cov,
           if (length(rec) && !is.na(rec) && nzchar(rec))
             paste0("; metaannot suggests emapper_id_transform: ", rec) else "")
    else gate("emapper coverage %.1f%%", 100 * cov)
  }
}
gp <- file.path(RD, "quant", "group_conflicts.tsv")
if (file.exists(gp)) {
  gcf <- read_tsv(gp, show_col_types = FALSE)
  # A zero-row file still has a header, and readr types empty columns as
  # character, so the logical operators below need explicit coercion.
  lg <- function(x) if (is.null(x)) logical(0) else as.logical(x)
  gate("%d groups with inconsistent member annotation, %d wholly unannotated",
       sum(lg(gcf$bin_conflict) | lg(gcf$ko_conflict), na.rm = TRUE),
       sum(lg(gcf$unannotated_group), na.rm = TRUE))
}
bc <- col_or_na(aq, "bin_conflict")
if (!all(is.na(bc))) {
  cat("\nconflict rate within bin:\n")
  print(aq %>% mutate(conflict = bc) %>% group_by(bin) %>%
          summarise(n = n(), pct_conflicting = round(100 * mean(conflict, na.rm = TRUE), 1),
                    .groups = "drop"))
}
```

A group whose members disagree carries the leading protein's annotation by
convention, not by evidence. Those rows stay in the analysis and are flagged
in every output table; drop them if a specific claim depends on the
annotation being right.

# Quantification

```{r prep}
X <- M
X[!is.finite(X) | X <= 0] <- NA
X <- log2(X)
# The per-sample log2 shift this normalisation applied, kept because the taxon
# diagnostics further down have to remove the SAME shift from the size factors.
# Fitted on raw factors they carried the loading difference between groups,
# which made every taxon look as if it had moved.
NORM_OFFSET <- setNames(rep(0, ncol(X)), colnames(X))
X <- switch(params$normalise,
            median   = {
              med <- apply(X, 2, median, na.rm = TRUE)
              NORM_OFFSET <- med - mean(med)
              sweep(X, 2, NORM_OFFSET)
            },
            quantile = {
              Xq <- limma::normalizeBetweenArrays(X, method = "quantile")
              # Quantile normalisation is not a shift; the per-sample median
              # difference is the closest constant to it, so the consistency
              # check below is approximate under this setting.
              NORM_OFFSET <- apply(X - Xq, 2, median, na.rm = TRUE)
              Xq
            },
            none     = X,
            stop("normalise must be median, quantile or none"))

fgrp <- meta[[params$group_col_for_filtering]]
if (is.null(fgrp))
  stop("group_col_for_filtering names a missing column: ",
       params$group_col_for_filtering, call. = FALSE)
fgrp <- droplevels(factor(fgrp))[match(colnames(X), meta$.col)]

# A level smaller than min_valid_per_group can never satisfy the filter, so
# every protein is dropped and the failure surfaces two sections later as a
# cryptic limma error on a 0-row matrix.
smallest <- table(fgrp)
if (params$min_valid_per_group > min(smallest))
  stop("min_valid_per_group is ", params$min_valid_per_group, " but level '",
       names(smallest)[which.min(smallest)], "' of ",
       params$group_col_for_filtering, " has only ", min(smallest),
       " sample(s).\n  No protein can pass the filter. Lower ",
       "min_valid_per_group, or filter on a\n  coarser column.", call. = FALSE)

# isTRUE(): tapply over a factor with an unused level returns NA, all(NA) is
# NA, and X[NA, ] silently produces a row of NAs rather than dropping it.
keep_valid <- apply(X, 1, function(r)
  isTRUE(all(tapply(is.finite(r), fgrp, sum) >= params$min_valid_per_group)))
keep_valid[is.na(keep_valid)] <- FALSE

# min_valid_per_group counts SAMPLES, and that is not enough for an isobaric
# run. Missingness there is structured by plex: a protein identified in one
# plex only is all-NA in every other, so "3 valid values in every group" can
# be satisfied entirely inside one batch, and the difference the model then
# reports is that batch. min_plexes counts PLEXES instead. The two ask
# different questions, so both are applied and each is reported on its own.
pbatch <- if (IS_ISOBARIC && !is.null(PLEX_OF)) PLEX_OF[colnames(X)] else NULL
if (params$min_plexes > 1 && is.null(pbatch))
  stop("min_plexes is ", params$min_plexes, " but no per-sample plex is ",
       "available: ", if (!IS_ISOBARIC)
         "this run is not isobaric (no quant/design_notes.txt), so there are no plexes"
       else paste0("neither the metadata nor ", design_path0, " gives every ",
                   "sample a plex"),
       ".\n  Set analysis.min_plexes: 1, or add a 'plex' column keyed on ",
       "sample.", call. = FALSE)
n_plex <- if (is.null(pbatch)) rep(NA_integer_, nrow(X)) else
  apply(X, 1, function(r) length(unique(pbatch[is.finite(r)])))
keep_plex <- if (is.null(pbatch)) rep(TRUE, nrow(X)) else
  n_plex >= params$min_plexes

# Separately, so it is visible which filter bit. Counted against the same
# starting set rather than in sequence: "min_plexes removed 40" has to mean
# 40 proteins, not "40 of whatever min_valid_per_group left".
cat(sprintf("min_valid_per_group >= %d in every level of %s: removes %d of %d\n",
            params$min_valid_per_group, params$group_col_for_filtering,
            sum(!keep_valid), nrow(X)))
if (is.null(pbatch)) {
  cat(sprintf("min_plexes: not applied (no per-sample plex%s)\n",
              if (IS_ISOBARIC) "" else "; this is not an isobaric run"))
} else {
  cat(sprintf("min_plexes >= %d: removes %d of %d, %d of which min_valid_per_group would have kept\n",
              params$min_plexes, sum(!keep_plex), nrow(X),
              sum(!keep_plex & keep_valid)))
  cat("\nproteins by number of plexes they are quantified in:\n")
  print(table(plexes = n_plex))
  # The number that says whether the filter is worth setting, printed whether
  # or not it is set: with min_plexes at 1 these proteins passed on a sample
  # count that one batch supplied on its own.
  confined <- sum(n_plex <= 1 & keep_valid & keep_plex)
  if (params$min_plexes < 2 && length(unique(pbatch)) > 1 && confined > 0)
    gate(paste("%d protein(s) pass min_valid_per_group but are quantified in a",
               "single plex — their group difference is inside one batch.",
               "Set analysis.min_plexes: 2 to drop them"), confined)
}
keep <- keep_valid & keep_plex
cat(sprintf("%d/%d groups retained by both filters\n", sum(keep), nrow(X)))
X   <- X[keep, , drop = FALSE]
aqk <- aq[keep, , drop = FALSE]

# A protein identical in every replicate has zero within-group variance. limma
# moderates it rather than failing, which turns a constant into an enormous
# t-statistic — usually a single shared peptide or an imputed value, not
# biology.
gv <- apply(X, 1, function(r) {
  v <- tapply(r, fgrp, function(z) stats::var(z, na.rm = TRUE))
  # A group with a single finite value has an UNDEFINED variance, not a zero
  # one. Counting NA as zero flagged sparse proteins as constants whenever
  # min_valid_per_group was 1.
  any(!is.na(v)) && all(v[!is.na(v)] == 0)
})
gv[is.na(gv)] <- FALSE
if (any(gv)) {
  gate("%d protein(s) have zero variance within every group", sum(gv))
  note(paste("they will carry an inflated moderated t. Usually a single shared",
             "peptide or an imputed constant.\n      Listed in",
             "zero_variance.tsv; consider dropping them."))
  write_tsv(tibble(group_id = rownames(X)[gv]),
            file.path(OUT, "zero_variance.tsv"))
  if (isTRUE(params$drop_zero_variance)) {
    X <- X[!gv, , drop = FALSE]; aqk <- aqk[!gv, , drop = FALSE]
    gate("dropped them (drop_zero_variance: true)")
  }
}

cat("\nretention by bin:\n")
print(tibble(bin = aq$bin, kept = keep) %>% group_by(bin) %>%
        summarise(n = n(), n_kept = sum(kept),
                  pct = round(100 * mean(kept), 1), .groups = "drop"))
```

No imputation. Missing values in DIA data are largely not missing at random,
and imputing pulls low-abundance proteins — which is where the KO-less
fraction lives — toward whatever the imputation assumes. `limma` fits each
row on the values it has. If MNAR modelling matters for a specific claim,
refit that subset with `proDA` or `msqrob2` rather than imputing here.

# Differential abundance

Two models. The second exists because a KO-less protein is far more likely to
be tracking the abundance of its source organism than a central-metabolism
enzyme is, and a protein riding its organism's abundance is not a regulatory
finding.

```{r fit-fn}
fit_de <- function(mat, label, amean = NULL) {
  cor_est <- NULL
  if (!is.null(blk)) {
    dc <- limma::duplicateCorrelation(mat, design, block = blk)
    cor_est <- dc$consensus.correlation
    gate("%s: consensus within-%s correlation = %.3f", label,
         params$block_col, cor_est)
    if (!is.na(cor_est) && cor_est < 0)
      note("negative consensus correlation; blocking is not helping here")
    f <- limma::lmFit(mat, design, block = blk, correlation = cor_est)
  } else {
    f <- limma::lmFit(mat, design)
  }
  f <- limma::contrasts.fit(f, ctr)
  # trend = TRUE models variance as a function of ABUNDANCE. For the ratio
  # model the response is a log-ratio centred near zero, so Amean carries no
  # abundance information at all and the trend is fitted against noise; the
  # caller passes the mean raw log2 intensity instead.
  if (!is.null(amean)) f$Amean <- unname(amean[rownames(f$coefficients)])
  f <- limma::eBayes(f, trend = TRUE, robust = TRUE)
  list(fit = f, correlation = cor_est)
}
tt_of <- function(f, coef) {
  limma::topTable(f, coef = coef, number = Inf, sort.by = "none") %>%
    rownames_to_column("group_id") %>% as_tibble()
}
sig_of <- function(d) !is.na(d$adj.P.Val) & d$adj.P.Val < params$fdr &
  abs(d$logFC) > params$min_lfc
```

## Model 1 — protein abundance

```{r naive}
USE_MS <- nzchar(params$msstats_comparison)
F1 <- NULL

if (USE_MS) {
  if (!file.exists(params$msstats_comparison))
    stop("msstats_comparison file not found: ", params$msstats_comparison,
         call. = FALSE)
  cr <- if (grepl("\\.csv$", params$msstats_comparison, ignore.case = TRUE))
    read_csv(params$msstats_comparison, show_col_types = FALSE, guess_max = 10000)
  else read_tsv(params$msstats_comparison, show_col_types = FALSE, guess_max = 10000)
  needed <- c("Protein", "Label", "log2FC", "pvalue", "adj.pvalue")
  miss <- setdiff(needed, names(cr))
  if (length(miss))
    stop("this does not look like MSstats groupComparison()$ComparisonResult; ",
         "missing: ", paste(miss, collapse = ", "),
         "\n  export it with write.table(cmp$ComparisonResult, sep = '\\t')",
         call. = FALSE)
  lab <- if (nzchar(params$msstats_label)) params$msstats_label else cr$Label[1]
  if (!lab %in% cr$Label)
    stop("msstats_label '", lab, "' not in the file; labels present: ",
         paste(unique(cr$Label), collapse = ", "), call. = FALSE)
  cr <- filter(cr, Label == lab)
  gate("using MSstats results for label '%s' (%d proteins)", lab, nrow(cr))

  if ("issue" %in% names(cr)) {
    iss <- cr %>% filter(!is.na(issue) & issue != "") %>% count(issue)
    if (nrow(iss)) {
      cat("\nMSstats flagged these fitting issues:\n"); print(iss)
      note(paste("oneConditionMissing / completeMissing rows carry an infinite",
                 "or NA log2FC.\n      They are kept and will simply fail the",
                 "significance filter."))
    }
  }
  res_naive <- cr %>%
    transmute(group_id = as.character(Protein), logFC = log2FC,
              P.Value = pvalue, adj.P.Val = adj.pvalue) %>%
    filter(group_id %in% aqk$group_id)
  gate("%d/%d MSstats proteins matched the annotated quant table",
       nrow(res_naive), nrow(cr))
  # Label the outputs with MSstats' own comparison name, but leave PRIMARY a
  # limma contrast name: the ratio model below still fits contrasts here, and
  # topTable(coef = "<MSstats label>") is a subscript-out-of-bounds error.
  PRIMARY_LABEL <- lab
  note(paste("the abundance model comes from MSstats label '%s' while the",
             "ratio model is fitted here for contrast '%s'.\n      Those two",
             "must describe the same comparison, or the verdict column is",
             "meaningless — check primary_contrast"), lab, PRIMARY)
} else {
  F1 <- fit_de(X, "abundance model")
  res_naive <- tt_of(F1$fit, PRIMARY)
  print(map_dfr(names(CTR), function(cn) {
    d <- tt_of(F1$fit, cn)
    tibble(contrast = cn, expression = CTR[[cn]], tested = nrow(d),
           significant = sum(sig_of(d)),
           up = sum(sig_of(d) & d$logFC > 0), down = sum(sig_of(d) & d$logFC < 0))
  }))
}
```

When MSstats supplies the abundance model, its normalisation, summarisation
and any imputation from `dataProcess()` are inherited wholesale. Check
`MBimpute` and `censoredInt` there: imputed low-abundance values are exactly
where the KO-less fraction sits, and imputing them manufactures the very
effects this document then tests. The ratio model below is still fitted here,
because MSstats has no notion of a protein's source organism.

## Model 2 — abundance relative to the source organism

Response becomes `log2(protein) − log2(that organism's size factor in the
same sample)`, where the size factor is a **median of ratios** across the
taxon's proteins, not a sum.

That distinction is not cosmetic. A summed reference is inflated by any
strongly changing member of the taxon, and every other member is then divided
by that inflated total, acquiring a bias in the opposite direction. On a
constructed test where one protein in eight was genuinely up 4× inside a taxon
that doubled, the summed reference read +1.46 instead of +1.00 and gave the
seven passengers a spurious −0.5; the median of ratios read +1.02 and left
them at −0.01. Taxa with fewer than `taxon_min_proteins_for_factor` proteins
fall back to the sum, and those groups are flagged, because a median over
three proteins is not robust either.

Any per-sample scaling cancels in the ratio, so this model uses raw
intensities and needs no normalisation of its own.

```{r taxon-adjusted}
sf_path <- file.path(RD, "quant", "taxon_size_factors.tsv")
tx_path <- file.path(RD, "quant", "taxon_intensity.tsv")
res_adj <- NULL

tax_of_all <- col_or_na(aq, "effective_taxid")
if (all(is.na(tax_of_all))) tax_of_all <- col_or_na(aq, "seed_taxid")
tax_of <- col_or_na(aqk, "effective_taxid")
if (all(is.na(tax_of)) || all(tax_of == "", na.rm = TRUE)) {
  tax_of <- col_or_na(aqk, "seed_taxid")
  note("no effective_taxid column; falling back to the eggNOG seed taxid")
}

if (file.exists(sf_path)) {
  # Typed as character on the way in: read numeric, a round taxid comes back
  # from as.character() as "1e+06" and never matches the quant table again.
  sf <- read_tsv_full(sf_path, c("effective_taxid", "seed_taxid"))
  id_col <- if ("effective_taxid" %in% names(sf)) "effective_taxid" else "seed_taxid"
  sf[[id_col]] <- as.character(sf[[id_col]])
  miss <- setdiff(int_cols, names(sf))
  if (length(miss)) {
    note("taxon_size_factors.tsv lacks %d sample column(s); skipping the ratio model",
         length(miss))
  } else {
    n_per_tax <- tibble(taxid = tax_of_all) %>%
      filter(!is.na(taxid), taxid != "") %>% count(taxid, name = "n")
    ok_tax <- n_per_tax$taxid[n_per_tax$n >= params$taxon_min_proteins]
    usable <- !is.na(tax_of) & tax_of != "" & tax_of %in% ok_tax &
      tax_of %in% sf[[id_col]]
    gate("%d/%d groups sit in a taxon with >= %d proteins and a size factor",
         sum(usable), nrow(aqk), params$taxon_min_proteins)
    fb <- sf$method[match(tax_of[usable], sf[[id_col]])] == "sum_fallback"
    if (any(fb, na.rm = TRUE))
      note(paste("%d group(s) belong to a taxon too small for a robust median,",
                 "so their reference is the plain sum and inherits its bias"),
           sum(fb, na.rm = TRUE))

    if (sum(usable) >= 20) {
      raw <- as.matrix(aqk[, int_cols]); rownames(raw) <- aqk$group_id
      raw[!is.finite(raw) | raw <= 0] <- NA
      fmat <- as.matrix(sf[match(tax_of[usable], sf[[id_col]]), int_cols])
      fmat[!is.finite(fmat) | fmat <= 0] <- NA
      R <- log2(raw[usable, , drop = FALSE]) - log2(fmat)
      rownames(R) <- aqk$group_id[usable]
      keep2 <- apply(R, 1, function(r)
        isTRUE(all(tapply(is.finite(r), fgrp, sum) >= params$min_valid_per_group)))
      keep2[is.na(keep2)] <- FALSE
      R <- R[keep2, , drop = FALSE]
      cat(sprintf("%d groups testable against the taxon size factor\n", nrow(R)))
      # The abundance covariate for eBayes(trend = TRUE) is the protein's mean
      # raw intensity, not the mean of its log-ratios.
      amean_raw <- rowMeans(log2(raw), na.rm = TRUE)
      if (nrow(R) >= 20)
        res_adj <- tt_of(fit_de(R, "ratio model", amean = amean_raw)$fit, PRIMARY)
    } else {
      note(paste("only %d group(s) have a usable taxon; at least 20 are needed",
                 "for the ratio model, so it is skipped"), sum(usable))
    }
  }
} else if (file.exists(tx_path)) {
  note(paste("only taxon_intensity.tsv is present. A summed reference is biased",
             "by any strongly changing member of the taxon;\n      rerun the",
             "join stage to produce taxon_size_factors.tsv."))
} else {
  note("no taxon reference available; skipping the ratio model")
}
adj_attempted <- !is.null(res_adj)
```

## Which findings survive the adjustment

```{r compare, fig.height=4.2}
de <- res_naive %>%
  transmute(group_id, logFC_naive = logFC, FDR_naive = adj.P.Val,
            sig_naive = sig_of(res_naive))

if (adj_attempted) {
  de <- de %>%
    left_join(res_adj %>% transmute(group_id, logFC_adj = logFC,
                                    FDR_adj = adj.P.Val, sig_adj = sig_of(res_adj)),
              by = "group_id") %>%
    mutate(verdict = case_when(
      is.na(sig_adj)        ~ "not testable (no usable taxon)",
      sig_naive & sig_adj   ~ "taxon-independent",
      # Not "tracks source organism": losing significance is also what a
      # noisier response does. The ratio adds the sampling noise of a median
      # over the taxon's proteins, and in a taxon of five the protein is a
      # fifth of its own reference, so a genuinely regulated protein in an
      # unchanged organism drops out for want of power alone.
      sig_naive & !sig_adj  ~ "attenuated by adjustment",
      !sig_naive & sig_adj  ~ "revealed by adjustment",
      TRUE                  ~ "not significant"))
  print(count(de, verdict))
  note(paste("'attenuated by adjustment' means the effect did not survive the",
             "ratio model. Organism tracking is the usual reason, but so is",
             "the\n      noise of a small reference: with taxon_min_proteins =",
             "%d a taxon can have as few as %d proteins behind its size",
             "factor.\n      Read it next to logFC_naive and logFC_adj, not on",
             "its own."), params$taxon_min_proteins, params$taxon_min_proteins)

  ggplot(filter(de, !is.na(logFC_adj)), aes(logFC_naive, logFC_adj, colour = verdict)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(alpha = 0.6, size = 1.1) +
    labs(x = "log2FC, protein abundance", y = "log2FC, relative to organism",
         colour = NULL, title = sprintf("Taxon adjustment: %s", PRIMARY_LABEL))
} else {
  de <- mutate(de, logFC_adj = NA_real_, FDR_adj = NA_real_,
               sig_adj = NA, verdict = "adjusted model not run")
}
```

## Is the normalisation assumption safe?

The size factor assumes most of a taxon's proteins are unchanged, so the
median tracks the organism rather than the biology. If a treatment genuinely
shifts the majority of one organism's proteome, the median follows the
majority, the shift is absorbed into the "abundance" term, and the minority
that did *not* change is reported as moving the other way.

One case is invisible on principle: if *every* protein in a taxon shifts
together, nothing is left unchanged to measure the shift against, and the size
factor absorbs all of it. The diagnostic below catches partial majority shifts,
which is the realistic case; a total one it cannot.

**This is not a solvable estimation problem.** From relative abundances alone,
"the organism doubled and its proteins held steady" and "the organism held
steady and its proteins doubled" are the same numbers. No within-taxon
normalisation can separate them; an orthogonal measurement of organism
abundance — metagenomic coverage, qPCR, 16S — is what settles it. What this
section can do is show you which taxa the question actually bites for, instead
of assuming it away.

```{r taxon-assumption, fig.height=4.6}
risk_tab <- NULL
if (!adj_attempted || !exists("sf") || is.null(res_adj)) {
  note("ratio model was not fitted, so there is no normalisation to diagnose")
} else {
  dev_lfc <- if (nzchar(as.character(params$deviation_lfc)))
    as.numeric(params$deviation_lfc) else params$min_lfc

  # The taxon factor's own effect, estimated under the SAME design and
  # contrast as the proteins, so the two are directly comparable and the
  # diagnostic cannot drift from the model it is checking.
  sfm <- as.matrix(sf[, int_cols])
  sfm[!is.finite(sfm) | sfm <= 0] <- NA
  rownames(sfm) <- sf[[id_col]]
  # The SAME per-sample constants that were removed from the protein matrix
  # are removed here. On raw factors this term carried the loading difference
  # between the groups: every taxon then showed a factor_log2FC equal to that
  # difference, the AT RISK rule's |factor| > min_lfc was satisfied by all of
  # them, and the consistency gate below fired on every run with unequal
  # loading. What is wanted is the organism's share of the community, which is
  # what the factors say once the loading is out of them.
  sfm_n <- log2(sfm) - matrix(NORM_OFFSET[int_cols], nrow(sfm), ncol(sfm),
                              byrow = TRUE)
  f_coef <- limma::contrasts.fit(limma::lmFit(sfm_n, design), ctr)$coefficients
  tax_fc <- setNames(as.numeric(f_coef[, PRIMARY]), rownames(sfm))

  memb <- tibble(group_id = aqk$group_id, taxid = as.character(tax_of)) %>%
    filter(!is.na(taxid), taxid != "") %>%
    left_join(res_naive %>% transmute(group_id, naive = logFC), by = "group_id") %>%
    left_join(res_adj %>% transmute(group_id, ratio = logFC,
                                    ratio_fdr = adj.P.Val), by = "group_id") %>%
    filter(!is.na(ratio))

  # Consistency: the ratio effect should equal the naive effect minus the
  # taxon factor's effect. A large discrepancy means the two models did not
  # see the same samples, and every number below would be suspect.
  chk <- memb %>% mutate(expected = naive - tax_fc[taxid]) %>%
    filter(is.finite(expected), is.finite(ratio))
  if (nrow(chk) > 2) {
    dd <- max(abs(chk$ratio - chk$expected), na.rm = TRUE)
    gate("ratio == naive - taxon factor, max deviation %.3g%s", dd,
         if (dd > 0.05) "   <-- the two models disagree, investigate" else "")
    if (dd > 0.05 && identical(params$normalise, "quantile"))
      note(paste("quantile normalisation is not a per-sample shift, so this",
                 "identity is only approximate under it;\n      a deviation of",
                 "this size need not mean the models saw different samples"))
  }

  # The signal is NOT "many proteins deviate". When the majority of a taxon
  # shifts, the median follows it and only the unchanged MINORITY deviates —
  # in the direction opposite to the factor. That opposition, against a factor
  # that actually moved, is the diagnostic.
  memb <- memb %>%
    mutate(fac = unname(tax_fc[taxid]),
           deviates = abs(ratio) > dev_lfc &
             !is.na(ratio_fdr) & ratio_fdr < params$fdr)
  risk_tab <- memb %>%
    group_by(taxid) %>%
    summarise(n_tested = n(),
              factor_log2FC = dplyr::first(fac),
              n_dev = sum(deviates, na.rm = TRUE),
              n_opposing = sum(deviates &
                                 sign(ratio) != sign(dplyr::first(fac)),
                               na.rm = TRUE),
              .groups = "drop") %>%
    mutate(frac_deviating = n_dev / n_tested,
           opposing = ifelse(n_dev > 0, n_opposing / n_dev, 0),
           # Opposition as a share of everything TESTED, not only of the few
           # that deviated: one opposing protein out of one deviating gives
           # opposing = 1 and used to flag the whole taxon.
           opposing_frac = n_opposing / n_tested,
           n_proteins = sf$n_proteins[match(taxid, sf[[id_col]])],
           method = sf$method[match(taxid, sf[[id_col]])],
           assumption = case_when(
             method == "sum_fallback" ~ "fragile (too few proteins)",
             abs(factor_log2FC) > params$min_lfc &
               n_dev >= params$assumption_min_deviating &
               opposing > params$assumption_opposing &
               opposing_frac >= params$assumption_opposing_frac ~ "AT RISK",
             frac_deviating > params$assumption_frac_deviating ~ "check (median in a sparse region)",
             TRUE ~ "ok"),
           # If the organism did not really change, add the factor back: that
           # is exactly the abundance model, so the two models bracket the
           # ambiguity rather than one of them being right.
           alt_interpretation = ifelse(assumption == "AT RISK",
                                       "use the abundance model", "")) %>%
    # Strong cases first: a taxon flagged on a single opposing protein is much
    # weaker evidence than one flagged on many.
    arrange(desc(abs(factor_log2FC) * opposing * n_opposing))

  print(count(risk_tab, assumption))
  print(head(risk_tab %>%
    select(taxid, n_proteins, n_tested, factor_log2FC, n_dev, frac_deviating,
           opposing, opposing_frac, assumption), 20))

  at_risk <- risk_tab$taxid[risk_tab$assumption == "AT RISK"]
  if (length(at_risk)) {
    # The message states the rule that was applied, not a stronger one.
    gate(paste("%d taxon(a) AT RISK: size factor moved by more than %.2f log2,",
               "at least %d protein(s) deviate significantly, more than %.0f%%",
               "of those\n      oppose the factor and they are at least %.0f%%",
               "of the taxon's tested proteins"),
         length(at_risk), params$min_lfc, params$assumption_min_deviating,
         100 * params$assumption_opposing,
         100 * params$assumption_opposing_frac)
    note(paste("For those taxa the ratio model may have absorbed real regulation",
               "into the abundance term.\n      The two models bracket the",
               "ambiguity: the ratio model assumes the organism moved, the",
               "abundance\n      model assumes it did not. Report both for",
               "these taxa, and settle the direction with metagenomic\n",
               "     coverage or qPCR for that organism."))
  } else {
    note("no taxon shows a coherent majority shift; the size factors behave like organism abundance")
  }
  if (!is.null(risk_tab) && any(risk_tab$assumption == "AT RISK" &
                                risk_tab$n_opposing == 1, na.rm = TRUE))
    note(paste("some taxa are flagged on a single opposing protein. That is",
               "the weakest form of this evidence;\n      the table is sorted",
               "so the better-supported ones come first."))

  # The picture that makes it judgeable: where does the size factor sit
  # relative to the bulk of its taxon's proteins?
  top <- risk_tab %>% slice_max(n_tested, n = 12, with_ties = FALSE)
  pd <- memb %>% filter(taxid %in% top$taxid) %>%
    left_join(top %>% select(taxid, assumption), by = "taxid") %>%
    mutate(taxid = factor(taxid, levels = top$taxid))
  if (nrow(pd) > 0) {
    print(ggplot(pd, aes(taxid, naive)) +
      geom_hline(yintercept = 0, colour = "grey80") +
      geom_boxplot(aes(fill = assumption), outlier.size = 0.4, alpha = 0.85) +
      geom_point(data = top %>% mutate(taxid = factor(taxid, levels = top$taxid)),
                 aes(taxid, factor_log2FC), shape = 23, size = 2.6,
                 fill = "white", colour = "black") +
      coord_flip() +
      labs(x = "taxon", y = sprintf("log2FC, protein abundance (%s)", PRIMARY_LABEL),
           fill = NULL,
           title = "Size factor (diamond) against its taxon's protein distribution",
           subtitle = "the diamond sits at the median by construction; what matters is whether the proteins that miss it fall on one side"))
  }
}
```

```{r merge}
res <- aqk %>%
  transmute(group_id, bin,
            has_ko              = col_or_na(aqk, "has_ko"),
            surface_or_secreted = col_or_na(aqk, "surface_or_secreted"),
            export_score        = col_or_na(aqk, "export_score"),
            sp_class            = col_or_na(aqk, "sp_class"),
            n_tmb               = col_or_na(aqk, "n_tmb"),
            lpxtg               = col_or_na(aqk, "lpxtg"),
            small_protein       = col_or_na(aqk, "small_protein"),
            description         = col_or_na(aqk, "description"),
            pfam_hits           = col_or_na(aqk, "pfam_hits"),
            ncbifam_hits        = col_or_na(aqk, "ncbifam_hits"),
            kofam_ko            = col_or_na(aqk, "kofam_ko"),
            interpro_sigs       = col_or_na(aqk, "interpro_sigs"),
            jackhmmer_hit       = col_or_na(aqk, "jackhmmer_hit"),
            interpro_ipr        = col_or_na(aqk, "interpro_ipr"),
            ko_source           = col_or_na(aqk, "ko_source"),
            hh_hit              = col_or_na(aqk, "hh_hit"),
            hh_prob             = col_or_na(aqk, "hh_prob"),
            fold_cluster        = col_or_na(aqk, "fold_cluster"),
            fold_cluster_size   = col_or_na(aqk, "fold_cluster_size"),
            context_pul         = col_or_na(aqk, "context_pul"),
            foldseek_desc       = col_or_na(aqk, "foldseek_desc"),
            toxin_fold          = col_or_na(aqk, "toxin_fold"),
            context_flags       = col_or_na(aqk, "context_flags"),
            family_id           = col_or_na(aqk, "family_id"),
            seed_taxid          = col_or_na(aqk, "seed_taxid"),
            effective_taxid     = col_or_na(aqk, "effective_taxid"),
            taxonomy_verdict    = col_or_na(aqk, "taxonomy_verdict"),
            unipept_name        = col_or_na(aqk, "unipept_name"),
            bin_conflict        = col_or_na(aqk, "bin_conflict"),
            kegg_visible        = col_or_na(aqk, "kegg_enrichment_visible")) %>%
  left_join(de, by = "group_id") %>%
  mutate(taxid_key = as.character(col_or_na(aqk, "effective_taxid"))) %>%
  left_join(if (!is.null(risk_tab))
              risk_tab %>% transmute(taxid_key = taxid,
                                     taxon_factor_log2FC = factor_log2FC,
                                     taxon_frac_deviating = frac_deviating,
                                     taxon_assumption = assumption)
            else tibble(taxid_key = character(), taxon_factor_log2FC = numeric(),
                        taxon_frac_deviating = numeric(),
                        taxon_assumption = character()),
            by = "taxid_key") %>%
  mutate(contrast   = PRIMARY_LABEL,
         model_used = if (!adj_attempted) "abundance (ratio model not run)"
                      else ifelse(is.na(sig_adj), "abundance (no usable taxon)",
                                  "ratio to source organism"),
         significant = ifelse(is.na(sig_adj), sig_naive, sig_adj),
         logFC       = ifelse(is.na(logFC_adj), logFC_naive, logFC_adj),
         FDR         = ifelse(is.na(FDR_adj), FDR_naive, FDR_adj))
```

# Where the signal sits

```{r volcano, fig.height=5}
ggplot(res, aes(logFC, -log10(FDR))) +
  geom_point(aes(colour = bin, alpha = significant %in% TRUE), size = 1.2) +
  geom_hline(yintercept = -log10(params$fdr), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = c(-1, 1) * params$min_lfc, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = BIN_COLS, name = NULL) +
  scale_alpha_manual(values = c(`FALSE` = 0.22, `TRUE` = 0.9), guide = "none") +
  labs(title = sprintf("%s  (%s)", PRIMARY_LABEL,
                       if (adj_attempted) "taxon-adjusted where possible"
                       else "unadjusted"))
```

```{r visibility}
de_only <- filter(res, significant %in% TRUE)
cat(sprintf("%d significant groups; %.1f%% of them carry no KEGG pathway and\nare invisible to pathway enrichment.\n",
            nrow(de_only), 100 * mean(!de_only$kegg_visible, na.rm = TRUE)))
count(de_only, bin, .drop = FALSE) %>% mutate(pct = round(100 * n / sum(n), 1))
```

Whether a bin is enriched among the significant proteins is asked once per
model, never on the mixed `significant` column. A protein without a usable
taxon is only ever tested by the abundance model, and that is exactly what
`4_dark`, `3s_structure_only` and `3p_profile_only` are full of: pooling the
two would let the dark bin look "more regulated" for no reason but never
having been adjusted.

```{r bin-enrichment}
# Testability is a property of the bin, so it is reported next to every
# enrichment: a bin that is 5% testable under model 2 cannot be compared with
# one that is 90% testable.
testable_by_bin <- res %>% group_by(bin) %>%
  summarise(n = n(),
            n_testable_ratio = if (adj_attempted) sum(!is.na(sig_adj)) else 0L,
            pct_testable_ratio = round(100 * n_testable_ratio / n, 1),
            .groups = "drop")
print(testable_by_bin)

bin_fisher_one <- function(d, sig, model) {
  map_dfr(BIN_LEVELS, function(b) {
    tt <- table(factor(d$bin == b, c(FALSE, TRUE)),
                factor(sig %in% TRUE, c(FALSE, TRUE)))
    if (any(dim(tt) < 2))
      return(tibble(model = model, bin = b, n_tested = nrow(d),
                    n_in_bin = sum(d$bin == b, na.rm = TRUE),
                    n_significant = NA_integer_, odds_ratio = NA_real_,
                    p = NA_real_))
    ft <- fisher.test(tt)
    tibble(model = model, bin = b, n_tested = nrow(d),
           n_in_bin = sum(d$bin == b, na.rm = TRUE), n_significant = tt[2, 2],
           odds_ratio = unname(ft$estimate), p = ft$p.value)
  }) %>% mutate(FDR = p.adjust(p, "BH")) %>% arrange(p)
}
bin_fisher <- bin_fisher_one(res, res$sig_naive, "abundance (all proteins)")
if (adj_attempted) {
  ra <- filter(res, !is.na(sig_adj))
  bin_fisher <- bind_rows(bin_fisher,
                          bin_fisher_one(ra, ra$sig_adj,
                                         "ratio (testable proteins only)"))
}
bin_fisher
```

# Enrichment without KEGG

```{r enrichment}
# One model per call. Mixing model 1's p-values for the untestable proteins
# with model 2's for the rest makes any term that is concentrated in the dark
# bins look enriched by construction.
run_enrich <- function(t2g_file, label, sig_ids, universe_ids) {
  if (!file.exists(t2g_file)) { note("absent: %s", t2g_file); return(NULL) }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    note("clusterProfiler not installed; skipping %s", label); return(NULL) }
  if (!length(sig_ids)) {
    note("no significant proteins under this model; skipping %s", label)
    return(NULL)
  }
  t2g <- read_tsv_full(t2g_file) %>% filter(gene %in% universe_ids)
  sizes <- count(t2g, term, name = "size")
  # %s, not %d: median() of an even number of integer sizes is x.5, and
  # sprintf refuses %d for a non-integer double — the knit died here, after
  # all the modelling and before anything was exported.
  cat(sprintf("\n%s: %d terms, %d with >= %d members (median size %s)\n",
              label, nrow(sizes), sum(sizes$size >= params$family_min_size),
              params$family_min_size, format(median(sizes$size))))
  t2g <- semi_join(t2g, filter(sizes, size >= params$family_min_size), by = "term")
  if (!nrow(t2g)) { note("no terms large enough for %s", label); return(NULL) }
  e <- clusterProfiler::enricher(
    gene = sig_ids, universe = universe_ids,
    TERM2GENE = t2g[, c("term", "gene")], pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = params$family_min_size, maxGSSize = 2000)
  if (is.null(e) || !nrow(as.data.frame(e))) return(NULL)
  as_tibble(as.data.frame(e)) %>% arrange(p.adjust) %>%
    mutate(model = label) %>% head(30)
}
enr_of <- function(t2g_file, what) {
  ab <- run_enrich(t2g_file, sprintf("%s, abundance model", what),
                   res$group_id[res$sig_naive %in% TRUE], res$group_id)
  rt <- NULL
  if (adj_attempted) {
    uni <- res$group_id[!is.na(res$sig_adj)]
    rt <- run_enrich(t2g_file, sprintf("%s, ratio model", what),
                     res$group_id[res$sig_adj %in% TRUE], uni)
  }
  list(abundance = ab, ratio = rt)
}
enr_family <- enr_of(file.path(RD, "quant", "term2gene_family.tsv"), "de novo families")
enr_pfam   <- enr_of(file.path(RD, "quant", "term2gene_pfam.tsv"), "Pfam")
enr_family$abundance
enr_family$ratio
enr_pfam$abundance
enr_pfam$ratio
```

At 50% identity most metagenome families are singletons, so the size
distribution printed above decides whether any family-level p-value is worth
reading.

# Effector candidates

Not the top of `export_score`: significant, no KO, and predicted to reach
the host. A protein that cannot be exported cannot act on the epithelium
whatever its fold suggests, so topology gates the list and the score only
orders it.

```{r shortlist}
short <- res %>%
  filter(significant %in% TRUE, has_ko %in% FALSE, surface_or_secreted %in% TRUE) %>%
  arrange(desc(export_score), FDR) %>%
  select(group_id, bin, contrast, logFC, FDR, model_used,
         # Both models' calls travel with the row: `significant` is the
         # adjusted call where there was one and the unadjusted call where
         # there was not, so on its own it cannot be compared across rows.
         sig_naive, sig_adj, logFC_naive, FDR_naive, logFC_adj, FDR_adj,
         export_score,
         sp_class, n_tmb, lpxtg, small_protein, toxin_fold, pfam_hits,
         ncbifam_hits, kofam_ko, interpro_ipr, interpro_sigs, hh_hit, hh_prob,
         jackhmmer_hit, fold_cluster,
         fold_cluster_size, context_pul, foldseek_desc, context_flags,
         seed_taxid, unipept_name, taxonomy_verdict,
         taxon_factor_log2FC, taxon_assumption, bin_conflict)
# metaannot no longer ingests external effector predictions. This survives so
# that a pred_* column a user joined in themselves still reaches the shortlist
# as a COLUMN - deliberately not as a term in the score, which is the correct
# status for somebody else's model.
pred_cols <- grep("^pred_", names(aqk), value = TRUE)
if (length(pred_cols)) {
  short <- bind_cols(short, aqk[match(short$group_id, aqk$group_id), pred_cols,
                               drop = FALSE])
  note("external prediction column(s) carried through: %s",
       paste(pred_cols, collapse = ", "))
}
cat(sprintf("%d candidates (significant, no KO, secreted or surface-exposed)\n", nrow(short)))
if (adj_attempted && nrow(short))
  note(paste("%d of them were called by the abundance model because they have",
             "no usable taxon (model_used says which).\n      They were never",
             "adjusted for organism abundance, so they are not on the same",
             "footing as the rest of the list."),
       sum(short$model_used == "abundance (no usable taxon)"))
head(short, 40)
```

```{r shortlist-plot, fig.height=4}
if (nrow(short)) {
  ggplot(short, aes(export_score, -log10(FDR), colour = bin)) +
    geom_point(aes(size = abs(logFC)), alpha = 0.8) +
    scale_colour_manual(values = BIN_COLS, name = NULL) +
    scale_size_continuous(name = "|log2FC|", range = c(1, 5)) +
    labs(x = "export score", title = "KO-less, secreted, differentially abundant")
}
```

# Export

```{r export}
# The label can be an MSstats comparison name with spaces or slashes in it,
# so it is sanitised before it becomes a filename.
FILE_TAG <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9._-]+", "_", PRIMARY_LABEL))
write_tsv(res,        file.path(OUT, sprintf("differential_abundance_%s.tsv", FILE_TAG)))
write_tsv(short,      file.path(OUT, sprintf("effector_candidates_%s.tsv", FILE_TAG)))
write_tsv(bin_tab,    file.path(OUT, "bin_composition.tsv"))
write_tsv(bin_fisher, file.path(OUT, "bin_enrichment.tsv"))
if (!is.null(risk_tab))
  write_tsv(risk_tab, file.path(OUT, "taxon_normalisation_risk.tsv"))
if (!is.null(F1))
  for (cn in names(CTR))
    write_tsv(tt_of(F1$fit, cn), file.path(OUT, sprintf("abundance_model_%s.tsv", cn)))
# enrichment_*.tsv keeps its name and holds the abundance model, which is the
# only one every protein was tested under; the ratio model gets its own file
# rather than being merged into it.
if (!is.null(enr_family$abundance))
  write_tsv(enr_family$abundance, file.path(OUT, "enrichment_family.tsv"))
if (!is.null(enr_family$ratio))
  write_tsv(enr_family$ratio, file.path(OUT, "enrichment_family_ratio.tsv"))
if (!is.null(enr_pfam$abundance))
  write_tsv(enr_pfam$abundance, file.path(OUT, "enrichment_pfam.tsv"))
if (!is.null(enr_pfam$ratio))
  write_tsv(enr_pfam$ratio, file.path(OUT, "enrichment_pfam_ratio.tsv"))

# Kept next to the numbers it produced: the log scrolls away, and a
# design_record.txt that does not say where the condition came from, or what
# was done with a reference channel, cannot be checked afterwards.
notes_path <- file.path(RD, "quant", "design_notes.txt")
writeLines(c(
  paste("design_formula:", params$design_formula),
  paste("coefficients:  ", paste(colnames(design), collapse = ", ")),
  if ("plex" %in% all.vars(as.formula(params$design_formula)))
    paste("plex:           in the model as a batch term;",
          "contrasts are over the condition") else NULL,
  if (file.exists(notes_path)) readLines(notes_path, warn = FALSE) else NULL,
  paste("block:         ", if (nzchar(params$block_col)) params$block_col else "none"),
  paste("abundance model:", if (USE_MS) paste("MSstats:", params$msstats_comparison)
        else "limma"),
  paste("correlation:   ", if (!is.null(F1) && !is.null(F1$correlation))
        round(F1$correlation, 4) else "NA"),
  paste("contrasts:     ", paste(sprintf("%s = %s", names(CTR), CTR), collapse = "; ")),
  paste("primary:       ", PRIMARY_LABEL),
  paste("primary contrast:", PRIMARY),
  paste("normalise:     ", params$normalise),
  paste("min_valid:     ", params$min_valid_per_group),
  paste("min_plexes:    ", params$min_plexes,
        if (!IS_ISOBARIC) "(not an isobaric run; not applied)"
        else if (is.null(PLEX_OF)) "(no per-sample plex; not applied)"
        else "(protein level, counted across plexes)"),
  paste("FDR / min_lfc: ", params$fdr, "/", params$min_lfc),
  paste("adjusted model:", if (adj_attempted) "fitted" else "not run")
), file.path(OUT, "design_record.txt"))
cat("written to ", OUT, "\n", sep = "")
```

# What this analysis does not settle

```{r caveat-gene-calling, results='asis'}
# Whether small ORFs were in the search space is a property of the database,
# not a constant: a database with a dedicated smORF tier (uhgpSM_, AMPSphere)
# puts them there, and printing the blanket caveat would then state something
# false about the very run being reported.
aa_len  <- suppressWarnings(as.numeric(col_or_na(aq, "length")))
n_small <- sum(aa_len <= 100, na.rm = TRUE)
cat("- **Gene calling.** Prodigal in meta mode has a hard floor at 90 nt and ",
    "calls\n  genes below ~100 aa with reduced sensitivity, rather than ",
    "discarding them.", sep = "")
if (n_small > 0) {
  cat(sprintf(paste0(" %d\n  quantified group(s) here are <= 100 aa, so small proteins were in the\n",
                     "  search space — but a plain Prodigal call under-samples them, and a\n",
                     "  smORF-augmented database is what makes this class countable.\n\n"),
              n_small))
} else {
  cat(" No quantified group here is <= 100 aa, so\n  bacteriocins, TA toxins",
      " and RiPPs were effectively never in the search space\n  and their",
      " absence is not evidence of absence.\n\n", sep = "")
}
```

```{r caveat-ratio-compression, results='asis'}
# Printed only for isobaric input, because it is false for label-free: this is
# a property of measuring several samples in ONE MS2 scan.
if (IS_ISOBARIC) {
  cat("- **Ratio compression, uncorrected.** Reporter ions are read from a\n",
      "  spectrum whose precursor window admitted more than one peptide, so\n",
      "  every channel carries some signal from co-isolated species that do\n",
      "  not share the true fold change. The measured ratio is therefore\n",
      "  pulled toward 1: log2 fold changes here are LOWER BOUNDS on the\n",
      "  real ones, and the shrinkage is not a constant — it is worse for\n",
      "  low-abundance proteins, in crowded windows, and in exactly the\n",
      "  strain-redundant regions of a metagenome database where one\n",
      "  peptide's neighbours are its own near-identical paralogues.\n",
      "  metaannot does **not** correct for it: there is no interference\n",
      "  model and no purity-weighted rescaling here, because every such\n",
      "  correction divides by an estimate of the contamination and turns a\n",
      "  known bias into an unknown variance. The consequences are that\n",
      "  direction and ranking are more trustworthy than magnitude, that\n",
      "  `min_lfc` is a stricter filter on these data than on label-free\n",
      "  data, and that an effect size read off this report should not be\n",
      "  compared with one from a label-free experiment. `tmt.min_purity`\n",
      "  limits how co-isolated the accepted spectra were; it does not undo\n",
      "  the compression in the ones that pass.\n\n", sep = "")
}
```

- **Group-level annotation.** Rows with `bin_conflict` inherit the leading
  protein's annotation, and KO-less proteins are disproportionately
  strain-specific.
- **The taxon adjustment is a ratio, not a mechanism.** It removes the most
  common confound; it does not show that a surviving change is regulatory
  rather than a shift in strain composition within the taxon.
- **Compositional non-identifiability.** Within a taxon, "the organism doubled"
  and "its proteome doubled" are the same numbers. `taxon_normalisation_risk.tsv`
  says which taxa that ambiguity actually bites for; only an orthogonal measure
  of organism abundance resolves it.
- **A Foldseek hit is a fold, not a function.** `toxin_fold` is a hypothesis
  for an assay, not an annotation.

```{r session}
sessionInfo()
```
"""


# ----------------------------------------------------------------------
# R literals
# ----------------------------------------------------------------------
# Reserved words cannot be names in R, so make.names() appends a dot to them.
_R_RESERVED = frozenset((
    "if", "else", "repeat", "while", "function", "for", "next", "break",
    "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_",
    "NA_character_"))


def r_make_names(name):
    """Python port of R's make.names(), applied to the same strings R will.

    model.matrix() pastes a factor level onto the variable name and the report
    then runs make.names() over the result, so the level "high fiber" in
    `group` is the coefficient `grouphigh.fiber` and "2wk" alone would be
    `X2wk`. A contrast written against the raw level names therefore names a
    coefficient that does not exist and the knit stops with "contrast term(s)
    not among the design coefficients".

    The character class is unicode-aware because R's is: in a UTF-8 locale
    make.names() keeps an accented letter and only mangles punctuation. Any
    residual mismatch is loud rather than silent - the report prints every
    available coefficient when a contrast term is not among them.
    """
    s = "".join(c if (c.isalnum() or c in "._") else "." for c in str(name))
    if not s:
        return "X"
    if s[0].isdigit() or s[0] == "_" or (s[0] == "." and s[1:2].isdigit()):
        s = "X" + s
    if s in _R_RESERVED:
        s += "."
    return s


def r_literal(v):
    """Render a Python value as R source text.

    The params block used to be plain YAML, which meant every value inherited
    YAML's escaping rules on the way into a document that reads them as R
    data: a Windows path lost its backslashes, an apostrophe closed the
    scalar, and a numeric the user left empty arrived as the character "" and
    broke the first comparison it reached. Each param is now emitted as
    `!r <expression>`, which rmarkdown evaluates as R, so the value the
    document sees is exactly the value written here - with its type.
    """
    if v is None:
        return "NULL"
    if isinstance(v, bool):                      # before int: bool IS an int
        return "TRUE" if v else "FALSE"
    if isinstance(v, int):
        # `3L`, not `3`: YAML handed R an integer here, and a param that
        # quietly became a double would change nothing today but would be a
        # different type for any downstream code that tests for one.
        return f"{v}L"
    if isinstance(v, float):
        if v != v:
            return "NA_real_"
        if v == float("inf"):
            return "Inf"
        if v == float("-inf"):
            return "-Inf"
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "c(" + ", ".join(r_literal(x) for x in v) + ")" if v else "c()"
    if isinstance(v, dict):
        die(f"a report parameter is a mapping ({v!r}); the Rmd's params are "
            "scalars, so write it as a string or a list")
    s = str(v)
    for a, b in (("\\", "\\\\"), ('"', '\\"'), ("\n", "\\n"), ("\r", "\\r"),
                 ("\t", "\\t")):
        s = s.replace(a, b)
    return '"' + s + '"'


def _yaml_r_param(key, value):
    """One `key: !r <expr>` line of the generated params block.

    The expression is wrapped in YAML single quotes (doubling any of its own):
    a double-quoted YAML scalar would apply YAML's backslash escapes to the R
    source and eat the very characters r_literal() just escaped.
    """
    return "  %s: !r '%s'" % (key, r_literal(value).replace("'", "''"))


def _formula_terms(formula):
    """Main-effect terms on the right of `~`, in the order they are written.

    Interactions and the intercept markers are dropped: neither names a single
    factor whose levels could be contrasted pairwise.
    """
    rhs = str(formula or "").split("~", 1)[-1]
    out = []
    for part in re.split(r"[+\-]", rhs):
        t = part.strip().strip("`")
        if not t or t in ("0", "1"):
            continue
        if not re.fullmatch(r"[A-Za-z._][A-Za-z0-9._]*", t):
            continue                    # an interaction, a call, an offset
        if t not in out:
            out.append(t)
    return out


def _has_intercept(formula):
    """`~ 0 + group` and `~ group - 1` drop the intercept; `~ group` keeps it."""
    rhs = str(formula or "").split("~", 1)[-1]
    return not re.search(r"(^|[+~])\s*0\s*(\+|$)|-\s*1\b", rhs)


def auto_contrasts(design_path, formula, factor_cols=""):
    """Every pairwise contrast between levels of the design's primary factor,
    so a run with a manifest needs no hand-written contrast to produce a
    result.

    The factor is the FIRST main-effect term of design_formula, not the last
    identifier in it: `~ 0 + group + sex` describes group adjusted for sex, and
    taking the last identifier silently produced contrasts over sex - a
    published number for the wrong comparison. When analysis.factor_cols
    disagrees with the formula about which term that is, the two statements of
    intent conflict and stopping beats guessing.
    """
    if not os.path.exists(design_path):
        return ""
    d = pd.read_csv(design_path, sep="\t", encoding="utf-8", encoding_errors="replace")
    terms = _formula_terms(formula)
    factors = [c.strip() for c in str(factor_cols or "").split(",") if c.strip()]
    cands = [t for t in terms if t in d.columns]
    if not cands and "group" in d.columns:
        # No usable term (an empty or unparsable formula): the manifest-derived
        # design only ever has `group`, so that is the honest fallback.
        cands = ["group"]
    if not cands:
        log(f"report: design_formula '{formula}' names no column of "
            f"{design_path} ({', '.join(map(str, d.columns))}), so no "
            "contrast could be derived", "WARN")
        return ""
    if factors:
        declared = [t for t in cands if t in factors]
        if declared and declared[0] != cands[0]:
            die(f"cannot tell which factor the contrasts are meant to "
                f"compare: design_formula '{formula}' puts '{cands[0]}' first, "
                f"but analysis.factor_cols declares '{declared[0]}' and not "
                f"'{cands[0]}'. Set analysis.contrasts explicitly, or put the "
                "factor of interest first in the formula.")
        if declared:
            cands = declared
        else:
            log(f"report: analysis.factor_cols ({factor_cols}) names no term "
                f"of design_formula '{formula}'; using '{cands[0]}'", "WARN")
    var = cands[0]
    levels = [str(x) for x in sorted(d[var].dropna().unique()) if str(x).strip()]
    if len(levels) < 2:
        log(f"report: '{var}' has {len(levels)} level(s) in {design_path}; "
            "a contrast needs two", "WARN")
        return ""
    ref = levels[0]
    log(f"report: contrasts over '{var}' ({len(levels)} levels, reference "
        f"'{ref}'), the first factor in design_formula '{formula}'")
    renamed = []

    def coef(level):
        raw = f"{var}{level}"
        name = r_make_names(raw)
        if name != raw:
            renamed.append((raw, name))
        return name

    intercept = _has_intercept(formula)
    if intercept:
        # With an intercept the reference level HAS no coefficient of its own,
        # so `a - ref` names something that is not in the design and the knit
        # stops. Each level's own coefficient already is the contrast against
        # the model's baseline.
        log(f"report: design_formula '{formula}' keeps the intercept, so each "
            "contrast is a single coefficient against the model's own "
            "reference level")
    out = []
    for lev in levels[1:]:
        lhs = re.sub(r"[^A-Za-z0-9_.]", ".", f"{lev}_vs_{ref}")
        rhs = coef(lev) if intercept else f"{coef(lev)} - {coef(ref)}"
        out.append(f"{lhs} = {rhs}")
    for raw, name in dict.fromkeys(renamed):
        log(f"report: level name gives the non-syntactic coefficient '{raw}'; "
            f"R will call it '{name}' (make.names), so the contrast uses that")
    return "; ".join(out)


def template_params():
    """Parameters the embedded report actually uses.

    Read from `params$x` references in the body, not from the YAML header:
    the header's parameter block is the placeholder this tool fills in, so
    there is nothing to parse there. Usage is the real requirement anyway —
    a param the code reads but the header omits is NULL at knit time.
    """
    return sorted(set(re.findall(r"params\$([A-Za-z_][A-Za-z0-9_.]*)",
                                 RMD_TEMPLATE)))


# Checked as soon as the template exists, not lazily: the bin vocabulary is
# duplicated between Python and the embedded R, and a drift is invisible in
# the rendered report.
_check_bin_vocabulary()


# Params that name a file or directory. rmarkdown::render() evaluates the
# document with the working directory set to the Rmd's own folder, so anything
# relative in the generated header resolves against results/analysis/ and is
# not found. The generated document therefore carries absolute paths.
REPORT_PATH_PARAMS = ("results_dir", "metadata", "msstats_comparison")

# Params that read like paths but are deliberately relative. out_subdir is
# joined onto results_dir inside the document (and passed to build_object.R as
# a bare name); absolutising it would break both.
REPORT_RELATIVE_PARAMS = ("out_subdir",)

# Anything named like a path is treated as one, so a path-valued param added
# later does not silently inherit the bug this list was written to fix.
_PATHISH_SUFFIXES = ("_dir", "_file", "_path", "_tsv", "_csv", "_fasta",
                     "_faa", "_fa", "_txt")


def _is_report_path_param(name):
    if name in REPORT_RELATIVE_PARAMS:
        return False
    return name in REPORT_PATH_PARAMS or str(name).endswith(_PATHISH_SUFFIXES)


def _render_params_block(a):
    """The generated `params:` block, and a check that it says what we mean.

    Every value is emitted as an R expression behind a `!r` tag, so the block
    is read back with a loader that keeps the tag: if YAML quoting mangled an
    expression on the way in, the round-trip mismatch is caught here rather
    than as an unparsable header three stages later.
    """
    lines = [_yaml_r_param(k, v) for k, v in a.items()]
    block = "\n".join(lines)

    class _KeepRTag(yaml.SafeLoader):
        pass

    _KeepRTag.add_constructor(
        "!r", lambda ldr, node: ("!r", ldr.construct_scalar(node)))
    try:
        back = yaml.load("params:\n" + block + "\n", Loader=_KeepRTag)
    except yaml.YAMLError as e:
        die(f"the generated report header is not valid YAML: {e}\n"
            "This is a bug in metaannot; the offending parameters are:\n  "
            + "\n  ".join(lines))
    got = (back or {}).get("params") or {}
    for k, v in a.items():
        want = ("!r", r_literal(v))
        if got.get(k) != want:
            die(f"the generated report header does not round-trip parameter "
                f"'{k}': wrote {want!r}, read back {got.get(k)!r}. This is a "
                "bug in metaannot.")
    return block


def _tmt_report_design(cfg, a, design_auto):
    """Point the report's model at the condition, with the plex as a batch.

    A TMT design has one structural difference from every label-free one this
    tool reads: the thing the input names per sample is the PLEX, and the plex
    is a batch. Left to the ordinary defaults the contrasts would come out
    plex-versus-plex — a batch effect presented as a hypothesis — so the
    default formula gains `+ plex`, and a condition that is nowhere to be
    found is refused rather than approximated.

    Only the DEFAULTS are changed. A formula the user wrote is their
    statement of the model and is left exactly as written.
    """
    d = None
    if os.path.exists(design_auto):
        d = pd.read_csv(design_auto, sep="\t", dtype=str,
                        encoding="utf-8", encoding_errors="replace")
    plexes = sorted(set(d["plex"].dropna())) if (
        d is not None and "plex" in d.columns) else []
    defaults = DEFAULT_CONFIG["analysis"]
    if len(plexes) > 1:
        if a.get("design_formula") == defaults["design_formula"]:
            a["design_formula"] = "~ 0 + group + plex"
            log(f"report: {len(plexes)} TMT plexes, so design_formula is "
                f"'{a['design_formula']}' — the condition is the hypothesis "
                "and the plex is a batch term that absorbs it. Set "
                "analysis.design_formula to override")
        # factor_cols follows the FORMULA, not the plex count. The two used to
        # be decided independently, so a user who wrote their own formula
        # without the plex ("~ 0 + group + sex") kept that formula and still
        # got factor_cols 'group,plex' — and the report's own check on
        # factor_cols then aborted the knit over a column their metadata had
        # no reason to carry, naming a key they never set. The whole right-hand
        # side is searched rather than the main-effect terms, because
        # "group * plex" models the plex too and it still has to be typed.
        rhs = str(a.get("design_formula", "")).split("~", 1)[-1]
        if (a.get("factor_cols") == defaults["factor_cols"]
                and re.search(r"(?<![\w.])plex(?![\w.])", rhs)):
            a["factor_cols"] = "group,plex"
            log(f"report: factor_cols is '{a['factor_cols']}', so the plex is "
                "typed as a factor; an untyped plex code like '1' would be "
                "fitted as a continuous slope through the batches")
    elif plexes:
        # model.matrix() cannot make a contrast for a one-level factor, so a
        # single-plex run must not carry the term at all.
        log(f"report: a single plex ({plexes[0]}), so plex is not added to "
            "design_formula — one level is not a batch effect", "WARN")

    terms = _formula_terms(a.get("design_formula", ""))
    meta_path = a.get("metadata") or design_auto
    same = os.path.abspath(meta_path) == os.path.abspath(design_auto)
    have = None
    if os.path.exists(meta_path):
        have = list(pd.read_csv(meta_path, sep="\t", nrows=0,
                                encoding="utf-8",
                                encoding_errors="replace").columns)
    if have is None or "group" not in terms:
        return
    rows = []
    if d is not None:
        cols = [c for c in ("sample", "plex", "channel") if c in d.columns]
        rows = [d[cols].iloc[i].tolist() for i in range(min(2, len(d)))]
        head = "\t".join(cols + ["group"])
    else:
        head = "sample\tplex\tgroup"
    if "group" not in have:
        die("the FragPipe TMT files do not carry the condition, and none was "
            f"derived from the sample names, so {meta_path} has no 'group' "
            "column and there is nothing to contrast. The plex is a TMT "
            "batch, not a condition, and metaannot will not use it as one — a "
            "plex-versus-plex contrast is a batch effect presented as a "
            "hypothesis.\n\n"
            "Write the conditions down. The recovered design is the head "
            "start, since it already names every sample and its plex:\n\n"
            f"    cp {design_auto} metadata.tsv\n\n"
            "then add a 'group' column, one condition per sample:\n\n"
            f"    {head}\n" +
            "".join(f"    {chr(9).join(map(str, r))}\t<condition>\n"
                    for r in rows) +
            f"    ...  ({len(d) if d is not None else 0} sample(s) in all)\n\n"
            "and set, in the config:\n\n"
            "    analysis:\n"
            "      metadata: metadata.tsv\n"
            f"      design_formula: \"{a.get('design_formula', '')}\"\n"
            f"      factor_cols: \"{a.get('factor_cols', '')}\"\n\n"
            "If the condition IS in the annotated sample names, set "
            "tmt.condition_from_name to a regular expression with one capture "
            "group and re-run the join stage instead.")
    if "plex" in terms and "plex" not in have:
        die(f"design_formula '{a.get('design_formula')}' models the plex, but "
            f"{meta_path} has no 'plex' column (it has: {have}). The plex is "
            "the batch a TMT experiment has to be adjusted for, and it is "
            f"recorded per sample in {design_auto}" +
            ("" if same else " — copy that column across, keyed on sample") +
            ". Or drop '+ plex' from analysis.design_formula, which reports "
            "the batch effect as biology.")


def write_report_rmd(cfg, p):
    a = dict(cfg.get("analysis") or {})
    missing = [k for k in template_params()
               if k not in a and k not in ("results_dir",)]
    if missing:
        die(f"the analysis config block is missing parameters the report "
            f"declares: {missing}. The generated header would omit them and "
            "the knit would fail on a NULL param.")
    a.setdefault("results_dir", cfg["results_dir"])
    a["results_dir"] = cfg["results_dir"]
    design_auto = f"{p.quant_dir}/design_from_input.tsv"
    if not a.get("metadata"):
        a["metadata"] = design_auto
    # Before the contrasts are derived, not after: for TMT the formula decides
    # which column they are taken over, and the default one would take them
    # over the plex.
    if str(cfg.get("quant_format", "")) == "fragpipe_tmt":
        _tmt_report_design(cfg, a, design_auto)
    if not a.get("contrasts"):
        src = design_auto
        # TMT only, and only because of where the condition lives: the
        # recovered design carries the plex, and the condition is often only
        # in the hand-written metadata, so deriving from the design would find
        # no group and silently produce no contrast. Label-free keeps reading
        # the recovered design, as it always has.
        if (str(cfg.get("quant_format", "")) == "fragpipe_tmt"
                and a.get("metadata") and os.path.exists(a["metadata"])
                and os.path.abspath(a["metadata"]) != os.path.abspath(design_auto)):
            src = a["metadata"]
            log(f"report: contrasts will be derived from analysis.metadata "
                f"({src}), which is where a TMT run's condition is written")
        a["contrasts"] = auto_contrasts(src, a.get("design_formula", ""),
                                        a.get("factor_cols", ""))
        if a["contrasts"]:
            log(f"report: contrasts derived from the design -> {a['contrasts']}")
        else:
            log("report: no contrasts given and none could be derived; set "
                "analysis.contrasts in the config", "WARN")
    for k in list(a):
        if _is_report_path_param(k) and isinstance(a[k], str) and a[k].strip():
            absolute = os.path.abspath(a[k])
            if absolute != a[k]:
                log(f"report: {k} -> {absolute} (the Rmd is knitted from its "
                    "own folder, so a relative path would not resolve)")
            a[k] = absolute
    if yaml is None:
        die("pyyaml is needed to write the report")
    rmd = RMD_TEMPLATE.replace("__METAANNOT_PARAMS__", _render_params_block(a))
    os.makedirs(f"{p.R}/{a.get('out_subdir', 'analysis')}", exist_ok=True)
    out = f"{p.R}/{a.get('out_subdir', 'analysis')}/analyse_metaannot.Rmd"
    # encoding="utf-8" explicitly: the template carries a real minus sign and
    # an em dash, and open()'s default is the locale's codec - cp1252 on
    # Windows, ASCII under a C locale - so the write died before writing
    # anything rather than producing a document.
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(rmd)
    log(f"report: wrote {out}")
    return out


def write_object_script(p, sub="analysis"):
    # Same subdirectory the report writes to, or the script would look for the
    # differential-abundance tables somewhere they were never written.
    os.makedirs(f"{p.R}/{sub}", exist_ok=True)
    path = f"{p.R}/{sub}/build_object.R"
    # utf-8 for the same reason the Rmd is: the locale's codec is not a
    # property of the script being written.
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(ROBJECT_SCRIPT)
    return path


def _run_rscript(cmd, what, hint=""):
    """Run Rscript, and make every way it can fail visible to the caller.

    run_cmd() discards stderr on success, but "Rscript exited 0 and wrote
    nothing" is the failure that used to look like success to a shell script,
    and the diagnosis for it is in the stderr of that successful run. The tail
    is returned so the caller can quote it when the expected output is missing.
    """
    log("$ " + " ".join(str(c) for c in cmd))
    # Launched by the path PATH resolves to, for the reason resolve_tool
    # states: on Windows CreateProcess ignores PATHEXT, so an Rscript.exe
    # later on PATH beat the Rscript.bat that have() and the logged line both
    # named. This is the one launcher that cannot go through run_cmd - the
    # stderr of a SUCCESSFUL run is what diagnoses "Rscript exited 0 and wrote
    # nothing", and run_cmd discards it - so it has to resolve the name
    # itself. The line above keeps the bare name, which is what a reader would
    # type.
    argv = [str(c) for c in cmd]
    argv[0] = resolve_tool(argv[0])
    proc = subprocess.run(argv, stdout=subprocess.DEVNULL,
                          stderr=subprocess.PIPE, text=True)
    tail = "\n".join((proc.stderr or "").strip().splitlines()[-20:])
    if proc.returncode != 0:
        die(f"{what} failed: Rscript exited {proc.returncode}"
            + (f"\n{hint}" if hint else "")
            + f"\n--- Rscript stderr (last 20 lines) ---\n{tail}")
    return tail


def cmd_object(args):
    """Assemble the pipeline's tables into one R object."""
    cfg = load_config(args.config)
    if args.results_dir:
        cfg["results_dir"] = args.results_dir
    p = Paths(cfg)
    p.mkdirs()
    if not os.path.exists(f"{p.quant_dir}/annotated_quant.tsv"):
        die(f"no {p.quant_dir}/annotated_quant.tsv; run the join stage first")
    sub = (cfg.get("analysis") or {}).get("out_subdir", "analysis") or "analysis"
    script = write_object_script(p, sub)
    out = args.out or p.robject
    if args.no_run:
        print(script)
        return 0
    if not have("Rscript"):
        log("Rscript not found; the script was written but not run. Run it "
            f"yourself:  Rscript {script} {p.R} {out} {sub}", "WARN")
        print(script)
        return 0
    hint = (f"The script is at {script}; run it interactively to see the "
            "error in context.")
    tail = _run_rscript(["Rscript", script, p.R, out, sub],
                        "building the R object", hint)
    # An R script that stops inside a tryCatch, or saves to a path it could
    # not create, can still exit 0. Without this the caller is told nothing
    # and a shell script carries on with an .rds that was never written.
    if not nonempty(out):
        die(f"building the R object reported success but wrote no {out}.\n"
            + hint
            + (f"\n--- Rscript stderr (last 20 lines) ---\n{tail}" if tail else ""))
    log(f"R object: {out} ({os.path.getsize(out)/2**20:.1f} MB)")
    return 0


def cmd_report(args):
    cfg = load_config(args.config)
    if args.results_dir:
        cfg["results_dir"] = args.results_dir
    p = Paths(cfg)
    p.mkdirs()
    rmd = write_report_rmd(cfg, p)
    if args.no_render:
        print(rmd)
        return 0
    if not have("Rscript"):
        log("Rscript not found; the Rmd was written but not rendered. Knit it "
            "in RStudio, or install R and rmarkdown.", "WARN")
        print(rmd)
        return 0
    rmd_abs = os.path.abspath(rmd)
    html = os.path.splitext(os.path.basename(rmd))[0] + ".html"
    # r_literal(), not an f-string: a results directory holding a backslash
    # (every Windows path) or a quote used to produce an R expression that
    # would not parse, or - worse - one that parsed as something else.
    expr = (f"rmarkdown::render({r_literal(rmd_abs)}, "
            f"output_file={r_literal(html)}, quiet=TRUE)")
    hint = (f"The Rmd is at {rmd}; knit it interactively to see the error in "
            "context.")
    tail = _run_rscript(["Rscript", "-e", expr], "rendering the report", hint)
    out_html = os.path.join(os.path.dirname(rmd_abs), html)
    if not nonempty(out_html):
        die(f"rendering reported success but wrote no {out_html}.\n" + hint
            + (f"\n--- Rscript stderr (last 20 lines) ---\n{tail}" if tail else ""))
    log(f"report: {out_html}")
    return 0


def cmd_all(args):
    rc = cmd_run(args)
    if rc or getattr(args, "dry_run", False):
        return rc
    cfg = load_config(args.config)
    if args.results_dir:
        cfg["results_dir"] = args.results_dir
    aq = f"{Paths(cfg).quant_dir}/annotated_quant.tsv"
    if not os.path.exists(aq):
        # `all --only pfam` runs the named stages and nothing else, so join
        # never ran and there is no quantified table for either consumer to
        # read. That is the run the user asked for, not a failure: say what is
        # being skipped instead of ending in a FATAL from cmd_object.
        if getattr(args, "only", None) or getattr(args, "from_stage", None):
            log(f"all: {aq} does not exist because the stage selection did "
                "not include join, so the report and the object are skipped. "
                "Rerun `all` without --only/--from, or run `run --only join` "
                "first.", "WARN")
            return 0
        die(f"the pipeline finished but wrote no {aq}, so there is nothing "
            "for the report to read; see the join stage in the log above.")
    # The object folds the report's differential-abundance tables into
    # rowData, so it runs second: the other order writes an object with no DE
    # columns on every first run.
    log("all: running the report first, then the object (the object reads the "
        "report's differential-abundance tables)")
    rc = cmd_report(args)
    if rc:
        return rc
    return cmd_object(args)

# ======================================================================
# requirements: one registry that both reports status and fixes it
# ======================================================================
# doctor renders its status lines and its install commands from the same
# entries. A separate list of fix commands would drift from the checks the
# moment either changed.
#
# Download URLs live in the config (`sources:`) rather than here, because a
# moved URL should be a config edit, not a patch to the tool. Anything that
# needs a licence or a click-through is marked `manual` and never attempted.

def _pkg_mgr():
    for m in ("mamba", "micromamba", "conda"):
        if have(m):
            return m
    return "conda"


def _exists(path):
    """A file with content, or a directory that is not empty.

    This used to be `os.path.exists(path) or glob(path + "*")`, which is the
    same as trusting the download: a curl killed halfway leaves the partial
    `.gz`, a gunzip that ran out of disk leaves the `.gz`, an unpacked-but-not-
    concatenated archive leaves the `.tgz`, and every one of those matched
    `path*` and counted as installed forever. Verifying the file the config
    actually names is what makes the post-`--fix` UNVERIFIED check mean
    something.
    """
    try:
        if not path:
            return False
        if os.path.isdir(path):
            return bool(os.listdir(path))
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def _prefix_exists(path):
    """For databases addressed by a prefix rather than a file.

    foldseek and hh-suite name a set of sibling files after one stem, so the
    stem itself may not exist. Only these two may use the glob.
    """
    if not path:
        return False
    return any(os.path.getsize(f) > 0 for f in glob.glob(str(path) + "*")
               if os.path.isfile(f))


def _pressed(hmm):
    """An HMM library that hmmsearch can actually open.

    hmmpress writes four binary siblings; without them `hmmsearch --cut_tc`
    stops with a format error. Requiring them also rejects the HTML landing
    page and the truncated download, which are otherwise indistinguishable
    from the real file by existence alone.
    """
    return _exists(hmm) and all(_exists(hmm + s)
                                for s in (".h3f", ".h3i", ".h3m", ".h3p"))


def q(v):
    """Shell-quote an interpolated value.

    Everything spliced into these commands comes from the config: paths and
    URLs the user wrote. Unquoted, a path containing a space silently breaks
    every command, and one containing `;` or `$(...)` executes.
    """
    return shlex.quote(str(v))


def _dbdir(path):
    return os.path.dirname(path) or "."


def _gb(x):
    """Format a size in GB.

    `:.0f` printed the 0.08 GB taxdump as "~0 GB to download", which is the
    number the docs tell the user to read before approving a download.
    """
    return f"{x:.1f}" if x < 10 else f"{x:.0f}"


def _reject_html(path):
    """A shell guard that deletes a download that is really a web page.

    `curl -fL` only fails on an HTTP error status. A provider that has moved a
    file usually answers 302 and then 200 with an HTML landing page, so curl
    exits 0 and the landing page lands on disk under the database's name. It
    then either fails in hmmpress with a baffling message or, worse, satisfies
    the existence check and counts as installed forever. Checking the first
    bytes costs nothing and turns a silent corruption into a clear error.
    """
    return (f"head -c 512 {q(path)} | grep -qiE '<!doctype html|<html' && "
            f"{{ echo 'ERROR: {path} is an HTML page, not the database. The "
            f"download URL has probably moved; fix it under sources: in "
            f"config.yaml.' >&2; rm -f {q(path)}; exit 1; }} || true")


def requirements(cfg, p):
    """Everything this configuration needs, with how to obtain it."""
    R, db = cfg.get("run") or {}, cfg.get("db") or {}
    src = cfg.get("sources") or {}
    pre = cfg.get("emapper_precomputed") or ""
    pre = [pre] if isinstance(pre, str) and pre else list(pre or [])
    M = _pkg_mgr()
    out = []

    def add(**kw):
        # size_gb is what comes down the wire; disk_gb is what is left after
        # unpacking, when the two differ enough to matter. The confirmation
        # prompt quotes both, because the number the user approves is the one
        # that decides whether the volume survives the night.
        out.append({"manual": None, "size_gb": 0.0, "disk_gb": 0.0, "cmds": [],
                    "kind": "tool", "note": "", **kw})

    # ---- tools -------------------------------------------------------
    if R.get("pfam") or R.get("dbcan") or R.get("ncbifam") or R.get("jackhmmer"):
        add(id="hmmer", label="hmmsearch/jackhmmer", ok=have("hmmsearch"),
            cmds=[f"{M} install -y -c bioconda hmmer"], size_gb=0.05)
    if R.get("diamond"):
        add(id="diamond", label="diamond", ok=have("diamond"),
            cmds=[f"{M} install -y -c bioconda diamond"], size_gb=0.05)
    if R.get("cluster"):
        add(id="mmseqs2", label="mmseqs", ok=have("mmseqs"),
            cmds=[f"{M} install -y -c bioconda mmseqs2"], size_gb=0.1)
    if R.get("structure"):
        add(id="foldseek", label="foldseek", ok=have("foldseek"),
            cmds=[f"{M} install -y -c bioconda foldseek"], size_gb=0.1)
        add(id="esmfold", label="fair-esm (ESMFold)", ok=_pyhas("esm"),
            cmds=['pip install "fair-esm[esmfold]"'], size_gb=3.0,
            note="GPU host only; pulls torch and openfold")
    if R.get("topology"):
        add(id="tmbed", label="tmbed", ok=have("tmbed"),
            cmds=["pip install tmbed", "tmbed download"], size_gb=2.5,
            note="GPU host; the second command fetches ProtT5 weights")
        add(id="signalp6", label="SignalP 6.0", ok=have("signalp6"),
            manual="academic licence required",
            note="register at services.healthtech.dtu.dk, then `pip install <tarball>`")
    if R.get("kofam"):
        add(id="kofamscan", label="KOfamScan", ok=have("exec_annotation"),
            cmds=[f"{M} install -y -c bioconda kofamscan"], size_gb=0.05)
    if R.get("hhblits"):
        add(id="hhsuite", label="hhblits", ok=have("hhblits"),
            cmds=[f"{M} install -y -c bioconda hhsuite"], size_gb=0.2)
    if R.get("interpro"):
        add(id="interproscan", label="InterProScan",
            ok=_exists(db.get("interproscan_sh")),
            manual="large Java distribution, version-specific",
            note=f"see {src.get('interproscan', '')} (that directory has "
                 "404'd before; check the current release path) -- untar, run "
                 "`interproscan.sh -i test_proteins.fasta`, then set "
                 "db.interproscan_sh")
    if R.get("smorf"):
        # SmORFinder installs its CLI as `smorf`, not `smorfinder`, so probing
        # the package name could never find it and this half was always
        # reported (and skipped) as absent.
        add(id="smorf", label="smorfinder (smorf) / macrel",
            ok=have("smorf") or have("macrel"),
            cmds=["pip install macrel", "pip install smorfinder"], size_gb=0.3,
            note="either one is enough to start; SmORFinder's command is "
                 "`smorf`, macrel's is `macrel`")
    if R.get("eggnog") and not pre:
        add(id="eggnog-mapper", label="eggnog-mapper", ok=have("emapper.py"),
            cmds=[f"{M} install -y -c bioconda eggnog-mapper"], size_gb=0.2,
            note="not needed if emapper_precomputed is set")

    # ---- databases ---------------------------------------------------
    def dbentry(rid, label, target, size_gb, build, note="", ok=None, disk_gb=0.0):
        """Register a database, or explain that its path is not configured.

        An unset path used to produce commands like `curl -o .gz <url>` and
        `hmmpress -f ` with no argument, which fail in confusing ways or, worse,
        write a file called `.gz` into the working directory.
        """
        if not target:
            add(id=rid, label=label, kind="db", ok=False,
                manual="no path configured",
                note=f"set db.{rid} to where this should live, then rerun doctor")
            return
        add(id=rid, label=label, kind="db",
            ok=_exists(target) if ok is None else ok,
            size_gb=size_gb, disk_gb=disk_gb, cmds=build(target), note=note)

    if R.get("pfam"):
        t = db.get("pfam_hmm", "")
        dbentry("pfam_hmm", "Pfam-A", t, 0.45, lambda t: [
            f"mkdir -p {q(_dbdir(t))}",
            f"curl -fL -o {q(t + '.gz')} {q(src.get('pfam',''))}",
            _reject_html(t + '.gz'),
            f"gunzip -f {q(t + '.gz')}", f"hmmpress -f {q(t)}"],
            disk_gb=1.7, ok=_pressed(t))
    if R.get("dbcan"):
        t = db.get("dbcan_hmm", "")
        dbentry("dbcan_hmm", "dbCAN HMMs", t, 0.2, lambda t: [
            f"mkdir -p {q(_dbdir(t))}",
            f"curl -fL -o {q(t)} {q(src.get('dbcan',''))}",
            _reject_html(t),
            f"hmmpress -f {q(t)}"],
            ok=_pressed(t),
            note="the dbCAN site has moved before and the old link answered "
                 "302 with an HTML page; confirm sources.dbcan against the "
                 "provider before a download")
    if R.get("ncbifam"):
        # hmm_PGAP.HMM.tgz unpacks to ~19,000 single-model .HMM files, which
        # hmmpress cannot take and hmmsearch cannot open; NCBI publishes the
        # concatenated, press-ready library as hmm_PGAP.LIB in the same
        # directory. The recipe now downloads that one file directly.
        t = db.get("ncbifam_hmm", "")
        dbentry("ncbifam_hmm", "NCBIfam/TIGRFAM HMMs", t, 2.7, lambda t: [
            f"mkdir -p {q(_dbdir(t))}",
            f"curl -fL -o {q(t)} {q(src.get('ncbifam',''))}",
            _reject_html(t),
            f"hmmpress -f {q(t)}"],
            disk_gb=4.0, ok=_pressed(t),
            note="sources.ncbifam must point at the concatenated hmm_PGAP.LIB, "
                 "not at hmm_PGAP.HMM.tgz (a directory of ~19,000 single "
                 "models) and not at the AMR-only AMRFinder library")
    if R.get("kofam"):
        prof, kl = db.get("kofam_profiles", ""), db.get("kofam_ko_list", "")
        if not (prof and kl):
            add(id="kofam_db", label="KOfam profiles + ko_list", kind="db",
                ok=False, manual="no path configured",
                note="set db.kofam_profiles and db.kofam_ko_list")
        else:
            # The tarball lands next to the database, not in /tmp: a fixed
            # /tmp name collides with whoever ran doctor last on a shared
            # server (curl -o then fails with permission denied) and /tmp is
            # often a tmpfs too small for 1.5 GB. It is removed afterwards.
            tar = os.path.join(_dbdir(prof), "kofam_profiles.tar.gz")
            add(id="kofam_db", label="KOfam profiles + ko_list", kind="db",
                ok=_exists(prof) and _exists(kl), size_gb=1.6, disk_gb=5.0,
                cmds=[
                    f"mkdir -p {q(_dbdir(prof))} {q(_dbdir(kl))}",
                    f"curl -fL -o {q(tar)} {q(src.get('kofam_profiles',''))}",
                    f"tar xzf {q(tar)} -C {q(_dbdir(prof))}",
                    f"rm -f {q(tar)}",
                    f"curl -fL -o {q(kl + '.gz')} {q(src.get('kofam_ko_list',''))}",
                    f"gunzip -f {q(kl + '.gz')}"])
    if R.get("jackhmmer"):
        t = db.get("jackhmmer_db", "")
        dbentry("jackhmmer_db", "UniRef50", t, 8.8, lambda t: [
            f"mkdir -p {q(_dbdir(t))}",
            f"curl -fL -o {q(t + '.gz')} {q(src.get('uniref50',''))}",
            _reject_html(t + '.gz'),
            f"gunzip -f {q(t + '.gz')}"],
            disk_gb=27.0,
            # The plain fasta, not the .gz: a gunzip that ran out of disk used
            # to satisfy the check through the old glob.
            ok=_exists(t) and not str(t).endswith(".gz"))
    if R.get("hhblits"):
        add(id="hhblits_db", label="HH-suite profile database", kind="db",
            ok=_prefix_exists(db.get("hhblits_db")),
            manual="hosting moves between releases",
            note="fetch a current hhsuite database (Pfam or UniRef30) and point "
                 "db.hhblits_db at its prefix")
    if R.get("structure"):
        t = db.get("foldseek_target", "")
        # afdb50.tar.gz is ~123 GB, not the ~1 TB this used to claim; the full
        # AFDB is ~490 GB. The old figure was the one number users read before
        # deciding whether to enable the structure stage at all.
        dbentry("foldseek_target", "Foldseek AFDB50", t, 123.0, lambda t: [
            f"mkdir -p {q(_dbdir(t))}",
            f"foldseek databases Alphafold/UniProt50 {q(t)} {q(_dbdir(t) + '/tmp')}"],
            disk_gb=200.0,
            ok=_exists(t + ".dbtype") and _exists(t + ".index"),
            note="hours, and the Ca index wants ~150 GB RAM; `foldseek "
                 "databases PDB <path> tmp` is ~2 GB and still finds most "
                 "classical toxin folds")
    if R.get("taxonomy"):
        t = db.get("ncbi_taxonomy", "")
        dbentry("ncbi_taxonomy", "NCBI taxdump", t, 0.08, lambda t: [
            f"mkdir -p {q(t)}",
            f"curl -fL -o {q(os.path.join(t, 'taxdump.tar.gz'))} "
            f"{q(src.get('taxdump',''))}",
            f"tar xzf {q(os.path.join(t, 'taxdump.tar.gz'))} -C {q(t)}"],
            ok=_exists(os.path.join(t, "nodes.dmp")) if t else False)
    if R.get("eggnog") and not pre:
        d = db.get("eggnog_data", "")
        dbentry("eggnog_data", "eggNOG data", d, 50.0, lambda d: [
            f"mkdir -p {q(d)}",
            f"download_eggnog_data.py -y --data_dir {q(d)}"],
            ok=_exists(os.path.join(d, "eggnog.db")) if d else False,
            note="not needed if emapper_precomputed is set")
    if R.get("diamond"):
        for tag, path in (db.get("diamond") or {}).items():
            url = (src.get("diamond") or {}).get(tag, "")
            if not path:
                add(id=f"diamond:{tag}", label=f"DIAMOND {tag}", kind="db",
                    ok=False, manual="no path configured",
                    note=f"set db.diamond.{tag}")
            elif url:
                stem = os.path.splitext(path)[0]
                # Staged next to the database, then removed: fixed /tmp names
                # collide with other users on a shared server and nothing ever
                # cleaned them up.
                src_f, fas = f"{stem}.src", f"{stem}.fas"
                add(id=f"diamond:{tag}", label=f"DIAMOND {tag}", kind="db",
                    ok=_exists(path), size_gb=0.3, cmds=[
                        f"mkdir -p {q(_dbdir(path))}",
                        f"curl -fL -o {q(src_f)} {q(url)}",
                        _reject_html(src_f),
                        f"gzip -dc {q(src_f)} > {q(fas)} || "
                        f"cp {q(src_f)} {q(fas)}",
                        f"diamond makedb --in {q(fas)} -d {q(stem)}",
                        f"rm -f {q(src_f)} {q(fas)}"])
            else:
                add(id=f"diamond:{tag}", label=f"DIAMOND {tag}", kind="db",
                    ok=_exists(path), manual="no URL in sources.diamond",
                    note="download the FASTA, then `diamond makedb --in <fas> "
                         f"-d {os.path.splitext(path)[0]}`")
    return out


def _pyhas(mod):
    import importlib.util as _u
    try:
        return _u.find_spec(mod) is not None
    except (ImportError, ValueError):
        return False


# ======================================================================
# commands
# ======================================================================
def cmd_init(args):
    if yaml is None:
        sys.exit("pyyaml is needed to write a config:  pip install pyyaml")
    if os.path.exists(args.out) and not args.force:
        sys.exit(f"{args.out} exists (use --force to overwrite)")
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("# metaannot config — every key is optional; defaults are "
                 "merged recursively.\n")
        yaml.safe_dump(DEFAULT_CONFIG, fh, sort_keys=False, default_flow_style=False)
    print(f"wrote {args.out}")


def cmd_doctor(args):
    cfg = load_config(args.config)
    p = Paths(cfg)
    ok = True

    if args.config and yaml is not None and os.path.exists(args.config):
        with open(args.config, encoding="utf-8") as fh:
            raw = yaml.safe_load(fh) or {}
        bad = unknown_keys(raw, DEFAULT_CONFIG) if isinstance(raw, dict) else []
        if bad:
            print("== config ==")
            for b in bad:
                if b in RETIRED_KEYS:
                    # A retired key is not a failure: the config predates a
                    # removal, the setting is inert, and doctor exiting
                    # non-zero over it would block a run that is otherwise
                    # correct. Say so and carry on.
                    print(f"  {'WARN':6s} '{b}' is no longer a setting - "
                          f"{RETIRED_KEYS[b]}")
                    continue
                ok = False
                print(f"  {'MISS':6s} unrecognised key '{b}'"
                      + (nearest_config_key(b) or ".")
                      + " It is being ignored, so this setting is NOT in "
                        "effect.")

    print("== inputs ==")
    for label, path in [("proteins_faa", cfg["proteins_faa"]),
                        ("quant_table", cfg["quant_table"]),
                        ("gff", cfg.get("gff") or "")]:
        if not path:
            print(f"  {'-':6s} {label:16s} (not set)")
            continue
        good = os.path.exists(path)
        ok &= good or label == "gff"
        print(f"  {'OK' if good else 'MISS':6s} {label:16s} {path}")

    pre = cfg.get("emapper_precomputed") or ""
    pre = [pre] if isinstance(pre, str) and pre else list(pre or [])
    if pre:
        print("== precomputed emapper ==")
        for path in pre:
            good = os.path.exists(path)
            ok &= good
            print(f"  {'OK' if good else 'MISS':6s} {path}")

    reqs = requirements(cfg, p)     # evaluated once; every view below uses it
    tools = [r for r in reqs if r["kind"] == "tool"]
    dbs = [r for r in reqs if r["kind"] == "db"]
    for title, group in (("tools", tools), ("databases", dbs)):
        if not group:
            continue
        print(f"== {title} ==")
        for r in group:
            if r["ok"]:
                mark = "OK"
            elif r["manual"]:
                mark = "MANUAL"
            else:
                mark = "MISS"
            # A MANUAL item used to count as satisfied, so doctor said "all
            # checks passed" with no InterProScan and no SignalP 6 and the run
            # died at that stage hours later. `manual` says doctor will not
            # fetch it, not that the stage can do without it.
            ok &= r["ok"]
            size = ""
            if not r["ok"] and r["size_gb"] >= 0.05:
                size = f"  ~{_gb(r['size_gb'])} GB"
                if r["disk_gb"] > r["size_gb"]:
                    size += f" ({_gb(r['disk_gb'])} GB on disk)"
            print(f"  {mark:6s} {r['label']}{size}")
            if not r["ok"]:
                if r["manual"]:
                    print(f"         manual: {r['manual']}")
                if r["note"]:
                    print(f"         {r['note']}")
                for c in r["cmds"]:
                    print(f"         $ {c}")

    if cfg["run"].get("diamond"):
        unweighted = set((cfg["db"].get("diamond") or {})) - set(cfg.get("diamond_weights") or {})
        if unweighted:
            print(f"  {'WARN':6s} diamond_weights missing for "
                  f"{sorted(unweighted)}; those hits still count as annotation "
                  "but score 0 in the effector ranking")
        # Existence is not usability. A failed `diamond makedb` leaves a
        # zero-byte .dmnd that passes every check above, and a database of
        # 15-residue peptides passes them all while being unable to reach the
        # configured e-value. Both return 0 hits, and 0 hits is also what a
        # real absence looks like, so doctor is the last place either can be
        # caught before hours of searching say nothing.
        for tag, path in sorted((cfg["db"].get("diamond") or {}).items()):
            if not path or not os.path.exists(path):
                continue          # already reported as MISS above
            try:
                bad, warn = diamond_db_check(cfg, tag, path)
            except StageError as e:      # a non-numeric diamond_evalues entry
                ok, bad, warn = False, f"{tag}: {e}", None
            if bad:
                ok = False
                print(f"  {'MISS':6s} {bad}")
            elif warn:
                print(f"  {'WARN':6s} {warn}")

    if cfg["run"].get("structure") or cfg["run"].get("topology"):
        # Without this, a machine with no usable GPU passed every check and the
        # user learned the truth hours later, when esmfold finally ran and
        # died. The GPU is the one requirement doctor could not see.
        print("== gpu ==")
        usable, why = cuda_probe()
        print(f"  {'OK' if usable else 'WARN':6s} {why}")
        if not usable:
            if cfg["run"].get("structure"):
                print(f"  {'MISS':6s} run.structure needs CUDA: stage_esmfold "
                      "exits rather than fold on CPU, which is impractical at "
                      "any real scale. Set run.structure: false, or fold on a "
                      "GPU host and copy results/structures/ back - foldseek "
                      "itself is CPU-only and will search whatever models are "
                      "there.")
                ok = False
            if cfg["run"].get("topology"):
                # Not a MISS: signalp is CPU-only and useful on its own, and
                # tmbed does run without a GPU - just not at this size.
                print(f"  {'WARN':6s} run.topology: SignalP 6 is CPU-only and "
                      "unaffected. tmbed will fall back to CPU, where it is "
                      "one to two orders of magnitude slower - fine for a few "
                      "thousand proteins, not for a few hundred thousand. Set "
                      "tmbed_use_gpu: false to say so deliberately, or true to "
                      "make a missing GPU fatal instead of slow.")
        elif cfg["run"].get("structure") and cfg["run"].get("topology"):
            print(f"  {'OK':6s} esmfold and tmbed share gpu_device "
                  f"{cfg['gpu_device']}, so at most "
                  f"{cfg.get('gpu_workers', 1)} of them runs at a time")

    if cfg["quant_format"] == "fragpipe_tmt":
        # Without this block doctor says nothing at all about a TMT run: the
        # only TMT line it printed lived inside `== manifest ==`, and a TMT
        # config normally sets no manifest. So the one layout this format is
        # most particular about - a run directory of per-plex folders, each
        # with its own annotation file - went unchecked until the run itself
        # died on it.
        print("== tmt ==")
        root = cfg.get("quant_table") or ""
        t = cfg.get("tmt") or {}
        lvl = str(t.get("level") or "ion").lower()
        fname = TMT_LEVEL_FILES.get(lvl, "ion.tsv")
        ref_name = str(t.get("reference_name") or "")
        ref_chan = str(t.get("reference_channel") or "")
        if ref_name and ref_chan:
            ok = False
            print(f"  {'MISS':6s} tmt.reference_name ('{ref_name}') and "
                  f"tmt.reference_channel ('{ref_chan}') are both set and can "
                  "disagree per plex; set one")
        if os.path.isfile(root):
            ok = False
            print(f"  {'MISS':6s} quant_table is a file: {root}. This format "
                  "reads the run DIRECTORY holding the per-plex folders, "
                  "because a tmt-report matrix has already collapsed the "
                  "peptides this tool needs")
        elif not os.path.isdir(root):
            ok = False
            print(f"  {'MISS':6s} quant_table not found: {root}")
        else:
            try:
                plexes = tmt_plex_dirs(root, cfg)
            except StageError as e:
                ok, plexes = False, []
                print(f"  {'MISS':6s} {e}")
            if plexes:
                print(f"  {'OK':6s} {len(plexes)} plex(es): "
                      f"{[n for n, _ in plexes]}")
            # Every rule below is read_fragpipe_tmt's own - which channels it
            # keeps, how it resolves the reference, which names it refuses -
            # because a verdict that disagrees with the reader is worse than
            # no verdict: the run doctor blessed still dies, on the same
            # config, hours later.
            drop_empty = bool(t.get("drop_empty_channels", True))
            sizes, seen_names, ref_hits = {}, {}, {}
            for plex, pdir in plexes:
                lvl_path = os.path.join(pdir, fname)
                if not os.path.exists(lvl_path):
                    ok = False
                    print(f"  {'MISS':6s} {plex}: no {fname} (tmt.level "
                          f"'{lvl}'); {pdir} holds "
                          f"{sorted(os.listdir(pdir))[:8]}")
                try:
                    apath = tmt_annotation_path(plex, pdir, cfg)
                    rows = read_tmt_annotation(apath, plex)
                except StageError as e:
                    ok = False
                    print(f"  {'MISS':6s} {plex}: {e}")
                    continue
                sizes[plex] = len(rows)
                # The reader's keep_samples: an unassigned <PLEX>_<CHANNEL>
                # placeholder is not a sample when drop_empty_channels is on.
                keep = [(c, s) for c, s in rows
                        if not (drop_empty and s == f"{plex}_{c}")]
                if ref_name:
                    hits = [s for _c, s in keep
                            if fnmatch.fnmatchcase(s, ref_name)]
                elif ref_chan:
                    # fnmatchcase on the CHANNEL, because the reader and the
                    # README both make reference_channel a glob ('131*'),
                    # while this compared it to the channel as a literal and
                    # so failed a config the run accepts. fnmatchcase, not
                    # fnmatch, because it is what the reader uses: fnmatch
                    # normalises through os.path.normcase, which lowercases on
                    # Windows only, so fnmatch would make doctor's verdict
                    # differ from the run's by platform. '131c' would pass
                    # here and be refused by the reader.
                    hits = [s for c, s in keep
                            if fnmatch.fnmatchcase(c, ref_chan)]
                else:
                    hits = []
                if ref_name or ref_chan:
                    ref_hits[plex] = hits
                    if len(hits) == 1:
                        # Under both treatments the reference stops being a
                        # sample column before the reader's collision check,
                        # so a bridge carrying one name in every plex is the
                        # design and must not be reported as a duplicate.
                        keep = [(c, s) for c, s in keep if s != hits[0]]
                for _c, s in keep:
                    seen_names.setdefault(s, []).append(plex)
            # "every plex" can only mean the ones doctor got as far as
            # reading: a plex whose annotation failed above is in none of
            # these counts and must not be summarised as though it had passed.
            scope = ("every plex" if len(sizes) == len(plexes)
                     else f"each of the {len(sizes)} plex(es) read")
            if sizes:
                dist = sorted(set(sizes.values()))
                if len(dist) > 1:
                    # Not fatal - the reader handles ragged plexes - but it is
                    # the kind of thing that is a typo far more often than it
                    # is the design.
                    print(f"  {'WARN':6s} plexes differ in channel count "
                          f"{dist}: "
                          f"{ {k: v for k, v in sorted(sizes.items())} }")
                else:
                    print(f"  {'OK':6s} {dist[0]} channels in {scope}, "
                          f"{sum(sizes.values())} in total")
            dupes = {n: pl for n, pl in seen_names.items() if len(pl) > 1}
            if dupes:
                # This was a WARN promising the plexes would be "treated as
                # one sample measured in each". They are not: the reader dies
                # on the second plex to claim a name, so doctor was exiting 0
                # on a config the run refuses outright.
                ok = False
                shown = dict(sorted(dupes.items())[:4])
                print(f"  {'MISS':6s} {len(dupes)} sample name(s) appear in "
                      f"more than one plex: {shown}. Sample names are the "
                      "columns of the joined matrix, so the reader refuses "
                      "two plexes that claim one; rename them per plex, or, "
                      "if this is a bridge, name it with tmt.reference_name "
                      "so it stops being a sample")
            if ref_name or ref_chan:
                which = f"reference_name '{ref_name}'" if ref_name                     else f"reference_channel '{ref_chan}'"
                missing = sorted(p for p, h in ref_hits.items() if not h)
                ambig = {p: h for p, h in sorted(ref_hits.items())
                         if len(h) > 1}
                if missing:
                    ok = False
                    print(f"  {'MISS':6s} tmt.{which} matches nothing in "
                          f"{len(missing)} plex(es): {missing[:6]}. A "
                          "reference absent from a plex leaves that plex "
                          "without a denominator")
                if ambig:
                    # The reader takes exactly one reference per plex and
                    # dies on anything else, so a pattern matching two
                    # channels is a failure to report, not a resolution.
                    ok = False
                    print(f"  {'MISS':6s} tmt.{which} matches more than one "
                          f"channel in {len(ambig)} plex(es): "
                          f"{dict(list(ambig.items())[:4])}. The reference is "
                          "one channel per plex; narrow the pattern until it "
                          "names it")
                if ref_hits and not missing and not ambig:
                    # Guarded on ref_hits: with no plex read there is nothing
                    # the pattern resolved in, and this line was printed
                    # anyway - an OK about a run doctor never opened.
                    print(f"  {'OK':6s} tmt.{which} resolves in {scope}")
            elif bool(t.get("use_reference_ratios", False)):
                ok = False
                print(f"  {'MISS':6s} tmt.use_reference_ratios is on but "
                      "neither tmt.reference_name nor tmt.reference_channel "
                      "is set")
            else:
                print(f"  {'WARN':6s} no reference channel named, so plexes "
                      "are compared on within-plex normalised intensity "
                      "alone; set tmt.reference_name if this design has a "
                      "bridge")

    if cfg.get("manifest"):
        print("== manifest ==")
        if not os.path.exists(cfg["manifest"]):
            ok = False
            print(f"  {'MISS':6s} {cfg['manifest']} (config key `manifest`)")
        else:
            # doctor is the one command that has to survive every other
            # failure: a manifest read_manifest rejects must not take the
            # tools, databases, resources and R blocks down with it.
            try:
                m = read_manifest(cfg["manifest"])
            except StageError as e:                         # noqa: PERF203
                ok, m = False, None
                print(f"  {'MISS':6s} {e}")
            if m is not None:
                print(f"  {'OK':6s} {len(m)} runs, groups: "
                      f"{sorted(m['experiment'].unique())}")
            if m is not None and cfg["quant_format"] == "fragpipe_tmt":
                print(f"  {'WARN':6s} quant_format is 'fragpipe_tmt', so the "
                      "manifest is not used: sample names come from each "
                      "plex's annotation file, and a TMT manifest's "
                      "experiment column is the plex, not a condition")
            elif m is not None and os.path.isfile(cfg["quant_table"]):
                try:
                    head = pd.read_csv(cfg["quant_table"], sep=None,
                                       engine="python", nrows=0, encoding="utf-8", encoding_errors="replace")
                    sfx = (" Intensity" if cfg["quant_format"].startswith("fragpipe")
                           else "")
                    cand = [c for c in head.columns
                            if (c.endswith(sfx) if sfx else True)
                            and not c.endswith(("MaxLFQ Intensity",
                                                "Spectral Count"))]
                    _, miss_rows, miss_cols = map_manifest_to_columns(
                        m, cand, sfx)
                    if miss_rows:
                        ok = False
                        print(f"  {'MISS':6s} {len(miss_rows)} manifest "
                              "run(s) match no column in "
                              f"{cfg['quant_table']}: {miss_rows[:5]}. Fix the "
                              "run names in the manifest, or unset `manifest:`")
                    else:
                        print(f"  {'OK':6s} every run maps to a quant column")
                    if miss_cols:
                        print(f"  {'WARN':6s} {len(miss_cols)} quant column(s) "
                              "absent from the manifest and will be dropped: "
                              f"{miss_cols[:4]}")
                except Exception as e:                      # noqa: BLE001
                    print(f"  {'WARN':6s} could not pre-check the mapping: {e}")

    if cfg["run"].get("unipept") or cfg["run"].get("taxonomy"):
        print("== taxonomy ==")
        u = cfg.get("unipept") or {}
        # The taxonomy stage consumes the unipept stage's output, so with
        # run.unipept off the honest answer is "turn it on", not "MISS unipept
        # cache".
        if cfg["run"].get("taxonomy") and not cfg["run"].get("unipept"):
            ok = False
            print(f"  {'MISS':6s} run.taxonomy needs run.unipept (it "
                  "compares the eggNOG lineage against Unipept's); set "
                  "run.unipept: true, or run.taxonomy: false")
        if u.get("result"):
            good = os.path.exists(u["result"])
            ok &= good
            print(f"  {'OK' if good else 'MISS':6s} unipept.result   "
                  f"{u['result']}")
        elif u.get("allow_http"):
            print(f"  {'WARN':6s} unipept.allow_http is on; the API version, "
                  "field names and rate limits are outside this tool's control")
        else:
            cached = os.path.exists(p.unipept_cache)
            ok &= cached
            print(f"  {'OK' if cached else 'MISS':6s} unipept cache    "
                  f"{p.unipept_cache}  (set unipept.result to an existing "
                  "pept2lca export, or unipept.allow_http: true)")
        if cfg["quant_format"] not in FEATURE_FORMATS and not u.get("result"):
            ok = False
            print(f"  {'MISS':6s} quant_format is '{cfg['quant_format']}', "
                  "which is protein level; Unipept needs peptides. Use a "
                  "peptide-level quant_format "
                  f"({', '.join(sorted(FEATURE_FORMATS))}), or point "
                  "unipept.result at an existing pept2lca export")
        src = cfg.get("taxonomy_source", "eggnog")
        if src != "eggnog" and not cfg["run"].get("taxonomy"):
            ok = False
            print(f"  {'MISS':6s} taxonomy_source is '{src}' but "
                  "run.taxonomy is off; set run.taxonomy: true, or "
                  "taxonomy_source: eggnog")

    w = max(1, int(cfg.get("stage_workers", 4)))
    # parse_ram, not int(): `ram_gb: 64G` is the form the --ram help text
    # advertises, and int() turned it into an unhandled ValueError traceback
    # halfway through doctor's output. A budget doctor cannot read is a
    # problem, but not one worth losing the R block over.
    ram, ram_err = 0, ""
    try:
        ram = parse_ram(getattr(args, "ram", None) or cfg.get("ram_gb"))
    except StageError as e:
        ok, ram_err = False, str(e)
    ram = ram or int(detect_ram_gb() * 0.8)
    print("== resources ==")
    if ram_err:
        print(f"  {'MISS':6s} ram_gb: {ram_err}")
    print(f"  {'OK':6s} {cfg['threads']} cpu, {ram or 'unknown'} GB budget, "
          f"up to {w} stage(s) at once")
    print(f"  {'OK' if ram else 'WARN':6s} per stage: "
          f"{max(1, int(cfg['threads']) // w)} cpu"
          + (f", {max(1, ram // w)} GB" if ram else ", memory budget unknown"))
    if ram:
        per = max(1, ram // w)
        need = int(cfg.get("emapper_dbmem_min_gb", 64))
        if cfg["run"].get("eggnog") and not (cfg.get("emapper_precomputed") or ""):
            print(f"  {'OK' if per >= need else 'WARN':6s} eggNOG --dbmem "
                  f"{'enabled' if per >= need else f'off ({per} < {need} GB)'}"
                  + ("" if per >= need else "; raise ram_gb or lower "
                     "stage_workers, or set emapper_dbmem_min_gb"))
        per = max(1, per)
        if per < 4:
            print(f"  {'WARN':6s} {per} GB per stage is tight; lower "
                  "stage_workers or raise ram_gb")
    else:
        # doctor has no --ram flag, so telling the user to pass one was advice
        # they could not follow; ram_gb in the config is what run reads too.
        print(f"  {'WARN':6s} set ram_gb: in the config (e.g. 64G) so "
              "memory-sensitive flags can be set (diamond -b, mmseqs/foldseek "
              "--split-memory-limit, hhblits -maxmem, InterProScan heap)")

    if have("Rscript"):
        # Ask for exactly what the generated Rmd and build_object.R load. The
        # old four-package probe let an environment pass doctor and then die
        # at the report's first `library(readr)`, after the whole pipeline had
        # run. Keep this list in step with RMD_TEMPLATE's setup chunk.
        RNEED = ["readr", "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
                 "purrr", "limma", "knitr", "rmarkdown", "SummarizedExperiment"]
        ROPT = ["QFeatures", "patchwork", "clusterProfiler"]
        # Bioconductor packages need BiocManager; everything else is CRAN.
        BIOC = {"limma", "SummarizedExperiment", "QFeatures", "clusterProfiler"}
        want = RNEED + ROPT
        chk = ('cat(paste(vapply(c(%s), function(p) paste0(p, "=", '
               'requireNamespace(p, quietly=TRUE)), character(1)), collapse=" "))'
               % ",".join(f'"{n}"' for n in want))
        print("== R ==")
        try:
            r = subprocess.run([resolve_tool("Rscript"), "-e", chk], capture_output=True,
                               text=True, timeout=180)
            seen = {}
            for tok in (r.stdout or "").split():
                name, _, val = tok.partition("=")
                seen[name] = val == "TRUE"
            if not seen:
                # An Rscript that answers nothing is not a green R stack.
                ok = False
                print(f"  {'MISS':6s} the R package query returned nothing "
                      f"(rc={r.returncode}): {(r.stderr or '').strip()[:200]}")
            for name in want:
                if name not in seen:
                    continue
                good, need = seen[name], name in RNEED
                how = (f'BiocManager::install("{name}")' if name in BIOC
                       else f'install.packages("{name}")')
                # A missing required package is a problem, not a note: the
                # verdict and the exit code used to ignore this whole block.
                if need and not good:
                    ok = False
                print(f"  {'OK' if good else ('MISS' if need else 'WARN'):6s} "
                      f"{name}" + ("" if good else "  " + how))
        except Exception as e:                              # noqa: BLE001
            ok = False
            print(f"  {'MISS':6s} could not query R packages: {str(e)[:160]}")
        if not have("pandoc"):
            # rmarkdown shells out to pandoc; RStudio bundles one, a bare
            # R install does not.
            print(f"  {'WARN':6s} pandoc not on PATH; rmarkdown::render "
                  "will fail unless R finds its own copy (RSTUDIO_PANDOC)")
    else:
        print("== R ==\n" + f"  {'WARN':6s} Rscript not found; the report "
              "and the R object cannot be built here (run.report/run.object "
              "still write the scripts, so they can be knitted elsewhere)")

    missing = [r for r in reqs if not r["ok"] and not r["manual"]]
    manual = [r for r in reqs if not r["ok"] and r["manual"]]
    total = sum(r["size_gb"] for r in missing)
    total_disk = sum(max(r["disk_gb"], r["size_gb"]) for r in missing)

    if args.install_plan or (args.fix and missing):
        # ASCII only: this file is read by a shell, and on a non-UTF-8 console
        # the em dash arrived as a replacement glyph inside the comment.
        lines = ["#!/usr/bin/env bash",
                 "# generated by metaannot doctor -- review before running",
                 "set -euo pipefail", ""]
        if not missing:
            lines.append("# nothing to install; every requirement of this "
                         "config is already satisfied")
        for r in missing:
            lines += [f"# {r['label']}"
                      + (f"  (~{_gb(r['size_gb'])} GB)" if r["size_gb"] >= 0.05
                         else "")]
            if r["note"]:
                lines.append(f"#   {r['note']}")
            lines += r["cmds"] + [""]
        script = "\n".join(lines) + "\n"
        if args.install_plan:
            # Always write the file that was asked for, even when it is empty:
            # `doctor --install-plan i.sh && bash i.sh` used to fail with "No
            # such file" on a machine that had everything.
            with open(args.install_plan, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(script)
            os.chmod(args.install_plan, 0o755)
            print(f"\ninstall plan -> {args.install_plan}  "
                  f"({len(missing)} item(s), ~{_gb(total)} GB)")
            print("review it, then run it, or rerun doctor with --fix"
                  if missing else "nothing to install")

    if args.fix and os.name == "nt":
        # The commands come out of `requirements()` as POSIX shell -- mkdir
        # -p, curl, tar, gunzip, hmmpress -- and are run through
        # subprocess(shell=True), which on Windows is cmd.exe. `mkdir -p
        # 'C:\\db'` there creates a directory called -p; curl and tar may or
        # may not exist; and every one of them can exit 0 having done nothing,
        # which is precisely the failure the post-install verification was
        # written to catch. Refuse rather than half-work.
        die("doctor --fix cannot run on Windows: the install commands it "
            "generates are POSIX shell, and cmd.exe silently mis-executes "
            "them (`mkdir -p C:\\db` makes a directory called -p).\n"
            "  Write the plan and run it where the tools live:\n"
            "    python metaannot.py doctor --config <cfg> "
            "--install-plan install.sh\n"
            "    wsl bash install.sh          # or Git Bash, or the Linux "
            "host that will do the run\n"
            "  Every other doctor check works here; only --fix is refused.")

    if args.fix:
        if not missing:
            print("\nnothing to install")
        else:
            print(f"\n{len(missing)} item(s) to install, ~{_gb(total)} GB to "
                  f"download (~{_gb(total_disk)} GB on disk):")
            for r in missing:
                print(f"  - {r['label']}"
                      + (f"  ~{_gb(r['size_gb'])} GB" if r["size_gb"] >= 0.05
                         else ""))
            if not args.yes:
                try:
                    reply = input("\nproceed? [y/N] ").strip().lower()
                except EOFError:
                    reply = ""
                if reply not in ("y", "yes"):
                    print("aborted; nothing was downloaded")
                    return 1
            done, failed = [], []
            for r in missing:
                log(f"installing {r['label']}")
                try:
                    for c in r["cmds"]:
                        log(f"$ {c}")
                        rc = subprocess.run(c, shell=True).returncode
                        if rc != 0:
                            raise RuntimeError(f"command exited {rc}: {c}")
                except Exception as e:                      # noqa: BLE001
                    failed.append((r["label"], str(e)[:160]))
                    continue
                done.append(r["label"])
            # Verify rather than trust the exit codes: a download that wrote a
            # 404 page exits 0.
            after = {r["id"]: r["ok"] for r in requirements(cfg, p)}
            unverified = [r["label"] for r in missing
                          if not after.get(r["id"]) and r["label"] in done]
            print(f"\ninstalled {len(done)}, failed {len(failed)}")
            for lab, err in failed:
                print(f"  FAILED  {lab}: {err}")
            for lab in unverified:
                print(f"  UNVERIFIED  {lab}: commands succeeded but the file is "
                      "still not what the config expects (missing, empty, or "
                      "not hmmpress'd) -- check the URL in `sources` and the "
                      "path in `db`")
            if failed or unverified:
                return 1
            print("rerun `doctor` to confirm")

    if manual:
        print(f"\n{len(manual)} item(s) need a manual download "
              f"({', '.join(r['label'] for r in manual)}); doctor will not "
              "attempt those. Install them or turn the stage off -- an "
              "enabled stage with a missing licence-gated tool dies at run "
              "time, hours in.")

    print("\n" + ("all checks passed" if ok else
                  "problems found -- fix these before a long run"))
    if missing and not (args.fix or args.install_plan):
        print(f"  `doctor --install-plan install.sh` writes the commands "
              f"(~{_gb(total)} GB), `doctor --fix` runs them")
    return 0 if ok else 1


def cmd_subset(args):
    if args.format in FEATURE_FORMATS:
        cfg = json.loads(json.dumps(DEFAULT_CONFIG))
        feats, _, _ = read_feature_table(args.quant, args.format, cfg)
        want = set()
        for cands in feats["candidates"]:
            want.update(c for c in cands if c)
        log(f"{len(want)} distinct protein ids across {len(feats)} features")
    elif args.format == "ids":
        with opener(args.quant) as fh:
            want = {l.strip() for l in fh if l.strip()}
    else:
        q = pd.read_csv(args.quant, sep="\t", low_memory=False, encoding="utf-8", encoding_errors="replace")
        if args.format == "diann":
            cands = ["Protein.Group", "Protein.Ids", "Protein.Names"]
        else:
            cands = ["Protein", "Protein ID", "Indistinguishable Proteins"]
        cols = [c for c in cands if c in q.columns]
        if not cols:
            die(f"none of {cands} found in {args.quant}")
        log(f"reading ids from: {cols}")
        want = set()
        for c in cols:
            for v in q[c].dropna().astype(str):
                for pid in re.split(r"[;,]", v):
                    pid = pid.strip()
                    if pid:
                        want.add(pid)
    if args.format not in FEATURE_FORMATS:
        log(f"{len(want)} distinct protein ids requested")

    found, keep = set(), False
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    # utf-8 explicitly: a description line carrying a non-ASCII character kills
    # the write under a cp1252/C locale, after the file has been truncated.
    with opener(args.db) as fh, open(args.out, "w", encoding="utf-8",
                                     newline="\n") as out:
        for line in fh:
            if line.startswith(">"):
                tok = line[1:].split()
                pid = tok[0] if tok else ""
                keep = pid in want and pid not in found
                if keep:
                    found.add(pid)
            if keep:
                out.write(line)
    log(f"wrote {len(found)} sequences -> {args.out}")
    missing = want - found
    if missing:
        log(f"{len(missing)} ids absent from the database fasta, e.g. "
            f"{sorted(missing)[:5]} -- usually a decoy/contaminant prefix "
            "(rev_, sp|, CON__) added by the search engine", "WARN")


def cmd_run(args):
    global _LOGFH
    cfg = load_config(args.config)
    if args.threads is not None:
        if int(args.threads) < 1:
            die("--threads must be at least 1")
        cfg["threads"] = args.threads
    if getattr(args, "ram", None) is not None:
        cfg["ram_gb"] = parse_ram(args.ram)
    try:
        set_progress_interval(cfg["progress_interval_s"])
    except (TypeError, ValueError):
        die(f"progress_interval_s must be a number of seconds (0 disables), "
            f"not {cfg['progress_interval_s']!r}")
    # Absolute, as load_config already makes the config's own paths: the cache
    # must not see 'results' and './results' as two different projects.
    if args.faa:
        cfg["proteins_faa"] = os.path.abspath(args.faa)
    if args.results_dir:
        cfg["results_dir"] = os.path.abspath(args.results_dir)

    p = Paths(cfg)
    if not args.dry_run:
        # Not before the dry-run branch: a plan check should not leave fifteen
        # new directories behind for the next person to wonder about.
        p.mkdirs()
        _LOGFH = open(p.logfile, "a", encoding="utf-8")
        log(f"metaannot {__version__} starting")
        lock = ResultsLock(p.lock, force=getattr(args, "force_unlock", False))
        lock.__enter__()
        # The BOUND METHOD, not a closure over `lock`. Both of these used to
        # be reachable from a name that this function reassigns further down,
        # and the signal handler duly called threading.Lock.__exit__ and died
        # with "release unlocked lock" - which took the run out with an
        # uncaught RuntimeError and exit 1 rather than releasing anything.
        release_results_lock = lock.__exit__
        atexit.register(release_results_lock)

        def _release_lock_on_signal(sig, _frame):
            """Release the results lock when the run is killed, not only when
            it exits.

            atexit does not run on SIGTERM or SIGHUP - Python's default
            handler terminates the process outright - so a run stopped by
            `kill`, by a scheduler hitting its time limit, or by a closing ssh
            session left a lock file behind naming a pid that no longer
            exists. On the same host the next run can prove that and reclaim
            it; from another node of a cluster it cannot, and the resume
            became a stale-lock refusal needing --force-unlock.

            What this does NOT do: stop the tools already running. TMbed,
            InterProScan and DIAMOND are separate processes that outlive us,
            and the state file is what records which stages were mid-flight.
            It releases the lock, flushes the log, and exits 128+N so a
            wrapper script still sees a killed process rather than a clean
            one.

            Windows delivers almost none of this: subprocess.terminate() is
            TerminateProcess, which runs no handler at all. SIGBREAK
            (Ctrl-Break) is the one that does arrive, so it is registered too.
            """
            # NOTHING in here may take a lock. A Python signal handler runs
            # IN THE MAIN THREAD, between two bytecodes of whatever that
            # thread was doing, so any lock the interrupted frame is holding
            # is still held while this runs and is NOT reentrant.
            #
            # The first version called log(), which goes through sys.stderr,
            # whose buffer lock is exactly such a lock. On a run that emits a
            # progress line per stage per minute the signal eventually lands
            # mid-write, the handler blocks forever on a lock its own frame
            # holds, and the process HANGS instead of releasing the lock --
            # strictly worse than the stale lock this exists to prevent. CI
            # caught it on one job of seven; it is a race, not a certainty,
            # which is the worst kind.
            #
            # So: os.remove (a raw syscall, inside ResultsLock.__exit__) and
            # os.write to fd 2, with the message pre-formatted and pre-encoded
            # at registration time. No formatting, no buffered I/O, no locks.
            # The cost is that the final line reaches stderr but not the log
            # FILE, whose buffer cannot be safely touched from here.
            release_results_lock()
            try:
                os.write(2, _sig_msgs.get(int(sig), b"\nstopping on a signal; "
                                          b"results lock released\n"))
            except OSError:
                pass
            # os._exit, not sys.exit: SystemExit here would unwind through the
            # stage pool's `with`, which WAITS for its workers, and a tmbed
            # chunk can be an hour. A kill has to mean now.
            os._exit(128 + int(sig))

        _sig_msgs = {}
        for _name in ("SIGTERM", "SIGHUP", "SIGBREAK"):
            _sig = getattr(signal, _name, None)
            if _sig is None:
                continue
            _sig_msgs[int(_sig)] = (
                f"\nWARN  stopping on {_name}: releasing the results lock "
                f"{p.lock}. Any tool already running is a separate process "
                "and is not stopped by this, so its output may be "
                f"incomplete; {p.state} records which stages were running.\n"
            ).encode("utf-8", "replace")
            # ValueError when this is not the main thread, OSError when the
            # platform refuses the signal. Neither is worth failing a run
            # over: the lock still comes off on a normal exit.
            with contextlib.suppress(ValueError, OSError, AttributeError):
                signal.signal(_sig, _release_lock_on_signal)

    for name in (args.only or []) + ([args.from_stage] if args.from_stage else []):
        if name not in STAGE_NAMES:
            die(f"unknown stage '{name}'; choose from {STAGE_NAMES}")
    if args.only and args.from_stage:
        die("--only and --from cannot be combined: --only would silently win "
            "over --from. Pick one.")

    if not os.path.exists(cfg["proteins_faa"]):
        die(f"proteins_faa not found: {cfg['proteins_faa']}")

    selected = STAGE_NAMES
    if args.from_stage:
        selected = STAGE_NAMES[STAGE_NAMES.index(args.from_stage):]
    if args.only:
        selected = [n for n in STAGE_NAMES if n in set(args.only)]
    only_set = set(args.only or [])
    by_name = {st["name"]: st for st in STAGES}

    # --force discards what the user asked to redo, and only that. Wiping the
    # whole file made every unselected stage look as if it had never been
    # recorded here, so the next run adopted its stale output without ever
    # comparing signatures — and the config change that prompted the redo was
    # then ignored for good.
    state = load_state(p.state)
    if args.force:
        for n in selected:
            state.pop(n, None)

    def unmet_deps(st):
        """Enabled dependencies that were excluded by --only/--from and have
        produced nothing. A disabled dependency is fine — its evidence is
        legitimately absent — but one that was merely not selected means the
        stage would run against inputs that do not exist yet, and several
        stages will happily produce a confident, empty answer from that."""
        # Only the dependency outputs this stage actually reads count. esmfold
        # depends on integrate for dark.faa alone, and demanding integrate's
        # other outputs (annotation_pass1.tsv, dark_all.faa) refused the
        # documented GPU hand-off, where only dark.faa is copied to the
        # laptop. hhblits and jackhmmer read dark_all.faa in the same way, so
        # a machine that runs only those needs that file copied instead.
        needed = set(x for x in st["inp"](cfg, p) if x)
        bad = []
        for d in st["deps"]:
            dep = by_name[d]
            if d in selected:
                continue
            if dep["enabled"] and not cfg["run"].get(dep["enabled"], False):
                continue                      # disabled on purpose
            outs = [o for o in dep["out"](p) if o in needed] or dep["out"](p)
            if not exists_all(outs):
                bad.append(d)
        return bad

    def decide(st):
        """Cached / adopted / RUN, for the stage as things stand right now."""
        name = st["name"]
        if name not in selected:
            return "not selected"
        if st["enabled"] and not cfg["run"].get(st["enabled"], False):
            # Naming a stage with --only is a clearer statement of intent than
            # a run flag left off for another machine: it is how the GPU box
            # runs esmfold/tmbed against the server's config.
            if name not in only_set:
                return f"disabled (run.{st['enabled']})"
            log(f"{name}: run.{st['enabled']} is false, but the stage was "
                "named with --only, so it is running anyway", "WARN")
        if args.force:
            return "RUN"
        outs = st["out"](p)
        prev = state.get(name)
        if prev and prev.get("status") == "running":
            # finish() never ran, so the writer was killed (OOM, Ctrl-C, a
            # dropped ssh). What is on disk is as likely to be half a file as
            # a whole one, and nothing here can tell the difference.
            log(f"{name}: the previous run was interrupted while this stage "
                "was writing, so its output may be truncated; recomputing",
                "WARN")
            return "RUN"
        if prev and prev.get("signature") == signature(st, cfg, p) \
                and exists_all(outs):
            if name in only_set:
                log(f"{name}: cached — nothing it reads has changed. Use "
                    f"--force --only {name} to recompute it anyway.")
            return "cached"
        if args.no_adopt or not exists_all(outs):
            return "RUN"
        # Outputs present but not recorded as produced by a finished run here:
        # made elsewhere — the GPU box folding structures, an hmmsearch run on
        # the cluster. Adopt them instead of recomputing, unless told
        # otherwise.
        if prev is not None and prev.get("status") != "failed":
            return "RUN"
        empty = [o for o in outs
                 if os.path.isfile(o) and os.path.getsize(o) == 0]
        if empty and not st.get("empty_ok"):
            log(f"{name}: not adopting {', '.join(empty)} — the file "
                "is empty and there is no record of producing it, which "
                "usually means an interrupted writer. Rerunning.", "WARN")
            return "RUN"
        if prev is None:
            if state.unreadable:
                log(f"{name}: outputs are here but the state file could not be "
                    "read, so there is no record of what made them; "
                    "recomputing instead of adopting them", "WARN")
                return "RUN"
            return "adopt"
        # A failed record used to block adoption for ever, which broke the
        # rsync-back half of the two-machine workflow: this box tried the
        # stage and died, the GPU box produced the real thing. Only files
        # written after that failure can be that real thing.
        when = _stamp_epoch(prev.get("finished"))
        if when is not None and all(_mtime(o) > when for o in outs):
            log(f"{name}: the last attempt here failed, but every output is "
                "newer than that failure, so it came from somewhere else; "
                "adopting it", "WARN")
            return "adopt"
        return "RUN"

    if args.dry_run:
        print(f"{'stage':12s} {'action':28s} outputs")
        for st in STAGES:
            action = decide(st)
            if action == "RUN":
                # The commonest refusal of a --only/--from run; showing RUN
                # here and dying at dispatch helps nobody.
                unmet = unmet_deps(st)
                if unmet:
                    action = f"refused: needs {' '.join(unmet)}"
            print(f"{st['name']:12s} {action:28s} {', '.join(st['out'](p))}")
        print("\nnote: a dry run evaluates every stage against the files as they "
              "are now,\nso a stage shown as cached may still rerun once an "
              "upstream stage rewrites its input.")
        return 0

    # ---- schedule over the DAG ------------------------------------------
    workers = 1 if args.serial else max(1, int(cfg.get("stage_workers", 4)))
    try:
        gpu_slots = max(1, int(cfg.get("gpu_workers", 1) or 1))
    except (TypeError, ValueError):
        die(f"gpu_workers must be a whole number of GPU stages, not "
            f"{cfg.get('gpu_workers')!r}")
    total_cpu = max(1, int(cfg["threads"]))
    total_ram = parse_ram(cfg.get("ram_gb"))   # the config may say '64G' too
    if not total_ram:
        det = detect_ram_gb()
        if det:
            # Leave headroom: the budget is what tools may take, not what the
            # machine has, and the rest of the pipeline needs to fit too.
            total_ram = max(1, int(det * 0.8))
            log(f"detected {det} GB of RAM, budgeting {total_ram} GB "
                "(set ram_gb or --ram to override)")
    if total_ram and total_ram < workers:
        log(f"{total_ram} GB across {workers} stages leaves under 1 GB each; "
            "using 1 GB per stage, but lower stage_workers or use --serial",
            "WARN")
    # Same "will it actually run" test decide() uses, so the announcement is
    # not made about a stage the run is not going to reach: enabled, or named
    # with --only, which is how the GPU box runs these against a server config.
    gpu_stages = [st["name"] for st in STAGES if st.get("gpu")
                  and st["name"] in selected
                  and (cfg["run"].get(st["enabled"], False)
                       or st["name"] in only_set)]
    if workers > 1 and len(gpu_stages) > gpu_slots:
        log(f"{', '.join(gpu_stages)} all use gpu_device {cfg['gpu_device']}, "
            f"so at most {gpu_slots} of them runs at a time; CPU-only stages "
            "keep running alongside. Raising gpu_workers puts them on the "
            "SAME card, not one per card.")
    if workers > 1:
        log(f"scheduling up to {workers} stages at a time, sharing "
            f"{total_cpu} cpu"
            + (f" and {total_ram} GB" if total_ram else "")
            + "; whatever is free goes to the stages that are starting")
        if not total_ram:
            log("no memory budget set and none detected; tools will use their "
                "own defaults and concurrent stages can oversubscribe RAM. "
                "Pass --ram.", "WARN")

    done = set()
    ran = adopted = skipped = 0
    failure = []
    # state_lock, not `lock`: the results lock taken at the top of this
    # function is also called lock, and a closure over the name (the signal
    # handler was one) got whichever had been assigned most recently.
    state_lock = threading.Lock()

    def finish(name, action, sig=None, err=None, secs=None):
        nonlocal ran, adopted, skipped
        with state_lock:
            if err is not None:
                state[name] = {"signature": None, "status": "failed",
                               "error": str(err)[:500],
                               "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
                failure.append((name, err))
                # Said here, not after the drain: a stage that fails in its
                # first second used to stay invisible until the longest
                # concurrent stage finished, which can be hours.
                log(f"stage '{name}' failed: {err}", "FATAL")
            elif action == "adopt":
                state[name] = {"signature": sig, "status": "adopted",
                               "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
                adopted += 1
            elif action == "RUN":
                state[name] = {"signature": sig, "status": "ok",
                               "seconds": round(secs, 1),
                               "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
                ran += 1
            else:
                skipped += 1
            save_state(p.state, state)
            done.add(name)

    def mark_running(name):
        """Recorded before the stage starts, so that a run killed mid-write
        leaves a trace. Without it the half-written file was the only evidence
        left, and the next run adopted it as a finished one."""
        with state_lock:
            state[name] = {"signature": None, "status": "running",
                           "started": time.strftime("%Y-%m-%dT%H:%M:%S")}
            save_state(p.state, state)

    def share(pending):
        """CPU and RAM for one stage about to start, out of what is free.

        Splitting the machine statically by stage_workers left the serial tail
        — InterProScan, jackhmmer, foldseek, and the emapper --dbmem decision —
        on a quarter of the box even when nothing else was running.
        """
        slots = max(1, min(workers - len(futures), pending))
        cpu = max(1, (total_cpu - sum(a[0] for a in alloc.values())) // slots)
        ram = 0
        if total_ram:
            ram = max(1, (total_ram - sum(a[1] for a in alloc.values()))
                      // slots)
        return cpu, ram

    def worker(st, cpu, ram):
        name = st["name"]
        set_log_context(name)
        t0 = time.time()
        try:
            st["fn"](stage_cfg(cfg, cpu, ram), p)
        except BaseException as e:                    # noqa: BLE001
            # BaseException on purpose: die() raises SystemExit, which in a
            # worker thread would otherwise end only that thread and let the
            # run continue as if the stage had succeeded.
            return name, None, e, 0.0
        return name, signature(st, cfg, p), None, time.time() - t0

    def needs_gpu(name):
        return bool(by_name[name].get("gpu"))

    remaining = [st["name"] for st in STAGES]
    alloc = {}                    # future -> (cpu, ram) committed to a stage
    gpu_waiting = set()           # said once per stage, not once per round
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {}
        while (remaining or futures) and not failure:
            ready = [n for n in remaining
                     if all(d in done for d in by_name[n]["deps"])]
            progressed = False
            run_now = []
            for name in ready:
                st = by_name[name]
                # Decided at dispatch, after every dependency has finished, so
                # a stage always sees its inputs in their final state.
                action = decide(st)
                if action != "RUN":
                    remaining.remove(name)
                    progressed = True
                    if action == "adopt":
                        log(f"--- {name}: adopting output this run did not "
                            "produce: "
                            + "; ".join(_file_note(o) for o in st["out"](p))
                            + ". Check that it is complete — nothing here can "
                              "tell a finished file from an interrupted one.",
                            "WARN")
                        finish(name, "adopt", sig=signature(st, cfg, p))
                    else:
                        log(f"--- {name}: {action}")
                        finish(name, action)
                    continue
                unmet = unmet_deps(st)
                if unmet:
                    remaining.remove(name)
                    progressed = True
                    finish(name, "RUN", err=RuntimeError(
                        # The stage names come after the phrase, not before
                        # it, so the sentence reads the same for one stage as
                        # for five and the troubleshooting table can quote it.
                        "cannot run: the stage(s) this one reads produced "
                        "nothing and were not selected: "
                        f"{', '.join(unmet)}. Each of them is enabled in "
                        "run:, so its evidence is expected here, not optional "
                        "— running without it would produce a confident, "
                        "empty answer. Rerun without --only/--from, or add "
                        f"{' '.join(unmet)} to the selection."))
                    continue
                run_now.append(name)
            # Longest first, so a long stage late in the table does not wait
            # behind a short one ahead of it for a worker. Python's sort is
            # stable, so stages of equal rank keep table order and the run
            # log reads the way it always did.
            run_now.sort(key=stage_priority, reverse=True)
            # The GPU is not divisible the way the CPU and RAM budgets are, so
            # it is leased rather than shared. A deferred stage stays in
            # `remaining` and is reconsidered next round; it never occupies a
            # worker while it waits.
            held = [n for n in futures.values() if needs_gpu(n)]
            run_now, waiting = gpu_lease(run_now, futures.values(), gpu_slots,
                                         needs_gpu)
            for name in waiting:
                # Once per stage, not once per round: an enabled stage that has
                # not started should be explained, not repeated at.
                if name not in gpu_waiting:
                    gpu_waiting.add(name)
                    log(f"--- {name}: waiting for the GPU — "
                        f"{', '.join(held) or 'another stage'} is using it "
                        f"and gpu_workers is {gpu_slots}. It starts when that "
                        "stage finishes; everything else carries on "
                        "meanwhile.")
            # Dispatched only once every ready stage has been decided, so the
            # share each one gets is measured against the stages that really
            # start alongside it.
            for i, name in enumerate(run_now):
                if len(futures) >= workers:
                    break
                st = by_name[name]
                cpu, ram = share(len(run_now) - i)
                remaining.remove(name)
                progressed = True
                log(f"=== {name}: running ({cpu} cpu"
                    + (f", {ram} GB" if ram else "") + ")")
                mark_running(name)
                fut = ex.submit(worker, st, cpu, ram)
                futures[fut] = name
                alloc[fut] = (cpu, ram)
            if not futures:
                # Skipping stages is progress: their dependents may have become
                # ready in this same round, so loop again before concluding
                # anything is stuck.
                if progressed:
                    continue
                if remaining:
                    stuck = [n for n in remaining
                             if not all(d in done for d in by_name[n]["deps"])]
                    die(f"deadlock: {stuck} can never become ready "
                        f"(unsatisfied deps: "
                        f"{ {n: [d for d in by_name[n]['deps'] if d not in done] for n in stuck} })")
                break
            for fut in concurrent.futures.as_completed(list(futures)):
                name, sig, err, secs = fut.result()
                del futures[fut]
                alloc.pop(fut, None)
                set_log_context(None)
                finish(name, "RUN", sig=sig, err=err, secs=secs)
                break

    # Drain anything still running so every failure is reported, not just the
    # one that happened to be noticed first.
    if failure and futures:
        log(f"stopping, but {len(futures)} stage(s) are still running and "
            "cannot be interrupted; waiting for them", "WARN")
    for fut in list(futures):
        if not fut.cancel():
            try:
                name, sig, err, secs = fut.result()
                finish(name, "RUN", sig=sig, err=err, secs=secs)
            except BaseException as e:                 # noqa: BLE001
                # finish() never saw this one, so say it here: the summary
                # below prints names only.
                log(f"stage '{futures[fut]}' failed: {e}", "FATAL")
                failure.append((futures[fut], e))
        futures.pop(fut, None)

    if failure:
        # Each one was already reported the moment it happened; this is the
        # summary at the end of the log.
        log("failed: " + ", ".join(n for n, _ in failure), "FATAL")
        if len(failure) > 1:
            log("more than one failed because they were running concurrently; "
                "fix them together, or use --serial to fail on the first",
                "FATAL")
        return 1

    log(f"done: {ran} run, {adopted} adopted, {skipped} skipped. "
        f"Results in {p.R}/")
    return 0


def main():
    # Before anything can log: a description carrying one U+FFFD used to be
    # able to kill a multi-hour run on a cp1252 console.
    configure_console_streams()
    ap = argparse.ArgumentParser(
        prog="metaannot",
        description="Single-file metaproteome functional annotation pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("init", help="write a config template")
    s.add_argument("--out", default="config.yaml")
    s.add_argument("--force", action="store_true")
    s.set_defaults(func=cmd_init)

    s = sub.add_parser("doctor",
                       help="check tools, databases and inputs, and offer to install them")
    s.add_argument("--config", default=None)
    s.add_argument("--install-plan", metavar="FILE",
                   help="write a reviewable shell script that installs what is missing")
    s.add_argument("--fix", action="store_true",
                   help="download and install what is missing, after confirmation")
    s.add_argument("--yes", "-y", action="store_true",
                   help="skip the confirmation prompt for --fix")
    s.set_defaults(func=cmd_doctor)

    s = sub.add_parser("subset", help="build the identified-protein fasta")
    s.add_argument("--db", required=True, help="full protein database fasta")
    s.add_argument("--quant", required=True)
    s.add_argument("--format", default="diann",
                   choices=sorted(ALL_FORMATS | {"ids"}))
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_subset)

    s = sub.add_parser("run", help="run the pipeline")
    s.add_argument("--config", default=None)
    s.add_argument("--faa", default=None, help="override proteins_faa")
    s.add_argument("--results-dir", default=None)
    s.add_argument("--threads", type=int, default=None)
    s.add_argument("--ram", default=None, metavar="GB",
                   help="total memory budget, split across concurrent stages "
                        "(e.g. 64, 64G, 512M). Default: 80%% of detected RAM")
    s.add_argument("--only", nargs="+", metavar="STAGE",
                   help=f"run only these: {' '.join(STAGE_NAMES)}")
    s.add_argument("--from", dest="from_stage", metavar="STAGE",
                   help="run from this stage onwards")
    s.add_argument("--force", action="store_true",
                   help="ignore cached stage state")
    s.add_argument("--no-adopt", action="store_true",
                   help="recompute stages whose outputs exist but were "
                        "produced outside this results directory")
    s.add_argument("--serial", action="store_true",
                   help="one stage at a time (lower peak memory, simpler logs)")
    s.add_argument("--force-unlock", action="store_true",
                   help="take over a results directory locked by a process "
                        "that is still alive")
    s.add_argument("--dry-run", action="store_true")
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("report", help="write and render the R report")
    s.add_argument("--config", default=None)
    s.add_argument("--results-dir", default=None)
    s.add_argument("--no-render", action="store_true",
                   help="write the Rmd but do not call Rscript")
    s.set_defaults(func=cmd_report)

    s = sub.add_parser("object", help="assemble the results into one R object")
    s.add_argument("--config", default=None)
    s.add_argument("--results-dir", default=None)
    s.add_argument("--out", default=None, help="output .rds path")
    s.add_argument("--no-run", action="store_true",
                   help="write the R script but do not call Rscript")
    s.set_defaults(func=cmd_object)

    s = sub.add_parser("all", help="run the pipeline, then the report")
    s.add_argument("--config", default=None)
    s.add_argument("--faa", default=None)
    s.add_argument("--results-dir", default=None)
    s.add_argument("--threads", type=int, default=None)
    s.add_argument("--ram", default=None, metavar="GB")
    s.add_argument("--force", action="store_true")
    s.add_argument("--no-adopt", action="store_true")
    s.add_argument("--serial", action="store_true")
    s.add_argument("--force-unlock", action="store_true")
    s.add_argument("--no-render", action="store_true")
    s.add_argument("--no-run", action="store_true")
    s.add_argument("--out", default=None)
    s.add_argument("--only", nargs="+", metavar="STAGE", default=None)
    s.add_argument("--from", dest="from_stage", metavar="STAGE", default=None)
    s.add_argument("--dry-run", action="store_true")
    s.set_defaults(func=cmd_all)

    args = ap.parse_args()
    try:
        rc = args.func(args) or 0
    except StageError as e:
        log(str(e), "FATAL")
        sys.exit(1)
    except BrokenPipeError:
        # Piping into head/less closes stdout early; exit quietly rather than
        # dumping a traceback over the user's terminal.
        try:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
        except OSError:
            pass
        sys.exit(0)
    except KeyboardInterrupt:
        log("interrupted. Stages that FINISHED are cached and skipped next "
            "time, but the stage that was running was never recorded, so its "
            "half-written output would be adopted as if it were complete: "
            "rerun that one with --only <stage> --force before trusting it.",
            "WARN")
        sys.exit(130)
    try:
        sys.stdout.flush()
    except BrokenPipeError:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    sys.exit(rc)


if __name__ == "__main__":
    main()
