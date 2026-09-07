# Changelog

## v0.1.0 — 2026-09-07

First public release. This is the version that has actually been run end to end
on a real label-free metaproteomics dataset, which no earlier state of the code
had been.

### Verified on real data

A 36-sample label-free FragPipe dataset against a matched-metagenome Prokka
database with a precomputed eggNOG table: all identified proteins resolved to
the database, 122,278 peptide features rolled up to 16,070 protein groups
across 832 taxa, the R Markdown report knitted to HTML with figures and tables,
and the QFeatures object built with both assays linked. Bin percentages fell
inside the ranges TUTORIAL.md gives as sane.

Not yet exercised on real data: the external search tools (hmmsearch, DIAMOND,
InterProScan, KOfamScan, HHblits, jackhmmer, ESMFold, Foldseek) and the Unipept
API. Their parsers and integration logic are tested, and the tools install
cleanly, but a full run with them enabled has not been completed.

### Fixed

An audit of the previous state found 211 defects. 250 fixes landed across two
passes; the counts differ because several findings needed changes in more than
one place. The ones that changed results rather than tidying code:

- **Razor identifiers never matched the FASTA.** The quantification reader
  preferred FragPipe's `Protein ID` column, which for a non-UniProt database
  holds `<id> <description>`. Identifiers are now taken from `Protein` and
  reduced to their first whitespace token, the same key the FASTA uses.
- **eggNOG's own domain calls were ignored when binning.** Only hmmsearch Pfam
  hits counted, so with the Pfam stage off, proteins eggNOG had already
  assigned a domain to were reported as having no evidence at all.
- **Fractionated manifests were rejected.** Repeated sample names are how
  FragPipe denotes fractions; they are now collapsed into one sample per group
  instead of aborting the run.
- **Taxonomy identifiers were corrupted on reuse.** A cached table was re-read
  without string dtypes, so a numeric taxid round-tripped through float and was
  written back as `821.0`.
- **A search database that prefixes its identifiers could not be joined** to an
  annotation table that does not carry the prefix, and the diagnostic then
  reported the proteins as genuinely unannotated. New `emapper_strip_id_prefix`.
- **FragPipe's zero means "not quantified"** and was being summed as a real
  zero, turning missingness into fold change. Now missing by default, under
  `zero_intensity_is_missing`.
- **Text files were written with the platform codec**, so on Windows the report
  died on a Unicode character and left an empty file. Every text file is now
  UTF-8.
- **Interrupted writers left truncated outputs** that the next run adopted as
  complete. Every stage now writes atomically.
- **The dark-protein work-list was capped for structure prediction**, and the
  profile-search stages reused that capped file, so a GPU budget silently
  truncated an unrelated search. The two lists are now separate.
- **Two configured download URLs were dead**, and one returned an HTML landing
  page that `curl -fL` accepted with exit 0, writing a web page to disk as a
  database. Downloads are now rejected if they begin with HTML.
- **Report and object plumbing**: contrasts were derived from the wrong term of
  the design formula, level names were not sanitised the way R sanitises them,
  parameters reached the document unquoted, and a failing `Rscript` did not set
  a non-zero exit status.

Silent behaviour is now reported: proteins dropped by the feature-count filter,
proteins withheld from folding by the structure cap, decoy and contaminant rows
removed, taxa whose size factor falls back to a plain sum, and which taxonomic
rank is in force.

### Added, all defaulting to existing behaviour

`rollup_method` (a log-space alternative to the plain sum), a
`taxon_or_family_unique` peptide-assignment mode, `taxon_rank`,
`ncbifam_uninformative_test`, `foldseek_target_priority`, `full_content_digest`,
`peptide_only_reader`, and `exclude_id_prefixes`.

### Known limitations

- **Isobaric labelling is not supported.** Reporter-ion channels are not read.
  FragPipe TMT output is refused with an explanatory error rather than being
  quantified from the pooled MS1 column, which is what earlier code did
  silently. See the limitation section in README.md.
- **Only the first contrast gets the full analysis.** Differential abundance is
  fitted for every contrast, but the ratio model, the effector shortlist and
  the R object cover the first one only.
- **SignalP 6.0 is licence-gated** and must be installed by hand; see
  `docs/signalp-6.md`. Without it, and without a GPU for tmbed, the topology
  stage is off and the report's effector shortlist is empty by construction
  rather than by result.
