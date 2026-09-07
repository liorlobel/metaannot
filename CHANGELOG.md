# Changelog

## v0.2.0 — 2026-09-07

A test suite, and the first three defects it and the first real run turned up.
`SIGNATURE_VERSION` is unchanged at 1: no stage's output means anything
different, so an existing results directory stays valid and nothing recomputes.

### Added

**`tests/` and CI, where there were none.** 419 tests, one named for each
defect found by hand during development or by the first real run, every one
carrying a one-line comment stating the symptom it protects against. The
suite runs offline with no external tool: where a stage shells out to
hmmsearch, DIAMOND or MMseqs2 the binary is a stub on PATH, so the stage's own
plumbing — atomic writes, adoption, signatures, the results lock — is
exercised without the tool. Fixtures are generated from fixed seeds; nothing
binary is committed. The default suite is about a hundred seconds and excludes
a `slow` mark carrying the resume sweep and parallel-versus-serial across
every `quant_format` x `peptide_assignment`. R tests skip cleanly where R is
absent and include a real knit where it is present. CI covers Python 3.9-3.13.

This also settles two claims the README made and nothing checked: parallel and
serial execution now provably produce identical output, and five repeated runs
produce one hash.

Eight tests are `xfail(strict)` rather than passing. Each names a guard that
is still missing — among them `parse_emapper` raising a bare `KeyError` on a
zero-row table, and a gap in a Unipept lineage truncating a protein's
consensus to the rank before it. They are marked so a fix turns them green
rather than being quietly forgotten.

**`quant/sample_columns.txt`.** `build_object.R` and the report both read it
as the authoritative statement of which columns carry intensities, and nothing
wrote it, so every run fell through to `design_from_input.tsv` — or, with no
manifest, to guessing from column types. The `join` stage now records what it
detected, renamed and filtered. Not one of the stage's declared outputs, so an
existing results directory keeps its old fallback until `join` reruns;
`--force --only join` produces the file without recomputing anything else.

### Fixed

**A non-UTF-8 byte in a tool's output is no longer fatal.** One `0xa0` — a
latin-1 non-breaking space — in a search result killed the `integrate` stage
of the first full run on real data, after InterProScan had already spent three
hours, with a message that named neither the file nor the stage. Every parser
reads through `opener()`, and all six died on it: `parse_diamond`,
`parse_interproscan`, `parse_kofam`, `parse_hmm_tblout`, `parse_hmm_lib_desc`
and `read_fasta`. VFDB subject titles, InterPro descriptions, HMM `DESC` lines
and FASTA headers all carry latin-1 in the wild. Now `errors="replace"`, the
same as `read_delim_table` beside it, so a bad byte costs one character of a
description rather than the stage.

The gzip branch of the same function passed no encoding at all, so it used the
locale's codec and was never UTF-8 by construction. `emapper_precomputed` is
routinely a `.gz`, which made it the input most likely to be read differently
on the server than on the laptop.

**The effector shortlist is not "empty by construction" without SignalP.**
`docs/signalp-6.md` and README said `surface_or_secreted` is False for
everything when the topology stage is off, so an empty shortlist was a
missing-tool artefact. Two of the gate's four terms need no topology tool: the
LPxTG sortase motif, computed from the sequence, and an anchor domain from the
`pfam` stage. On the first real run 604 proteins passed the gate with topology
off, 155 of them KO-less; the shortlist was empty because nothing was
significant — none of the 212 groups that reached the model passed FDR 0.05.
The docs now say what SignalP genuinely adds, which is every secreted protein
carrying neither motif, and that a shortlist without it is not empty but
biased toward cell-wall-anchored surface proteins.

### Known limitations

Unchanged from v0.1.0, and still the ones that decide how a result is read:
isobaric labelling is refused rather than quantified, only the first contrast
gets the full analysis, and the InterProScan memory budget is inherited by
child JVMs in distributed mode. See the v0.1.0 entry below.

Issue #5 records a further one found by the same run: the report prints
retention by bin, including a bin that contributes zero proteins to the model,
but nothing escalates it — and on that run neither `4_dark` nor `3d_duf_only`
reached the statistics at all.

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
- **The InterProScan memory budget is passed through `_JAVA_OPTIONS`**, which
  the JVM applies at higher precedence than a command-line `-Xmx` and which
  every child JVM inherits. That is correct for the default standalone mode,
  where InterProScan's workers are threads inside one JVM, and it is how the
  stage budget actually reaches it. If you run InterProScan in distributed
  mode, its worker JVMs are configured for 9 GB each in
  `interproscan.properties` and would instead inherit the whole stage budget,
  so divide `--ram` by the worker count there.
- **SignalP 6.0 is licence-gated** and must be installed by hand; see
  `docs/signalp-6.md`. Without it, and without a GPU for tmbed, the topology
  stage is off, so the signal-peptide and beta-barrel terms of
  `surface_or_secreted` are dead. That narrows the effector shortlist's gate
  rather than closing it: an LPxTG sortase motif or an anchor domain still
  qualifies a protein, and neither needs a topology tool.
  *(Corrected after v0.2.0. This bullet originally said the shortlist is then
  "empty by construction rather than by result", which is wrong — see the
  v0.2.0 entry above. Nothing else in the v0.1.0 entry has been changed.)*
