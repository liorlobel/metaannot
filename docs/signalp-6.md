# SignalP 6.0 — still to install (Lior said he would download it)

SignalP 6.0 is the only tool metaannot needs that cannot be fetched
automatically. It is licence-gated: DTU requires you to register and accept the
academic terms yourself, so `doctor` marks it MANUAL and never attempts it.

## Get it

1. Go to https://services.healthtech.dtu.dk/services/SignalP-6.0/ and request
   the **"fast"** model for academic use. They email a download link.
2. Install into the environment that holds the other search tools:

```bash
wsl -d FedoraLinux-43
source ~/miniconda3/etc/profile.d/conda.sh && conda activate smorf
pip install signalp-6-package/          # the unpacked tarball from DTU
```

3. Check it is visible:

```bash
signalp6 --version && which signalp6
```

## Then turn the stage on

`signalp` and `tmbed` are both driven by one config flag. In `config.yaml`:

```yaml
run:
  topology: true
signalp_mode: "fast"      # fast | slow | slow-sequential
```

Then confirm before a long run:

```bash
python metaannot.py doctor --config config.yaml
```

The SignalP row should read `OK` instead of `MANUAL`.

## What is missing until then

`run.topology` also covers **tmbed**, which needs a GPU and is not installed
here either, so the topology stage is off entirely. Without it:

- `sp_class` (SEC/SPI, LIPO/SPII, TAT, PILIN) is empty for every protein, so
  the signal-peptide weights contribute nothing to `effector_score`.
- `n_tmb` and the beta-barrel weight are empty, so outer-membrane proteins are
  not scored.
- `surface_or_secreted` — the **gate** on the report's effector shortlist —
  is narrowed, **not** closed. It is the OR of four terms, and two of them do
  not come from topology at all: the C-terminal LPxTG sortase motif, computed
  from the sequence, and an anchor domain from `anchor_pfams`, which comes
  from the `pfam` stage. Those two keep working.

So an empty shortlist here does **not** explain itself. On the 38,204-protein
UC test run, with topology off, 604 proteins still passed the gate — 573 by
anchor domain, 43 by LPxTG — and 155 of those were KO-less. The shortlist was
empty for an unrelated reason: nothing was significant. Not one of the 212
protein groups that reached the model in the first contrast passed FDR 0.05
(the best was 0.083), so the list would have been empty with SignalP installed
too.

Read the differential-abundance table before blaming the missing tool. Check
the shortlist's other two conditions first — significance, then `has_ko` —
because those are what emptied it there, and only then the gate.

What SignalP genuinely adds is every secreted protein carrying neither an
LPxTG motif nor an anchor domain, which is most of them: a Sec/SPI substrate
with an ordinary N-terminal signal peptide is invisible to both surviving
terms. The shortlist without it is not empty, it is biased toward
cell-wall-anchored surface proteins.

## Note on the version

metaannot calls the executable `signalp6` with `--fastafile / --organism other
/ --output_dir / --format none / --mode <signalp_mode>` and reads
`prediction_results.txt`. That is the SignalP **6.0** command line. SignalP 5
takes different flags and writes a different file, so 6.0 specifically is what
is needed.
