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

- `sp_class` (SEC/SPI, LIPO/SPII, TAT) is empty for every protein, so the
  signal-peptide weights contribute nothing to `effector_score`.
- `n_tmb` and the beta-barrel weight are empty, so outer-membrane proteins are
  not scored.
- `surface_or_secreted` is False for everything, which is the **gate** on the
  report's effector shortlist. The shortlist is therefore empty by
  construction, not because no candidate exists.

That last point is the one that matters: an empty effector shortlist in the
current report is a missing-tool artefact, not a scientific result. Anything
depending on secretion prediction has to wait for SignalP.

## Note on the version

metaannot calls the executable `signalp6` with `--fastafile / --organism other
/ --output_dir / --format none / --mode <signalp_mode>` and reads
`prediction_results.txt`. That is the SignalP **6.0** command line. SignalP 5
takes different flags and writes a different file, so 6.0 specifically is what
is needed.
