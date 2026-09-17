# Cross-Polarity Acquisition Grouping: Synthetic Spectrum Batches Labelled By Which Acquisitions Share A Compound

## Overview

**What this dataset is.** Each record is one acquisition: a list of peaks from a tandem mass spectrometry run, in either positive or negative ionization mode, belonging to a batch. The label is not a pairing between the two polarities and not a property of any single spectrum. It is a **partition of each batch**: which acquisitions in that batch came from the same compound. A compound may contribute several acquisitions in one polarity, acquisitions in both, or acquisitions in only one, so groups have no fixed size and no one-to-one structure.

This is an original, fully synthetic dataset of tandem mass spectrometry batches generated from a hidden fragmentation grammar. It contains no record of any spectral library, no real compound and no measurement of any real instrument. Every spectrum is produced by the generation procedure described below from a secret that is not published, so nothing here can be found, matched or regenerated from any external source.

Concretely, each compound in a batch is acquired zero to three times in positive mode and zero to three times in negative mode, with at least one acquisition overall. The correspondence between the two polarities is therefore many-to-many.

Acquisitions are degraded in three ways that are never marked on the individual spectrum: co-isolation of a neighbouring compound, detector saturation, and low injection.

## Release At A Glance

- Raw files: 6
- Compound families: 900, each a set of mass-matched, spectrally confusable compounds
- Batches: 2,700, three per family drawn with independent acquisitions
- Compounds per batch: 12.6 on average
- Acquisitions: 103,849, 38.5 per batch on average, 29.2 peaks each
- Peak rows: 3,033,843
- Acquisitions per compound: 1.56 in positive mode and 1.50 in negative mode on average; 43.7 percent of compounds appear in one polarity only, and 75.0 percent have two or more acquisitions in at least one polarity
- Degradation rates: 15 percent co-isolation, 15 percent saturation, 15 percent low injection, applied independently per acquisition
- Hidden vocabulary: 64 substructure units; compounds are chains of 3 to 6 units
- Prepared release: 585 training families (1,755 batches, all three repeats each) and 315 test families (315 batches, one each); no compound is shared across the split
- Data origin: creator-generated synthetic data

## Raw File Structure

The uploaded ZIP is flat and contains exactly these six files at its root: `acquisitions.csv`, `peaks.csv`, `grouping.csv`, `LICENSE`, `DATASET_DESCRIPTION.md` and `PACKAGE_MANIFEST.sha256`. Every column of every data file is documented below.

There is no `spectra.csv` in this release, because a spectrum is not one row. A spectrum is the set of `peaks.csv` rows sharing a `case_id` and a `slot`, and its metadata is the matching row of `acquisitions.csv`.

**`acquisitions.csv`** — 103,849 rows, one per acquisition. This is the index of the release.

- `case_id` (string): the batch this acquisition belongs to. A batch is the scoring unit.
- `block_id` (string): the compound family the batch was drawn from. Three batches share each `block_id`. This is the unit of the train/test split; splitting on `case_id` instead will give optimistic validation.
- `slot` (string): the acquisition's label within its batch, `P00, P01, ...` for positive mode and `N00, N01, ...` for negative. Unique within a `case_id`. Slot order is randomised and carries no information about the grouping.
- `polarity` (string): `POSITIVE` or `NEGATIVE`. Redundant with the first character of `slot`, and kept because solvers group on it.

**`peaks.csv`** — 3,033,843 rows, one per peak, 29.2 peaks per acquisition on average.

- `case_id` (string), `slot` (string): join keys to `acquisitions.csv`.
- `mz` (float): the peak's mass-to-charge position.
- `log_intensity` (float): `log1p` of a relative intensity scaled so each spectrum's base peak is 999. Apply `expm1` to recover the linear intensity.

**`grouping.csv`** — 103,849 rows, one per acquisition. **This file is the answer key.**

- `case_id` (string), `slot` (string): join keys to `acquisitions.csv`, one row for every acquisition in the release.
- `compound_group` (string): an opaque identifier for the compound that produced this acquisition. Two acquisitions came from the same compound exactly when this value matches, and only within the same `case_id`; the identifier is not comparable between batches and encodes nothing about the compound itself.

`prepare.py` reads `grouping.csv` to build `train_labels.csv` for the training split, and to build the private answer file for the test split. The test split's groupings are never written to any public prepared file. It exposes no field beyond the three columns above, and no compound identity, formula, mass, collision energy, injection amount or degradation flag exists anywhere in this release.

**`LICENSE`** — the full CC BY 4.0 legal text.

**`DATASET_DESCRIPTION.md`** — this document.

**`PACKAGE_MANIFEST.sha256`** — SHA-256 checksum of every other file in the package.

## How The Data Was Generated

Every draw and every identifier comes from HMAC-SHA256 keyed to a withheld 256-bit secret; the generator code needs that secret to produce anything and is published without it.

1. **Grammar.** Draw 64 substructure units, each with a mass between 44 and 176. For each unit and each polarity, draw 2 to 5 fragment ions with base intensities, a propensity for each of the polarity's 8 neutral losses, and an end-of-chain intensity boost; with probability 0.15 the unit is silent in that polarity. Draw the two loss tables and a sparse table of neighbour-context intensity effects for ordered unit pairs.
2. **Families.** Each family starts from a seed chain of 3 to 6 units and adds up to 17 variants by substituting one or two units for units of similar mass, so members are mass-matched and share most of their structure.
3. **Batches.** A batch draws 10 to 15 compounds from one family. Each compound is acquired 0 to 3 times per polarity, with at least one acquisition overall; the draw is repaired if needed so that every batch has at least two compounds present in both polarities and at least one compound acquired twice in some polarity.
4. **Clean spectra.** For each acquisition: every non-silent unit emits its fragments with position, context and log-normal intensity jitter, each retained with probability 0.85, plus neutral-loss peaks by propensity; sequence ions run along the chain in the polarity's own direction with probability 0.90 and in the opposite direction with probability 0.60; the precursor ion and its losses appear with fixed probabilities; 3 to 8 low-intensity noise peaks are added; every m/z receives Gaussian jitter of 0.003.
5. **Degradation.** With probability 0.15 the acquisition is co-isolated: the clean spectrum of another compound of the same batch is added at 10 to 40 percent of the base intensity. With probability 0.15 it is saturated: all intensities are clipped at 25 to 55 percent of the base peak. With probability 0.15 it is low-injection: only the 5 to 12 most intense peaks are retained. None of these is recorded anywhere in the release.
6. **Anonymisation.** Peaks below 1/999 of the base peak are dropped and at most 80 are kept. Acquisitions are shuffled within each polarity and relabelled `P00, P01, ...` and `N00, N01, ...`, so slot order carries no information about the grouping.

## Intended Use And Limitations

- Intended use: research and benchmarking of many-to-many cross-modal correspondence recovered as a constrained partition, under degradations that break the assumption that one spectrum represents one source.
- Out of scope: any claim about real fragmentation chemistry, real compounds or real instruments. The grammar is a stylised stand-in, not a model of mass spectrometry, and nothing learned here transfers to real spectra.
- No personal data, no real laboratory data and no third-party material is included.
