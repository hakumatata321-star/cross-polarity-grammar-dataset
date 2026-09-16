# Cross-Polarity Fragmentation Grammar Dataset: Paired Positive- and Negative-Mode Spectra From A Withheld Rule Set

## Overview

This is an original, fully synthetic dataset of tandem mass spectra generated from a hidden fragmentation grammar. It contains no record of any public spectral library, no real compound and no measurement of any real instrument. Every spectrum is produced by the generation procedure described below from a secret that is not published, so no spectrum in this release can be found, matched or regenerated from any external source.

A compound is a short chain of substructure units drawn from a hidden vocabulary. Each unit carries two different fragment signatures, one for positive-mode ionization and one for negative-mode, plus mode-specific neutral losses; some units produce no fragments at all in one of the two modes. The positive-mode and negative-mode spectra of one compound are therefore generated from the same chain through two different tables, share almost none of their peaks, and correspond to each other only through the hidden grammar.

The labels are the correspondence: which negative-mode acquisition, if any, was produced by the same compound as each positive-mode acquisition.

## Release At A Glance

- Raw files: 5
- Compounds: 11,200, in 700 families of up to 16 mass-matched, spectrally confusable variants
- Compounds acquired in both modes: 6,983; the rest in one mode only, so that batches contain single-polarity decoys
- Spectra: 36,366, two replicate acquisitions per compound and mode, 29.3 peaks on average
- Peak rows: 1,064,775
- Hidden vocabulary: 64 substructure units; chains of 3 to 6 units
- Prepared release: 337 training batches, each drawn three times with different replicate spectra (1,011 training cases), and 159 test batches; test compounds appear in no training batch
- Data origin: creator-generated synthetic data

## Raw File Structure

The uploaded ZIP is flat and contains exactly these five files at its root:

- `spectra.csv`: one record per acquisition: `spectrum_id`, `compound_id`, `ion_mode` (`POSITIVE` or `NEGATIVE`), `precursor_mz`.
- `peaks.csv`: one record per peak: `spectrum_id`, `mz`, `log_intensity`. Intensities are `log1p` of a relative intensity scaled so that each spectrum's base peak is 999.
- `LICENSE`: CC BY 4.0.
- `DATASET_DESCRIPTION.md`: this document.
- `PACKAGE_MANIFEST.sha256`: SHA-256 checksum of every other file in the package.

`compound_id` is the creator-side key that `prepare.py` uses to build batches and labels. It never appears in the public prepared files, where spectra are identified only by anonymous slots within a batch.

## How The Data Was Generated

Every draw and every identifier comes from HMAC-SHA256 keyed to a withheld 256-bit secret; the generator code and the secret are not released.

1. **Grammar.** Draw 64 units, each with a mass between 44 and 176. For each unit and each mode, draw 2 to 5 fragment ions (a fixed sub-mass of the unit plus a mode-specific adduct shift) with base intensities, a propensity for each of the mode's 8 neutral losses, and an end-of-chain intensity boost; with probability 0.15 the unit is silent in that mode. Draw the two loss tables and a sparse table of neighbour-context intensity effects for ordered unit pairs.
2. **Families.** Each family starts from a seed chain of 3 to 6 units and adds up to 15 variants by substituting one or two units for units of similar mass, so members are mass-matched and share most of their structure.
3. **Modes.** Each compound is acquired in both modes with probability 0.62, otherwise in one mode only.
4. **Spectra.** For each acquisition: every non-silent unit emits its fragments with position, context and log-normal intensity jitter (sd 0.35), each retained with probability 0.85, plus neutral-loss peaks by propensity; sequence ions are prefix sums in positive mode and suffix sums in negative mode; the precursor ion and its losses appear with fixed probabilities; 3 to 8 low-intensity noise peaks are added; every m/z receives Gaussian jitter of 0.003. Peaks below 1/999 of the base peak are dropped and at most 80 are kept.
5. **Batches.** The published `prepare.py` groups compounds acquired in both modes into batches of 7 to 10 by mass window and negative-mode spectral similarity, adds mass-matched single-mode decoys on each side, anonymises every spectrum to a slot, and holds out about a third of the batches with their compounds excluded from every training batch.

## Intended Use And Limitations

- Intended use: research and benchmarking of cross-modal correspondence under a one-to-one constraint with abstention, where the mapping must be learned from paired examples rather than measured by a distance.
- Out of scope: any claim about real fragmentation chemistry, real compounds or real instruments. The grammar is a stylised stand-in, not a model of mass spectrometry.
- No personal data, no real laboratory data and no third-party material is included.

## Licence

Creative Commons Attribution 4.0 International (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/
