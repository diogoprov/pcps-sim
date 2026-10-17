# pcps-sim

Research compendium for a simulation study on the validity of inference based on
phylogeny-weighted species composition (phylogenetic fuzzy weighting, matrix **P**,
and Principal Coordinates of Phylogenetic Structure, PCPS), and on how tree shape
affects it. A secondary aim compares PCPS with phylogenetic community composition
(Peres-Neto & Kembel 2015), evoPCA (Pavoine 2016) and phylofactorization.

Author: Diogo B. Provete.

Status: planning. See `docs/PLANO_SIMULACAO.md` (in Portuguese).

## Structure

```
R/              functions
analysis/       targets pipeline and run scripts
data/raw/       immutable inputs (e.g., empirical phylogenies)
data/derived/   simulated trees and simulation outputs
manuscript/     Quarto manuscript
docs/           project notes and simulation plan
legacy/         material from the 2012–2016 version of the project (read-only)
references/     literature (not tracked)
```
