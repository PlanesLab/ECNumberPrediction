# Oxidoreductase (EC class 1) no-cofactor study

Rhea Stratified split (Case2), reaction SMILES with metabolic cofactors
stripped (see `methods/SIMMER/SIMMER_scripts/cofactors.py` /
`strip_cofactors.py`), evaluated with EC class 1 called out specifically.
Two methods run: SIMMER (`methods/SIMMER/run_SIMMER_rhea_nocofactor.sh`)
and BEC-Pred (`methods/BEC-Pred/run_BECPred_rhea_stratified_nocofactor.sh`).

## Headline finding: the two methods disagree

- **SIMMER**: stripping cofactors *helped* class 1 (MCC 0.694 -> 0.794,
  +0.100) and 4 of the other 6 classes; only class 6 (ligases) got worse.
- **BEC-Pred**: stripping cofactors *hurt* class 1 (MCC 0.855 -> 0.795,
  -0.06).