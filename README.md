## IsoSim for flux estimation in a small network

In this work, we tried to replicate the use of IsoSim, specifically the model "prenylpyrophosphate pathway" in https://github.com/MetaSys-LISBP/IsoSim with a small network with only three metabolites in *Synechocystis sp.* PCC 6803 with 6 different strains. Methods for cultivation and labeling can be found here https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.70412.
Follow the authors' guide to install and run IsoSim first with their example models to guarantee that the environment is ready to be tested with other files.
These strains differ in the functionality or lack thereof of the protein CP12, which downregulates the Calvin-Benson Cycle under dark conditions:

1. Wild Type
2. Δcp12
3. Δcp12::cp12 (complement, wild type-like)
4. Δcp12::cp12ΔCysC (complement, with mutation in the Cysteine pair in C-terminal)
5. Δcp12::cp12ΔCysN (complement, with mutation in the Cysteine pair in N-terminal)
6. Δcp12::cp12ΔCysNC (complement, with mutation in both Cysteine pairs in C-terminal and N-terminal)

The three metabolites for which we had measurements of relative enrichment and metabolite concentration were 3PGA, 2PGA and PEP under dark conditions only with labeled <sup>13</sup>CO<sub>2</sub> as carbon source.
Therefore, we created scenarios to compare results in flux estimation, either by reducing the number of timepoints, utilizing or not the actual metabolite concentrations and considering or not reversibility. 
A combination of these factors were used in the scenarios created:

*Scenario 1: All 7 timepoints (0, 5, 10, 15, 30, 60, 90 mins), no metabolite concentrations
  + Full 3 metabolite system: 3PGA <-> 2PGA <-> PEP -> ∅
  + Partial system 3PGA <-> 2PGA -> ∅
  + Partial system 2PGA <-> PEP -> ∅

*Scenario 2: All 7 timepoints (0, 5, 10, 15, 30, 60, 90 mins), with metabolite concentrations
  + Full 3 metabolite system: 3PGA <-> 2PGA <-> PEP -> ∅
  + Partial system 3PGA <-> 2PGA -> ∅
  + Partial system 2PGA <-> PEP -> ∅

*Scenario 3: All 7 timepoints except initial used (5, 10, 15, 30, 60, 90 mins), with metabolite concentrations
  + Full 3 metabolite system: 3PGA <-> 2PGA <-> PEP -> ∅
  + Partial system 3PGA <-> 2PGA -> ∅
  + Partial system 2PGA <-> PEP -> ∅

*Scenario 4: All 7 timepoints (0, 5, 10, 15, 30, 60, 90 mins), with metabolite concentrations, irreversible
  + Full 3 metabolite system: 3PGA -> 2PGA -> PEP -> ∅
  + Partial system 3PGA -> 2PGA -> ∅
  + Partial system 2PGA -> PEP -> ∅

Unfortunately, the equations for determined fluxes were not possible to implement with running codes as in the models "example_network.R" from the authors (https://github.com/MetaSys-LISBP/IsoSim).
That is why we conducted the experiments under the assumption of a global flux through the reactions.

The independent runs with added noise can be obtained through "Repeat and source 10 runs.R", just be careful to change the names of folders and files accordingly.

"General Plots.R" helps to visualize both the label propagation simulation created by IsoSim through the different strains and scenarios, with their respective cases and some other important figures, namely, plots containing the flux estimation and confidence interval through the 10 runs.
