\# Quantum Depletion in the BEC with LHY Correction – A Variational Approach and Numerical Analysis



This repository contains two implementations that compute and visualize the energy and density profiles of an isotropically trapped Bose–Einstein condensate (BEC) including \*\*quantum depletion\*\* through the \*\*Lee–Huang–Yang (LHY)\*\* correction with the variational approach.



Both codes explore the effects of depletion on physical quantities such as condensate fraction, density and total energy. While one implementation examined it under the constraint of fixed total particle number, the other investigated the behaviors for varying interaction strengths. 



---



\## Contents



\- `julia/bec\_gp\_lda\_variational.jl` – Variational approach in \*\*Julia\*\* using `Optimization.jl` and `Integrals.jl`.

\- `python/Energy\_and\_Density\_with\_depletion.ipynb` – Numerical computation and plotting in \*\*Python\*\* using SciPy and Matplotlib.



---



\## Requirements



\### Julia

```bash

using Pkg

Pkg.add(\["Plots", "Roots", "Polynomials", "LaTeXStrings", "Integrals", "Optimization", "OptimizationOptimJL"])

julia julia/bec\_gp\_lda\_variational.jl



