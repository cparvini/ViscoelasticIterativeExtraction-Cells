# "Viscoelastic Parameterization of Human Skin Cells to Characterize Material Behavior at Multiple Timescales"
## Accompanying Code Distribution
This repository accompanies the manuscript "Viscoelastic Parameterization of Human Skin Cells to Characterize Material Behavior at Multiple Timescales" by Parvini, Cartagena-Rivera, and Solares. The results from analyzing several adherent skin cell lines (Human Primary Melanocytes, 2D Adherent Human Melanoma, and Human Fibroblasts) are presented and discussed in the manuscript.

As laid out in the main document, this code uses AFM Force Curve information to parameterize viscoelastic models of varying complexity. The optimal parameter sets are then analyzed, and conclusions can be drawn from how the viscoelastic harmonic quantities (Storage Modulus, Loss Modulus, and Loss angle) vary as a function of frequency.

The Viscoelastic Parameter Extraction methodology has been previously outlined in the literature ([1](https://www.beilstein-journals.org/bjnano/articles/11/77),[2](https://onlinelibrary.wiley.com/doi/abs/10.1002/polb.24327)). However, the addition of iterative term introduction both increases optimization stability and the overall likelihood of acquiring an optimal parameter set. In the main document, the performance of the "open search" (legacy) methodology is compared to the new iterative one, and the increased stability is shown to incur negligible computational costs.

## Software Requirements
* Matlab (>= 9.7), including:
  * Signal Processing Toolbox
  * Optimization Toolbox
  * Statistics and Machine Learning Toolbox
  * Symbolic Math Toolbox
  * DSP System Toolbox
  * Parallel Computing Toolbox
  * Matlab Parallel Server
  * Curve Fitting Toolbox
  * Fixed-Point Designer
  * Matlab Coder
  * (Optional) NSMatlabUtilities (to read Bruker files)

## Citations
[[1]](https://www.beilstein-journals.org/bjnano/articles/11/77) Parvini, C. H.; Saadi, M. A. S. R.; Solares, S. D. Beilstein J. Nanotechnol. 2020, 11, 922–937. doi:10.3762/bjnano.11.77

[[2]](https://onlinelibrary.wiley.com/doi/abs/10.1002/polb.24327) López‐Guerra, E.A., Eslami, B. and Solares, S.D. (2017), Calculation of standard viscoelastic responses with multiple retardation times through analysis of static force spectroscopy AFM data. J. Polym. Sci. Part B: Polym. Phys., 55: 804-813. https://doi.org/10.1002/polb.24327
