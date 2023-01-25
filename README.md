
# PMSCCA: Preference Matrix Guided Sparse Canonical Correlation Analysis

Canonical Correlation Analysis (CCA) is a statistical technique used to identify and quantify the relationship between two datasets for linear combinations of their features (loading weight vectors) that maximize the datasets' correlation. 

Additionally, Preference Matrix Guided Sparse Canonical Correlation Analysis (PMSCCA) takes priors encoded as a preference matrix $E$ and enforces sparsity on both loading vector weight vectors, which has been shown to successfully increase the interpretability and the generalizability of the loading weight vectors.

Please refer to the following paper for the technical details of the model.

***J. Sha et al., "Preference Matrix Guided Sparse Canonical Correlation Analysis for Genetic Study of Quantitative Traits in Alzheimer’s Disease," 2022 IEEE International Conference on Bioinformatics and Biomedicine (BIBM), Las Vegas, NV, USA, 2022, pp. 541-548, doi: 10.1109/BIBM55620.2022.9995342.***

This repository contains the implementation of PMSCCA in R with the files necessary for reproducing the simulation studies presented in the paper.


## How should I use the code?

*"/Simulation.R"* is the file to execute. I recommend using RStudio as it helps install package dependencies recursively.

*"/Simulated datasets/"* contains 
- two datasets $X,Y$, whose canonical correlation is to be maximized (target canonical correlation $\rho=0.3$);
- two ground truth loading weight vectors $u,v$;
- and an additional file in the folder, $C$, as the target joint covariance matrix, which is for reference only and should not contribute to the computation.

*"/simulation/"* contains a sample set of results that this implementation should be expected to produce.
