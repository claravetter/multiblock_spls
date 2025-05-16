---
title: 'MB-SPLS: A multiblock sparse partial least squares toolbox for multimodal association analysis in translational psychiatry'
tags:
  - machine learning
  - multiblock analysis
  - neuroimaging
  - psychiatry
  - dimensionality reduction
  - MATLAB
authors:
  - name: Clara S. Vetter
    orcid: 0000-0000-0000-0000  # replace with your actual ORCID
    affiliation: 1
  - name: Florian Eichin
    affiliation: 2
  - name: David Popovic
    affiliation: 1
  - name: Nikolaos Koutsouleris
    affiliation: 1
affiliations:
  - name: Department of Psychiatry, LMU Munich & Max Planck Institute of Psychiatry, Munich, Germany
    index: 1
  - name: Munich Center for Machine Learning (MCML), Munich, Germany
    index: 2
date: 2025-05-15
---

# Summary

**MB-SPLS** is a multiblock sparse partial least squares toolbox developed to enable joint modeling of multiple high-dimensional data blocks. It allows researchers to uncover shared latent structures across modalities such as neuroimaging, genetics, clinical assessments, and environmental data—especially in psychiatry and related biomedical domains.

Unlike standard sparse PLS implementations, MB-SPLS supports:
- **Simultaneous integration of more than two data blocks**
- **View-specific sparsity constraints** for interpretability
- **Frobenius-norm based association optimization**
- **Nested cross-validation** for hyperparameter tuning
- **Permutation testing** for statistical significance
- **Bootstrap-based stability analysis**
- **Projection-based matrix deflation** to extract successive associative effects

The toolbox is implemented in MATLAB and distributed as a compiled app, with extensive documentation and example workflows.

# Statement of Need

Many biomedical machine learning applications require integration across diverse data modalities. Classical methods such as Canonical Correlation Analysis (CCA) or Partial Least Squares (PLS) only operate on pairs of datasets and lack built-in support for sparsity, model validation, or multiblock structure.

While sparse PLS (SPLS) improves interpretability, its two-block limitation makes it inadequate for more complex settings. MB-SPLS addresses this gap by providing:
- **True multiblock integration** for more than two data views
- **Nested cross-validation** for reliable model selection
- **Permutation testing and bootstrap resampling** to assess statistical reliability and feature importance

This enables researchers to move beyond predictive modeling and toward interpretable, reproducible multimodal signature discovery—critical in translational psychiatry and other data-rich fields.

# Implementation and Features

The MB-SPLS software consists of:
- A core implementation of the multiblock SPLS algorithm, based on Frobenius norm maximization of between-block covariances
- A nested cross-validation framework for hyperparameter selection (`c_i` sparsity constraints per view)
- Significance testing via label permutation of one block, corrected for multiple comparisons
- Feature weight stability estimation using bootstrap sampling
- Iterative matrix projection deflation for extracting successive latent variables
- Output visualizations of identified latent components

The software is implemented in MATLAB and packaged as a compiled application. 


# Repository

The software, example data, and usage instructions are available at:  
**https://github.com/claravetter/multiblock_spls/**  
(Replace with actual URL)

# License

MB-SPLS is released under the MIT License.

# Acknowledgements

This work was developed within the context of the PRONIA and OMLS projects. We thank the NeuroMiner development team for foundational components. Support was provided by the Munich Center for Machine Learning (MCML), the Helmholtz Association, and the Max Planck Institute of Psychiatry.

# References

- Witten, D. M., Tibshirani, R., & Hastie, T. (2009). A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. *Biostatistics*, 10(3), 515–534.
- Monteiro, J. M. (2016). *Multiple sparse partial least squares: A multiblock framework for the integration of omics data sets*. PhD Thesis.
- Tenenhaus, A., Philippe, C., Guillemot, V., Le Cao, K. A., Grill, J., & Frouin, V. (2014). Variable selection for generalized canonical correlation analysis. *Biostatistics*, 15(3), 569–583.
- Le Cao, K. A., Gonzalez, I., & Dejean, S. (2009). integrOmics: An R package to unravel relationships between two omics datasets. *Bioinformatics*, 25(21), 2855–2856.

