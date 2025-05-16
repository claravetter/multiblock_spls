---
title: 'MB-SPLS: A multiblock sparse partial least squares toolbox for multimodal association analysis in translational psychiatry'
tags:
  - machine learning
  - psychiatry
  - multiblock analysis
  - neuroimaging
  - sparse PLS
  - MATLAB
  - data integration
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

MB-SPLS is a multiblock sparse partial least squares toolbox designed to jointly analyze multiple high-dimensional datasets (e.g., neuroimaging, genetics, clinical) and identify latent factors that capture cross-modal associations. Built on an extension of the SPLS method to the multiblock case, it employs Frobenius norm maximization, view-specific sparsity, and projection-based deflation to extract interpretable associative effects.

The toolbox includes a nested cross-validation framework for hyperparameter optimization, permutation testing for statistical significance, and bootstrap resampling to assess weight stability. It is tailored to translational psychiatry but applicable to any multimodal setting with structured data blocks.

# Statement of Need

Existing implementations of PLS and SPLS do not support multiblock extensions or integrated model validation frameworks. MB-SPLS addresses this gap by enabling simultaneous modeling of multiple data blocks with view-specific sparsity and robust evaluation procedures. It allows users to uncover multimodal signatures that are statistically validated, interpretable, and generalizable—key requirements for applications in psychiatric research and beyond.

# Implementation

The MB-SPLS toolbox is implemented in MATLAB and packaged as a compiled app with detailed documentation. The toolbox allows:
- Exploratory or predictive modeling with multiple data blocks
- Frobenius-norm based optimization
- Cross-validation (nested or split-half)
- Permutation testing
- Bootstrap-based stability analysis
- Projection deflation for iterative LV extraction

# Repository

The software is available at: https://github.com/YOURNAME/mbspls  
(Replace with actual GitHub URL)

# Acknowledgements

This work was developed within the context of the PRONIA and OMLS projects. We thank the NeuroMiner team for foundational infrastructure and the Munich Center for Machine Learning (MCML) for support.

# References

References are automatically pulled from your `paper.bib` or inline citations.
