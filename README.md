<a href ="https://www.mrcctu.ucl.ac.uk/"><img src="logo_ukri-mrc-ctu_transparent-background.png" width="50%" /></a>

# metafloat
 v0.15  30aug2023

# A Stata package to estimate covariate interactions and subgroup-specific treatment effects in meta-analysis

Meta-analysis is a statistical technique for combining results from multiple independent studies, with the aim of estimating a single overall effect. However, it can also be used for assessing how treatment effects vary across participant subgroups. Ideally, this assessment will be based on treatment-covariate interaction effects derived _within each trial separately_, so that the pooled effects are free from aggregation bias. In particular, `metafloat` provides a simple way of estimating pooled interactions and "floating" subgroup-specific treatment effects from a set of observed (published, or otherwise pre-aggregated) treatment effects by trial and by subgroup.

See [Godolphin et al (2022)](https://doi.org/10.1002/jrsm.1590) for more information on the underlying methodology.

# Installation

To install the package directly from GitHub, type from within Stata:

    . net describe metafloat, from("https://raw.githubusercontent.com/UCL/metafloat/master/src/")

# Usage and documentation

Currently, documentation on usage and options may be found in the documentation files within Stata.  After installation, type in Stata:

    . help metafloat
