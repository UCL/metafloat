<a href ="https://www.ctu.mrc.ac.uk/"><img src="MRCCTU_at_UCL_Logo.png" width="50%" /></a>

# metafloat
 v0.15  30aug2023

# A Stata package to estimate covariate interactions and subgroup-specific treatment effects in meta-analysis

Meta-analysis is a statistical technique for combining results from multiple independent studies, with the aim of estimating a single overall effect. However, it can also be used for assessing how treatment effects vary across participant subgroups. Ideally, this assessment will be based on treatment-covariate interaction effects derived _within each trial separately_, so that the pooled effects are free from aggregation bias. In particular, **metafloat** provides a simple way of estimating pooled interactions and "floating" subgroup-specific treatment effects from a set of observed (published, or otherwise pre-aggregated) treatment effects by trial and by subgroup.

See Godolphin et al (Research Synthesis Methods, 2022) for more information on the underlying methodology.
