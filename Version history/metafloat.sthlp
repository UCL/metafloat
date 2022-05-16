{smcl}
{* *! version 0.13 beta 05apr2022}{...}
{cmd:help metafloat}
{hline}

{title:Title}

{phang}
{hi:metafloat} {hline 2} Routine for estimating covariate interactions and subgroup-specific treatment effects in aggregate data meta-analysis

{title:Syntax}

{p 8 17 2}
{cmd:metafloat} {it:ES seES} {ifin} {cmd:,} {opt study(varname)} {opt subgroup(varname)} [{it:{help metafloat##options:options}}]

{pstd}
where {it:ES seES} are variables containing effect sizes and standard errors within subgroups within studies.
Effect sizes must be based on a Normal distribution; for example, log odds-ratios rather than odds ratios.

{marker options}{...}
{synoptset 24 tabbed}{...}
{synopthdr :options}
{synoptline}
{syntab :Required options}
{synopt :{opt study(varname)}}specifies the variable containing the study identifier{p_end}
{synopt :{opt subgroup(varname)}}specifies the variable containing the subgroup identifier{p_end}

{syntab :Heterogeneity covariance structures}
{synopt :{opt unstr:uctured}}unstructured random-effects for both SigmaGamma and SigmaBeta (default){p_end}
{synopt :{opt fixed}}all fixed (common) effects{p_end}
{synopt :{opt exch:angeable}}exchangeable structures for both SigmaGamma and SigmaBeta{p_end}
{synopt :{opt randomb:eta}}special case of {opt exchangeable} with common-effects on SigmaGamma{p_end}
{synopt :{opt wscorrzero}}specia -case of {opt exchangeable} with zero within-study covariances for SigmaBeta{p_end}

{syntab :Other options}
{synopt :{opt notrend}}don't test for linear trend across subgroups{p_end}
{synopt :{opt augvar:iance(string)}}specify the augmentation variance for missing/imprecise observations{p_end} 
{synopt :{opt des:ign}}additional parameters in final model describing the available subgroups per trial (e.g. "single-subgroup" trials){p_end} 
{synopt :{opt naive}}unstructured random-effects for both SigmaGamma and SigmaBeta{p_end}
{synopt :{opt show:models}}all {bf:{help mvmeta}} models are displayed with {bf:{help noisily}}{p_end} 
{synopt :{opt nounc:ertainv}}{bf:{help mvmeta}} option; see mvmeta documentation{p_end} 
{synopt :{opt listc:onstraints}}list the constraints applied to the {bf:{help mvmeta}} model{p_end}
{synopt :{opt noaugtr:end}}don't use variance augmentation for estimating trend{p_end} 

{syntab :Saved datasets}
{synopt :{opt saving(saving_option)}}save data in the form of a "forestplot results set" to {it:filename}{p_end} 
{synopt :{opt clear}}replace the data in memory with the "results set", instead of saving to a separate file{p_end} 
{synopt :{opt keepvars}}carry forward variables from the original dataset into the saved dataset{p_end} 

{title:Description}

{pstd}{cmd:metafloat} estimates "floating" subgroup-specific treatment effects in an aggregate data meta-analysis. These floating 
subgroup effects are consistent with the pooled treatment-covariate within-trial interaction(s), with these pooled 
within-trial interactions free from aggregation bias. {cmd:metafloat} requires a set of observed (published, or otherwise pre-aggregated) 
treatment effects and standard errors by trial and by subgroup.

Requires {cmd:metan} version 4+ and {cmd:mvmeta} version 4+.

{pstd}For more information on the underlying methodology utilised in {cmd:metafloat}, see {help metafloat##refs:Godolphin et al (2022)}.

{marker options}{...}
{title:Options}

{dlgtab:Required options}

{phang}
{opt study(varname)} specifies a variable containing the study identifier, which must be either integer-valued or string.

{phang}
{opt subgroup(varname)} specifies a variable containing the subgroup identifier, which must be either interger-valued or string.

{dlgtab:Heterogeneity covariance structures}

{phang}
{opt unstructured} unstructured random-effects for both SigmaGamma and SigmaBeta (default).

{phang}
{opt fixed} all fixed (common) effects.

{phang}
{opt exchangeable} exchangeable structures for both SigmaGamma and SigmaBeta.

{phang}
{opt randombeta} special-case of exchangeable with common-effects on SigmaGamma, implying constant within and between-subgroup heterogeneity for SigmaBeta. 

{phang}
{opt wscorrzero} special-case of exchangeable with zero within-study covariances for SigmaBeta.

{dlgtab:Other options}

{phang}
{opt notrend} doesn't test for a linear trend across subgroups. This option is only applicable when the number of subgroups is greater then two. N.B. for a binary subgroup a test of trend is not carried out. 

{phang}
{opt design} additional parameters in final model describing the available subgroups per trial (e.g. "single-subgroup" trials)

{phang}
{opt naive} unstructured random-effects for both SigmaGamma and SigmaBeta, removing the constraint that they be related via the {bf:M} matrix.

{phang}
{opt showmodels} all {bf:{help mvmeta}} models are displayed with {bf:{help noisily}}.

{phang}
{opt nouncertainv} {bf:{help mvmeta}} option; see mvmeta documentation.

{phang}
{opt noaugtrend} don't use variance augmentation for estimating trend.

{dlgtab:Saved datasets}

{phang}
{opt saving(saving_option)} save dataset of results to {it:filename}, structured for use with {bf:{help forestplot}} (i.e. matched subgroups and interactions).

{phang}
{opt clear} leave dataset of results in memory, rather than saving to {it:filename}, replacing the original dataset.

{phang}
{opt keepvars} carry forward variables from the original dataset into the saved dataset, including design (if specified).


{title:Examples}

{phang}{to do}{p_end}

{title:Notes for us}
Have not included text on the following options: forcetrend, constraint, listconstraints


{title:Authors}

{phang}David Fisher, Peter Godolphin, Ian White{break}
MRC Clinical Trials Unit at UCL

{marker refs}{...}
{title:References}

{phang}
Godolphin PJ, White IR, Tierney JF, Fisher DJ. 2022.
Estimating interactions and subgroup-specific treatment effects in meta-analysis without aggregation bias: A within-trial framework.
Research Synthesis Methods 

