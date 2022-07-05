{smcl}
{* *! version 0.13 beta 05apr2022}{...}
{cmd:help metafloat}
{hline}

{vieweralsosee "metan" "help metan"}{...}
{vieweralsosee "mvmeta" "help mvmeta"}{...}
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
{synopt :{opt unstr:uctured}}unstructured random effects for both SigmaGamma and SigmaBeta (default){p_end}
{synopt :{opt fixed}}all fixed (common) effects{p_end}
{synopt :{opt exch:angeable}}exchangeable structures for both SigmaGamma and SigmaBeta{p_end}
{synopt :{opt randomb:eta}}special case of {opt exchangeable} with common effect on Gamma (i.e. SigmaGamma = 0){p_end}
{synopt :{opt wscorrzero}}special case of {opt exchangeable} with zero within-study covariances for SigmaBeta{p_end}

{syntab :Other options}
{synopt :{opt augvar:iance(string)}}specify the augmentation variance for missing/imprecise observations{p_end} 
{synopt :{opt des:ign}}additional parameters in final model describing the available subgroups per trial (e.g. "single-subgroup" trials){p_end} 
{synopt :{opt eform}}report exponentiated effect sizes{p_end} 
{synopt :{opt naive}}unstructured random-effects for both SigmaGamma and SigmaBeta{p_end}
{synopt :{opt show:models}}display all intermediate {bf:{help mvmeta}} models{p_end} 
{synopt :{opt nounc:ertainv}}{bf:{help mvmeta}} option; see mvmeta documentation{p_end}
{synopt :{opt notrend}}suppress test for linear trend across subgroups{p_end}

{syntab :Saved datasets}
{synopt :{opt saving(saving_option)}}save data in the form of a "{help metan:forestplot results set}" to {it:filename}{p_end} 
{synopt :{opt clear}}replace the data in memory with the "results set", instead of saving to a separate file{p_end} 
{synopt :{opt keepvars(varlist)}}carry forward variables from the original dataset into the saved dataset{p_end} 

{title:Description}

{pstd}{cmd:metafloat} estimates "floating" subgroup-specific treatment effects in an aggregate data meta-analysis. These floating 
subgroup effects are consistent with the pooled treatment-covariate within-trial interaction(s), with these pooled 
within-trial interactions free from aggregation bias. {cmd:metafloat} requires a set of observed (published, or otherwise pre-aggregated) 
treatment effects and standard errors by trial and by subgroup.

{pstd}Requires {bf:{help metan}} version 4+ and {bf:{help mvmeta}} version 4+.

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

{pmore}The mathematics of an exchangeable structure suggests two "special cases". It is currently unclear what use these might have.

{pmore}{opt randombeta} places a random effect on the Betas only, with a common effect for the Gammas, implying constant within and between-subgroup heterogeneity for SigmaBeta;

{pmore}{opt wscorrzero} forces zero within-study covariances for SigmaBeta.

{dlgtab:Other options}

{phang}
{opt augvariance(string)} allows the augmentation variance for missing/imprecise observations to be changed; the default is 100,000

{phang}
{opt design} includes additional parameters in final model to identify subgroups of trials on the basis of their available subgroups (e.g. "single-subgroup" trials)

{phang}
{opt eform} reports exponentiated effect sizes, if the supplied effect data represents log-odds ratios or similar.

{phang}
{opt naive} unstructured random-effects for both SigmaGamma and SigmaBeta, removing the constraint that they be related via the {bf:M} matrix.

{phang}
{opt showmodels} all intermediate {bf:{help mvmeta}} models are displayed with {bf:{help noisily}}.

{phang}
{opt nouncertainv} {bf:{help mvmeta}} option; see mvmeta documentation.

{phang}
{opt notrend} suppresses the test for a linear trend across subgroups.
Note that this option is only applicable when the number of subgroups is greater then two; for a binary subgroup, a test of trend is not carried out.

{dlgtab:Saved datasets}

{phang}
{opt saving(saving_option)} save dataset of results to {it:filename}, structured for use with {bf:{help forestplot}} (i.e. matched subgroups and interactions).

{phang}
{opt clear} leave dataset of results in memory, rather than saving to {it:filename}, replacing the original dataset.

{phang}
{opt keepvars} carry forward variables from the original dataset into the saved dataset, including design (if specified).


{title:Examples}

{pstd}
All examples use aggregate data from the WHO REACT Working Group Prospective meta-analysis of interleukin-6 antagonists for hospitalised COVID-19 patients (WHO REACT 2021).
This data is also provided in Appendix 2 of Godolphin et al.

{pstd}
Load the data

{cmd}{...}
{phang2}. clear{p_end}
{phang2}. input str17 TrialName str3 Subgroup int(n0 e0 n1 e1){p_end}
{phang2}"ARCHITECTS"        "No"  	0   0   1   0{p_end}
{phang2}"ARCHITECTS"        "Yes"   	11  2   9   0{p_end}
{phang2}"BACC-Bay"          "No"  	81  4   158 9{p_end}
{phang2}"BACC-Bay"          "Yes"   	1   0   3   0{p_end}
{phang2}"CORIMUNO-TOCI-1"   "No"  	55  5   53  6{p_end}
{phang2}"CORIMUNO-TOCI-1"   "Yes"   	12  3   10  1{p_end}
{phang2}"CORIMUNO-TOCI-ICU" "No"  	39  8   41  4{p_end}
{phang2}"CORIMUNO-TOCI-ICU" "Yes"   	4   2   8   4{p_end}
{phang2}"COV-AID"           "No"  	30  4   33  3{p_end}
{phang2}"COV-AID"           "Yes"   	42  3   48  6{p_end}
{phang2}"COVACTA"           "No"  	103 16  238 44{p_end}
{phang2}"COVACTA"           "Yes"   	41  12  56  14{p_end}
{phang2}"COVIDOSE2-SS-A"    "No"  	6   1   13  0{p_end}
{phang2}"COVIDOSE2-SS-A"    "Yes"   	2   1   6   0{p_end}
{phang2}"EMPACTA"           "No"  	16  0   49  2{p_end}
{phang2}"EMPACTA"           "Yes"   	112 11  200 24{p_end}
{phang2}"HMO-020-0224"      "No"  	2   0   6   1{p_end}
{phang2}"HMO-020-0224"      "Yes"   	15  8   31  10{p_end}
{phang2}"ImmCoVA"           "No"  	1   0   1   0{p_end}
{phang2}"ImmCoVA"           "Yes"   	26  2   21  2{p_end}
{phang2}"PreToVid"          "No"  	11  2   11  1{p_end}
{phang2}"PreToVid"          "Yes"   	167 32  159 20{p_end}
{phang2}"RECOVERY"          "No"  	367 127 357 139{p_end}
{phang2}"RECOVERY"          "Yes"   	1721 600 1664 482{p_end}
{phang2}"REMAP-CAP"         "No"  	129 41  127 30{p_end}
{phang2}"REMAP-CAP"         "Yes"   	217 73  214 53{p_end}
{phang2}"REMDACTA"          "No"  	29  2   72  9{p_end}
{phang2}"REMDACTA"          "Yes"   	181 39  358 69{p_end}
{phang2}"TOCIBRAS"          "No"  	30  1   34  6{p_end}
{phang2}"TOCIBRAS"          "Yes"   	34  5   31  8{p_end}
{phang2}. end{p_end}
{txt}{...}

{pstd}
Use {bf:{help metan}} to derive log-odds ratios and standard errors from the raw counts

{cmd}{...}
{phang2}. generate long f1 = n1 - e1{p_end}
{phang2}. generate long f0 = n0 - e0{p_end}
{phang2}. quietly metan e1 f1 e0 f0, or nogr study(Trial) by(Subgroup) nooverall iv{p_end}
{txt}{...}

{pstd}
Finally, run {cmd:metafloat}. By default, unstructured random-effects will be used.

{cmd}{...}
{phang2}. metafloat _ES _seES, study(Trial) subgroup(Subgroup) 
{txt}{...}


{title:Authors}

{phang}David Fisher, Peter Godolphin, Ian White{break}
MRC Clinical Trials Unit at UCL

{marker refs}{...}
{title:References}

{phang}
Godolphin PJ, White IR, Tierney JF, Fisher DJ. 2022.
Estimating interactions and subgroup-specific treatment effects in meta-analysis without aggregation bias: A within-trial framework.
Research Synthesis Methods 

{phang}
WHO Rapid Evidence Appraisal for COVID-19 Therapies Working Group. 
Association Between Administration of IL-6 Antagonists and Mortality Among Patients Hospitalized for COVID-19: A Meta-analysis. 
JAMA 2021; 326: 499-518.

