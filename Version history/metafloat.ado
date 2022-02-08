*  v0.1 beta  David Fisher 10jan2021
*  v0.2 beta  David Fisher 30mar2021
*  v0.3 beta  David Fisher 21apr2021
*  v0.4 beta  David Fisher 07may2021
*  v0.5 beta  David Fisher 01jul2021
*  v0.6 beta  David Fisher 12jul2021
*  v0.7 beta  David Fisher 04aug2021
*  v0.8 beta  David Fisher 20sep2021
*  v0.9 beta  David Fisher 14dec2021
*! v0.10 beta  David Fisher 11jan2022


*** METAFLOAT ***
* Routine for fitting the models described in Godolphin et al (2021)
// by David Fisher, Peter Godolphin, Ian White
// (MRC Clinical Trials Unit at UCL)


* Release notes:
// v0.5: following discusion between IW, PG, corrected code so that results remain invariant to change of reference
//         under random-effects, with or without `nouncertainv'
// v0.6: addition of `design' option; analogous to inconsistency parameters in NMA
//       observations with extreme variances now set to be augmented
// v0.7: dataset constructed containing floating subgroups and interactions side-by-side to facilitate "two-panel" forest plots
//       options saving() and `clear' added
// v0.8: finalizing random-effects options, particularly "randombeta" and "exchangeable"
// v0.9: moved sections into subroutines
//       added keepvars() option; improvements to trend options
// v0.10: changes to matrix setup to aid correspondence with -ipdfloat-
//       specifically: study-level variance matrices are derived from mvmeta_make rather than "manually" from the raw data
//       but since the parameterisation used by mvmeta_make is different, these matrices then need to be used in a slightly different way

// TO DO:
// debugging/test script
// for the future: are augmentations *always* necessary?  might there be scenarios in which they introduce inaccuracy??
// for the future: user-defined structures for covariance and/or trend (removed here to keep things simple)

* Syntax:
// metafloat ES seES [if] [in], study(varname) subgroup(varname) [ options ]
// where ES, seES, study() and subgroup() are all required

* Heterogeneity covariance structures:
// fixed = all fixed (common) effects
// default (no options) = unstructured random-effects for both SigmaGamma and SigmaBeta
// exchangeable = exchangeable structures for both SigmaGamma and SigmaBeta
// randombeta = special-case of exchangeable with common-effects on SigmaGamma
//   (implies constant within and between-subgroup heterogeneity for SigmaBeta)
// wscorrzero = special-case of exchangeable with zero within-study covariances for SigmaBeta

* Trend across subgroups:
* [Note: if k > 2, *test* for trend is reported automatically, but is *not* imposed upon the estimates]
// notrend = don't test for trend
// forcetrend = impose linear trend upon subgroup estimates

* Other options:
// augvariance = specify the augmentation variance for missing/imprecise observations
// design = additional parameters in final model describing the available subgroups per trial (e.g. "single-subgroup" trials)
// naive = unstructured random-effects for both SigmaGamma and SigmaBeta, removing the constraint that they be related via the `M' matrix
// showmodels = all -mvmeta- models are displayed with -noisily- (by default some intermediate steps are done -quietly-)
// listconstraints = list the constraints applied to the -mvmeta- model for estimating SigmaBeta, for debugging purposes
// nouncertainv = -mvmeta- option; see mvmeta documentation
// `options' = eform options, plus other options to pass to -mvmeta- e.g. tolerance()

* Saved datasets:
// saving() = save dataset of results, structured for use with -forestplot-  (i.e. matched subgroups and constrasts)
// clear = leave dataset of results in memory, replacing original dataset
// keepvars = carry forward variables from the original dataset into the saved dataset (includes "_Design" if specified)


program define metafloat, eclass

	version 11.0
	local version : di "version " string(_caller()) ":"
	* NOTE: mata requires v9.x
	* factor variable syntax requires 11.0

	// Check that -metan- v4.0+ is installed
	cap metan
	if "`r(metan_version)'"=="" {
		nois disp as err "This program requires {bf:metan} version 4.04 or higher"
		exit 499
	}
	else {
		local current_version = 4.04
		if `r(metan_version)' < `current_version' {
			nois disp as err "This program requires {bf:metan} version 4.04 or higher"
			exit 499
		}
	}
	// NOTE: if random-effects, also requires -mvmeta- version 3.5.1 # Ian White # 21dec2021
	// but currently no easy way to test for this

	syntax varlist(numeric min=2 max=2) [if] [in], STUDY(varname) SUBGROUP(varname) ///
		[ SHOWmodels DESign AUGVARiance(real 1e5) FORCEtrend ///
		  LISTConstraints noUNCertainv PRINT(string) ///		/* -mvmeta- options */
		  SAVING(passthru) CLEAR KEEPVars(varlist) * ]
		
	_get_eformopts, eformopts(`options') allowed(__all__) soptions
	local eformopt eform(`"`s(str)'"')
	local soptions `"`s(options)'"'
	marksample touse

	
	** Preserve data, prior to editing
	preserve
	qui keep if `touse' & !missing(`study', `subgroup')
	keep `varlist' `study' `subgroup' `keepvars'
		
	tokenize `varlist'
	args eff stderr

	cap assert `stderr' >= 0
	if _rc {
		nois disp as err "Standard error cannot be negative"
		exit 198
	}
	
	// check that study-subgroup combinations are unique
	tempvar obs Ngroup
	qui gen int `obs' = _n
	qui bysort `study' `subgroup' (`obs') : gen int `Ngroup' = _N
	summ `Ngroup', meanonly
	if r(max) > 1 {
		nois disp as err "The following combination(s) of study and subgroup are not unique:"
		nois list `study' `subgroup' `eff' `stderr' if `Ngroup' > 1
		exit 198
	}
	drop `obs' `Ngroup'
	
	tempvar V
	qui gen double `V' = `stderr'^2		// Generate estimate of the variance (s-squared)
	drop `stderr'

	capture {
		if `"`eff'"'!=`"b"' {			// use {b, V} here, and {y, S} later with mvmeta_make
			qui rename `eff' b
		}
		if `"`V'"'!=`"V"' {
			qui rename `V' V
		}
	
		// augment if missing subgroup data
		qui fillin `study' `subgroup'
		assert _fillin == missing(b)
	}
	if _rc {
		nois disp as err `"Something that should be true of your data is not"'
		nois disp as err `"In particular, variables {it:ES}, {it:seES}, {bf:study()} and {bf:subgroup()} must all be distinct"'
		if `"`keepvars'"'!=`""' {
			nois disp as err `"as should any variables in the {bf:keepvars()} option)"'
		}
		exit 459
	}

	// if any observed variances are larger than `augvariance'
	//  (e.g. if estimated from very small sample)
	// set to _fillin==2
	qui replace _fillin = 2 if V > `augvariance'

	summ _fillin, meanonly
	if r(max) {
		nois disp as text "The following subgroups are unobserved (or contain missing or extremely imprecise data)"
		nois disp as text " and will be given augmented variances of " as res `augvariance' as text ":"
		list `study' `subgroup' if _fillin
	}
	qui replace V = `augvariance' if _fillin
	qui replace b = 0 if _fillin
	
	// obtain "design", i.e. which subgroups are present
	if "`design'"!="" {
		local desvar _Design
		qui gen `desvar' = ""
		GetDesign `desvar', study(`study') subgroup(`subgroup')
		local showmodels showmodels
	}
	if trim(`"`desvar'`keepvars'"')!=`""' {
		local collapse_opt `"collapse((firstnm) `desvar' `keepvars')"'
	}

	if `"`saving'"'!=`""' | `"`clear'"'!=`""' {

		* Save subgroup value labels
		// Note: from Stata documentation:
		// " -encode- looks at a string variable and makes an internal table of all the values it takes on...
		// ... It then alphabetizes that list and assigns numeric codes "
		// so ordering of -encode- ought to match with that of -levelsof-
		cap confirm string variable `subgroup'
		if !_rc {
			tempvar subgroup2
			tempname sglab
			qui encode `subgroup', gen(`subgroup2') label(`sglab')
			qui drop `subgroup'
			qui rename `subgroup2' `subgroup'
		}
		else local sglab : value label `subgroup'

		// later, we will need to merge on `subgroup'
		// but it appears -merge- cannot keep its value label
		// therefore, attach it to y instead, and retrieve it later
		if `"`sglab'"'!=`""' {
			label values b
			label values b `sglab'
		}

		local trunclen = 9
		local sgtrunc = substr(`"`subgroup'"', 1, `trunclen')
		if `"`sgtrunc'"'!=`"`subgroup'"' {
			qui rename `subgroup' `sgtrunc'
			local subgroup `sgtrunc'
		}
		
		tempfile main_effects
		qui save `main_effects'
	}

	sreturn clear	// in advance of running s-class subroutines TransMatrix and CovStruct
	

	// Generate transformation matrix `T' based on choice of reference subgroup (= identity matrix if ref = 1st subgroup)
	// also define trend vector `d' if appropriate
	tempname T d
	TransMatrix `T' `d', subgroup(`subgroup') `soptions'
	local sglist `"`s(sglist)'"'
	local k `s(k)'
	local ref `s(ref)'
	if `k' < 3 {
		local trend notrend
		local forcetrend
	}
	
	// Setup covariance structures
	tempname sgmat propU
	CovStruct `sgmat' `propU', k(`k') trans(`T') `s(options)'
	local structure  `s(structure)'		// must be one of: fixed, unstructured, exchangeable, randombeta, wscorrzero
	local sigmagamma `s(sigmagamma)'
	local sigmabeta  `s(sigmabeta)'
	local naive      `s(naive)'			// "naive" = undocumented option; don't constrain SigmaBeta and SigmaGamma to be linked
	local soptions `"`s(options)'"'
	
	if `"`forcetrend'"'!=`""' {
		cap assert inlist(`"`structure'"', "fixed", "randombeta") | `"`s(default)'"'!=`""'		// either "fixed" or "randombeta" specified, or no structure was specified
		if _rc {
			nois disp as err `"cannot specify {bf:bscovariance(`structure')} with {bf:forcetrend}; random-effects structure is implied by use of trend"'
			exit 184
		}
	}
	
	
		
	// Generate contrasts, using prefix _J
	// Note: y_cons is the reference subgroup (i.e NOT a contrast)
	qui xi, noomit prefix(_J) i.`subgroup'
	local OldVarsPrefix : char _dta[__xi__Vars__Prefix__]
	local OldVarsToDrop: char _dta[__xi__Vars__To__Drop__]
	local i = 0
	foreach prefix of local OldVarsPrefix {
		local ++i
		if `"`prefix'"'==`"_J"' {
			local el : word `i' of `OldVarsToDrop'
			local VarsToDrop `VarsToDrop' `el'
			local varlab : variable label `el'
			tokenize `"`varlab'"', parse("=")
			args sgname eq sgval
			cap {
				assert `"`sgname'"'==`"`subgroup'"'
				assert `"`eq'"'==`"=="'
				confirm number `sgval'
				assert `"`4'"'==`""'
			}
			if _rc {
				nois disp as err `"Something has gone wrong involving -xi-"'
				exit 198
			}
			local Iyvar = subinstr(`"y`el'"', `"_J"', `"_I"', 1)
			local Iyvars `Iyvars' `Iyvar'
			if `sgval'!=`ref' {
				local usevars `usevars' `el'
				local Jyvars `Jyvars' y`el'
			}
		}
	}
	
	// Extract abbreviated version of `subgroup' used in xi-derived variable names
	// and check that it can be successfully unabbreviated
	gettoken sgabbrev : usevars
	local sgabbrev = reverse(substr(`"`sgabbrev'"', 3, .))	// remove prefix "_J", and reverse...
	gettoken number rest : sgabbrev, parse("_")				// ...so that we can remove the final "_#" stub
	confirm number `number'
	gettoken us rest : rest, parse("_")
	assert `"`us'"'==`"_"'
	local sgabbrev = reverse(`"`rest'"')
	cap {
		unab sgunab : `sgabbrev', max(1)
		assert `"`sgunab'"'==`"`subgroup'"'
	}
	if _rc {
		nois disp as err `"Variable name conflict"'
		nois disp as err `"variable name {bf:`subgroup'} in required option {bf:subgroup()} is not sufficiently different"'
		nois disp as err `"from another variable name supplied to this command; please check"'
		if `"`keepvars'"'!=`""' {
			nois disp as err `"(this includes any variables in the {bf:keepvars()} option)"'
		}
		exit _rc
	}
	
	local matexist : all matrices "S*"
	if `"`matexist'"'!=`""' {
		nois disp `"{error}Note: the following matrices in memory will be cleared as a result of running this command:"'
		nois disp `"{error}`matexist'"'
		matrix drop `matexist'
	}
			

	** Generate -mvmeta- dataset containing contrasts plus constant (reference subgroup)	
	qui mvmeta_make regress b `usevars' [iw=V^-1], mse1 by(`study') keepmat `collapse_opt' clear nodetails useconstant
	cap confirm numeric variable `Jyvars' y_cons
	if _rc {
		nois disp as err "Error: inconsistency in matrix stripe elements"
		exit 198
	}	

	label variable `study' `"`studylab'"'
	label variable y_cons `"`subgroup'"'

	// N.B. regress ... [iw], mse1  is equivalent to using -vwls- 
	//  but -regress, mse1- appears to be faster for some reason


	** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
	if "`showmodels'"=="" local quietly quietly
	if "`print'"=="" local print bscov

	// fit as standard model (assuming no trend)
	tempname VarGammaHat GammaHat SigmaGamma Chol
	nois disp as text _n "Test of interaction(s):"
	`quietly' mvmeta y S, vars(`Jyvars') print(`print') `sigmagamma' `soptions' `eformopt'
	cap assert `"`e(yvars)'"'==`"`Jyvars'"'
	if _rc {
		nois disp as err "Error: inconsistency in matrix stripe elements"
		exit 198
	}	
	nois test `Jyvars'
	local chi2_int = r(chi2)
	local df_int = r(df)
	local p_int = r(p)

	local k1 = `k' - 1
	matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
	matrix define `GammaHat'    = e(b)[1,       1..`k1']
	matrix define `SigmaGamma'  = e(Sigma)	
	
	// extract cholesky matrix of Sigma
	if `"`structure'"'=="unstructured" {
		// special case: if `unstructured', -mvmeta- outputs the elements of `Chol' as part of e(b)
		// so don't bother packing them up into a matrix; we're only going to unpack them again later		
		matrix define `Chol' = e(b)[1, `k'...]
	}
	else if inlist(`"`structure'"', "exchangeable", "wscorrzero") {
		local tau = e(b)[1, `k']
		cholesky2 `SigmaGamma' `Chol'
		matrix define `Chol' = sign(`tau') * `Chol'
	}
	else if `"`structure'"'=="randombeta" {
		matrix define `Chol' = J(`k1', `k1', 0)
	}
	
	local rownames_eb : rownames e(b)
	matrix rownames `GammaHat' = `: word 1 of `rownames_eb''	// REVISIT
	matrix coleq `GammaHat' = :
	matrix coleq `VarGammaHat' = :
	matrix roweq `VarGammaHat' = :
	
	// if `k' > 2 (and unless `notrend'), also test for trend
	local eqlist
	if `"`trend'"'==`""' {
		local i = 0		
		tokenize `Jyvars'
		while `"`1'"'!=`""' {
			local ++i
			qui gen _Trend_`i' = `d'[`i', 1]
			local eqlist `eqlist' `1':_Trend_`i',
			macro shift
		}
		
		local sgammatrend fixed
		if !inlist(`"`structure'"', "fixed", "randombeta") {
			tempname D
			matrix define `D' = `d'*`d''
			local sgammatrend bscov(proportional `D')
		}

		nois disp as text _n "Test for trend:"
		`quietly' mvmeta y S, vars(`Jyvars') print(`print') `sgammatrend' commonparm nocons equations(`eqlist') `soptions' `eformopt'
		nois test _Trend_1
		local chi2_trend = r(chi2)
		local df_trend = r(df)
		local p_trend = r(p)
		
		if `"`forcetrend'"'!=`""' {
			matrix define `VarGammaHat' = e(V)[1, 1]		// Gamma is a constant if `trend'
			matrix define `GammaHat'    = e(b)[1, 1]
			matrix define `SigmaGamma'  = e(Sigma)

			matrix define `VarGammaHat' = `VarGammaHat'*`d'*`d''
			matrix define `GammaHat'    = `GammaHat'*`d''
			
			if `"`structure'"'=="randombeta" {
				matrix define `Chol' = J(`k1', `k1', 0)
			}
			else if `"`structure'"'!="fixed" {
				local tau = e(b)[1, 2]
				cholesky2 `SigmaGamma' `Chol'
				matrix define `Chol' = sign(`tau') * `Chol'
			}
			local trendtext `" (with trend imposed)"'		// for later
		}
	}	
	cap assert `"`e(yvars)'"'==`"`Jyvars'"'
	if _rc {
		nois disp as err "Error: inconsistency in matrix stripe elements"
		exit 198
	}


	** STEP 2: ESTIMATE TREATMENT EFFECT IN REF GROUP (MEAN, BS VARIANCE AND BS CORRELATION WITH INTERACTIONS)
	// subtract within-trial interactions from observed interaction data
	// set up equations that make shifted interactions have mean zero

	// Note: use permanent names so that they display nicely in the -mvmeta- output
	qui gen byte _Zero = 0
	qui gen byte _One  = 1

	// optionally adjust for "design" : use pre-defined variable `desvar' to derive set of indicators
	if "`design'"!="" {
		qui levelsof `desvar', local(dlist)
		local designvars
		foreach d of local dlist {
			local destxt = subinstr(trim(`"`d'"'), `" "', `"_"', .)
			qui gen _Des_`destxt' = (`desvar'==`"`d'"')
			local designvars `designvars' _Des_`destxt'
		}
	}		

	local eqlist
	local i = 0
	tokenize `Jyvars'
	while `"`1'"'!=`""' {
	    local ++i
		qui replace `1' = `1' - `GammaHat'[1, `i']
		local eqlist `eqlist' `1':_Zero `designvars',
		macro shift
	}
	local eqlist equations(`eqlist' y_cons:_One `designvars')
	
	// set up constraints on the BS variance of the interactions
	// using the Cholesky matrix from STEP 1
	if `"`structure'"'!="fixed" & `"`naive'"'==`""' {
		SetupConstraints `Chol' `propU', k(`k') structure(`structure') `forcetrend'
		local constr_list `s(constr_list)'
		
		if `"`constr_list'"'!=`""' {
			if `"`listconstraints'"'!=`""' {
				nois disp as text _n "Constraints on elements of Cholesky factor for SigmaU:"
				nois constraint dir `constr_list'
			}
			local constr_opt constraints(`constr_list')
		}
	}
	
	// fit model: make sure we are using "NEW" mvmeta
	if "`showmodels'"!="" nois disp as text _n `"Estimation of floating subgroup for reference category`trendtext':"'
	`quietly' mvmeta y S, print(`print') nocons commonparm `eqlist' `constr_opt' `sigmabeta' `soptions' `eformopt'
	// Note: the coefficient of _cons refers to the reference subgroup
	
	// Define covariate data matrix Z (arranged such that the reference subgroup comes first)
	// Also define matrix L, which transforms "contrast + reference" parameterisation
	// into "subgroup-specific" parameterisation
	tempname ones Z L
	matrix define `ones' = J(`k', 1, 1)			// column vector of ones, of length k
	matrix define `Z' = J(`k', `k1', 0)
	matrix define `Z'[2, 1] = I(`k1')	
	matrix define `L' = `Z', `ones'
	
	// Compute implied means for beta
	tempname ThetaHat BetaHat
	matrix define `ThetaHat' = e(b)[1,1]
	matrix define `BetaHat' = `ones'*`ThetaHat' + `Z'*`GammaHat''
	
	// First, derive `VarThetaHat' ignoring the uncertainty in the heterogeneity matrix
	// because it's needed for deriving the matrix A...
	tempname VarThetaHat
	matrix define `VarThetaHat' = e(V)
	matrix define `VarThetaHat' = invsym(`VarThetaHat')
	matrix define `VarThetaHat' = `VarThetaHat'[1,1]
	matrix define `VarThetaHat' = invsym(`VarThetaHat')	

	// Define matrix A, the independent multiplier of GammaHat
	// In order to do this, we need to extract matrices U_i, where L * U_i * L' = S_i = diagonal matrix of subgroup-specific variances
	// then L*(U`i' + SigmaU)*L' = S`i' + SigmaBeta ... and we need to sum these over all studies `i'
	// The required matrices are available as _S* via the "keepmat" option to -mvmeta_make-
	local n = e(N)							// number of trials
	local Ulist : all matrices "S*"
	assert `n'==`: word count `Ulist''
	tempname SigmaU A W
	matrix define `SigmaU' = e(Sigma)
	matrix define `A' = J(1, `k1', 0)
	foreach U of local Ulist {
		matrix define `W' = `L'*(`U' + `SigmaU')*`L''
		matrix define `A' = `A' + `ones''*invsym(`W')*`Z'
	}
	matrix define `A' = `Z' - `ones'*`VarThetaHat'*`A'
	matrix drop `Ulist'
		
	// ... and now correct VarThetaHat for uncertainty in the heterogeneity matrix if requested
	if `"`uncertainv'"'==`""' {
		matrix define `VarThetaHat' = e(V)[1,1]
	}
		
	// correct variance for beta
	tempname VarBetaHat
	matrix define `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A'*`VarGammaHat'*`A''
	
	// reverse transformation
	matrix define `BetaHat'    = `T''*`BetaHat'
	matrix define `VarBetaHat' = `T''*`VarBetaHat'*`T'

	// ...also applies to SigmaBeta
	tempname SigmaBeta
	matrix define `SigmaBeta' = `T''*`L'*`SigmaU'*`L''*`T'
	matrix coleq `SigmaBeta' = ""
	matname `SigmaBeta' `Iyvars', explicit
	
	// collect and display under "ereturn"
	matrix define `BetaHat' = `BetaHat''
	matrix coleq `BetaHat' = :
	matname `BetaHat'    `Iyvars', columns(.) explicit
	matname `VarBetaHat' `Iyvars', explicit

	
	*** Construct saved dataset
	if `"`saving'"'!=`""' | `"`clear'"'!=`""' {

		// reverse previous operation
		local i = 0
		tokenize `Jyvars'
		while `"`1'"'!=`""' {
			local ++i
			qui replace `1' = `1' + `GammaHat'[1, `i']
			macro shift
		}
		drop _Zero _One
		if `"`trend'"'==`""' {
			drop _Trend*
		}
		if `"`design'"'!=`""' {
			drop _Des*
		}

		// generate matrices containing pooled estimates and variances
		tempname BetaTable GammaTable
		matrix define `GammaTable' = `GammaHat'', vecdiag(`VarGammaHat')'
		local sgnoref : list sglist - ref
		matrix rownames `GammaTable' = `sgnoref'	// guaranteed to be numeric
		matrix colnames `GammaTable' = _yInt _yInt_S
	
		matrix define `BetaTable' = `BetaHat'', vecdiag(`VarBetaHat')'
		matrix rownames `BetaTable' = `sglist'		// guaranteed to be numeric
		matrix colnames `BetaTable' = _y _y_S
		
		ProcessSavedData `GammaTable' `BetaTable', saving(`saving') ///
			study(`study') subgroup(`sgabbrev') ///
			mainfile(`main_effects') augvariance(`augvariance')
	
		if `"`clear'"'!=`""' {
			restore, not
		}
	}
	
	
	*** Finally, present the results
	ereturn post `BetaHat' `VarBetaHat', depname(Subgroup)
	ereturn matrix GammaHat    = `GammaHat'
	ereturn matrix VarGammaHat = `VarGammaHat'
	ereturn matrix SigmaGamma  = `SigmaGamma'
	ereturn matrix SigmaBeta   = `SigmaBeta'

	ereturn scalar ThetaHat    = `ThetaHat'[1,1]
	ereturn scalar VarThetaHat = `VarThetaHat'[1,1]
	ereturn scalar n = `n'		// number of studies
	ereturn scalar k = `k'		// number of subgroups

	ereturn scalar chi2_int = `chi2_int'
	ereturn scalar df_int = `df_int'
	ereturn scalar p_int = `p_int'
	
	if `"`trend'"'==`""' {
		ereturn scalar chi2_trend = `chi2_trend'
		ereturn scalar df_trend = `df_trend'
		ereturn scalar p_trend = `p_trend'
	}
	
	ereturn local bscov_Gamma `"`sigmagamma'"' 
	ereturn local bscov_Beta `"`sigmabeta'"' 
	
	nois disp as text _n `"Floating subgroups`trendtext':"'
	ereturn display, `eformopt'

	if `"`constr_list'"'!=`""' & `"`structure'"'!=`"user"' {	// only drop internally-defined constraints
	    constraint drop `constr_list'
	}
	
end



// process covariance structures
program define CovStruct, sclass
	syntax anything, K(integer) TRANS(name) [ BSCOVariance(string) NAIVE ///
		FIXED UNStructured EXCHangeable RANDOMBeta WSCORRZero WSC0 * ]
	
	local soptions : copy local options
	if `"`wsc0'"'!=`""' local wscorrzero wscorrzero
	opts_exclusive `"`fixed' `unstructured' `exchangeable' `randombeta' `wscorrzero'"' `""' 184
	local structure `fixed'`unstructured'`exchangeable'`randombeta'`wscorrzero'

	if `"`structure'"'!=`""' {
		if `"`bscovariance'"'!=`""' {
			nois disp as err `"Cannot specify both {bf:`structure'} and {bf:bscovariance()}"'
			exit 198
		}
	}
	else {
		local 0 `bscovariance'
		syntax [ , FIXED UNStructured EXCHangeable RANDOMBeta WSCORRZero WSC0 ]
		if `"`wsc0'"'!=`""' local wscorrzero wscorrzero
		opts_exclusive `"`fixed' `unstructured' `exchangeable' `randombeta' `wscorrzero'"' `""' 184
		local structure `fixed'`unstructured'`exchangeable'`randombeta'`wscorrzero'
	}
	
	// parse pre-defined tempname matrices
	tokenize `anything'
	args sgmat propU
	local T : copy local trans
	
	if `"`structure'"'=="fixed" {
		local sigmagamma fixed
		local sigmabeta fixed
	}
	else if `"`structure'"'=="randombeta" {
		local sigmagamma fixed
		matrix define `propU' = J(`k', `k', 0)
		matrix define `propU'[`k', `k'] = 1
		matrix define `propU' = `T''*`propU'*`T'
		local sigmabeta bscov(proportional `propU')
	}		
	else {
		local sigmagamma bscov(unstructured)
		local sigmabeta bscov(unstructured)

		if inlist(`"`structure'"', "exchangeable", "wscorrzero") {
			local sigmagamma bscov(exchangeable .5)
			
			if `"`structure'"'=="wscorrzero" {
				local sigmabeta bscov(equals `propU')	// matrix yet to be defined -- depends on SigmaGamma
			}
		}
		else if `"`structure'"'==`""' {
			sreturn local default default
			local structure unstructured
		}
	}

	if `"`naive'"'!=`""' & `"`structure'"'!=`"unstructured"' {
		nois disp as err `"Option {bf:naive} is only valid with {bf:unstructured} covariance structure"'
		exit 198
	}	
	
	sreturn local structure    `structure'
	sreturn local sigmagamma `"`sigmagamma'"'
	sreturn local sigmabeta  `"`sigmabeta'"'
	sreturn local naive      `"`naive'"'
	sreturn local options    `"`soptions'"'	

end



program define GetDesign
	syntax varname(string), study(varname) subgroup(varname)

	cap confirm variable _fillin
	if _rc {
		nois disp as err "variable {bf:_fillin} not found"
		exit _rc
	}
	qui levelsof `study', local(slist)
	
	// `subgroup' must either be numeric, or single-character string
	cap confirm string variable `subgroup'
	if !_rc local string string

	if `"`string'"'!=`""' {
		tempvar lensub
		gen `lensub' = length(`subgroup')
		summ `lensub', meanonly
		if r(max) > 1 {
			nois disp as err "With option {bf:design}, variable {bf:subgroup()} must be either numeric, or use single string characters only"
			exit 198
		}
	}
	
	foreach s of local slist {
		if `"`string'"'!=`""' {
			qui levelsof `subgroup' if `study'==`"`s'"' & !_fillin, clean
			qui replace `varlist' = `"`r(levels)'"' if `study'==`"`s'"'
		}
		else {
			qui levelsof `subgroup' if `study'==`s' & !_fillin
			qui replace `varlist' = `"`r(levels)'"' if `study'==`s'			
		}
	}
end
	
	
// Generate transformation matrix based on choice of reference subgroup
// ...and set up trend if appropriate
program define TransMatrix, sclass

	syntax namelist, subgroup(varname) [ noTRend * ]
	tokenize `namelist'
	args T d
	
	qui tab `subgroup'
	local k = r(r)
	local k1 = `k' - 1
	cap confirm string variable `subgroup'
	if !_rc local string string
	qui levelsof `subgroup', local(sglist)

	if `"`: char `subgroup'[omit]'"'==`""' {
		local r = 1
		gettoken ref sgnoref : sglist
		matrix define `T' = I(`k')
	}
	else {
		local ref : char `subgroup'[omit]
		local r : list posof `"`ref'"' in sglist
		cap assert `r' > 0
		if _rc {
			nois disp as err `"Reference value `ref' not taken on by variable {bf:`subgroup'} within the analysis sample"'
			exit 198
		}
		if `r'==1 matrix define `T' = I(`k')
		else {
			local kr = `k' - `r'
			local kr1 = `kr' + 1
			matrix define `T' = J(1, `k', 0)
			matrix define `T'[1, `r'] = 1
			if `r'==`k' {
				matrix define `T' = `T' \ ( I(`k1'), J(`k1', 1, 0))
			}
			else {
				tempname Id
				matrix define `Id' = I(`k1')
				matrix define `T' = `T' \ ( `Id'[1..`k1', 1..`kr'], J(`k1', 1, 0), `Id'[1..`k1', `kr1'..`k1'] )
			}
		}
	}
	if `"`string'"'!=`""' {
		qui numlist "1(1)`k'"
		local sglist `"`r(numlist)'"'
		local ref : copy local r
	}
	
	// Note:
	// `sglist' is an **ordered** (according to -levelsof- ) list of subgroup identifiers
	// `r' is the index w.r.t. `sglist'
	// `ref' is the relevant element of `sglist' (e.g. for labelling the matrix stripe elements)
	// (and, unless `string', `ref' will also be a numeric value, which may have a value label)
	
	* Setup test for trend
	if `"`trend'"'==`""' {
		forvalues i = 1 / `k' {
			matrix define `d' = nullmat(`d') \ `i'	// linear trend (could in theory allow more general forms for `d')
		}

		// Convert `d' into a matrix of contrasts
		// ... and apply transformation if reference is not first subgroup
		tempname M
		matrix define `M' = J(`k1', 1, -1), I(`k1')
		matrix define `d' = `M' * `T' * `d'
	}
	
	sreturn local sglist `"`sglist'"'
	sreturn local k `"`k'"'
	sreturn local ref  `"`ref'"'	
	sreturn local options `"`options'"'
	
end

	
program define cholesky2
* produce an approximate cholesky decomposition, even if matrix is not positive definite
* (Note: taken from -mvmeta- by Ian White)
args M C
cap matrix `C' = cholesky(`M')
if _rc {
    local eps = trace(`M')/1000
    if `eps'<=0 {
        di as error "cholesky2: matrix `M' has non-positive trace"
        matrix list `M'
        exit 498
    }
    local p = rowsof(`M')
    cap matrix `C' = cholesky(`M'+`eps'*I(`p'))
    if _rc {
        di as error "cholesky2: can't decompose matrix `M'"
        matrix list `M'
        exit 498
    }
}
end



program define SetupConstraints, sclass

	syntax namelist(min=2 max=2), k(integer) structure(string) [forcetrend]
	tokenize `namelist'
	args Chol propU
	
	local k1 = `k' - 1
	if `"`structure'"'!="wscorrzero" {
		local c = 0
		forvalues i = 1 / `k1' {
			forvalues j = 1 / `i' {
				local ++c
				
				if `"`structure'"'=="unstructured" & `"`forcetrend'"'==`""' {
					local el = `Chol'[1, `c']	// special case: see above
				}
				else local el = `Chol'[`i', `j']
				
				constraint free
				constraint `r(free)' [chol`j'`i']_b[_cons] = `el'
				local constr_list `constr_list' `r(free)'
			}
		}
	}
	
	// if exchangeable, need to do a little more work
	if `"`structure'"'=="exchangeable" {
		local tau = `Chol'[1, 1]
		tempname Chol_k
		matrix define `Chol_k' = inv(`Chol') * J(`k1', 1, -.5*(`tau')^2)
		forvalues j = 1 / `k1' {
			constraint free
			constraint `r(free)' [chol`j'`k']_b[_cons] = `=`Chol_k'[`j', 1]'
			local constr_list `constr_list' `r(free)'
		}
	}
	
	// similarly if trend
	else if `"`forcetrend'"'!=`""' {
		forvalues j = 1 / `k1' {
			constraint free
			constraint `r(free)' [chol`j'`k']_b[_cons] = 0
			local constr_list `constr_list' `r(free)'
		}
	}
	
	// finally, wscorrzero is a special case: SigmaBeta is completely defined by SigmaGamma
	else if `"`structure'"'=="wscorrzero" {
		local tau = `Chol'[1, 1]
		
		tempname LinvT		// = Linv * Linv-transpose, where `L' is as defined in the main routine
		matrix define `LinvT' = I(`k') + J(`k', `k', 1)
		matrix define `LinvT'[1, `k'] = J(`k1', 1, -1) 
		matrix define `LinvT'[`k', 1] = J(1, `k1', -1) 
		matrix define `LinvT'[`k', `k'] = 1

		matrix define `propU' = ((`tau')^2)*`LinvT'
	}
	
	// constraint 4 y_Isubgroup_2=0
	// constraint 5 y_Isubgroup_3=0

	sreturn local constr_list `"`constr_list'"'
	
end



program define ProcessSavedData

	syntax anything, ///
		study(name) subgroup(name) mainfile(string) ///
		[AUGVARiance(real 1e5) SAVING(string asis) ] 
	
	tokenize `anything'
	args GammaTable BetaTable
	
	local sgnoref : rownames `GammaTable'
	local sglist : rownames `BetaTable'
	local ref : list sglist - sgnoref
	local r : list posof `"`ref'"' in sglist	
	local k : word count `sglist'
	local k1 = `k' - 1
	
	// insert missing reference row into `GammaTable' so that its elements match with those of `BetaTable'
	if `r'==1 {
		matrix define `GammaTable' = (., .) \ `GammaTable'[1..`k1', 1..2]
	}
	else if `r'==`k' {
		matrix define `GammaTable' = `GammaTable'[1..`k1', 1..2] \ (., .)
	}
	else {
		local r1 = `r' - 1
		matrix define `GammaTable' = `GammaTable'[1..`r1', 1..2] \ (., .) \ `GammaTable'[`r'..`k1', 1..2]
	}
	matrix rownames `GammaTable' = `sglist'
	matrix colnames `GammaTable' = _yInt _yInt_S

	// identify y-vars and S-vars relating to the reference subgroup
	qui ds *_cons*
	local vlist `"`r(varlist)'"'
	local sgunab : variable label y_cons		// retrieve this, stored earlier

	// prepare for reshape long, by replacing "_cons" with reference value
	foreach x of local sglist {
		if `x'==`ref' {
			foreach v of local vlist {
				local newname = subinstr(`"`v'"', `"_cons"', `"_J`subgroup'_`x'"', .)
				qui rename `v' `newname'
			}
		}
		local rsvlist `rsvlist' S_J`subgroup'_@_J`subgroup'_`x'
	}
	qui reshape long y_J`subgroup'_ `rsvlist', i(`study') j(`subgroup')
	
	qui gen double S_J`subgroup'_ = .
	foreach x of local sglist {
		if `x'!=`ref' {
			qui replace S_J`subgroup'_ = S_J`subgroup'__J`subgroup'_`x' if `subgroup'==`x'
		}
		drop S_J`subgroup'__J`subgroup'_`x'
	}
	qui replace y_J`subgroup'_ = .  if `subgroup'==`ref'
	qui replace S_J`subgroup'_ = .  if `subgroup'==`ref'
	
	cap rename `subgroup' `sgunab'
	qui merge 1:1 `study' `sgunab' using `mainfile', assert(match) nogen
	
	local sglab : value label b
	if `"`sglab'"'!=`""' {
		label values `subgroup' `sglab'
		label values b
	}
	
	// add in pooled estimates
	tempvar pooled obs
	qui gen byte `pooled' = 1
	local oldN = _N
	local newN = `oldN' + `k'
	qui set obs `newN'
	qui replace `pooled' = 0 if _n > `oldN'
	
	foreach x of local sglist {
		local ++oldN
		qui replace `subgroup' = `x' in `oldN'
	}
	
	cap confirm numeric variable `study'
	if !_rc {
		summ `study', meanonly
		qui replace `study' = `r(max)' + 1 if !`pooled'
	}
	else {
		qui replace `study' = "_Pooled" if !`pooled'
	}
	
	// -svmat- places results starting at the first observation
	// so need to re-sort so that new observations are at the top, then re-sort again afterwards
	qui gen `obs' = _n
	sort `pooled' `obs'
	
	svmat `BetaTable', names(col)
	qui replace b = _y if !`pooled'
	qui replace V = _y_S if !`pooled'
	drop _y _y_S
	
	svmat `GammaTable', names(col)
	qui replace y_J`subgroup'_ = _yInt if !`pooled'
	qui replace S_J`subgroup'_ = _yInt_S if !`pooled'
	drop _yInt _yInt_S

	qui replace `pooled' = 1 - `pooled'
	sort `pooled' `obs'
	drop `obs'
	
	// use -metan- to restructure
	qui rename y_J`subgroup'_ _yInt_ES
	qui replace _yInt_ES = . if S_J`subgroup'_ >= `augvariance'
	qui gen double _yInt_seES = sqrt(S_J`subgroup'_) if S_J`subgroup'_ < `augvariance'
	drop S_J`subgroup'*
	
	qui replace b = . if V >= `augvariance'
	qui gen double stderr = sqrt(V) if V < `augvariance'
	
	qui ds `study' `subgroup' b V stderr, not
	local lcolvars `"`r(varlist)'"'
	
	qui metan b stderr, nograph nohet nosubgroup nooverall keeporder ///
		study(`subgroup') by(`study') rcols(`lcolvars') ///
		clear prefix(_y)

	qui drop *_EFFECT
	qui rename _y_BY _BY
	qui rename _y_STUDY _STUDY
	qui rename _y_LABELS _LABELS
	
	summ _BY if `pooled'==1, meanonly
	qui replace `pooled' = (_BY==`r(min)')
	
	qui replace _BY = . if `pooled'
	qui replace _LABELS = "Pooled estimates" if _y_USE==0 & `pooled'
	qui drop in L
	
	qui order _y_USE, before(_y_ES)
	qui replace _y_USE = 5 if inlist(_y_USE, 1, 2) & `pooled'

	summ _y_WT if !`pooled', meanonly
	qui replace _y_WT = 100*_y_WT / r(sum) if !`pooled'
	summ _y_WT if `pooled', meanonly
	qui replace _y_WT = 100*_y_WT / r(sum) if `pooled'
	
	qui gen double _yInt_LCI = _yInt_ES - invnormal(.975)*_yInt_seES
	qui gen double _yInt_UCI = _yInt_ES + invnormal(.975)*_yInt_seES

	qui gen double _yInt_WT = 1/_yInt_seES^2 if _y_USE!=5
	summ _yInt_WT, meanonly
	qui replace _yInt_WT = 100*_yInt_WT / r(sum)
	qui replace _yInt_WT = 100 if _y_USE==5 & _STUDY!=`ref'
	label variable _yInt_WT "% Weight"
	format _yInt_WT %6.2f
	
	qui gen byte _yInt_USE = _y_USE
	qui replace _yInt_USE = 9 if _STUDY==`ref'
	qui replace _yInt_USE = 2 if _yInt_USE==1 & missing(_yInt_seES)
	
	order _yInt_USE _yInt_ES _yInt_seES _yInt_LCI _yInt_UCI _yInt_WT, after(_y_WT)
	qui compress	
	
	if `"`saving'"'!=`""' {
		_prefix_saving `saving'
		qui save `s(filename)', `s(replace)'
	}
	
end
