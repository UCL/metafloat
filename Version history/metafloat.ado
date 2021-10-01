*  v0.1 beta  David Fisher 10jan2021
*  v0.2 beta  David Fisher 30mar2021
*  v0.3 beta  David Fisher 21apr2021
*  v0.4 beta  David Fisher 07may2021
*  v0.5 beta  David Fisher 01jul2021
*! v0.6 beta  David Fisher 12jul2021


*** METAFLOAT ***
* Routine for fitting the models described in Godolphin et al (2021)
// by David Fisher, Peter Godolphin, Ian White
// (MRC Clinical Trials Unit at UCL)

* Release notes:
// v0.5: following discusion between IW, PG, corrected code so that results remain invariant to change of reference
//         under random-effects, with or without `nouncertainv'
// v0.6: addition of `design' option; analogous to inconsistency parameters in NMA
//       observations with extreme variances now set to be augmented

* Syntax:
// varlist (required) = effect size and std err
// study, subgroup (required) = study & subgroup identifiers

// default (no options) = unstructured random-effects
// fixed = all fixed (common) effects
// sigmagamma, sigmabeta = specify structures
// trend = test for linear trend across subgroups (if applicable);
//         note: presented subgroup estimates do *not* assume linear trend (i.e. trend is only *tested for*, not implemented)
//         note: actual internal option name is `trend2'
// trend(mat [, testonly]) = fit specific trend defined by matrix "mat"
//         testonly = only *test* for this trend, do not present subgroup estimates assuming trend
// noaugtrend = don't use variance augmentation for estimating trend
// design = additional parameters in final model describing the available subgroups per trial (e.g. "single-subgroup" trials)

// nouncertainv = -mvmeta- option; see mvmeta documentation
// `options' = eform options, plus other options to pass to -mvmeta- e.g. tolerance()


program define metafloat, eclass sortpreserve

	syntax varlist(numeric min=2 max=2) [if] [in], ///
		STUDY(varname) SUBGROUP(varname) ///
		[ FIXED RANDOMBeta EXCHangeable UNStructured SIGMAGamma(string) SIGMABeta(string) ///
		 TREND2 TREND(string) AUGVARiance(real 1e5) noAUGTRend SHOWmodels noUNCertainv DESign * ]

	marksample touse
	
	tokenize `varlist'
	args y stderr
	
	cap assert `stderr' >= 0 if `touse'
	if _rc {
		nois disp as err "Standard error cannot be negative"
		exit 198
	}
	
	_get_eformopts, eformopts(`options') allowed(__all__) soptions
	local eformopt eform(`"`s(str)'"')
	local soptions `"`s(options)'"'

	
	** Check for potential conflicts in varnames:
	// check that varnames do not have prefix _I*
	foreach x in study subgroup {
		cap assert `"`=substr(`"``x''"', 1, 2)'"'!=`"_I"'
		if _rc {
			nois disp as err `"Cannot have {bf:`x'()} variable name beginning with {bf:_I}; please rename this variable"'
		}
		if "`design'"!="" {
			cap assert `"`=substr(`"``x''"', 1, 8)'"'!=`"_Design"'
			if _rc {
				nois disp as err `"Cannot have {bf:`x'()} variable name beginning with {bf:_Design}; please rename this variable"'
			}
		}
	}
	
	// check that varnames are not called _fillin, _Trend, _One, _Zero
	foreach v1 in study subgroup {
		foreach v2 in fillin Trend One Zero {
			cap assert `"``v1''"'!=`"_`v2'"'
			if _rc {
				nois disp as err `"Cannot have {bf:`v1'()} variable name equal to {bf:_`v2'}; please rename this variable"'
			}
		}
	}
	
	// find number of subgroups
	qui tab `subgroup' if `touse'
	local k = r(r)	
	local k1 = `k' - 1

	// find number of studies
	qui tab `study' if `touse'
	local t = r(r)
	
	// check that study-subgroup combinations are unique
	tempvar Ngroup
	qui bysort `touse' `study' `subgroup' (`obs') : gen int `Ngroup' = _N
	summ `Ngroup' if `touse', meanonly
	if r(max) > 1 {
		nois disp as err "The following combination(s) of study and subgroup are not unique:"
		nois list `study' `subgroup' `y' `stderr' if `Ngroup' > 1 & `touse'
		exit 198
	}
	
	// setup fixed/random and covariance structures
	opts_exclusive `"`fixed' `randombeta' `exchangeable' `unstructured'"' `""' 184
	local structure `fixed'`randombeta'`exchangeable'`unstructured'
	if "`structure'"!="" {
	    if "`sigmagamma'`sigmabeta'"!="" {
		    nois disp as err `"Covariance structure {bf:`structure'} specified; options {bf:sigmagamma()} and {bf:sigmabeta()} are invalid"'
			exit 198
		}
	}
	
	
	** Preserve data, prior to editing
	preserve

	qui keep if `touse' & !missing(`study', `subgroup')
	keep `y' `stderr' `study' `subgroup'
	
	tempvar V
	qui gen double `V' = `stderr'^2				// Generate estimate of the variance (s-squared)
	drop `stderr'
	
	// augment if missing subgroup data
	qui fillin `study' `subgroup'
	assert _fillin == missing(`y')
	
	// if any observed variances are larger than `augvariance'
	//  (e.g. if estimated from very small sample)
	// also set to _fillin
	qui replace _fillin = 1 if `V' > `augvariance'

	summ _fillin, meanonly
	if r(sum) {
		nois disp as text "The following subgroups are unobserved (or contain missing or extremely imprecise data)"
		nois disp as text " and will be given augmented variances of " as res `augvariance' as text ":"
		list `study' `subgroup' if _fillin
	}
	qui replace `V' = `augvariance' if _fillin
	qui replace `y' = 0 if _fillin

	// obtain "design", i.e. which subgroups are present
	if "`design'"!="" {
	    
		// `subgroup' must either be numeric, or single-character string
		cap confirm numeric variable `subgroup'
		if _rc {
		    tempvar lensub
			gen `lensub' = length("`subgroup'")
			summ `lensub', meanonly
			if r(max) > 1 {
			    nois disp as err "With option {bf:design}, variable {bf:subgroup()} must be either numeric, or use single string characters only"
				exit 198
			}
		}
		
		qui gen _Design = ""
		qui levelsof `study', local(slist)
		foreach s of local slist {
			cap confirm string variable `study'
			if _rc {
				qui levelsof `subgroup' if `study'==`s' & !_fillin
				qui replace _Design = `"`r(levels)'"' if `study'==`s'			
			}
			else {
				qui levelsof `subgroup' if `study'==`"`s'"' & !_fillin
				qui replace _Design = `"`r(levels)'"' if `study'==`"`s'"'
			}
		}
		local design_opt collapse((firstnm) _Design)
		local showmodels showmodels
	}
	drop _fillin
	
	// obtain subgroup names (for attaching to final matrices)
	// NOTE: already ensured that no vars `y', `V', `study', `subgroup' have prefix _I*]
	cap drop _I*
	xi, noomit i.`subgroup'
	qui ds _I*
	local colnames `"`r(varlist)'"'
	
	qui levelsof `subgroup', local(sgvals)
	cap assert `: word count `sgvals''==`: word count `colnames''
	if _rc {
		nois disp as err `"Error in identification of subgroups"'
		exit 198
	}	
	
	// identify reference subgroup
	tempname T
	if `"`: char `subgroup'[omit]'"'==`""' {
	    local ref : word 1 of `sgvals'
		local r = 1
		matrix define `T' = I(`k')
	}
	
	// ... and generate transformation matrix if reference is not first subgroup
	else {
	    local ref : char `subgroup'[omit]
	    qui ds _I*_`ref'
		assert `: word count `r(varlist)''==1
		local r : list posof `"`ref'"' in sgvals

		if `r'==1 {
			matrix define `T' = I(`k')
		}
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
	
	// store Wi matrices, for later
	// (with transformation applied)
	qui levelsof `study', local(slist)
	local i = 0
	foreach s of local slist {
		local ++i
		tempname W`i'
		
		cap confirm string variable `study'
		if _rc {	
			mkmat `V' if `study'==`s', matrix(`W`i'')
		}
		else {
			mkmat `V' if `study'==`"`s'"', matrix(`W`i'')			
		}
		
		cap assert rowsof(`W`i'')==`k'
		if _rc {
			nois disp as err `"Conformability error: study `s' has `=rowsof(`W`i'')' subgroups where it should have `k'"'
			exit 503
		}
		matrix define `W`i'' = diag(`T'*`W`i'')
	}
	
	
	// setup tests for trend
	if `"`trend'"'!=`""' & `"`trend2'"'!=`""' {
		nois disp as err `"only one of {bf:trend} or {bf:trend()} is allowed"'
		exit 184
	}
	if "`trend2'"!="" {
		if `k' < 3 {
			nois disp `"{error}Cannot test for trend with fewer than 3 subgroups"'
			local trend2
		}
		else {
			tempname d
			forvalues i = 1 / `k' {
				matrix define `d' = nullmat(`d') \ `i'	// linear trend
			}
			local testonly testonly
			
			// Convert `d' into a matrix of contrasts
			// ... and apply transformation if reference is not first subgroup
			tempname M
			matrix define `M' = J(`k1', 1, -1), I(`k1')
			matrix define `d' = `M' * `T' * `d'
		}
	}
	else if "`trend'"!="" {
		if `k' < 3 {
			nois disp `"{error}Cannot test for trend with fewer than 3 subgroups"'
			local trend
		}
		else {
			local 0 `"`trend'"'
			cap syntax anything [, TESTONLY ]
			
			* If -syntax- didn't work, assume due to commas but no brackets
			if _rc {
				disp as err "invalid syntax for option {bf:trend()}"
				disp as err "Valid syntaxes are:"
				disp as err `" the name of an existing matrix; or"'
				disp as err `" matrix-style input {bf:(} {it:a}{bf:,} {it:b} {bf:\} {it:c}{bf:,} {it:d} {bf:)}; see {help matrix define}"'
				exit _rc
			}
			tempname d
			cap {
				matrix define `d' = `anything'
				assert rowsof(`d')==`k'
				assert colsof(`d')==1
			}
			if _rc {
				nois disp as err "matrix (`anything') not found, or has invalid dimensions"
				nois disp as err "(should be a column vector of length equal to the number of subgroups)"
				exit 198
			}
			
			// Convert `d' into a matrix of contrasts
			// ... and apply transformation if reference is not first subgroup
			tempname M
			matrix define `M' = J(`k1', 1, -1), I(`k1')
			matrix define `d' = `M' * `T' * `d'
		}
	}
	
	// generate contrasts (use prefix "_Int" = "interaction")
	// Note: _y_cons is thereference subgroup (i.e NOT a contrast)
	qui xi, prefix(_Int) : mvmeta_make regress `y' i.`subgroup' [iw=`V'^-1], mse1 by(`study') `design_opt' clear nodetails useconstant
	// N.B. regress ... [iw], mse1  is equivalent to using -vwls- 
	//  but currently -mvmeta_make- and -vwls- don't like each other for some reason
	
	
	** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
	if "`showmodels'"=="" local quietly quietly

	tempname VarGammaHat GammaHat SigmaGamma
	if `"`trend'`trend2'"'!=`""' {
		local eq
		local i = 0
		local yvars : colnames e(b) 
		
		tokenize `yvars'
		while `"`1'"'!=`""' {
			local ++i
			qui gen _Trend`i' = `d'[`i', 1]
			local eq `eq' y`1':_Trend`i',
			
			if "`augtrend'"!="" {
				tempvar y`1'
				qui clonevar `y`1'' = y`1'
			}
			
			macro shift
		}

		// optional request to use "exact" data (no variance augmentation) for estimating trend
		if "`augtrend'"!="" {
			foreach v1 of local yvars {
				foreach v2 of local yvars {
					cap confirm variable S`v1'`v2'
					if !_rc {
						qui replace y`v1' = . if abs(S`v1'`v1') >= `augvariance'
						qui replace y`v2' = . if abs(S`v2'`v2') >= `augvariance'

						tempvar S`v1'`v2'
						qui clonevar `S`v1'`v2'' = S`v1'`v2'
						qui replace S`v1'`v2' = . if abs(S`v1'`v2') >= `augvariance'
					}
				}
			}
		}
		
		nois disp as text _n "Test for trend:"
		`quietly' mvmeta y S, vars(y_Int*) `fixed' `sigmagamma' commonparm nocons eq(`eq') `soptions'
		test _Trend1

		if "`augtrend'"!="" {
			foreach v1 of local yvars {
				qui replace y`v1' = `y`v1''
				drop `y`v1''
				
				foreach v2 of local yvars {
					cap confirm variable S`v1'`v2'
					if !_rc {
						qui replace S`v1'`v2' = `S`v1'`v2''
						drop `S`v1'`v2''
					}
				}
			}
		}		
		
		if `"`trend'"'!=`""' {
			matrix define `VarGammaHat' = e(V)*`d''*`d'
			matrix define `GammaHat'    = e(b)*`d'
			matrix define `SigmaGamma'  = e(Sigma)
		}
		
		// if `trend2' , *test* for (linear) trend only;  now re-fit as standard model assuming no trend
		else {
			`quietly' mvmeta y S, vars(y_Int*) `fixed' `sigmagamma' `soptions'
			matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
			matrix define `GammaHat'    = e(b)[1,       1..`k1']'
			matrix define `SigmaGamma'  = e(Sigma)
		}
	}
	
	// default: test for individual contrasts
	else {
		nois disp as text _n "Test of interaction(s):"
		`quietly' mvmeta y S, vars(y_Int*) `fixed' `sigmagamma' `soptions'
		nois testparm y*

		matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		matrix define `GammaHat'    = e(b)[1,       1..`k1']'		
		matrix define `SigmaGamma'  = e(Sigma)
	}
	
	// DF: revisit this bit; what to do if we *do* have `trend' ??
	if `"`trend'"'==`""' {		
		if `"`fixed'"'==`""' {
			tempname Chol
			matrix define `Chol' = e(b)[1, `k'...]
		}
	
		local colnames_eV : colnames e(V)
		local colnames_eb : colnames e(b)
		local rownames_eb : rownames e(b)
		
		forvalues i = 1 / `k1' {
			local colnames_eV2 `colnames_eV2' `: word `i' of `colnames_eV''
			local colnames_eb2 `colnames_eb2' `: word `i' of `colnames_eb''
		}
		
		matname `VarGammaHat' `colnames_eV2', explicit
		matname `SigmaGamma'  `colnames_eV2', explicit
		matrix colnames `GammaHat' = `: word 1 of `rownames_eb''
		matrix rownames `GammaHat' = `colnames_eb2'
	}

	
	** STEP 2: ESTIMATE TREATMENT EFFECT IN REF GROUP (MEAN, BS VARIANCE AND BS CORRELATION WITH INTERACTIONS)
	// subtract within-trial interactions from observed interaction data
	// set up equations that make shifted interactions have mean zero
	
	qui gen byte _Zero = 0
	qui gen byte _One  = 1

	// DF July 2021: optionally adjust for "design"
	if "`design'"!="" {
	    tempvar obs
		gen long `obs' = _n
		local designvars
		qui tab _Design, gen(_Design) sort
		forvalues i = 2 / `r(r)' {
		    summ `obs' if _Design`i'==1, meanonly
			local destxt = _Design[`r(min)']
			local destxt = subinstr(trim(`"`destxt'"'), " ", "_", .)
			rename _Design`i' _Des_`destxt'
			local designvars `designvars' _Des_`destxt'
		}
		drop _Design1
	}		
		
	local eq
	local i = 0
	tokenize `e(yvars)'
	while `"`1'"'!=`""' {
	    local ++i
		qui replace `1' = `1' - `GammaHat'[`i', 1]
		local eq `eq' `1':_Zero `designvars',
		// local eq `eq' `1':_Zero,
		macro shift
	}	
	local eq eq(`eq' y_cons:_One `designvars')
	// local eq eq(`eq' y_cons:_One `designvars')
	
	cap assert `i'==`k'-1
	if _rc {
		nois disp as err "Error detected when subtracting contrasts from non-reference subgroup"
		exit 198
	}
	
	// set up constraints on the BS variance of the interactions
	// DF Note: Need to look at this again
	if "`fixed'"=="" {
		constraint 1 [chol11]_b[_cons] = `=`Chol'[1,1]'
		constraint 2 [chol12]_b[_cons] = `=`Chol'[1,2]'
		constraint 3 [chol22]_b[_cons] = `=`Chol'[1,3]'
		// constraint 4 y_Isubgroup_2=0
		// constraint 5 y_Isubgroup_3=0
		// constraint dir
		local constr_opt constraints(1 2 3)
	}
	
	// fit model: make sure we are using "NEW" mvmeta
	if "`showmodels'"!="" nois disp as text _n "Estimation of floating subgroup for reference category:"
	`quietly' mvmeta_new y S, `eq' nocons commonparm print(bscov) `constr_opt' `fixed' `sigmabeta' `soptions'
	// Note: the coefficient of _cons refers to the reference subgroup
	
	// extract fitted BS variance and derive SigmaBeta
	tempname L SigmaU SigmaBeta
	matrix define `L' = J(`k', `k1', 0), J(`k', 1, 1)
	matrix define `L'[2, 1] = I(`k1')
	matrix define `SigmaU' = e(Sigma)
	matrix define `SigmaBeta' = `L'*`SigmaU'*`L''
	
	// compute implied means for beta
	tempname ones ThetaHat BetaHat
	matrix define `ones' = J(`k', 1, 1)			// column vector of ones, of length k
	matrix define `ThetaHat' = e(b)[1,1]
	matrix define `BetaHat' = `ThetaHat'*`ones' + (0 \ `GammaHat')
	// ^^ (0 \ `GammaHat') reflects the fact that the reference comes first
	
	// IW June 2021: First derive `VarThetaHat' ignoring the uncertainty in the heterogeneity matrix
	// because it's needed for deriving the gradient vector
	tempname VarThetaHat
	matrix define `VarThetaHat' = e(V)
	matrix define `VarThetaHat' = invsym(`VarThetaHat')
	matrix define `VarThetaHat' = `VarThetaHat'[1,1]
	matrix define `VarThetaHat' = invsym(`VarThetaHat')	

	tempname D GradSum Grad
	matrix define `D' = J(`k', `k1', 0)
	matrix define `D'[2, 1] = I(`k1')
	matrix define `GradSum' = J(1, `k1', 0)
	forvalues i = 1 / `t' {
		matrix define `GradSum' = `GradSum' + `ones''*invsym(`W`i'' + `SigmaBeta')*`D'
	}
	matrix define `Grad' = -`VarThetaHat'*`GradSum'
	matrix rownames `Grad' = `: word 1 of `rownames_eb''
	matrix colnames `Grad' = `colnames_eb2'
	
	// now correct VarThetaHat for uncertainty in the heterogeneity matrix, if requested
	if `"`uncertainv'"'==`""' {
		matrix define `VarThetaHat' = e(V)[1,1]
	}
		
	// correct variance for beta
	tempname A VarBetaHat
	matrix define `A' = `ones'*`Grad' + `D'
	matrix define `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A'*`VarGammaHat'*`A''
	
	// reverse transformation
	matrix define `BetaHat'    = `T''*`BetaHat'
	matrix define `VarBetaHat' = `T''*`VarBetaHat'*`T'

	// ...also applies to SigmaBeta
	matrix define `SigmaBeta' = `T''*`SigmaBeta'*`T'
	matrix coleq `SigmaBeta' = ""
	matname `SigmaBeta' `colnames', explicit
	
	// collect and display under "ereturn"
	matrix define `BetaHat' = `BetaHat''
	matrix coleq `BetaHat' = ""
	matname `BetaHat'    `colnames', columns(.) explicit
	matname `VarBetaHat' `colnames', explicit

	ereturn post `BetaHat' `VarBetaHat', depname(Subgroup)
	ereturn matrix GammaHat    = `GammaHat'
	ereturn matrix VarGammaHat = `VarGammaHat'
	ereturn matrix ThetaHat    = `ThetaHat'
	ereturn matrix VarThetaHat = `VarThetaHat'
	ereturn matrix Grad        = `Grad'
	ereturn matrix SigmaGamma  = `SigmaGamma'
	ereturn matrix SigmaBeta   = `SigmaBeta'
	ereturn scalar n = `t'						// number of studies
	ereturn scalar k = `k'						// number of subgroups
	
	nois disp as text _n "Floating subgroups:"
	ereturn display, `eformopt'

end
