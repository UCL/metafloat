*  v0.1 beta  David Fisher 10jan2021
*  v0.2 beta  David Fisher 30mar2021
*! v0.3 beta  David Fisher 21apr2021

// varlist (required) = effect size and std err
// study, subgroup (required) = study & subgroup identifiers
// default (no options) = unstructured random-effects
// fixed = all fixed (common) effects
// sigmagamma, sigmabeta = specify structures
// trend = fit linear trend across subgroups, if applicable


program define metafloat, eclass sortpreserve

	syntax varlist(numeric min=2 max=2) [if] [in], ///
		STUDY(varname) SUBGROUP(varname) ///
		[ FIXED RANDOMBeta EXCHangeable UNStructured SIGMAGamma(string) SIGMABeta(string) ///
		 TREND2 TREND(string) TREND AUGVARiance(real 1e5) SHOWmodels * ]

	marksample touse
	
	tokenize `varlist'
	args y stderr

	_get_eformopts, eformopts(`options')
	local eformopt eform(`"`s(str)'"')

	// check that varnames do not have prefix _I*
	cap assert `"`=substr(`"`study'"', 1, 2)'"'!=`"_I"'
	if _rc {
	    nois disp as err `"Cannot have {bf:study()} variable name beginning with {bf:_I}; please rename this variable"'
	}
	cap assert `"`=substr(`"`subgroup'"', 1, 2)'"'!=`"_I"'
	if _rc {
	    nois disp as err `"Cannot have {bf:subgroup()} variable name beginning with {bf:_I}; please rename this variable"'
	}
	
	// find number of subgroups
	qui tab `subgroup' if `touse'
	local k = r(r)	
	if "`trend'"!="" {
		if `k'<3 {
			nois disp `"{error}Cannot test for trend with fewer than 3 subgroups"'
			local trend
		}
		else local trend commonparm
	}

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
	
	// setup tests for trend
	if `"`trend'"'!=`""' & `"`trend2'"'!=`""' {
		nois disp as err `"only one of {bf:trend} or {bf:trend()} is allowed"'
		exit 184
	}
	if "`trend2'"!="" {
		if `k'<3 {
			nois disp `"{error}Cannot test for trend with fewer than 3 subgroups"'
			local trend
		}
		else {
			tempname d
			forvalues i = 1 / `k' {
				mat `d' = nullmat(`d') \ `i'	// linear trend
			}
			local trend commonparm noconstant
			local testonly testonly
		}
	}
	else if "`trend'"!="" {
		if `k'<3 {
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
				assert rowsof(`d') = `k'
				assert colsof(`d') = 1
			}
			if _rc {
				nois disp as err "matrix (`anything') not found, or has invalid dimensions"
				exit 198
			}
			local trend commonparm noconstant
		}
	}
	else local trend
	
	// preserve data, prior to editing
	preserve

	qui keep if `touse' & !missing(`study', `subgroup')
	keep `y' `stderr' `study' `subgroup'
	
	tempvar V
	qui gen double `V' = `stderr'^2				// Generate estimate of the variance (s-squared)
	drop `stderr'

	// augment if missing subgroup data
	qui fillin `study' `subgroup'
	summ _fillin, meanonly
	if r(sum) {
		nois disp as text "The following subgroups are missing, and will be given augmented variances of " as res `augvariance' as text ":"
		list `study' `subgroup' if _fillin
	}
	drop _fillin
	qui replace `V' = `augvariance' if missing(`y')
	qui replace `y' = 0   if missing(`y')
	
	// store Wi matrices, for later
	qui levelsof `study', local(slist)
	local i = 0
	foreach s of local slist {
		local ++i
		tempname W`i'
		mkmat `V' if `study'==`s', matrix(`W`i'')
		cap assert rowsof(`W`i'')==`k'
		if _rc {
			nois disp as err `"Conformability error: study `s' has `=rowsof(`W`i'')' subgroups where it should have `k'"'
			exit 503
		}
		matrix `W`i'' = diag(`W`i'')
	}

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
	if `"`: char `subgroup'[omit]'"'==`""' {
	    local refgroup : word 1 of `colnames'
		local i = 1
		local c = substr(`"`refgroup'"', -`i', 1)
		while `"`c'"'!=`"_"' {
			local ++i
			local c = substr(`"`refgroup'"', -`i', 1)
		}
		local --i
		local ref = substr(`"`refgroup'"', -`i', `i')
	}
	else {
		// DF NOTE: temporary ifstmt
	    if `"`trend'"'!=`""' {
		    nois disp as err "Trend option does not currently work properly with non-standard reference category"
			exit 198
		}
		
	    local ref : char `subgroup'[omit]
	    qui ds _I*_`ref'
		local refgroup `"`r(varlist)'"'
		assert `: word count `refgroup''==1
	}
	
	// re-order `colnames', placing reference first
	local colnamesnotref : list colnames - refgroup
	local colnames `refgroup' `colnamesnotref'
	
	// generate contrasts (use prefix "_Int" = "interaction")
	// Note: _y_cons is the reference subgroup (i.e NOT a contrast)
	qui xi, prefix(_Int) : mvmeta_make regress `y' i.`subgroup' [iw=`V'^-1], mse1 by(`study') clear nodetails useconstant
	// N.B. regress ... [iw], mse1  is equivalent to using -vwls- 
	//  but currently -mvmeta_make- and -vwls- don't like each other for some reason

	
	** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
	if "`showmodels'"=="" local quietly quietly
	`quietly' mvmeta y S, vars(y_Int*) `fixed' `sigmagamma' `trend'
		
	tempname VarGammaHat GammaHat SigmaGamma
	local k1 = `k' - 1
	mat `SigmaGamma'  = e(Sigma)
	
	if "`fixed'"=="" {
		tempname Chol
		mat `Chol' = e(b)[1, `k'...]
	}
	if "`trend'"!="" {		// DF Note: Need to look again at this
		nois disp as text _n "Test for trend:"
		test _cons
			
		tempname integers
		forvalues i = 1 / `k1' {
			mat `integers' = nullmat(`integers') \ `i'
		}
		mat `VarGammaHat' = e(V)*`integers''*`integers'
		mat `GammaHat'    = e(b)*`integers'
	}
	else {
		nois disp as text _n "Test of interactions(s):"
		testparm y*
		mat `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		mat `GammaHat'    = e(b)[1,       1..`k1']'			
	}

	
	** STEP 2: ESTIMATE TREATMENT EFFECT IN REF GROUP (MEAN, BS VARIANCE AND BS CORRELATION WITH INTERACTIONS)
	// subtract within-trial interactions from observed interaction data
	// set up equations that make shifted interactions have mean zero
	tempname AugGammaHat gMat
	
	qui gen byte zero = 0
	qui gen byte one  = 1
	
	local i = 0
	tokenize `e(yvars)'
	while `"`1'"'!=`""' {
	    local ++i
		local g = `GammaHat'[`i', 1]
		qui replace `1' = `1' - `g'
		local eq `eq' `1':zero,
		
		// AugGammaHat (for later)
		if `: word `i' of `sgvals''==`ref' mat `gMat' = 0 \ `g'
		else mat `gMat' = `g'
		mat `AugGammaHat' = nullmat(`AugGammaHat') \ `gMat'

		macro shift
	}
	if `: word `k' of `sgvals''==`ref' {
		mat `AugGammaHat' = `AugGammaHat' \ 0
	}
	
	local eq eq(`eq' y_cons:one)
	
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
	`quietly' mvmeta_new y S, `eq' nocons commonparm print(bscov) `constr_opt' `fixed' `sigmabeta'
	// Note: the coefficient of _cons refers to the reference subgroup

	// extract fitted BS variance and derive SigmaBeta
	tempname L SigmaU SigmaBeta
	mat `L' = J(`k', `k1', 0), J(`k', 1, 1)
	mat `L'[2, 1] = I(`k1')
	mat `SigmaU' = e(Sigma)
	mat `SigmaBeta' = `L' * `SigmaU' * `L''

	// compute implied means for beta
	tempname ones ThetaHat VarThetaHat BetaHat
	mat `ones' = J(`k', 1, 1)			// column vector of ones, of length k
	mat `ThetaHat'    = e(b)[1,1]
	mat `VarThetaHat' = e(V)[1,1]
	// mat `BetaHat' = `ThetaHat'*`ones' + (0 \ `GammaHat')
	// ^^ (0 \ `GammaHat') reflects the fact that the reference comes first
	mat `BetaHat' = `ThetaHat'*`ones' + `AugGammaHat'
	
	// correct variance for beta
	tempname D GradSum Grad
	mat `D' = J(`k', `k1', 0)
	mat `D'[2, 1] = I(`k1'')
	mat `GradSum' = J(1, `k1', 0)
	forvalues i = 1 / `t' {
		mat `GradSum' = `GradSum' + `ones''*invsym(`W`i'' + `SigmaBeta')*`D'
	}
	mat `Grad' = -`VarThetaHat'*`GradSum'
	
	tempname A VarBetaHat
	mat `A' = `ones'*`Grad' + `D'
	mat `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A' * `VarGammaHat' * `A''
	
	// collect and display under "ereturn"
	matrix `BetaHat' = `BetaHat''
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
