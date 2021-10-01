*  v0.1 beta  David Fisher 10jan2021
*  v0.2 beta  David Fisher 30mar2021
*  v0.3 beta  David Fisher 21apr2021
*! v0.4 beta  David Fisher 07may2021

// varlist (required) = effect size and std err
// study, subgroup (required) = study & subgroup identifiers
// default (no options) = unstructured random-effects
// fixed = all fixed (common) effects
// sigmagamma, sigmabeta = specify structures
// trend = fit linear trend across subgroups, if applicable


// Note: under fixed-effects, all works fine
//  but with unstructured random-effects, changing reference has an effect -- need to investigate


program define metafloat, eclass sortpreserve

	syntax varlist(numeric min=2 max=2) [if] [in], ///
		STUDY(varname) SUBGROUP(varname) ///
		[ FIXED RANDOMBeta EXCHangeable UNStructured SIGMAGamma(string) SIGMABeta(string) ///
		 TREND2 TREND(string) AUGVARiance(real 1e5) SHOWmodels noUNCertainv * ]

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
			matrix `T'[1, `r'] = 1
			if `r'==`k' {
				matrix `T' = `T' \ ( I(`k1'), J(`k1', 1, 0))
			}
			else {
				tempname Id
				matrix define `Id' = I(`k1')
				matrix `T' = `T' \ ( `Id'[1..`k1', 1..`kr'], J(`k1', 1, 0), `Id'[1..`k1', `kr1'..`k1'] )
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
		mkmat `V' if `study'==`s', matrix(`W`i'')
		cap assert rowsof(`W`i'')==`k'
		if _rc {
			nois disp as err `"Conformability error: study `s' has `=rowsof(`W`i'')' subgroups where it should have `k'"'
			exit 503
		}
		matrix `W`i'' = diag(`T'*`W`i'')
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
			forvalues i = 1 / `k1' {
				matrix `d' = nullmat(`d') \ `i'	// linear trend
			}
			// local trend2 commonparm /*noconstant*/
			local testonly testonly
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
				assert rowsof(`d')==`k1'
				assert colsof(`d')==1
			}
			if _rc {
				nois disp as err "matrix (`anything') not found, or has invalid dimensions"
				nois disp as err "(should be a column vector of length equal to one fewer than the number of subgroups)"
				exit 198
			}			
			// local trend commonparm /*noconstant*/

		}
		// matrix `d' = `T' * `d'		// NO!!  Revisit this.
	}
	/*
	else {
		local trend
		local trend2
	}
	*/
	
	// generate contrasts (use prefix "_Int" = "interaction")
	// Note: _y_cons is thereference subgroup (i.e NOT a contrast)
	qui xi, prefix(_Int) : mvmeta_make regress `y' i.`subgroup' [iw=`V'^-1], mse1 by(`study') clear nodetails useconstant
	// N.B. regress ... [iw], mse1  is equivalent to using -vwls- 
	//  but currently -mvmeta_make- and -vwls- don't like each other for some reason

	
	** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
	if "`showmodels'"=="" local quietly quietly
	quietly mvmeta y S, vars(y_Int*) `fixed' `sigmagamma' /* `trend' `trend2' */
	
	/*
	tempname integers
	forvalues i = 1 / `k1' {
		matrix `integers' = nullmat(`integers') \ `i'
	}
	matrix `VarGammaHat' = e(V)*`integers''*`integers'
	matrix `GammaHat'    = e(b)*`integers'
	*/

	tempname VarGammaHat GammaHat SigmaGamma
	if `"`trend'`trend2'"'!=`""' {
		local eq
		local i = 0
		tokenize `e(yvars)'
		while `"`1'"'!=`""' {
			local ++i
			tempvar x`i'
			qui gen `x`i'' = `d'[`i', 1]
			local eq `eq' `1':`x`i'',
			macro shift
		}	
		`quietly' mvmeta y S, vars(y_Int*) `fixed' `sigmagamma' commonparm nocons eq(`eq')
	
		if `"`trend'"'!=`""' {
			matrix `VarGammaHat' = e(V)*`d''*`d'
			matrix `GammaHat'    = e(b)*`d'
			matrix `SigmaGamma'  = e(Sigma)
		}
		else {		// if `trend2' , *test* for (linear) trend only;  now re-fit as standard model assuming no trend
			quietly mvmeta y S, vars(y_Int*) `fixed' `sigmagamma'
		}
	
		nois disp as text _n "Test for trend:"
		// test _cons
		test `x1'
	}
	else {
		matrix `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		matrix `GammaHat'    = e(b)[1,       1..`k1']'		
		matrix `SigmaGamma'  = e(Sigma)

		nois disp as text _n "Test of interactions(s):"
		testparm y*
	}
	
	// DF: revisit this bit; what to do if we *do* have `trend' ??
	if `"`trend'"'==`""' {		
		if `"`fixed'"'==`""' {
			tempname Chol
			matrix `Chol' = e(b)[1, `k'...]
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
	qui gen byte zero = 0
	qui gen byte one  = 1
	
	local eq
	local i = 0
	tokenize `e(yvars)'
	while `"`1'"'!=`""' {
	    local ++i
		qui replace `1' = `1' - `GammaHat'[`i', 1]
		local eq `eq' `1':zero,
		macro shift
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
	matrix define `L' = J(`k', `k1', 0), J(`k', 1, 1)
	matrix `L'[2, 1] = I(`k1')
	matrix define `SigmaU' = e(Sigma)
	matrix define `SigmaBeta' = `L'*`SigmaU'*`L''
	
	// compute implied means for beta
	tempname ones ThetaHat VarThetaHat BetaHat
	matrix define `ones' = J(`k', 1, 1)			// column vector of ones, of length k
	matrix define `ThetaHat'    = e(b)[1,1]
	// matrix define `VarThetaHat' = e(V)[1,1]
	matrix define `BetaHat' = `ThetaHat'*`ones' + (0 \ `GammaHat')
	// ^^ (0 \ `GammaHat') reflects the fact that the reference comes first
	
	matrix define `VarThetaHat' = e(V)

	// implement -nouncertainv- option
	if `"`uncertainv'"'!=`""' {
		matrix `VarThetaHat' = invsym(`VarThetaHat')
		matrix `VarThetaHat' = `VarThetaHat'[1,1]
		matrix `VarThetaHat' = invsym(`VarThetaHat')
	}
	else {
		matrix `VarThetaHat' = `VarThetaHat'[1,1]
	}
	
	// correct variance for beta
	tempname D GradSum Grad
	matrix define `D' = J(`k', `k1', 0)
	matrix `D'[2, 1] = I(`k1')
	matrix define `GradSum' = J(1, `k1', 0)
	forvalues i = 1 / `t' {
		matrix `GradSum' = `GradSum' + `ones''*invsym(`W`i'' + `SigmaBeta')*`D'
	}
	matrix define `Grad' = -`VarThetaHat'*`GradSum'
	matrix rownames `Grad' = `: word 1 of `rownames_eb''
	matrix colnames `Grad' = `colnames_eb2'
	
	tempname A VarBetaHat
	matrix define `A' = `ones'*`Grad' + `D'
	matrix define `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A'*`VarGammaHat'*`A''
	
	// reverse transformation
	matrix `BetaHat' = `T''*`BetaHat'
	matrix `VarBetaHat' = `T''*`VarBetaHat'*`T'

	// ...also applies to SigmaBeta
	matrix `SigmaBeta' = `T''*`SigmaBeta'*`T'
	matrix coleq `SigmaBeta' = ""
	matname `SigmaBeta' `colnames', explicit
	
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
