*  v0.1 beta  David Fisher 10jan2021
*! v0.2 beta  David Fisher 30mar2021

// varlist (required) = effect size and std err
// study, subgroup (required) = study & subgroup identifiers
// default (no options) = unstructured random-effects
// fixed = all fixed (common) effects
// sigmagamma, sigmabeta = specify structures
// trend = fit linear trend across subgroups, if applicable


program define metafloat, eclass sortpreserve

	syntax varlist(numeric min=2 max=2) [if] [in], ///
		STUDY(varname) SUBGROUP(varname) [ Base(string) IBase(string) REFerence(string) ///
		FIXED RANDOMBeta EXCHangeable UNStructured SIGMAGamma(string) SIGMABeta(string) ///
		TREND2 TREND(string) AUGVARiance(real 1e5) * ]

	marksample touse
	
	tokenize `varlist'
	args y stderr

	_get_eformopts, eformopts(`options')
	local eformopt eform(`"`s(str)'"')
	
	// convert study variable into "canonical" form, i.e. with consecutive values 1, 2, ...
	tempname s
	tempname slab
	qui egen `s' = group(`study') if `touse', label(`slab')
	summ `s', meanonly
	local t = r(max)

	// convert subgroup variable into "canonical" form, i.e. with consecutive values 1, 2, ...
	tempvar g
	tempname glab
	qui egen `g' = group(`subgroup') if `touse', label(`glab')
	summ `g', meanonly
	local k = r(max)
	local k1 = `k' - 1

	local b = 1
	if "`base'"!=""      local baseopt base(`base')
	if "`ibase'"!=""     local baseopt `baseopt' ibase(`ibase')
	if "`reference'"!="" local baseopt `baseopt' reference(`reference')
	opts_exclusive `"`baseopt'"' `""' 184
	
	local base = trim(itrim(`"`base'`ibase'`reference'"'))
	if `"`base'"'==`""' local base first	// default
		
	tempvar obs
	qui gen int `obs' = _n

	// base can be:
	// (1) string = value of variable (including numeric)
	// (1) ## = #th ordered value of variable
	// (3) first | last | freq
	if substr("`base'", 1, 1)=="#" {
		local b = substr("`base'", 2, .)
		cap numlist "`b'", min(1) max(1) integer range(>=0)
	}
	
	// else, try to match `base' with a value of `subgroup'			
	else if !inlist("`base'", "first", "last", "freq") {
		cap confirm string variable `subgroup'
		if !_rc {
			summ `obs' if `subgroup'==`"`base'"' & `touse', meanonly
		}
		else {
			summ `obs' if `subgroup'==`base' & `touse', meanonly
		}
		if r(N) {
			local b = `g'[`r(min)']
		}
		else {
			nois disp as err `"Error in option {bf:`baseopt'}"'
			exit 198
		}
	}
	
	// else:  first | last | freq
	else {			
		summ `g' if `touse', meanonly
		if "`base'"=="first" local b = r(min)
		else if "`base'"=="last" local b = r(max)
		else {		// freq
			tempvar c ismode
			qui bysort `touse' `g' : gen int `c' = _N if `touse'
			qui bysort `touse' (`c') : gen byte `ismode' = (`c' == `c'[_N]) if `touse'
			summ `g' if `ismode' & `touse', meanonly
			local b = r(min)
			drop `c' `ismode'
		}
	}
	summ `g' if `touse', meanonly
	if `b' == r(min) local bfirst bfirst		// marker of default behaviour

	/*
	// now sort so that chosen reference comes first
	// (`bvar' is "opposite" so that sorting puts the zeros before the ones)
	tempvar bvar
	qui gen byte `bvar' = (`g'!=`b')
	*/
	
	// check that study-subgroup combinations are unique
	tempvar Ngroup
	qui bysort `touse' `s' /*`bvar'*/ `g' (`obs') : gen int `Ngroup' = _N
	summ `Ngroup' if `touse', meanonly
	if r(max) > 1 {
		nois disp as err "The following combination(s) of study and subgroup are not unique:"
		list `study' `subgroup' `y' `stderr' if `Ngroup' > 1 & `touse'
		exit 198
	}
	
	/*
	// re-generate `g' using new sort order with alternative reference first
	if "`bfirst'"=="" {
		drop `g'
		lab drop `glab'
		qui egen `g' = group(`subgroup') if `touse', label(`glab')
	}
	*/
	
	// setup fixed/random and covariance structures
	opts_exclusive `"`fixed' `randombeta' `exchangeable' `unstructured'"' `""' 184
	local structure `fixed'`randombeta'`exchangeable'`unstructured'
	if "`structure'"!="" {
	    if "`sigmagamma'`sigmabeta'"!="" {
		    nois disp as err `"Covariance structure {bf:`structure'} specified; options {bf:sigmagamma()} and {bf:sigmabeta()} are invalid"'
			exit 198
		}
	}
	/*
	else {
	    if "`sigmagamma'"!="" local sigmagamma "bscov(`sigmagamma')"
	    if "`sigmabeta'" !="" local sigmabeta  "bscov(`sigmabeta')"
	}
	*/

	
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

	qui keep if `touse' & !missing(`s', `g')
	keep `y' `stderr' `s' `g' /*`bvar'*/
	cap rename `y' y
	tempvar V
	qui gen double `V' = `stderr'^2				// Generate estimate of the variance (s-squared)
	drop `stderr'
	rename (`V' `s' `g') (V study subgroup)
	
	// augment if missing subgroup data
	qui isid study subgroup
	cap drop _fillin
	qui fillin study subgroup
	summ _fillin, meanonly
	if r(sum) {
		nois disp as text "The following subgroups are missing, and will be given augmented variances of " as res `augvariance' as text ":"
		list study subgroup if _fillin
	}
	// qui replace `bvar' = (subgroup!=`b')			// fill in `bvar' after -fillin- (new obs will be missing)
	qui isid study /*`bvar'*/ subgroup, sort miss	
	qui replace V = `augvariance' if missing(y)
	qui replace y = 0   if missing(y)

	// store Wi matrices, for later
	forvalues i = 1 / `t' {
		tempname W`i'
		mkmat V if study==`i', matrix(`W`i'')
		cap assert rowsof(`W`i'')==`k'
		if _rc {
			nois disp as err `"Conformability error: study `i' has `=rowsof(`W`i'')' subgroups where it should have `k'"'
			exit 503
		}
		matrix `W`i'' = diag(`W`i'')
	}

	// generate contrasts, in the (default) form:
	//   _y_Isubgroup2 _y_Isubgroup3 ... _y_Isubgroup[k] _y_cons
	// where _y_cons is the reference subgroup (i.e NOT a contrast)
	
	// set alternative reference category
	char subgroup[omit] `b'
	
	// ... if `b' != 1 then instead contrasts are:
	//   _y_Isubgroup1 ... _y_Isubgroup[k or k-1] _y_cons
	// where _y_cons is the reference subgroup (i.e NOT a contrast)
	nois list study subgroup y V, clean noobs
	qui xi: mvmeta_make regress y i.subgroup [iw=V^-1], mse1 by(study) clear nodetails useconstant
	// N.B. regress ... [iw], mse1  is equivalent to using -vwls- 
	//  but currently -mvmeta_make- and -vwls- don't like each other for some reason

	

	** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
	tempname D
	mat `D' = J(`k', `k1', 0)
	mat `D'[2, 1] = I(`k1')
	
	local eq
	if "`trend'"!="" {
		
		// take contrasts of `d'
		tempname dd
		mat `dd' = `D'*`d'

		/*
			// take contrasts of `d'
		tempname dd
		if `b'==1 {					// default
			mat `dd' = `D'*`d'
		}
		else {
			tempname p
			mat `p' = J(`k', 1, .)
			local j = 1
			forvalues i = 1 / `k' {
				if `i'==`b' mat `p'[`b', 1]==1
				else {
					local ++j
					mat `p'[`i', 1]==`j'
				}
			}			
			
			// `b'==2 ==> 2, 1, 3, 4
			// `b'==3 ==> 2, 3, 1, 4
			// `b'==4 ==> 2, 3, 4, 1
			mat `dd' st_matrix()
			
			
		}
		*/
		
		forvalues i = 1 / `k1' {
			local val = `dd'[`i', 1]
			if `i'==1 {
			    gen byte Trend = `val'		// label first equation appropriately to appear in the model and for testing
				local val1 Trend			// (we are using "commonparm" if "trend" so it doesn't really matter which equation we choose)
			}
			else {
				tempvar val`i'
				qui gen byte `val`i'' = `val'
			}
			local ii = `i' + 1
			local eq `eq' y_Isubgroup_`ii' : `val`i'',
			local todrop `todrop' `val`i''
		}
		local eq eq(`eq')
	}

	/*
	A = (1 \ 2 \ 3)
	B2 = (-1, 0, 0 \ -1, 1, 0 \ -1, 0, 1)
	B3 = (0, -1, 0 \ 1, -1, 0 \ 0, 1, -1)
	B4 = (0, 0, -1 \ 1, 0, -1 \ 0, 1, -1)
	P3 = (0, 1, 0 \ 1, 0, 0 \ 0, 0, 1)
	P4 = (0, 0, 1 \ 1, 0, 0 \ 0, 1, 0)
	
	B2*A == B2*A
	B3*A == B2*P3*A
	B4*A == B2*P4*A
	*/
	
	qui mvmeta y S, vars(y_I*) `eq' `fixed' `sigmagamma' `trend'
	tempname VarGammaHat GammaHat SigmaGamma
	
	if "`trend'"!="" {
		nois disp as text _n "Test for trend:"
		test Trend
		mat `VarGammaHat' = e(V)*`dd'*`dd''
		mat `GammaHat'    = e(b)*`dd'
	    qui drop `todrop'
		
		// re-estimate, now we've tested for trend
		if "`testonly'"!="" {
			qui mvmeta y S, vars(y_I*) `fixed' `sigmagamma'
		}
	}
	if "`trend'"=="" | "`testonly'"!="" {
		nois disp as text _n "Test of interactions(s):"
		testparm y*
		mat `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		mat `GammaHat'    = e(b)[1,       1..`k1']'			
	}
	if "`fixed'"=="" {
		tempname Chol
		mat `Chol' = e(b)[1, `k'...]
	}
	mat `SigmaGamma'  = e(Sigma)

	
	** STEP 2: ESTIMATE TREATMENT EFFECT IN REF GROUP (MEAN, BS VARIANCE AND BS CORRELATION WITH INTERACTIONS)
	// subtract within-trial interactions from observed interaction data
	// set up equations that make shifted interactions have mean zero
	qui gen byte zero = 0
	qui gen byte one  = 1
	local eq
	local j = 0
	forvalues i = 1 / `k' {
		cap confirm var y_Isubgroup_`i'
		if !_rc {
			local ++j
			qui replace y_Isubgroup_`i' = y_Isubgroup_`i' - `GammaHat'[`j', 1]
			local eq `eq' y_Isubgroup_`i':zero,
		}
	}
	local eq eq(`eq' y_cons:one)
	cap assert `j'==`k'-1
	if _rc {
		nois disp as err "Error detected when subtracting contrasts from non-reference subgroup"
		exit 198
	}
	
	// set up constraints on the BS variance of the interactions
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
	qui mvmeta_new y S, `eq' nocons commonparm showall print(bscov) `constr_opt' `fixed' `sigmabeta'
	
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
	mat `BetaHat' = `ThetaHat'*`ones' + (0 \ `GammaHat')

	// correct variance for beta
	tempname /*D*/ GradSum Grad
	// mat `D' = J(`k', `k1', 0)
	// mat `D'[2, 1] = I(`k1')
	mat `GradSum' = J(1, `k1', 0)
	forvalues i = 1 / `t' {
		mat `GradSum' = `GradSum' + `ones''*invsym(`W`i'' + `SigmaBeta')*`D'
	}
	mat `Grad' = -`VarThetaHat'*`GradSum'
	
	tempname A VarBetaHat
	mat `A' = `ones'*`Grad' + `D'
	mat `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A' * `VarGammaHat' * `A''

	/*
	// permute if non-standard reference
	if `b'!=1 {
		tempname BetaHatTemp		
		local j = 1
		forvalues i = 1 / `k' {
			if `i'==`b' mat `BetaHatTemp' = nullmat(`BetaHatTemp') \ `BetaHat'[1, 1]
			else {
				local ++j
				mat `BetaHatTemp' = nullmat(`BetaHatTemp') \ `BetaHat'[`j', 1]
			}
		}
		mat `BetaHat' = `BetaHatTemp'
	}
	*/
	matrix `BetaHat' = `BetaHat''	
	
	// set colnames
	matrix coleq `BetaHat' = ""
	local colnames group`b'_ref
	forvalues i = 1 / `k' {
		if `i'!=`b'	local colnames `colnames' group`i'
	}
	
	matname `BetaHat'    `colnames', columns(.) explicit
	matname `VarBetaHat' `colnames', explicit

	// collect and display under "ereturn"
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


