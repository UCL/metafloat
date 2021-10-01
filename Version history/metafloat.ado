*! v0.1 beta  David Fisher 10jan2021

// varlist (required) = effect size and std err
// study, subgroup (required) = study & subgroup identifiers
// default (no options) = unstructured random-effects
// fixed = all fixed (common) effects
// sigmagamma, sigmabeta = specify structures
// trend = fit linear trend across subgroups, if applicable


program define metafloat, eclass

	syntax varlist(numeric min=2 max=2) [if] [in], ///
		STUDY(varname) SUBGROUP(varname) ///
		[ FIXED SIGMAGamma(string) SIGMABeta(string) ///
		 TREND AUGVARiance(real 1e5) * ]

	marksample touse
	
	tokenize `varlist'
	args y stderr

	_get_eformopts, eformopts(`options')
	local eformopt eform(`"`s(str)'"')
	
	// convert subgroup variable into "canonical" form, i.e. with consecutive values 1, 2, ...
	tempvar g
	tempname glab
	qui egen `g' = group(`subgroup') if `touse', label(`glab')

	summ `g', meanonly
	local k = r(max)
	
	if "`trend'"!="" {
		if `k'<3 {
			nois disp `"{error}Cannot test for trend with fewer than 3 subgroups"'
			local trend
		}
		else local trend commonparm
	}

	// convert study variable into "canonical" form, i.e. with consecutive values 1, 2, ...
	tempname s
	tempname slab
	qui egen `s' = group(`study') if `touse', label(`slab')
	summ `s', meanonly
	local t = r(max)
	
	// preserve data, prior to editing
	preserve

	qui keep if `touse' & !missing(`s', `g')
	keep `y' `stderr' `s' `g'
	
	cap rename `y' y
	tempvar V
	qui gen double `V' = `stderr'^2				// Generate estimate of the variance (s-squared)
	drop `stderr'
	rename (`V' `s' `g') (V study subgroup)

	// augment if missing subgroup data
	qui fillin study subgroup
	summ _fillin, meanonly
	if r(sum) {
		nois disp as text "The following subgroups are missing, and will be given augmented variances of " as res `augvariance' as text ":"
		list study subgroup if _fillin
	}
	drop _fillin
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

	// generate contrasts (use prefix "_Int" = "interaction")
	qui xi, prefix(_Int) : mvmeta_make reg y i.subgroup [iw=V^-1], mse1 by(study) clear nodetails useconstant

	
	** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
	qui mvmeta y S, vars(y_Int*) `fixed' `sigmagamma' `trend'
		
	tempname VarGammaHat GammaHat SigmaGamma
	local k1 = `k' - 1
	mat `SigmaGamma'  = e(Sigma)
	
	if "`fixed'"=="" {
		tempname Chol
		mat `Chol' = e(b)[1, `k'...]
	}
	if "`trend'"!="" {
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
	qui gen byte zero = 0
	qui gen byte one  = 1
	forvalues i = 2 / `k' {
		qui replace y_Isubgroup_`i' = y_Isubgroup_`i' - `GammaHat'[`i'-1, 1]
		local eq `eq' y_Isubgroup_`i':zero,
	}
	local eq eq(`eq' y_cons:one)
	
	// local eq eq(y_Isubgroup_2:zero, y_Isubgroup_3:zero, y_cons:one) 
	
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
	
	forvalues i = 1 / `k' {
		local colnames `colnames' group`i'
	}
	
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