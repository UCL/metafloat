*  v0.1 beta  David Fisher 10jan2021
*  v0.2 beta  David Fisher 30mar2021
*  v0.3 beta  David Fisher 21apr2021
*  v0.4 beta  David Fisher 07may2021
*  v0.5 beta  David Fisher 01jul2021
*  v0.6 beta  David Fisher 12jul2021
*  v0.7 beta  David Fisher 04aug2021
*  v0.8 beta  David Fisher 20sep2021
*! v0.8.1 beta  David Fisher 08oct2021


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
// v0.9: improvements to trend options, and to handling of `subgroup' as numeric/string
//       saving() code put into subroutine

// TO DO:
// are augmentations for interactions necessary (`noaugtrend')??
// sort into subroutines
// collapse() option to mvmeta_make
// is it necessary to have _Trend and/or _Design as permanent varnames not temp?
// debugging/test script


* Syntax:
// varlist (required) = effect size and std err
// study, subgroup (required) = study & subgroup identifiers

// fixed = all fixed (common) effects
// randombeta = common-effects on SigmaGamma; constant within and between-subgroup heterogeneity for SigmaBeta
// exchangeable = exchangeable (single-parameter) structures for both SigmaGamma and SigmaBeta
// default (no options) = unstructured random-effects for both SigmaGamma and SigmaBeta
// naive = unstructured random-effects for both SigmaGamma and SigmaBeta, removing the constraint that they be related via the `M' matrix
// sigmagamma, sigmabeta = manually specify structures

// notrend = don't test for trend [actual internal option name is `trend1']
// trend = impose linear trend upon subgroup estimates [actual internal option name is `trend2']
// trend(mat) = fit specific trend defined by matrix "mat"
// noaugtrend = don't use variance augmentation for estimating trend

//  [Note: if k > 2, *test* for trend is reported automatically, but is *not* imposed upon the estimates]

// design = additional parameters in final model describing the available subgroups per trial (e.g. "single-subgroup" trials)

// nouncertainv = -mvmeta- option; see mvmeta documentation
// `options' = eform options, plus other options to pass to -mvmeta- e.g. tolerance()


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
	
	syntax varlist(numeric min=2 max=2) [if] [in], ///
		STUDY(varname) SUBGROUP(varname) ///
		[ FIXED RANDOMBeta EXCHangeable UNStructured SIGMAGamma(string) SIGMABeta(string) NAIVE LISTConstraints CONSTRaints(numlist integer) ///
		 TREND(string) TREND2 noTREND1 AUGVARiance(real 1e5) noAUGTRend SHOWmodels DESign ///
		 noUNCertainv PRINT(string) ///		/* -mvmeta- options */
		 SAVING(passthru) CLEAR * ]

	marksample touse	
	tokenize `varlist'
	args y stderr
		
	_get_eformopts, eformopts(`options') allowed(__all__) soptions
	local eformopt eform(`"`s(str)'"')
	local soptions `"`s(options)'"'

	
	** Check for potential conflicts in varnames:
	// check that varnames do not conflict with each other
	
	// Note: -xi- with no interactions and a single-letter prefix (as above) will result in a character limit of 9 for the derived names
	//  so later on, we will need to truncate `subgroup' to 9 characters.
	// Therefore, check here that this does not cause a conflict with other relevant varnames
	local subgroup9 = substr(`"`subgroup'"', 1, 9)

	foreach v1 in y stderr study subgroup9 {
		foreach v2 in y stderr study subgroup9 {
			if `"`v1'"'==`"`v2'"' continue
			cap assert `"``v1''"' !=  `"``v2''"'
			if _rc {
				nois disp as err `"Variable name conflict"'
				nois disp as err `"Variables {it:y}, {it:stderr}, {bf:study()} and {bf:subgroup()} must all be distinct"'
				exit _rc
			}
		}
	}
	
	// check that varnames are not called y or V
	foreach v1 in stderr study subgroup {
		foreach v2 in y V {
			cap assert `"``v1''"'!=`"`v2'"'
			if _rc {
				nois disp as err `"Variable name conflict"'
				nois disp as err `"Variables {it:stderr}, {bf:study()} and {bf:subgroup()} must not be named {bf:`v2'}; please rename"'
				exit _rc
			}
		}
	}
	cap assert `"`y'"'!="V"
	if _rc {
		nois disp as err `"Variable name conflict"'
		nois disp as err `"Variables {it:y}, {it:stderr}, {bf:study()} and {bf:subgroup()} must not be named {bf:V}; please rename"'
		exit _rc
	}
	
	// check that varnames do not have prefix _I* or _J* (or _Design if applicable)
	foreach x in y stderr study subgroup {
	    foreach prefix in _I _J {
			cap assert `"`=substr(`"``x''"', 1, 2)'"'!=`"`prefix'"'
			if _rc {
				nois disp as err `"Variable name conflict"'
				nois disp as err `"Variables {it:y}, {it:stderr}, {bf:study()} and {bf:subgroup()} must not begin with {bf:`prefix'}; please rename"'
				exit _rc
			}
		}
		if "`design'"!="" {
			cap assert `"`=substr(`"``x''"', 1, 8)'"'!=`"_Design"'
			if _rc {
				nois disp as err `"Variable name conflict"'
				nois disp as err `"Variables {it:y}, {it:stderr}, {bf:study()} and {bf:subgroup()} must not begin with {bf:_Design}; please rename"'
				exit _rc
			}
		}
	}
	
	// check that varnames are not called _fillin, _Trend, _One, _Zero	
	foreach v1 in y stderr study subgroup {
		foreach v2 in _fillin _Trend _One _Zero {
			cap assert `"``v1''"'!=`"`v2'"'
			if _rc {
				nois disp as err `"Variable name conflict"'
				nois disp as err `"Variables {it:y}, {it:stderr}, {bf:study()} and {bf:subgroup()} must not be named {bf:`v2'}; please rename"'
				exit _rc
			}
		}
	}
	
	// covariance structure options
	opts_exclusive `"`fixed' `randombeta' `exchangeable' `unstructured'"' `""' 184
	local structure `fixed'`randombeta'`exchangeable'`unstructured'
	if `"`sigmagamma'`sigmabeta'"'!=`""' {
	    if `"`structure'"'!=`""' {
			nois disp as err `"Covariance structure {bf:`structure'} specified; options {bf:sigmagamma()} and {bf:sigmabeta()} are invalid"'
			exit 198
		}
		local structure user
		
		foreach x in sigmagamma sigmabeta {
			if `"``x''"'!=`""' & `"``x''"'!="fixed" {
				local 0 `", ``x''"'
				syntax [, UNStructured PROPortional EXCHangeable EQuals CORRelation * ]
				opts_exclusive `"`unstructured' `proportional' `exchangeable' `equals' `correlation'"' `""' 184
				if `"`equals'"'!=`""' | `"`correlation'"'!=`""' {
					nois disp as err `"Invalid option {bf:`x'(`equals'`correlation'} {it:matexp}{bf:)}"'
					exit 198
				}
				if `"`unstructured'"'!=`""' & `"`options'"'!=`""' {
					nois disp as err `"Invalid option {bf:`x'(`unstructured' `options')}"'
					exit 198
				}
				local sp = cond("`x'"=="sigmagamma", "sg", "sb")
				local `sp'struct `unstructured'`proportional'`exchangeable'`equals'`correlation'
				local `sp'opt `options'
				local `x' `"bscov(``sp'struct' ``sp'opt')"'
			}
		}
	}
	else if "`structure'"=="" & `"`trend'`trend2'"'==`""' local structure unstructured

	if "`structure'"!="unstructured" & "`naive'"!="" {
		nois disp as err `"Option {bf:naive} is only valid with {bf:unstructured} covariance structure"'
		exit 198
	}
	if "`sigmabeta'"!="" & "`constraints'"!="" {
		nois disp as err `"Option {bf:contraints()} is only valid with option {bf:sigmabeta()}"'
		exit 198
	}

		
	** Preserve data, prior to editing
	preserve

	qui keep if `touse' & !missing(`study', `subgroup')
	keep `y' `stderr' `study' `subgroup'
		
	cap assert `stderr' >= 0
	if _rc {
		nois disp as err "Standard error cannot be negative"
		exit 198
	}

	// find number of subgroups
	qui tab `subgroup'
	local k = r(r)	
	local k1 = `k' - 1
	cap confirm string variable `subgroup'
	if !_rc local string string
	qui levelsof `subgroup', local(sglist)
	
	// find number of studies
	qui tab `study'
	local t = r(r)
	
	// check that study-subgroup combinations are unique
	tempvar Ngroup
	qui bysort `study' `subgroup' (`obs') : gen int `Ngroup' = _N
	summ `Ngroup', meanonly
	if r(max) > 1 {
		nois disp as err "The following combination(s) of study and subgroup are not unique:"
		nois list `study' `subgroup' `y' `stderr' if `Ngroup' > 1
		exit 198
	}
	
	tempvar V
	qui gen double `V' = `stderr'^2				// Generate estimate of the variance (s-squared)
	drop `stderr'
	if `"`y'"'!=`"y"' {
		qui rename `y' y		// There should not be conflicts here, due to earlier checking
	}
	if `"`V'"'!=`"V"' {
		qui rename `V' V
	}
	
	// augment if missing subgroup data
	qui fillin `study' `subgroup'
	assert _fillin == missing(y)
	
	// if any observed variances are larger than `augvariance'
	//  (e.g. if estimated from very small sample)
	// also set to _fillin
	qui replace _fillin = 1 if V > `augvariance'

	summ _fillin, meanonly
	if r(sum) {
		nois disp as text "The following subgroups are unobserved (or contain missing or extremely imprecise data)"
		nois disp as text " and will be given augmented variances of " as res `augvariance' as text ":"
		list `study' `subgroup' if _fillin
	}
	qui replace V = `augvariance' if _fillin
	qui replace y = 0 if _fillin

	// obtain "design", i.e. which subgroups are present
	qui levelsof `study', local(slist)
	if "`design'"!="" {
	    
		// `subgroup' must either be numeric, or single-character string
		if `"`string'"'!=`""' {
			tempvar lensub
			gen `lensub' = length(`subgroup')
			summ `lensub', meanonly
			if r(max) > 1 {
			    nois disp as err "With option {bf:design}, variable {bf:subgroup()} must be either numeric, or use single string characters only"
				exit 198
			}
		}
		
		qui gen _Design = ""
		foreach s of local slist {
			if `"`string'"'!=`""' {
				qui levelsof `subgroup' if `study'==`"`s'"' & !_fillin
				qui replace _Design = `"`r(levels)'"' if `study'==`"`s'"'
			}
			else {
				qui levelsof `subgroup' if `study'==`s' & !_fillin
				qui replace _Design = `"`r(levels)'"' if `study'==`s'			
			}
		}
		local design_opt collapse((firstnm) _Design)
		local showmodels showmodels
	}
	drop _fillin
	

	** Generate transformation matrix based on choice of reference subgroup
	tempname T
	if `"`: char `subgroup'[omit]'"'==`""' {
	    local r = 1
	    gettoken ref sgnoref : sglist
		matrix define `T' = I(`k')
	}
	else {
	    local ref : char `subgroup'[omit]
		local r : list posof `"`ref'"' in sglist
		local sgnoref : list sglist - ref
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
	if `"`string'"'!=`""' {
		qui numlist "1(1)`k'"
		local sglist `"`r(numlist)'"'
		local ref : copy local r
		local sgnoref : list sglist - ref
	}
	// Note: `r' is the index; `ref' is the value (e.g. for labelling the matrix stripe elements)

	// store Wi matrices, for later
	// (with transformation applied)
	qui levelsof `study', local(slist)
	local i = 0
	foreach s of local slist {
		local ++i
		tempname W`i'
		
		cap confirm string variable `study'
		if _rc {	
			mkmat V if `study'==`s', matrix(`W`i'')
		}
		else {
			mkmat V if `study'==`"`s'"', matrix(`W`i'')			
		}
		
		cap assert rowsof(`W`i'')==`k'
		if _rc {
			nois disp as err `"Conformability error: study `s' has `=rowsof(`W`i'')' subgroups where it should have `k'"'
			exit 503
		}
		matrix define `W`i'' = diag(`T'*`W`i'')
	}
	
	// define useful matrices L and M
	tempname L M
	matrix define `L' = J(`k', `k1', 0), J(`k', 1, 1)
	matrix define `L'[2, 1] = I(`k1')
	matrix define `M' = J(`k1', 1, -1), I(`k1')	

	// setup tests for trend
	if `"`trend1'"'==`""' | `"`trend'`trend2'"'!=`""' {
		if `"`trend1'"'!=`""' {
			nois disp as err `"cannot specify both {bf:trend} and {bf:notrend}"'
			exit 184
		}
		if `"`trend'"'!=`""' & `"`trend2'"'!=`""' {
			nois disp as err `"only one of {bf:trend} or {bf:trend()} is allowed"'
			exit 184
		}
		
		if `k' < 3 {
			local trend1 notrend
			if `"`trend'`trend2'"'!=`""' {
				nois disp `"{error}Cannot test for trend with fewer than 3 subgroups"'
				local trend
				local trend2
			}
		}
		
		else if `"`trend'"'==`""' {
			tempname d
			forvalues i = 1 / `k' {
				matrix define `d' = nullmat(`d') \ `i'	// linear trend
			}
						
			// Convert `d' into a matrix of contrasts
			// ... and apply transformation if reference is not first subgroup
			matrix define `d' = `M' * `T' * `d'
		}
		
		else {
			* If -syntax- didn't work, assume due to commas but no brackets
			if _rc {
				disp as err "invalid syntax for option {bf:trend()}"
				disp as err "Valid syntaxes are:"
				disp as err `" the name of an existing matrix; or"'
				disp as err `" matrix-style input {bf:(} {it:a}{bf:,} {it:b} {bf:,} {it:c} {bf:)}; see {help matrix define}"'
				exit _rc
			}
			tempname d
			cap {
				matrix define `d' = `trend'
				if colsof(`d')==`k' & rowsof(`d')==1 {
					matrix define `d' = `d''
				}
				assert rowsof(`d')==`k'
				assert colsof(`d')==1
			}
			if _rc {
				nois disp as err `"matrix {bf:`trend'} not found, or has invalid dimensions"'
				nois disp as err `"(should be a vector of length equal to the number of subgroups)"'
				exit 198
			}
			local trend trend
			
			// Convert `d' into a matrix of contrasts
			// ... and apply transformation if reference is not first subgroup
			tempname M
			matrix define `M' = J(`k1', 1, -1), I(`k1')
			matrix define `d' = `M' * `T' * `d'
		}
		
		if `"`trend'`trend2'"'!=`""' {
		    local brackets = cond(`"`trend'"'!=`""', "()", "")
			local trend trend
			if !inlist(`"`structure'"', `"fixed"', `""') {
				nois disp as err `"cannot specify {bf:`structure'} with {bf:`trend'`brackets'}; random-effects structure is implied by use of trend"'
				exit 184
			}
			local trendtext `" (with trend imposed)"'		// for later
		}
	}

	if `"`saving'"'!=`""' | `"`clear'"'!=`""' {

		// save subgroup value labels
		tempfile sglabfile
		if `"`string'"'!=`""' {
			tempvar subgroup2
			tempname sglab
			qui encode `subgroup', gen(`subgroup2') label(`sglab')
			qui drop `subgroup'
			qui rename `subgroup2' `subgroup'
		}
		else local sglab : value label `subgroup'
		qui label save `sglab' using `sglabfile'

		cap rename `subgroup' `subgroup9'
		
		tempfile main_effects
		qui save `main_effects'
	}

	
	// generate contrasts, using prefix _J
	// Note: _y_cons is the reference subgroup (i.e NOT a contrast)
	qui xi, prefix(_J) : mvmeta_make regress y i.`subgroup' [iw=V^-1], mse1 by(`study') `design_opt' clear nodetails useconstant
	local yvars : colnames e(b)
	
	// N.B. regress ... [iw], mse1  is equivalent to using -vwls- 
	//  but currently -mvmeta_make- and -vwls- don't like each other for some reason

	// Note: -xi- with no interactions and a single-letter prefix (as above) will result in a character limit of 9 for the derived names
	// so in what follows, we will refer to `subgroup' truncated at 9 characters.
	// (we have already determined that this does not cause a conflict with other relevant varnames)	

	// form baseline expression
	local i = 0
	tokenize `yvars'
	local Ibase = subinstr(`"y`1'"', `"_J"', `"_I"', 1)
	local l = length(`"`Ibase'"') - 2
	local Ibase = substr(`"`Ibase'"', 1, `l')
	local Ibase `Ibase'_`ref'
	
	// insert baseline expression at appropriate point
	while `"`1'"'!=`""' {
		local ++i
		if `i' > `k' continue, break
		else if `i'==`r' local Iyvars `Iyvars' `Ibase'
		else {
		    local Iyvar = subinstr(`"y`1'"', `"_J"', `"_I"', 1)
			local Iyvars `Iyvars' `Iyvar'
			local Jyvars `Jyvars' y`1'
			macro shift
		}
	}	
	
	if `"`saving'"'!=`""' | `"`clear'"'!=`""' {
		tempfile interactions
		qui save `interactions'
	}

	
	** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
	if "`showmodels'"=="" local quietly quietly
	if "`print'"=="" local print bscov

	// user-specified matrices
	if "`structure'"=="user" {
		if `"`sigmabeta'"'!=`""' & `"`sigmabeta'"'!="fixed" & `"`sigmagamma'"'==`""' & `"`naive'"'==`""' {
			if `"`sbstruct'"'=="proportional" {
				tempname sgmat
				cap matrix define `sgmat' = `M'*`sbopt'*`M''
				if _rc {
					nois disp as err `"{bf:sigmabeta()} option: `sbopt' is not `k'x`k'"'
					exit 198
				}
				local sigmagamma `"bscov(proportional `sgmat')"'
			}
			else if `"`sbstruct'"'=="exchangeable" {
				tempname sgmat
				cap matrix define `sgmat' = `M'*( (1-`sbopt')*I(`k') + `sbopt'*J(`k', `k', 1) )*`M''
				if _rc {
					nois disp as err `"{bf:sigmabeta()} option: syntax is {bf:exchangeable} {it:#}{bf:)} with -1<=#<=1"'
					exit 198
				}
				local sigmagamma `"bscov(proportional `sgmat')"'
			}
		}
	}
	else if inlist("`structure'", "fixed", "randombeta") local sigmagamma fixed
	else if "`structure'"=="exchangeable" local sigmagamma bscov(exchangeable .5)
	else if `"`sigmagamma'"'==`""' local sigmagamma bscov(unstructured)
	
	
	// fit as standard model (assuming no trend)
	tempname VarGammaHat GammaHat SigmaGamma Chol
	nois disp as text _n "Test of interaction(s):"
	`quietly' mvmeta y S, vars(y_J*) print(`print') `sigmagamma' `soptions' `eformopt'
	nois testparm y*
	
	matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
	matrix define `GammaHat'    = e(b)[1,       1..`k1']
	matrix define `SigmaGamma'  = e(Sigma)
		
	if `"`sigmagamma'"'==`"fixed"' & `"`fixed'"'==`""' {
		matrix define `Chol' = J(`k1', `k1', 0)
	}
	
	else if "`structure'"=="exchangeable" {
		local tau = e(b)[1, `k']
		cholesky2 `SigmaGamma' `Chol'
		matrix define `Chol' = sign(`tau') * `Chol'
	}
	
	// special case: if `unstructured', -mvmeta- outputs the elements of `Chol' as part of e(b)
	// so don't bother packing them up into a matrix; we're only going to unpack them again later
	else if "`structure'"=="unstructured" {
		matrix define `Chol' = e(b)[1, `k'...]
	}
	
	local rownames_eb : rownames e(b)
	matrix rownames `GammaHat' = `: word 1 of `rownames_eb''
	matrix coleq `GammaHat' = :
	matrix coleq `VarGammaHat' = :
	matrix roweq `VarGammaHat' = :

	local eqlist
	if `"`trend1'"'==`""' {		// test for trend if k>3, unless `notrend'
		local i = 0		
		tokenize `yvars'
		while `"`2'"'!=`""' {	// this will ignore the final element, which is "_cons"
			local ++i
			qui gen _Trend`i' = `d'[`i', 1]
			local eqlist `eqlist' y`1':_Trend`i',
			
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
		
		local sgammatrend : copy local sigmagamma
		if "`sigmagamma'"!="fixed" {
			tempname D
			matrix define `D' = `d'*`d''
			local sgammatrend bscov(proportional `D')
		}
		
		nois disp as text _n "Test for trend:"
		`quietly' mvmeta y S, vars(y_J*) print(`print') `sgammatrend' commonparm nocons equations(`eqlist') `soptions' `eformopt'
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
			matrix define `VarGammaHat' = e(V)[1, 1]		// Gamma is a constant if `trend'
			matrix define `GammaHat'    = e(b)[1, 1]
			matrix define `SigmaGamma'  = e(Sigma)

			matrix define `VarGammaHat' = `VarGammaHat'*`d'*`d''
			matrix define `GammaHat'    = `GammaHat'*`d''
			
			if "`sigmagamma'"=="fixed" {
				matrix define `Chol' = J(`k1', `k1', 0)
			}
			else  {
				local tau = e(b)[1, 2]
				cholesky2 `SigmaGamma' `Chol'
				matrix define `Chol' = sign(`tau') * `Chol'
			}
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
	qui gen byte _Zero = 0
	qui gen byte _One  = 1

	// DF July 2021: optionally adjust for "design"
	if "`design'"!="" {
	    tempvar obs
		gen long `obs' = _n
		local designvars
		qui tab `_Design gen(_Design) sort
		forvalues i = 2 / `r(r)' {
		    summ `obs' if _Design`i'==1, meanonly
			local destxt = _Design[`r(min)']
			local destxt = subinstr(trim(`"`destxt'"'), " ", "_", .)
			rename _Design`i' _Des_`destxt'
			local designvars `designvars' _Des_`destxt'
		}
		drop _Design1
	}		

	local eqlist
	local i = 0
	tokenize `e(yvars)'
	while `"`1'"'!=`""' {
	    local ++i
		qui replace `1' = `1' - `GammaHat'[1, `i']
		local eqlist `eqlist' `1':_Zero `designvars',
		macro shift
	}	
	local eqlist equations(`eqlist' y_cons:_One `designvars')
	
	cap assert `i'==`k'-1
	if _rc {
		nois disp as err "Error detected when subtracting contrasts from non-reference subgroup"
		exit 198
	}

	// set up structure of sigmabeta: user-defined
	if "`structure'"=="user" {
		if `"`sbstruct'"'=="proportional" {
			tempname propU
			matrix define `propU' = inv(`L')*`sbopt'*inv(`L')'
			local sigmabeta bscov(proportional `propU')
		}
		else if `"`sbstruct'"'=="exchangeable" {
			tempname propU
			matrix define `propU' = inv(`L')*( (1-`sbopt')*I(`k') + `sbopt'*J(`k', `k', 1) )*inv(`L')'
			local sigmabeta bscov(proportional `propU')
		}
	}
	
	// set up structure of sigmabeta: options
	else if "`structure'"=="fixed" local sigmabeta fixed
	else if "`structure'"=="randombeta" {
		tempname propU
		matrix define `propU' = J(`k', `k', 0)
		// matrix `propU'[1, 1] = 1
		matrix `propU'[`k', `k'] = 1
		matrix `propU' = `T''*`propU'*`T'
		local sigmabeta bscov(proportional `propU')
	}
	
	// set up constraints on the BS variance of the interactions
	// DF Note: Need to look at this again
	else if "`naive'"=="" {
	    local sigmabeta bscov(unstructured)
		
		local c = 0
		forvalues i = 1 / `k1' {
			forvalues j = 1 / `i' {
				local ++c
				
				if "`structure'"=="unstructured" local el = `Chol'[1, `c']	// special case: see above
				else local el = `Chol'[`i', `j']
				
				constraint free
				constraint `r(free)' [chol`j'`i']_b[_cons] = `el'
				local constr_list `constr_list' `r(free)'
			}
		}
		
		// if exchangeable, need to do a little more work
		if "`structure'"=="exchangeable" {
			tempname Chol_k
			matrix define `Chol_k' = inv(`Chol') * J(`k1', 1, -.5*(`tau')^2)
			forvalues j = 1 / `k1' {
				constraint free
				constraint `r(free)' [chol`j'`k']_b[_cons] = `=`Chol_k'[`j', 1]'
				local constr_list `constr_list' `r(free)'
			}
		}
		
		// similarly if trend
		else if "`trend'"!="" {
			forvalues j = 1 / `k1' {
				constraint free
				constraint `r(free)' [chol`j'`k']_b[_cons] = 0
				local constr_list `constr_list' `r(free)'
			}
		}
		
		// constraint 4 y_Isubgroup_2=0
		// constraint 5 y_Isubgroup_3=0
		
		if "`listconstraints'"!="" {
		    nois disp as text _n "Constraints on elements of Cholesky factor for SigmaU:"
			nois constraint dir `constr_list'
		}
		local constr_opt constraints(`constr_list')
	}
	else if "`structure'"=="user" & "`listconstraints'"!="" {
		nois constraint dir `constraints'
		local constr_opt constraints(`constraints')
	}

	
	// fit model: make sure we are using "NEW" mvmeta
	if "`showmodels'"!="" nois disp as text _n `"Estimation of floating subgroup for reference category`trendtext':"'
	`quietly' mvmeta_new y S, print(`print') nocons commonparm `eqlist' `constr_opt' `sigmabeta' `soptions' `eformopt'
	// Note: the coefficient of _cons refers to the reference subgroup
	
	// extract fitted BS variance and derive SigmaBeta
	tempname SigmaU SigmaBeta
	matrix define `SigmaU' = e(Sigma)
	matrix define `SigmaBeta' = `L'*`SigmaU'*`L''
	
	// define covariate data matrix Z (arranged such that the reference subgroup comes first)
	tempname ones Z
	matrix define `ones' = J(`k', 1, 1)			// column vector of ones, of length k
	matrix define `Z' = J(`k', `k1', 0)
	matrix define `Z'[2, 1] = I(`k1')	
	
	// compute implied means for beta
	tempname ThetaHat BetaHat
	matrix define `ThetaHat' = e(b)[1,1]
	matrix define `BetaHat' = `ones'*`ThetaHat' + `Z'*`GammaHat''
	
	// IW June 2021: First derive `VarThetaHat' ignoring the uncertainty in the heterogeneity matrix
	// because it's needed for deriving the gradient vector
	tempname VarThetaHat
	matrix define `VarThetaHat' = e(V)
	matrix define `VarThetaHat' = invsym(`VarThetaHat')
	matrix define `VarThetaHat' = `VarThetaHat'[1,1]
	matrix define `VarThetaHat' = invsym(`VarThetaHat')	

	// Now define matrix A, the independent multiplier of GammaHat
	tempname A
	matrix define `A' = J(1, `k1', 0)
	forvalues i = 1 / `t' {
		matrix define `A' = `A' + `ones''*invsym(`W`i'' + `SigmaBeta')*`Z'
	}
	mat `A' = `Z' - `ones'*`VarThetaHat'*`A'

	// now correct VarThetaHat for uncertainty in the heterogeneity matrix, if requested
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
	matrix define `SigmaBeta' = `T''*`SigmaBeta'*`T'
	matrix coleq `SigmaBeta' = ""
	matname `SigmaBeta' `Iyvars', explicit
	
	// collect and display under "ereturn"
	matrix define `BetaHat' = `BetaHat''
	matrix coleq `BetaHat' = :
	matname `BetaHat'    `Iyvars', columns(.) explicit
	matname `VarBetaHat' `Iyvars', explicit

	
	*** Construct saved dataset
	if `"`saving'"'!=`""' | `"`clear'"'!=`""' {

		tempname BetaTable GammaTable
		mat `GammaTable' = `GammaHat'', vecdiag(`VarGammaHat')'
		matrix rownames `GammaTable' = `sgnoref'	
	
		mat `BetaTable' = `BetaHat'', vecdiag(`VarBetaHat')'
		matrix rownames `BetaTable' = `sglist'

		ProcessSavedData `GammaTable' `BetaTable', saving(`saving') ///
			study(`study') subgroup(`subgroup9') sglabfile(`sglabfile') ///
			mainfile(`main_effects') intfile(`interactions') augvariance(`augvariance')

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
	ereturn scalar n = `t'		// number of studies
	ereturn scalar k = `k'		// number of subgroups

	ereturn local bscov_Gamma `"`sigmagamma'"' 
	ereturn local bscov_Beta `"`sigmabeta'"' 
	
	nois disp as text _n `"Floating subgroups`trendtext':"'
	ereturn display, `eformopt'

	if `"`constr_list'"'!=`""' {
	    constraint drop `constr_list'
	}
	
end




program define cholesky2
* produce an approximate cholesky decomposition, even if matrix is not positive definite
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




program define ProcessSavedData

	syntax anything, ///
		study(name) subgroup(name) sglabfile(string) mainfile(string) intfile(string) ///
		[AUGVARiance(real 1e5) SAVING(string asis) ] 
	
	tokenize `anything'
	args GammaTable BetaTable
	
	local sgnoref : rownames `GammaTable'
	local sglist : rownames `BetaTable'
	local ref : list sglist - sgnoref
	
	// store pooled estimates:
	tempfile pooled_ests
	qui metani `GammaTable', variances nograph nooverall rownames rowtitle(Pooled estimates) ///
		saving(`pooled_ests', stacklabel) prefix(_yInt)

	qui metani `BetaTable', variances nograph nooverall rownames rowtitle(Pooled estimates) ///
		clearstack prefix(_y)

	rename _y_LABELS _yInt_LABELS
	
	tempvar sortorder
	qui gen int `sortorder' = _n		// -merge- will sort the data on _yInt_LABELS, so need to re-sort afterwards
	qui merge 1:1 _yInt_LABELS using `pooled_ests', assert(match master)
	cap {
		assert _merge==1 if _yInt_LABELS=="`ref'"
		assert _merge==3 if _yInt_LABELS!="`ref'"
	}
	if _rc {
		nois disp as err "Error when merging interaction data with main-effect data"
		exit _rc
	}

	drop _merge *_EFFECT _yInt_USE _yInt_STUDY
	rename _y_STUDY _STUDY
	rename _yInt_LABELS _LABELS
	qui replace _y_USE = 5 if inlist(_y_USE, 1, 2)
	sort `sortorder'
	qui save `pooled_ests', replace

	// merge subgroup estimates and interaction estimates into a single dataset
	qui use `intfile', clear
	qui ds *_cons*						// y-vars and S-vars relating to the reference subgroup
	local vlist `"`r(varlist)'"'

	// prepare for reshape long
	foreach v of local vlist {
		local newname = subinstr(`"`v'"', `"_cons"', `"_J`subgroup'_`ref'"', .)
		qui rename `v' `newname'
	}
	foreach x of local sglist {
		local rsvlist `rsvlist' S_J`subgroup'_@_J`subgroup'_`x'
	}
	qui reshape long y_J`subgroup'_ `rsvlist', i(`study') j(`subgroup')
	qui merge 1:1 `study' `subgroup' using `mainfile', assert(match) nogen
	qui do `sglabfile'
	label values `subgroup' `sglab'
	
	qui gen double S_J`subgroup'_ = .
	order S_J`subgroup'_, after(y_J`subgroup'_)
	foreach x of local sglist {
		if `x'==`ref' continue
		qui replace S_J`subgroup'_ = S_J`subgroup'__J`subgroup'_`x' if `subgroup'==`x'
		drop S_J`subgroup'__J`subgroup'_`x'
	}
	qui replace y_J`subgroup'_ = .  if `subgroup'==`ref'
	qui replace S_J`subgroup'_ = .  if `subgroup'==`ref'

	qui rename y_J`subgroup'_ _yInt_ES
	qui replace _yInt_ES = . if S_J`subgroup'_ >= `augvariance'
	qui gen double _yInt_seES = sqrt(S_J`subgroup'_) if S_J`subgroup'_ < `augvariance'
	drop S_J`subgroup'_
	
	// use -metan- to restructure the main effects
	qui replace y = . if V >= `augvariance'
	qui gen double stderr = sqrt(V) if V < `augvariance'
	qui metan y stderr, nograph nohet nosubgroup nooverall keeporder ///
		study(`subgroup') by(`study') rcols(_yInt_ES _yInt_seES) ///
		clear prefix(_y)
	qui drop *_EFFECT
	qui rename _y_BY _BY
	qui rename _y_STUDY _STUDY
	qui rename _y_LABELS _LABELS
	qui order _y_USE, before(_y_ES)

	qui gen double _yInt_LCI = _yInt_ES - invnormal(.975)*_yInt_seES
	qui gen double _yInt_UCI = _yInt_ES + invnormal(.975)*_yInt_seES
	qui gen double _yInt_WT = 1/_yInt_seES^2 if _y_USE!=5
	summ _yInt_WT, meanonly
	qui replace _yInt_WT = 100*_yInt_WT / r(sum)
	label variable _yInt_WT "% Weight"
	format _yInt_WT %6.2f
			
	// append subgroup estimates
	// ... and "correct" the contents of _STUDY (see above)
	tempvar pooled obs
	qui append using `pooled_ests', gen(`pooled')
	qui gen int `obs' = _n
	local i = 0
	foreach x of local sglist {
		local ++i
		summ `obs' if _STUDY==`x' & !`pooled', meanonly
		local labi = _LABELS[`r(min)']
		qui replace _LABELS = `"`labi'"' if `pooled' & _STUDY==`i'
		qui replace _STUDY = `x' if `pooled' & _STUDY==`i'
	}
	drop `pooled' `obs'
	
	qui gen byte _yInt_USE = _y_USE
	qui replace _yInt_USE = 9 if _STUDY==`ref'
	qui replace _yInt_USE = 2 if _yInt_USE==1 & missing(_yInt_seES)
	order _yInt_USE, before(_yInt_ES)	

	qui compress
	
	if `"`saving'"'!=`""' {
		_prefix_saving `saving'
		qui save `s(filename)', `s(replace)'
	}
	
end






