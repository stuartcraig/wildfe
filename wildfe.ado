/*---------------------------------------------------------------------wildfe.do

Stuart Craig
20190401
*/




cap prog drop wildfe
prog define wildfe, rclass
	syntax varlist [if] [in] [aweight] , absorb(varlist) [cluster(varname)] [weights(string)] [niter(real 1000)]

/*
-------------------------------------------

Preliminaries

-------------------------------------------
*/	
		
	// Sample
		marksample touse
		
	// Parse the variable list
		loc depvar: word 1 of `varlist'
		loc indvars ""
		local ncol=wordcount("`varlist'")
		forval i=2/`ncol' {
			loc in: word `i' of `varlist'
			loc indvars "`indvars' `in'"
		}

		
	// Demand well-behaved weights	
		if "`weights'"=="" loc weights "rademacher"
		cap assert inlist("`weights'","mammen","normal","rademacher")
		if _rc!=0 {
			di as error "`weights' weighting not supported"
			exit 12
		}
		
/*
-------------------------------------------

Implement the actual BS procedure

-------------------------------------------
*/				
		tempfile results
		cap postclose BS
		qui postfile BS b str40 ivar iter using `results', replace
	

	// Run baseline model and generate residuals and fitted values
		tempvar r
		tempvar xb
		di in yellow "+==============================================================+"
		di in yellow "+ Baseline Estimates                                           +"
		di in yellow "+==============================================================+"
		reghdfe `depvar' `indvars' if `touse' [`weight' `exp'], absorb(`absorb') resid(`r')
		qui gen `xb' = `depvar' - `r'
		foreach ivar of varlist `indvars' {
			return scalar b_`ivar'=_b[`ivar']
		}
		
	// Generate the Mammen weights up front if necessary 
	// (reduces computation)
		if inlist("`weights'","mammen") {
			tempvar W1
			tempvar W2
			qui gen `W1' = -(sqrt(5) - 1)/2 
			qui gen `W2'=(sqrt(5) + 1)/2
			loc cut = (sqrt(5) - 1)/(2*sqrt(5))
		}
	
	// Iterate:
	di in yellow "+==============================================================+"
	di in yellow "+ Running Bootstrap Iterations                                 +"
	di in yellow "+==============================================================+"
	
	sort `cluster'
	cap timer clear
	timer on 1
		forval iter=1/`niter' {
		qui {
			tempvar W
			if "`weights'"=="rademacher" {
				qui gen `W' = -1
				qui replace `W' = 1 if runiform()>0.5
				
			}
			if "`weights'"=="mammen" {
				qui gen `W' = `W1'
				qui replace `W' = `W2' if runiform()< `cut'
			}
			if "`weights'"=="normal" qui gen `W' = rnormal()
			
			// Cluster the errors if required
			if "`cluster'"!="" {
				qby `cluster': replace `W' = `W'[1]
			}
			tempvar ystar
			cap drop `ystar'
			qui gen `ystar' = `xb' + `r'*`W'
			* reg `ystar' `indvars'
			qui reghdfe `ystar' `indvars' if `touse' [`weight' `exp'] , absorb(`absorb') 
			
			loc cc=0
			foreach v of varlist `indvars' {
				loc ++cc
				post BS (_b[`v']) ("`v'") (`iter')
			}
			drop `ystar' `W'
			macro drop ystar W 
		}	
		di . _c
		if mod(`iter',50)==0 {
			qui timer off 1
			qui timer list
			di `iter', "in `r(t1)' seconds"
			timer clear
			timer on 1
		}
		}
		postclose BS
	
	// Display the results (crudely):
	preserve
		use `results', clear
		
		cap drop p_lower
		cap drop p_upper
		qui gen p_lower=.
		qui gen p_upper=.
		qui levelsof ivar, local(vars)
		foreach var of local vars {
			qui _pctile b if ivar=="`var'", p(2.5, 97.5)
			qui replace p_lower = r(r1) if ivar=="`var'"
			qui replace p_upper = r(r2) if ivar=="`var'"
		}
		collapse (mean) Coeff=b (sd) SE=b (max) p_ll=p_lower p_ul=p_upper, by(ivar) fast
		
		list ivar SE p_*, noobs
		// Return the SEs
		cap drop n
		qui gen n=_n
		qui levelsof ivar, local(ivars)
		foreach ivar of local ivars {
			qui summ n if ivar=="`ivar'", mean
			return scalar se_`ivar' = SE[r(mean)]
			return scalar p_ll_`ivar' = p_ll[r(mean)]
			return scalar p_ul_`ivar'= p_ul[r(mean)]
		}
	restore
	
end



exit
