/*******************************************************************************\
* Based on myfun_combined_v2.do in the version 3 folder
* Modified in 4 June, 2016                        
* version 3.0                                    
* Use STATA command qreg to do all quantile regressions  
\*******************************************************************************/

capture program drop eqreg
program eqreg,eclass 
	//version 13
	//syntax varlist(min=2 numeric) [if] [in] , HOMoskedastic(numlist >=0 <=1 integer) [Boots(real 600)]
	syntax varlist(min=2 numeric) [if] [in] , PHI(string) [Boots(real 600)]
	marksample touse
	tokenize `varlist'
	
	// Assume right now users pass -Y into the varlist!!!
	* PARAMETERS
	
	local G = 40
	local B = 300
	matrix l = (0.65, 0.85, 1.15, 1.45)
	local J = colsof(l)
	local d: word count `varlist'
	qui count
	local N = r(N) 
	local lower = min(80/`boots',0.1)
	local upper = 0.3
	local step = (`upper' - `lower')/`G'
	
	mat ss = J(1, colsof(`phi'),1) * (`phi')'
	local dbb = ss[1,1]
	mat par_bootb = J(`dbb',`B',0)
	mat chi_bootbb = J(`B',1,0)
	scalar disbb = 100
	
	di in gr "Begin"
	di "$S_TIME"
	
	* SELECT OPTIMAL TAU
	forvalues gg = 1/`G'{
		
		local tau = `lower' + `step' * `gg'
		qui qreg_simplified12 `varlist',quantile(`tau') cformat(%10.0g)
		mat theta = e(b)'
		mat thetatemp = J(`d',`B',0)
		
		* BOOSTRAP THE FIRST STAGE ESTIMATOR TO COMPUTE OMEGA_0
		
		forvalues bb = 1/`B'{
			preserve
			qui bsample 
			qui qreg_simplified12 `varlist',quantile(`tau') cformat(%10.0g)
			mat thetatemp[1,`bb'] = e(b)'
			restore
		}
		
		* Adjust the order of coefficients of independent variables
		mat thetahat = (theta[`d',1]\theta[1..`d'-1,1])
		mat thetadagger = (thetatemp[`d',1..`B']\thetatemp[1..`d'-1,1..`B'])
		mat Sigma=(thetadagger-thetahat*J(1,`B',1))*(thetadagger-thetahat*J(1,`B',1))'/`B'
		
		* CALL MYFUN_HOM_NEW()
			
		mat pf = J(`d',`J'+1,0)
		qui qreg_simplified12 `varlist',quantile(`tau') cformat(%10.0g)
		mat pf[1,1] = e(b)'
		forvalues j = 1/`J'{
			local tau1 = `tau'*l[1,`j']
			qui qreg_simplified12 `varlist',quantile(`tau1') cformat(%10.0g)
			mat pf[1,`j'+1] = e(b)'		
		} 
		
		local sigma Sigma
		local ll l
		* The order of pf is adjusted in myfun_hom_new()
	    local PF  pf
		mata: myfun_hom_new( "`phi'","`sigma'", "`ll'","`PF'", `d')
		mat beta = e(beta)
		local dis_b = e(dis_b)

		* Use subsample to compute the point estimator and J-test statistic
		forvalues bb = 1/`B'{
			preserve
			qui sample `boots',count
			
			* CALL MYFUN_HOM_NEW()
			
			mat pf = J(`d',`J'+1,0)
			qui qreg_simplified12 `varlist',quantile(`tau') cformat(%10.0g)
			mat pf[1,1] = e(b)'
			forvalues j = 1/`J'{
				local tau1 = `tau'*l[1,`j']
				qui qreg_simplified12 `varlist',quantile(`tau1') cformat(%10.0g)
				mat pf[1,`j'+1] = e(b)'		
			} 

			mata: myfun_hom_new("`phi'","`sigma'", "`ll'","`PF'", `d')
			mat par_bootb[1,`bb'] = e(beta)
			mat chi_bootbb[`bb',1] = e(dis_b)
			restore
		}
			
		local Chi_bootbb  chi_bootbb
		local Par_bootb  par_bootb
		mata: IC("`Chi_bootbb'", "`Par_bootb'",`J',`dbb',`B',`boots',`N',`tau')
		local disbbtemp = e(disbbtemp)
		if(`disbbtemp' < disbb){
			scalar disbb = `disbbtemp'
			mat Sigmahat = Sigma
			mat beta_hom = beta
			local chibb = `dis_b'
			local tau0 = `tau'
			
		}
 	}

	
		
	* BOOTSTRAP THE CONFIDENCE INTERVAL
	
	local sigmahat Sigmahat
	mat par_bootstraphomb = J(`dbb',`B',0)
	forvalues bb = 1/`B'{
			preserve
			qui bsample 
			
			
			* CALL MYFUN_HOM_NEW()
			
			mat pf = J(`d',`J'+1,0)
			qui qreg_simplified12 `varlist',quantile(`tau0') cformat(%10.0g)
			mat pf[1,1] = e(b)'
			forvalues j = 1/`J'{
				local tau1 = `tau0'*l[1,`j']
				qui qreg_simplified12 `varlist',quantile(`tau1') cformat(%10.0g)
				mat pf[1,`j'+1] = e(b)'		
			} 
			
			local PF  pf
			mata: myfun_hom_new("`phi'","`sigmahat'", "`ll'","`PF'", `d')
			mat par_bootstraphomb[1,`bb'] = e(beta)
			
			restore
	}
	local Par_bootstraphomb par_bootstraphomb
	mata: stat("`Par_bootstraphomb'", `J', `dbb',`chibb')
	mat std_b = e(std_b)
	mat V = diag(std_b)
	local specificationtest e(specificationtest)
	
	di in gr "End"
	di "$S_TIME"

	
	
	* RETURNS IN ECLASS
	ereturn scalar tau0 = `tau0'
	ereturn scalar specificationtest =`specificationtest'
	* DISPLAY
	di in gr "Optimal quantile index = " %10.0g e(tau0)
	//di in gr "Bootstrapped standard deviation = " %10.0g e(std_b)
	di in gr "Specification test = " %10.0g e(specificationtest)
	
	* Display the results in a table
	mat Phi = `phi'
	tokenize `varlist'
	local i=1
	mat colnames beta_hom = "`1'"
	mac shift
	while "`1'" != "" {
		scalar index =  Phi[1,`i']
		if (index == 1){
			local names `names' `1'
		}
		mac shift
		local i=`i'+1
	}
	
	mat rownames V = `names'
	mat colnames V = `names'
	mat rownames beta_hom = `names'
	mat b = beta_hom'
	_coef_table , bmatrix(b) vmatrix(V)
	
	* RETURNS IN ECLASS
	* MATRIX
	ereturn matrix beta_hom = beta_hom
	ereturn matrix std_b = std_b
	ereturn matrix v = V
		
	* MATRIX
	//matlist e(beta_hom)
	//matlist e(std_b)
	
end

//-----------------------------------------------------------------------
//qreg_simplified12
capture program drop qreg_simplified12
program qreg_simplified12, eclass byable(recall) sort prop(sw mi)
	version 6, missing
	local options "Level(cilevel)"
	if !replay() {
		local cmdline : copy local 0
		syntax varlist [aw fw] [if] [in] [, `options' /*
			*/ noConstant Quantile(real 0.5) /*
			*/ WLSiter(integer 1) noLOg * ]
		_get_diopts diopts options , `options'
		if "`constan'"!="" {
			di in red "nocons invalid"
			exit 198
		}
		if `quantil' >= 1 {
			local quant = `quantil'/100
		}
		else	local quant "`quantil'"
		if `quant' <= 0 | `quant' >= 1 {
			di in red "quantiles(`quantil') out of range"
			exit 198
		}
		if `wlsiter'<1 { error 198 }
		marksample touse
		qui count if `touse'
		if r(N)<2 { error cond(_n,2001,2000) }
		gettoken dep indep : varlist
		_rmcoll `indep' [`weight'`exp'] if `touse'
		local varlist `dep' `r(varlist)'
		tempvar r s2 p
		gen long `s2' = _n

		/* initial estimates via weighted least squares 	*/
		_qregwls `varlist' [`weight' `exp'] if `touse',	///
			iterate(`wlsiter') quant(`quant') r(`r') `log'

		quietly {
			if "`log'"=="" { local log "noisily" }
			else 	local log

			sort `r' `s2'
			drop `r'

			`log' _qreg `varlist' if `touse' [`weight'`exp'], /*
				*/ quant(`quant') `options'

end

//------------------------------------------------------------------------------
//Mata part: main function
mata:
// version 13
mata clear
mata set matastrict on

//------------------------------------------------------------------------------
//myfun_hom_new.m
void myfun_hom_new( string matrix Phi, string matrix sigma,string matrix ll, ///
					  string matrix PF, real scalar dd)
{

/*******************************************************************************\
*** Input                                                                       *
* myfun_hom_new computes the point estimator and the J-test statistic for a     *
  given quantile index tau.                                                     *
* tau: Quantile                                                                 *
* phi: User supplied vector to indicate which variable (in our application, it  *
  is "Black") is assumed to be homoskedastic. (The code allows for more than one*
  covariates to be homoskedastic.)                                              *
* X: Dependent variables                                                        *
* Y: Independent variable                                                       *
* l: Equations we will explore by minimum distance estimations are indexed by l.*
 In general, we will use quantile level tau, tau*l to estimate beta and delta   *
* Sigma: An estimator for \omega_0 in the paper                                 *
*** Output                                                                      *
* beta: The point estimator                                                     *
* dis_b: J-test statistic for this specification.                               *
\*******************************************************************************/

	real scalar  JJ, dbb, i, j
	real matrix  l, phi, Sigma, pf
	l = st_matrix(ll)
	JJ = cols(l)
	phi = st_matrix(Phi)
	dbb = sum(phi)
	Sigma = st_matrix(sigma)
	pf = st_matrix(PF)

    //Convert vector phi into a matrix which picks out the covariates that are 
	//homoskedastic.
	real scalar Count, counta
	real matrix phitemp
	phitemp = J(dbb,dd-1,0)
	Count = 1
	for(i=1;i<=dbb;i++){
		counta = 1
		for(j=Count;j<=dd-1;j++){
			if (phi[j] == 1  && counta == 1) {
				phitemp[i,j] = 1
				counta = counta + 1 
				Count = j + 1
			}
		}
	}
	
	
	//Compute matrix L, Gamma_2, Gamma_3 in the paper.
	real matrix l1, L, Gamma2, Gamma3
	l1 = (1,l)
	L = J(JJ+1,JJ+1,0)
	for(i = 1; i <= JJ+1; i++){
		for(j = 1; j <= JJ+1; j++){
			L[i,j] = min((l1[i],l1[j]))/sqrt(l1[i] * l1[j])
		}
	}
	
	Gamma2 = (J(dd-1,1,0),I(dd-1))
	Gamma3 = diag(1:/sqrt(l1))

	// The first step: extremal quantile regression
	real matrix tempy2
	tempy2 = J(dbb*(JJ+1),1,0)
	// Be careful aboout the subscript of pf
	tempy2[|1,1\dbb,1|] = phitemp * pf[|1,1\dd-1,1|]
	

	
	for(j=1;j<=JJ;j++){
		// Be careful aboout the subscript of pf
		tempy2[dbb*j+1..dbb*(j+1)] = phitemp * pf[|1,j+1\dd-1,j+1|]
	}
	// The second step: minimum distance estimation 
	real matrix omega_0, W2, mom2,  GpG, beta
	real scalar dis_b
	omega_0 = Sigma
	GpG = Gamma3 # (phitemp * Gamma2)

	W2 = invsym(GpG * (L # omega_0) * GpG' )                                     //W2 is the optimal weighting matrix for homo beta

	beta = - invsym(J(JJ+1,1,I(dbb))' * W2 * J(JJ+1,1,I(dbb))) * J(JJ+1,1,I(dbb))' * W2 * tempy2
	mom2 = J(JJ*dbb, 1, 0)
	for(j = 1; j <= JJ; ++j){
		// Be careful aboout the subscript of pf
		mom2[dbb*(j-1)+1..dbb*j,1] = phitemp * pf[|1,j\dd-1,j|] + beta             //mom2 is the value of moments for homo beta
	}
	// Be careful aboout the subscript of pf
	mom2 = mom2\(phitemp * pf[|1,JJ+1\dd-1,JJ+1|] + beta )
	dis_b =  mom2' * W2 * mom2                                                  //here we compute the distance of homo beta evaluated at extremal quantile estimator of delta
	
	st_matrix("e(beta)",beta)
	st_numscalar("e(dis_b)",dis_b)
	
}

//------------------------------------------------------------------------------
//
void IC(string matrix Chi_bootbb, string matrix Par_bootb, real scalar J, ///
			real scalar dbb, real scalar B, real scalar boots, real scalar N, ///
			real scalar tau){
	real scalar median_b1, disbbbias, disbbvar, disbbtemp
	real matrix tempmean, chi_bootbb, par_bootb
	chi_bootbb = st_matrix(Chi_bootbb)
	par_bootb = st_matrix(Par_bootb)
	
	median_b1 = mm_median(rchi2(100000,1,J*dbb))
	disbbbias=abs(mm_median(chi_bootbb)*boots/N-median_b1)'
	tempmean = J(1,B,mean(par_bootb')')
	disbbvar = mean(colsum((par_bootb-tempmean):^2)')
	disbbtemp = disbbvar*boots/N + disbbbias/sqrt(tau*boots)
	st_numscalar("e(disbbtemp)",disbbtemp)
}

void stat(string matrix Par_bootstraphomb, real scalar J, real scalar dbb, ///
		real scalar chibb	){
	real scalar  specificationtest
	real matrix std_b,par_bootstraphomb
	par_bootstraphomb = st_matrix(Par_bootstraphomb)
	std_b = mm_colvar(par_bootstraphomb')'                                      //Compute the standard deviation.
	std_b = sqrt(std_b)
	specificationtest = 1 - chi2(J*dbb,chibb)
	st_matrix("e(std_b)",std_b)
	st_numscalar("e(specificationtest)",specificationtest)
}

 
end
