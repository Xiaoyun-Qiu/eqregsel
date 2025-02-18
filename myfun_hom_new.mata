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

	W2 = luinv(GpG * (L # omega_0) * GpG' )                                     //W2 is the optimal weighting matrix for homo beta

	beta = - luinv(J(JJ+1,1,I(dbb))' * W2 * J(JJ+1,1,I(dbb))) * J(JJ+1,1,I(dbb))' * W2 * tempy2
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
end
mata: mata mosave myfun_hom_new(),dir(PERSONAL) replace
