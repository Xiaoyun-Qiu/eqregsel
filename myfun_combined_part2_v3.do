/************************************************\
* program myfun_combined Xiaoyun Qiu 5Aug2015    *
* Modified in 28, May, 2016                      *
* version 3.0                                    *
* Part 2: Compute the standard deviation         *
\************************************************/
//------------------------------------------------------------------------------
//Stata program: creates a new command
capture program drop std_b
program std_b,eclass
//	version 13
	syntax varlist(min=2 numeric) [if] [in] , TAU(real) PHI(string) SIGMAhat(string) 
	marksample touse
	tokenize `varlist'
	
	mata:std_b("`varlist'","`touse'","`phi'","`sigmahat'",`tau')
	
	ereturn list

end


//------------------------------------------------------------------------------
//Mata part: main function
mata:
// version 13
mata clear
mata set matastrict on

void std_b(string scalar varlist,string scalar touse, string matrix Phi, ///
			string scalar Sigmahat,	real scalar tau)
{

/*******************************************************************************\
*** Input                                                                       *
* myfun_combined computes the homoskedastic beta and the J-test statistic.      * 
* X: Dependent variables                                                        *
* Y: Independent variable                                                       *
* phi: User supplied vector to indicate which variable (in our application, it  *
  is "Black") is assumed to be homoskedastic. (The code allows for more than one*
  covariates to be homoskedastic.)                                              *
* tau0
* sigmahat                                              *
*** Output                                                                      *
* std_b: Standard deviation for the point estimator.                            *
\******************************************************************************/
	
	real matrix Z, X, Y, phi, sigmahat
	st_view(Z=.,.,tokens(varlist),touse)
	Y = Z[.,1]
	X = Z[|1,2\.,.|]
	phi = st_matrix(Phi)
	sigmahat = st_matrix(Sigmahat)
	
    real scalar N, dd, JJ, dbb, B, ss
	real matrix  l
	N = rows(X)
	X = (J(N,1,1) , X)                                                          //Note that users would not specify the constant term when they call the function, so we need to add a constant term here.
	dd = cols(X)
	B = 3                                                                  //B is the number of replications in our subsampling and bootstrap method.
	l = (0.65,0.85,1.15,1.45)                                                   //l is corresponding to the moment we use.
	JJ = cols(l)
	dbb = sum(phi)
   
	Z = (Y,X)                                                                   //Now there is a constant term in X, and thus in Z
	//bootstrap the confidence interval
	real matrix par_bootstraphomb,temp, Zz, tempidz
	real scalar tempv
	par_bootstraphomb = J(dbb,B,0)
	
	for(ss=1;ss<=B;ss++){                                                       //Bootstrap the first stage estimator, use it to compute omega_0.
		tempidz = rdiscrete(N,1,J(N,1,1/N))
		Zz = Z[tempidz,.]
		temp=(0)
		tempv=0
		myfun_hom_new(tau,Zz[|1,2\.,.|],Zz[|1,1\.,1|],l,phi,sigmahat,temp,tempv)
		par_bootstraphomb[|1,ss\.,ss|] = temp
	}

	real scalar std_b, specificationtest
	std_b = mm_colvar(par_bootstraphomb')'                                      //Compute the standard deviation.
	
	st_numscalar("e(std_b)",std_b)

}

//------------------------------------------------------------------------------
//myfun_hom_new.m
void myfun_hom_new(real scalar tau, real matrix X, real matrix Y, real matrix l, ///
                   real matrix phi, real matrix Sigma, real matrix beta, ///
				   real scalar dis_b)
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

	real scalar T, dd, JJ, dbb, i, j
	T = rows(X)
	dd = cols(X)
	JJ = cols(l)
	dbb = sum(phi)

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
	real matrix tempy2, pf
	tempy2 = J(dbb*(JJ+1),1,0)
	pf = J(dd,JJ+1,0)
	pf[.,1] = rq_fnm(X,-Y,tau)
	tempy2[|1,1\dbb,1|] = phitemp * pf[|2,1\.,1|]
	
	for(j=1;j<=JJ;j++){
		pf[.,j+1] = rq_fnm(X,-Y,tau*l[j])
		tempy2[dbb*j+1..dbb*(j+1)] = phitemp * pf[|2,j+1\.,j+1|]
	}

	// The second step: minimum distance estimation 
	real matrix omega_0, W2, mom2,  GpG
	omega_0 = Sigma
	GpG = Gamma3 # (phitemp * Gamma2)
	W2 = luinv(GpG * (L # omega_0) * GpG' )                                     //W2 is the optimal weighting matrix for homo beta
	beta = - luinv(J(JJ+1,1,I(dbb))' * W2 * J(JJ+1,1,I(dbb))) * J(JJ+1,1,I(dbb))' * W2 * tempy2
	mom2 = J(JJ*dbb, 1, 0)
	for(j = 1; j <= JJ; ++j){
		mom2[dbb*(j-1)+1..dbb*j,1] = phitemp * pf[|2,j\.,j|] + beta             //mom2 is the value of moments for homo beta
	}
	mom2 = mom2\(phitemp * pf[|2,JJ+1\.,JJ+1|] + beta )
	dis_b =  mom2' * W2 * mom2                                                  //here we compute the distance of homo beta evaluated at extremal quantile estimator of delta

}

//------------------------------------------------------------------------------
//rq_fnm

real matrix rq_fnm(real matrix X, real matrix y, real scalar q)
{
	real scalar m
	m = rows(X)
		
	real matrix u, a, b
	u = J(m,1,1)
	a = (1 - q) :* u
	b = -lq_fnm(X', -y', X' * a, u, a)'
	
	return (b)
}
//------------------------------------------------------------------------------
//
real matrix lq_fnm(A,c,b,u,x )
{
	// Set some constants
	real scalar beta, small, max_it, m, n
	beta = 0.9995
	small = 1e-5
	max_it = 50
	m = rows(A)
	n = cols(A)
	
	//generate initial feasible point
	real matrix s, y, r, z, w
	real scalar gap
	s =  u - x
	y = (invsym(A * A') * A * c')'
	r = c - y * A
	r = r + 0.001 * (r :== 0)
	z = r :* (r :> 0)
	w = z - r
	gap = c * x - y * b + w * u
	
	//Start iterations
	real scalar it
	real matrix q, Q, AQ, rhs, dy, dx, dss, dz, dw, fx, fs, fxfs, fw, fz, fwfz, ///
				fpp, fd, mu, g, dxdz, dsdw, xinv, sinv, xii
	it = 0
	while((gap > small) & (it < max_it)){
		it = it + 1
		
		//Compute affine step
		q = 1 :/ (z' :/ x + w' :/ s)
		r = z - w
		Q = diag(sqrt(q))
		AQ = A * Q
		rhs = Q * r'
		dy = (invsym(AQ * AQ') * AQ * rhs)'
		dx = q :* (dy * A - r)'
		dss = -dx
		dz = -z :* (1 :+ dx :/ x)'
		dw = -w :* (1 :+ dss :/ s)'
		
		//Compute maximum allowable step lengths
		fx = bound(x,dx)
		fs = bound(s,dss)
		fw = bound(w,dw)
		fz = bound(z,dz)
		fxfs = (fx, fs)
		fpp = rowmin(fxfs)
		fwfz = (fw, fz)
		fd = rowmin(fwfz)
		fpp = min((min(beta :* fpp),1))
		fd = min((min(beta :* fd),1))

	
	//If full step is feasile, take it. Otherwise modeify it
	if (min((fpp,fd)) < 1){
		//Update mu
		mu = z * x + w * s
		g = (z + fd * dz) * (x + fpp *dx) + (w + fd * dw) * (s + fpp * dss)
		mu = mu * (g/mu)^3 /(2*n)
		
		//Compute modified step
		dxdz = dx :* dz'
		dsdw = dss :* dw'
		xinv = 1 :/ x
		sinv = 1 :/ s
		xii = mu * (xinv - sinv)
		rhs = rhs + Q * (dxdz - dsdw - xii)
		dy = (invsym(AQ * AQ') * AQ * rhs)'
		dx = q :* (A' * dy' + xii - r' - dxdz +dsdw)
		dss = -dx
		dz = mu * xinv' - z - xinv' :* z :* dx' - dxdz'
		dw = mu * sinv' - w - sinv' :* w :* dss' - dsdw'
		
		//Compute maximum allowable step lengths
		fx = bound(x,dx)
		fs = bound(s,dss)
		fw = bound(w,dw)
		fz = bound(z,dz)
		fpp = bound(fx,fs)
		fxfs = (fx, fs)
		fpp = rowmin(fxfs)
		fwfz = (fw, fz)
		fd = rowmin(fwfz)
		fpp = min((min(beta :* fpp),1))
		fd = min((min(beta :* fd),1))
		
		}
	
	//Take the step
	x = x + fpp * dx
	s = s + fpp * dss
	y = y + fd * dy
	w = w + fd * dw
	z = z + fd * dz
	gap = c * x - y * b + w * u
	
	}
	return (y)
}
//------------------------------------------------------------------------------
//
real matrix bound(x,dx)
{
	real matrix b, f
	b = 1e20 :+ (0 :* x)
	f = (dx :< 0)
	
	real matrix ind
	ind = selectindex(f)
	b[ind] = - x[ind] :/ dx[ind]
	
	return (b)
}
end
