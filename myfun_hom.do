*! Program myfun_hom XQ 3Aug2015

capture program drop myfun_hom
program myfun_hom,eclass
	version 13
	syntax varlist(min=2 numeric) [if] [in], Tau(real) M(real) L(string) Phi(string) [ MSE_delta(string)]
	marksample touse
	tokenize `varlist'
	
	local indevar: list varlist - 1
	local d: word count `varlist'
	qui count if `touse'
	scalar T = r(N)
	
	tempvar J count ppl 
	
	matrix l_col = colsof(`l')
	local J = l_col[1,1]
	matrix psi = J(1,`d'-1,1) - `phi'
	
	mata: getpsiphi("`phi'","psi","`d'")
	
	ereturn list
	
	matrix phitemp = r(phitemp)
	matrix psitemp = r(psitemp)
	scalar dd = r(dd)
	scalar db = r(db)
	
	//------------------
	//the first case
	if dd > 0 {
	mat tempy1 = J(dd*`J',1,0)
	mat tempy2 = J(db*(`J'+1),1,0)
	
	qui qreg  `varlist',quantile(`tau')
	mat beta1 = e(b)'
	mat pf = beta1
	
	tempvar quant tau1 tau2 
	local tau1 = `tau'	
	local quant = `tau' * `m'
	
	qui qreg `varlist',quantile(`quant')
	mat beta2 = e(b)'
	mat pm = beta2
	
	mat tempy2[1,1] = phitemp * pf[1..`d'-1,1]
	mat tempx1 = J(`J',1,0)
		
	forvalues j = 1/`J'{
		local tau1 = `tau'*l[1,`j']
		qui qreg `varlist',quantile(`tau1')
		mat beta1`j' = e(b)'
		mat pf = (pf, beta1`j')			
			
		tempvar ind1 ind2
		local ind1 = dd * (`j' - 1) + 1
		mat tempy1[`ind1',1] = psitemp * (pf[1..`d'-1,`j'+1] - pf[1..`d'-1,1])
		local ind2 = db * `j' +1
		mat tempy2[`ind2',1] = phitemp * pf[1..`d'-1,`j'+1]
		mat tempx1[`j',1] = pf[`d',`j'+1] - pf[`d',1]
			
		local tau2 = `quant'* l[1,`j']
		qui qreg `varlist',quantile(`tau2')
		mat beta2`j' = e(b)'
		mat pm = (pm, beta2`j')
	} 
	if "`mse_delta'" == "" {
	mata: hom_pos("`indevar'", "`touse'", "tempx1", "tempy1","tempy2","phitemp","psitemp","l", "pf","pm", "`tau'","`m'", "T", "`d'", "dd", "db")
	ereturn list
	}
	else{
	mata: hom_pos("`indevar'", "`touse'", "tempx1", "tempy1","tempy2","phitemp","psitemp","l", "pf","pm", "`tau'","`m'", "T", "`d'", "dd", "db", "`mse_delta'")
	}
	}
	//-----------------------------------
	//another case
	else if dd == 0 {
	
	mat tempy2 = J(db*(`J'+1),1,0)
	
	qui qreg  `varlist',quantile(`tau')
	mat beta1 = e(b)'
	mat pf = beta1
	
	tempvar quant tau1 tau2 
	local tau1 = `tau'	
	local quant = `tau' * `m'
	
	qui qreg `varlist',quantile(`quant')
	mat beta2 = e(b)'
	mat pm = beta2
	
	mat tempy2[1,1] = phitemp * pf[1..`d'-1,1]
	mat tempx1 = J(`J',1,0)
		
	forvalues j = 1/`J'{
		local tau1 = `tau'*l[1,`j']
		qui qreg `varlist',quantile(`tau1')
		mat beta1`j' = e(b)'
		mat pf = (pf, beta1`j')			
			
		tempvar  ind2
		local ind2 = db * `j' +1
		mat tempy2[`ind2',1] = phitemp * pf[1..`d'-1,`j'+1]
		mat tempx1[`j',1] = pf[`d',`j'+1] - pf[`d',1]
			
		local tau2 = `quant'* l[1,`j']
		qui qreg `varlist',quantile(`tau2')
		mat beta2`j' = e(b)'
		mat pm = (pm, beta2`j')
	} 
	mata: hom_zero("`indevar'", "`touse'", "tempx1","tempy2", "phitemp" "l", "pf", "pm", "`tau'","`m'", "T", "`d'", "dd", "db")
	
	}
	
end
	
	
capture mata mata drop getpsiphi()
version 13
mata:
mata clear
mata set matastrict on	

//------------------------------------------------------------------------------
//hom_pos compute the estimation results when dd > 0
void hom_pos(string scalar indevar, string scalar touse, string scalar mtempx1, ///
		string scalar mtempy1, string scalar mtempy2, string scalar pphitemp, /// 
		string scalar ppsitemp, ///
		string scalar ll, string scalar ppf, string scalar ppm, string scalar tau, ///
		string scalar m, string scalar TT, string scalar d, string scalar ddd, ///
		string scalar ddb, | string scalar delta_mse)
{
	real matrix X, l, pf, pm, tempx1, tempy1, tempx, l1, Qx, L, psitemp, Gamma, ///
				Gamma2, Gamma3, Gd, Gb, tempy2, dtemp, dtemp2, Tri, Xd, tempXd, ///
				idx, QH, omega_0, W1, delta, psiTri, Vd, W2, beta, Vb, mom1, ///
				mom2, par, GpG, invQH, phitemp
	real scalar ttau, dv, dd, dbb, T, mm, JJ, i, j, Nb, dis_d, dis_b
	
	st_view(X=., .,tokens(indevar),touse)
	l = st_matrix(ll)
	pf = st_matrix(ppf)
	pm = st_matrix(ppm)
	tempx1 = st_matrix(mtempx1)
	tempy1 = st_matrix(mtempy1)
	tempy2 = st_matrix(mtempy2)
	psitemp = st_matrix(ppsitemp)
	phitemp = st_matrix(pphitemp)
	ttau = strtoreal(st_local("tau"))
	dv = strtoreal(st_local("d"))
	dd = st_numscalar(ddd)
	dbb = st_numscalar(ddb)
	T = st_numscalar(TT)
	mm = strtoreal(st_local("m"))
	
	X = (J(T,1,1),X)
	Qx = X'*X/T
	JJ = cols(l)
	
	l1 = (1,l)
	L = J(JJ+1,JJ+1,0)
	for(i = 1; i <= JJ+1; i++){
		for(j = 1; j <= JJ+1; j++){
			L[i,j] = min((l1[i],l1[j]))/sqrt(l1[i] * l1[j])
		}
	}
	Gamma = (-J(JJ,1,I(dv)),diag(1:/sqrt(l))#I(dv))
	Gamma2 = (J(dv-1,1,0),I(dv-1))
	Gamma3 = diag(1:/sqrt(l1))
	
	Gd = log(l') # I(dd)
	Gb = - J(JJ+1,1,1) # I(dbb)
	
	tempx = tempx1 # I(dd)	
	if (args() < 17){
		dtemp = luinv(tempx' * tempx) * tempx' * tempy1
	}
	else{
		dtemp = st_matrix("delta_mse")	
	}
	dtemp2 = dtemp' * psitemp
	
	for(i = 1; i <=2; ++i){
		Tri = (-dtemp2', I(dv-1))
		Xd = X:/(J(1,dv,abs(X*(1, dtemp2)')):^.5)
		tempXd = Xd[.,dv]
	
		idx = abs(tempXd :- mm_median(tempXd)) :< 3*(mm_quantile(tempXd,1,.75) - mm_quantile(tempXd,1,.25))
		Xd = select(Xd,idx)
		QH = Xd' * Xd/rows(idx)
		invQH = invsym(QH)
		omega_0 = invQH * Qx * invQH
		
		psiTri = psitemp * Tri
		W1 = (I(JJ) # psiTri) * Gamma * (L # omega_0) * Gamma' * (I(JJ) # psiTri')
		W1 = luinv(W1)
		delta = luinv(tempx'*W1*tempx) * tempx' * W1 * tempy1
		dtemp = delta
		dtemp2 = dtemp' * psitemp
	}		
	
	Vd = luinv(Gd' * W1 * Gd)
	GpG = Gamma3 # (phitemp * Gamma2)
	W2 = luinv(GpG * (L # omega_0) * GpG' )
	beta = - luinv(J(JJ+1,1,I(dbb))' * W2 * J(JJ+1,1,I(dbb))) * J(JJ+1,1,I(dbb))' * W2 * tempy2
	Vb = luinv(Gb' * W2 * Gb)
	
	mom1 = J(JJ*dd, 1, 0)
	mom2 = J(JJ*dbb, 1, 0)
	for(j = 1; j <= JJ; j++){
		mom1[dd*(j-1)+1..dd*j,1] = psitemp * (pf[1..dv-1,j+1] - pf[1..dv-1,1]) :- (pf[dv,j+1] - pf[dv,1]) * delta
		mom2[dbb*(j-1)+1..dbb*j,1] = phitemp * pf[1..dv-1,j] + beta
	}
	mom2 = mom2\(phitemp * pf[1..dv-1,JJ+1] + beta)
	Nb = log(mm) * sqrt(ttau*T)/mm_median(abs(pm[dv,1...] - pf[dv,1...])')
	
	dis_d = Nb^2 * mom1' * W1 * mom1
	dis_b = Nb^2 * mom2' * W2 * mom2
	par = delta
	par = par\beta
	
	st_matrix("e(par)",par)
	st_matrix("e(Vb)",Vb)
	st_numscalar("e(Vd)",Vd)
	st_numscalar("e(Nb)",Nb)
	st_numscalar("e(dis_d)",dis_d)
	st_numscalar("e(dis_b)",dis_b)
	
}

//------------------------------------------------------------------------------
//hom_zero compute the estimation results when dd == 0
void hom_zero(string scalar indevar, string scalar touse, string scalar mtempx1, ///
		 string scalar mtempy2, string scalar pphitemp, string scalar ll, ///
		 string scalar ppf, string scalar ppm, string scalar tau, ///
		string scalar m, string scalar TT, string scalar d, string scalar ddd, ///
		string scalar ddb)
{
	real matrix X, l, pf, pm, tempx1, l1, Qx, L,  phitemp, Gamma, Gamma2, ///
				Gamma3, Gd, Gb, tempy2, Vd, W2, beta, Vb, mom2, par, GpG
	real scalar ttau, dv, dd, dbb, T, mm, JJ, i, j, Nb, dis_d, dis_b
	
	st_view(X=., .,tokens(indevar),touse)
	l = st_matrix(ll)
	pf = st_matrix(ppf)
	pm = st_matrix(ppm)
	tempx1 = st_matrix(mtempx1)
	tempy2 = st_matrix(mtempy2)
	phitemp = st_matrix(pphitemp)
	ttau = strtoreal(st_local("tau"))
	dv = strtoreal(st_local("d"))
	dd = st_numscalar(ddd)
	dbb = st_numscalar(ddb)
	T = st_numscalar(TT)
	mm = strtoreal(st_local("m"))
	
	X = (J(T,1,1),X)
	Qx = X'*X/T
	JJ = cols(l)
	
	l1 = (1,l)
	L = J(JJ+1,JJ+1,0)
	for(i = 1; i <= JJ+1; i++){
		for(j = 1; j <= JJ+1; j++){
			L[i,j] = min((l1[i],l1[j]))/sqrt(l1[i] * l1[j])
		}
	}
	
	Gamma = (-J(JJ,1,I(dv)),diag(1:/sqrt(l))#I(dv))
	Gamma2 = (J(dv-1,1,0),I(dv-1))
	Gamma3 = diag(1:/sqrt(l1))
	
	Gb = - J(JJ+1,1,1) # I(dbb)
	
	GpG = Gamma3 # (phitemp * Gamma2)
	W2 = luinv(GpG * (L # invsym(Qx)) * GpG' ) 
	beta = - luinv(J(JJ+1,1,I(dbb))' * W2 * J(JJ+1,1,I(dbb))) * J(JJ+1,1,I(dbb))' * W2 * tempy2
	Vb = luinv(Gb' * W2 * Gb)
	
	mom2 = J(JJ*dbb, 1, 0)
	for(j = 1; j <= JJ; ++j){
		mom2[dbb*(j-1)+1..dbb*j,1] = phitemp * pf[1..dv-1] + beta
	}
	mom2 = mom2\(phitemp * pf[1..dv-1,JJ+1] + beta )
	Nb = log(mm) * sqrt(ttau*T)/mm_median(abs(pm[dv,1...] - pf[dv,1...])')
	
	dis_d = 0
	dis_b = Nb^2 * mom2' * W2 * mom2
	
	par = beta
	Vd = 0
	
	st_matrix("e(par)",par)
	st_matrix("e(Vb)",Vb)
	st_numscalar("e(Nb)",Nb)
	st_numscalar("e(Vd)",Vd)
	st_numscalar("e(dis_d)",dis_d)
	st_numscalar("e(dis_b)",dis_b)

}
//------------------------------------------------------------------------------
//return psi and phi
void getpsiphi(string scalar phi, string scalar psi, string scalar d)
{

	real matrix pphi, ppsi, phitemp, psitemp
	real scalar	dv, db, dd, count, i, j
	
	pphi = st_matrix(st_local("phi"))
	ppsi = st_matrix(psi)
	dv = strtoreal(st_local("d"))
	
	db = rowsum(pphi)
	dd = rowsum(ppsi)
	count = 1
	
	phitemp = J(db,dv-1,0)
	psitemp = J(dd,dv-1,0)
	
	for(i=1;i<=db;i++){
		for(j=count;j<=dv-1;j++){
			if (pphi[j] == 1) {
				phitemp[i,j] = 1
				count = j + 1
				break
			}
		}
	}
	count = 1
	for(i=1;i<=dd;i++){
		for(j=count;j<=dv-1;j++){
			if (ppsi[j] == 1) {
				psitemp[i,j] = 1
				count = j + 1
				break
			}
		}
	}
	
	st_matrix("r(psitemp)",psitemp)
	st_matrix("r(phitemp)",phitemp)
	st_numscalar("r(db)",db)
	st_numscalar("r(dd)",dd)
	st_numscalar("e(d)",dv)
	
}
end
