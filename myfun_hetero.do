*!Program myfun_hetero, Xiaoyun Qiu 30july2015

//------------------------------------------------------------------------------
//compute tempx1 for the mata function
capture program drop myfun_hetero
program myfun_hetero,eclass 
	version 13
	syntax varlist(min=2 numeric) [if] [in] , Tau(real) 
	marksample touse
	tokenize `varlist'
	
	local indevar: list varlist - 1
	
	scalar m = 1.2
	scalar b0 = 0
	scalar d0 = 1
	matrix l = (0.65, 0.85, 1.15, 1.45)
	matrix l_col = colsof(l)
	local J = l_col[1,1]
	qui count if `touse'
	scalar T = r(N)
	
	local d: word count `varlist'

	tempvar trow
	local trow = (`d'-1)*`J'
	mat tempy = J(`trow',1,0)
	mat tempx1 = J(`J',1,0)
		
	qui qreg  `varlist',quantile(`tau')
	mat beta1 = e(b)'
	mat pf = beta1
		
	tempvar quant tau1 tau2 
	local tau1 = `tau'	
	local quant = `tau'*m
		
	qui qreg `varlist',quantile(`quant')
	mat beta2 = e(b)'
	mat pm = beta2
		
	forvalues j = 1/`J'{
		local tau1 = `tau'*l[1,`j']
		qui qreg `varlist',quantile(`tau1')
		mat beta1`j' = e(b)'
		mat pf = (pf, beta1`j')			
			
		tempvar ind
		local ind = (`d'-1)*(`j'-1)+1
		mat tempy[`ind',1] = (pf[1..`d'-1,`j'+1] - pf[1..`d'-1,1])*d0
		mat tempx1[`j',1] = pf[`d',`j'+1] - pf[`d',1]
			
		local tau2 = `quant'* l[1,`j']
		qui qreg `varlist',quantile(`tau2')
		mat beta2`j' = e(b)'
		mat pm = (pm, beta2`j')
	} 
		
	mat pf[`d',1] = pf[`d',1...] + J(1,`J'+1,b0)
	
	mata: hetero("`d'","d0","T","m","`tau'","`indevar'","`touse'", "tempx1", "tempy","pf", "pm", "l")
	
	ereturn list
end



//Mata--------------------------------------------------------------------------

capture mata mata drop hetero()
mata:
version 13 
void hetero(string scalar d, string scalar d0, string scalar T, string scalar m, /// 
			string scalar tau, string scalar indevar, string scalar touse, ///
			string scalar tempx1, string scalar tempy, string scalar pf, ///
			string scalar pm, string scalar l)
{	
	real matrix  l1, L, G, tempx, Tri, Xd, tempXd, idx, QH, Qx, omega_0, W, delta, /// 
			Out, mom, Nb, Nb1, Dis, beta, Var, Gamma, ll, mtempx1, mtempy, ///
			ppf, ppm, invQH, kronTri, invkronTri, dtemp, X
	real scalar  JJ, ttau, dd, dd0, TT, mm, s, i, j
	
	st_view(X=.,.,tokens(indevar),touse)
	ll = st_matrix(l)
	ppf = st_matrix(pf)
	ppm = st_matrix(pm)
	mtempx1 = st_matrix(tempx1)
	mtempy = st_matrix(tempy)
	ttau = strtoreal(st_local("tau"))
	dd = strtoreal(st_local("d"))
	dd0 = st_numscalar(d0)
	TT = st_numscalar(T)
	mm = st_numscalar(m)
	
	tempx = mtempx1 # I(dd-1)
	X = (J(TT,1,1),X)
	Qx = X'*X/TT
	JJ = cols(ll)
	Gamma = (-J(JJ,1,I(dd)),diag(1:/sqrt(ll))#I(dd))
	l1 = (1,ll)
	L = J(JJ+1,JJ+1,0)
	for(i = 1; i <= JJ+1; ++i){
		for(j = 1; j <= JJ+1; ++j){
			L[i,j] = min((l1[i],l1[j]))/sqrt(l1[i] * l1[j])
		}
	}
	G = log(ll') # I(dd-1)
		
	dtemp = invsym(tempx' * tempx) * tempx' * mtempy
	
	for(i = 1; i <=2; ++i){
		Tri = (-dtemp:/dd0, I(dd-1))
		Xd = X:/(J(1,dd,abs(X*(dd0, dtemp')')):^.5)
		tempXd = Xd * (dd0\dtemp)
	
		idx = abs(tempXd :- mm_median(tempXd)) :< 3*(mm_quantile(tempXd,1,.75) - mm_quantile(tempXd,1,.25))
		Xd = select(Xd,idx)
		QH = Xd' * Xd/rows(idx)
		invQH = invsym(QH)
		omega_0 = invQH * Qx * invQH
	
		W = (I(JJ) # Tri) * Gamma * (L # omega_0) * Gamma' * (I(JJ) # Tri')
		W = luinv(W)
		dtemp = luinv(tempx'*W*tempx) * tempx' * W * mtempy	
	}
	
	delta = dtemp
	Var = luinv(G'*W*G)
	Out = delta
	
	mom = J(JJ*(dd-1), 1, 0)
	for(j = 1; j <= JJ; ++j){
		mom[(dd-1)*(j-1)+1..(dd-1)*j,1] = (ppf[1..dd-1,j+1] - ppf[1..dd-1,1]) :- (ppf[dd,j+1] - ppf[dd,1]) * delta/dd0
	}
	
	Nb = log(mm) * sqrt(ttau*TT)/mm_median(abs(ppm[dd,1...] - ppf[dd,1...])')
	beta = - ppf[|1,1\dd-1,1|] + delta/dd0 * ppf[dd,1]
	Out = (Out\beta)
	Nb1 = abs(sqrt(ttau*TT))/mm_median(ppf[|dd,1\dd,.|]')
	Dis = Nb:^2 * mom' * W * mom
	
	st_numscalar("e(Nb1)",Nb1)
	st_numscalar("e(Dis)",Dis)
	st_matrix("e(Out)",Out)
	st_matrix("e(Var)",Var)	
	
}
end





