*! program myfun_combined Xiaoyun Qiu 5Aug2015


//------------------------------------------------------------------------------
//Stata program: creates a new command
capture program drop myfun_combined_ori
program myfun_combined_ori,eclass
	version 13
	syntax varlist(min=2 numeric) [if] [in] 
	marksample touse
	tokenize `varlist'
	
	mata: myfun_combined("`varlist'","`touse'")
	ereturn list

end


//------------------------------------------------------------------------------
//Mata part: main function
mata:
version 13
mata clear
mata set matastrict on

//void myfun_combined(string scalar varlist, string scalar touse)
void myfun_combined(string scalar varlist,string scalar touse)
{
	real matrix Z, X, Y
	real scalar dv, N
	st_view(Z=.,.,tokens(varlist),touse)
	//st_view(Z=.,.,tokens(varlist))
	X = Z[|1,2\.,.|]
	N = rows(X)
	X = (J(N,1,1) , X)
	Y = Z[.,1]
	dv = cols(X)
	
	//--------------------------------------------------------------------------
	//defining parameters
	real scalar G, JJ, B, boots, mm, b0, d0, lower, upper, step, dd, dbb, i, j, b, gg
	real matrix l, ciboot1, phi, nb_mse, disd, disbd, disbb, disdbias, disbbbias, disbdbias
				
	G = 2
	l = (0.65,0.85,1.15,1.45)
	JJ = cols(l)
	B = 10
	boots = 500
	mm = 1.2
	b0 = 0
	d0 = 1
	
	ciboot1 = J(2*dv-2,2,0)
	
	lower = min((80/boots, 0.1))
	//lower = 0.05
	upper = 0.3
	step = (upper - lower)/G
	
	phi = (1,1,0,1,1)
	dbb = sum(phi)
	dd = dv - 1 - dbb
	
	nb_mse = J(G,1,0)
	disd = J(G,1,0)
	disbd = disd
	disbb = disd
	disdbias = J(G,1,0)
	disbbbias = disdbias
	disbdbias = disdbias
	
	pointer(pointer(real colvector)) colvector par_boot, par_bootb
	pointer(real matrix) rowvector V_homb, V_homd, V_mse
	par_boot = J(B,1,NULL)
	par_bootb = J(B,1,NULL)
	V_homb = J(1,G,NULL)
	V_homd = J(1,G,NULL)
	V_mse = J(1,G,NULL)
	
	for(gg=1;gg<=G;gg++){
		V_homb[gg] = &(0)
		V_homd[gg] = &(0)
		V_mse[gg] = &(0)
	}
	
	real matrix par_mse, par_hom, Nb_hom, chid, chibb
	real scalar median_d1, median_d2, median_b1
	par_mse = J(2*(dv-1),G,0)
	par_hom = J(dv-1,G,0)
	Nb_hom = J(G,1,0)
	chid = J(G,1,0)
	chibb = J(G,1,0)
	nb_mse = J(G,1,0)
	median_d1 = mm_median(rchi2(100000,1,(JJ-1)*(dv-1)))
	median_d2 = mm_median(rchi2(100000,1,(JJ-1)*dd))
	median_b1 = mm_median(rchi2(100000,1,JJ*dbb))
	
	real scalar nbtemp, tau, tempdis_d
	real matrix tempmse
	for(gg=1;gg<=G;gg++){
		tau = lower + step * gg
		myfun_hetero(tau,mm,b0,d0, X,Y,l,par_mse[.,gg],chid[gg],*V_mse[gg],nb_mse[gg])
		tempmse = par_mse[|1,gg\dv-1,gg|]
		myfun_hom(tau,mm,phi,X,Y,l,par_hom[.,gg],tempdis_d,chibb[gg],*V_homd[gg],*V_homb[gg],Nb_hom[gg])
	}
	
	//--------------------------------------------------------------------------
	//line 124
	real matrix tempchi_bootd, tempchi_bootbd, tempchi_bootbb
	tempchi_bootd = J(B,G,0)
	tempchi_bootbd = J(B,G,0)
	tempchi_bootbb = J(B,G,0)
	
	real matrix Zz, idz, tempidz, seltempidz, tempboot, tempvmse, ///
				 tempvhomb, temppar, tempout
	real scalar  tempnbmse, tempvhomd,tempnbhom
	temppar = (0)
	tempboot = J(3*(dv-1),G,0)
	for(b=1;b<=B;b++){
		idz = 1::N
		tempidz = jumble(idz)
		seltempidz = tempidz[|1,1\boots,1|]
		Zz = Z[seltempidz,.]
		X = Zz[|1,2\.,.|]
		X = (J(boots,1,1), X)
		Y = Zz[|1,1\.,1|]
		par_boot[b] = &(J(G,1,NULL))
		par_bootb[b] = &(J(G,1,NULL))
		for(gg=1;gg<=G;gg++){	
			tau = lower + step * gg
			myfun_hetero(tau,mm,b0,d0,X,Y,l,tempboot[|1,gg\2*(dv-1),gg|],tempchi_bootd[b,gg])
			(*par_boot[b])[gg] = &(tempboot[|1,gg\dv-1,gg|])
		    myfun_hom(tau,mm,phi,X,Y,l,tempboot[|2*(dv-1)+1,gg\3*(dv-1),gg|],tempchi_bootbd[b,gg],tempchi_bootbb[b,gg])
			(*par_bootb[b])[gg] = &(tempboot[|2*(dv-1)+1,gg\.,gg|])
		}

	}
	
	disdbias[|1,1\.,1|] = abs(mm_median(tempchi_bootd) :- median_d1)'
	disbdbias[|1,1\.,1|] = abs(mm_median(tempchi_bootbd) :- median_d2)'
	disbbbias[|1,1\.,1|] = abs(mm_median(tempchi_bootbb) :- median_b1)'	
	
	//--------------------------------------------------------------------------
	//line 153
	real matrix disdvar, disbdvar, disbbvar
	disdvar = J(G,1,0)
	disbdvar = disdvar
	disbbvar = disdvar
	
	real matrix tempboot1, tempboot2, tempbootd, tempbootb, temphom, ///
				temphomd, temphomb, w, beta_mse, ///
				delta_mse, delta_hom, beta_hom
	real scalar gstard, gstarbd, gstarbb
	w = (0,0)
	tempboot1 = J(dv-1,G,0)
	tempboot2 = J(dv-1,G,0)
	for(gg=1;gg<=G;gg++){
			tempboot1 = *(*par_boot[1])[gg]
			tempboot2 = *(*par_bootb[1])[gg]
		for(b=2;b<=B;b++){
			tempboot1 = (tempboot1,*(*par_boot[b])[gg])
			tempboot2 = (tempboot2,*(*par_bootb[b])[gg])
		}
		
		tempmse = par_mse[|1,gg\dv-1,gg|]
		tempbootd = tempboot2[|1,1\dd,.|]
		tempbootb = tempboot2[|dd+1,1\.,.|]
		temphom = par_hom[|1,gg\.,gg|]
		temphomd = temphom[|1,1\dd,.|]
		temphomb = temphom[|dd+1,1\.,.|]
		
		disdvar[gg] = sum((tempboot1 - J(1,B,mean(tempboot1')')):^2)
		disbdvar[gg] = sum((tempbootd - J(1,B,mean(tempbootd')')):^2)
		disbbvar[gg] = sum((tempbootb - J(1,B,mean(tempbootb')')):^2)
		
		tau = lower + step * gg
		
		disd[gg] = disdvar[gg] * boots/N + disdbias[gg]/sqrt(tau*boots)
		disbd[gg] = disbdvar[gg] * boots/N + disbdbias[gg]/sqrt(tau*boots)
		disbb[gg] = disbbvar[gg] * boots/N + disbbbias[gg]/sqrt(tau*boots)
		
		//-----------------------------
		//Optional,tested
		/*
		real matrix disdmse, disbdmse, disbbmse
		disdmse = J(G,1,0)
		disbdmse = disdmse
		disbbmse = disdmse
		disdmse[gg] = sum((tempboot1 - J(1,B,tempmse)):^2)
		disbdmse[gg] = sum((tempbootd - J(1,B,temphomd)):^2)
		disbbmse[gg] = sum((tempbootb - J(1,B,temphomb)):^2)
		
		tau = lower + step * gg
		disd[gg] = disdmse[gg] * boots/N + sqrt(disdbias[gg]/(tau*boots))
		disbd[gg] = disbdmse[gg] * boots/N + sqrt(disbdbias[gg]/(tau*boots))
		disbb[gg] = disbbmse[gg] * boots/N + sqrt(disbbbias[gg]/(tau*boots))
		*/
		//-----------------------------
	}	
	
	
	minindex(disd,1,gstard,w)
	minindex(disbd,1,gstarbd,w)
	minindex(disbb,1,gstarbb,w)
	
	delta_mse = par_mse[|1,gstard\dv-1,gstard|]
	beta_mse = par_mse[|dv,gstard\2*(dv-1),gstard|]
	delta_hom = par_hom[|1,gstarbd\dd,gstarbd|]
	beta_hom = par_hom[|dd+1,gstarbd\dv-1,gstarbd|]
	
	//--------------------------------------------------------------------------
	//bootstrap
	real matrix temp_bootstraphomb, temp_bootstraphomd, par_bootstrap, ///
				par_bootstraphomb, par_bootstraphomd
	real scalar tempdis, tempdis_b
	par_bootstrap = J(2*(dv-1),B,0)
	par_bootstraphomb = J(dv-dd-1,B,0)
	par_bootstraphomd = J(dd,B,0)
	tempboot = (0)
	temp_bootstraphomb = (0)
	temp_bootstraphomd = (0)
	
	for(b=1;b<=B;b++){
		idz = 1::N
		tempidz = jumble(idz)
		Zz = Z[tempidz,.]
		X = Zz[|1,2\.,.|]
	    X = (J(N,1,1) , X)
		Y = Zz[|1,1\.,1|]
		
		tau = lower + step * gstard
		myfun_hetero(tau,mm,b0,d0,X,Y,l,tempboot)
	
		tau = lower + step * gstarbb
		myfun_hom(tau,mm,phi,X,Y,l,temp_bootstraphomb)
		
		tau = lower + step * gstarbd
		myfun_hom(tau,mm,phi,X,Y,l,temp_bootstraphomd)
		
		par_bootstrap[|1,b\.,b|]= tempboot
		par_bootstraphomb[|1,b\.,b|] = temp_bootstraphomb[|dd+1,1\dv-1,1|]
		par_bootstraphomd[|1,b\.,b|] = temp_bootstraphomd[|1,1\dd,1|]
				
	}
	
	//--------------------------------------------------------------------------
	real matrix cibootd95t, cibootb95t, cibootdlogt, cibootd975t, cibootb975t, ///
				cibootbb95t, cibootbb975t
				
	tau = lower + step * gstard		
	cibootd95t = (delta_mse :- 1.96*sqrt(diag((*V_mse[gstard])):/(N*tau)), delta_mse :+ 1.96*sqrt(diag((*V_mse[gstard])):/(N*tau))                        )
	cibootb95t = (beta_mse :- 1.96*sqrt(diag((*V_mse[gstard])):/nb_mse[gstard]),beta_mse :+ 1.96*sqrt(diag((*V_mse[gstard])):/nb_mse[gstard]))
	cibootdlogt = (delta_mse :- 7*sqrt(log(N))*sqrt(diag((*V_mse[gstard])):/(N*tau)),delta_mse :+ 7*sqrt(log(N))*sqrt(diag((*V_mse[gstard])):/(N*tau)))
	
	cibootbb95t = (beta_hom :- 1.96*sqrt(diag((*V_homb[gstarbb]))):/Nb_hom[gstarbb],beta_hom :+ 1.96*sqrt(diag((*V_homb[gstarbb]))):/Nb_hom[gstarbb])
	
	cibootd975t = (delta_mse :- 2.24*sqrt(diag((*V_mse[gstard])):/(N*tau)),delta_mse :+ 2.24*sqrt(diag((*V_mse[gstard])):/(N*tau)))
	cibootb975t = (beta_mse :- 2.24*sqrt(diag((*V_mse[gstard])):/nb_mse[gstard]),beta_mse :+ 2.24*sqrt(diag((*V_mse[gstard])):/nb_mse[gstard]))
	
	cibootbb975t = (beta_hom :- 2.24*sqrt(diag((*V_homb[gstarbb]))):/Nb_hom[gstarbb],beta_hom :+ 2.24*sqrt(diag((*V_homb[gstarbb]))):/Nb_hom[gstarbb])
	
	//--------------------------------------------------------------------------
	// bootstrap CI
	real matrix ciboot95p, ciboot975p,cibootbd95p,cibootbd975p, cibootbb95p, cibootbb975p
	
	ciboot95p = (mm_quantile(par_bootstraphomd',1,0.025)',mm_quantile(par_bootstraphomd',1,0.975)')
	ciboot975p = (mm_quantile(par_bootstraphomd',1,0.0125)',mm_quantile(par_bootstraphomd',1,0.9875)')
	
	cibootbd95p = (mm_quantile(par_bootstraphomd',1,0.025)',mm_quantile(par_bootstraphomd',1,0.975)')
	cibootbd975p = (mm_quantile(par_bootstraphomd',1,0.0125)',mm_quantile(par_bootstraphomd',1,0.9875)')
	
	cibootbb95p = (mm_quantile(par_bootstraphomb',1,0.025)',mm_quantile(par_bootstraphomb',1,0.975)')
	cibootbb975p = (mm_quantile(par_bootstraphomb',1,0.0125)',mm_quantile(par_bootstraphomb',1,0.9875)')
	
	//--------------------------------------------------------------------------
	// bootstrap std
	real matrix std_msedboot, cibootd95ps, cibootd975ps, std_msebboot, cibootb95ps, ///
				cibootb975ps, std_homdboot, cibootbd95ps, cibootbd975ps, ///
				std_hombboot, cibootbb95ps, cibootbb975ps
	//std_msedboot = (mm_quantile(par_bootstrap[|1,1\dv-1,.|]',1,0.75)' - mm_quantile(par_bootstrap[|1,1\dv-1,.|]',1,0.25)') :/ (invnormal(0.75) - invnormal(0.25))
	std_msedboot = sqrt(mm_colvar(par_bootstrap[|1,1\dv-1,.|]')')
	
	cibootd95ps = (delta_mse - 1.96*std_msedboot,delta_mse + 1.96*std_msedboot)
	cibootd975ps = (delta_mse - 2.24*std_msedboot,delta_mse + 2.24*std_msedboot)
	
	std_msebboot = sqrt(mm_colvar(par_bootstrap[|dv,1\.,.|]')')
	
	cibootb95ps = (beta_mse - 1.96*std_msebboot,beta_mse + 1.96*std_msebboot)
	cibootb975ps = (beta_mse - 2.24*std_msebboot,beta_mse + 2.24*std_msebboot)
	
	std_homdboot = sqrt(mm_colvar(par_bootstraphomd')')
	
	cibootbd95ps = (delta_hom - 1.96*std_homdboot,delta_hom + 1.96*std_homdboot)
	cibootbd975ps = (delta_hom - 2.24*std_homdboot,delta_hom + 2.24*std_homdboot)
	
	std_hombboot = sqrt(mm_colvar(par_bootstraphomb')')
	
	cibootbb95ps = (beta_hom - 1.96*std_hombboot,beta_hom + 1.96*std_hombboot)
	cibootbb975ps = (beta_hom - 2.24*std_hombboot,beta_hom + 2.24*std_hombboot)
	
	//--------------------------------------------------------------------------
	
	
	//what is specificationtest
	
	real matrix beta_median, beta_median_boot, ci_med95, ci_med975
	
	beta_median = rq_fnm(X,Y,0.5)
	beta_median_boot = J(dv,500,0)
	
		
	for(j=1;j<=500;j++){
		idz = 1::N
		tempidz = jumble(idz)
		Zz = Z[tempidz,.]
		X= Zz[|1,2\.,.|]
		X = (J(N,1,1) , X)
		Y = Zz[|1,1\.,1|]
		beta_median_boot[.,j] = rq_fnm(X,Y,0.5)
	}
	
	ci_med95 = (mm_quantile(beta_median_boot',1,0.025)',mm_quantile(beta_median_boot',1,0.975)')
	ci_med975 = (mm_quantile(beta_median_boot',1,0.0125)',mm_quantile(beta_median_boot',1,0.9875)')
	
	st_matrix("e(cibootd95ps)",cibootd95ps)
	st_matrix("e(cibootd975ps)",cibootd975ps)
	st_matrix("e(cibootbd95ps)",cibootbd95ps)
	st_matrix("e(cibootbd975ps)",cibootbd975ps)
	st_matrix("e(cibootbb95ps)",cibootbb95ps)
	st_matrix("e(cibootbb975ps)",cibootbb975ps)
	st_matrix("e(par_bootstrap)",par_bootstrap)
	st_matrix("e(par_bootstraphomb)",par_bootstraphomb)
	st_matrix("e(par_bootstraphomd)",par_bootstraphomd)
	
	
}

//------------------------------------------------------------------------------
//myfun_hetero
void myfun_hetero(real scalar tau, real scalar mm, real scalar b0, real scalar d0, ///
				real matrix X, real matrix Y, real matrix l, real matrix Out, ///
				| real scalar Dis, real matrix Var,real scalar Nb1,real matrix delta_MSE)
{
	//struct myres scalar v
	real scalar T, dd, JJ,i, j
	real matrix Gamma, l1, L, G, Qx
	T = rows(X)
	dd = cols(X)
	JJ = cols(l)
	
	Gamma = (-J(JJ,1,I(dd)),diag(1:/sqrt(l))#I(dd))
	l1 = (1,l)
	L = J(JJ+1,JJ+1,0)
	for(i = 1; i <= JJ+1; ++i){
		for(j = 1; j <= JJ+1; ++j){
			L[i,j] = min((l1[i],l1[j]))/sqrt(l1[i] * l1[j])
		}
	}
	G = log(l') # I(dd-1)
	Qx = X' * X / T
	
	real matrix tempy, tempx1, tempx, dtemp, pf, pm
	tempy = J((dd-1)*JJ,1,0)
	tempx1 = J(JJ,1,0) 
	pf = J(dd,JJ+1,0)
	pm = J(dd,JJ+1,0)
	pf[.,1] = rq_fnm(X,-Y,tau)
	pm[.,1] = rq_fnm(X,-Y,mm*tau)
	
	for(j=1;j<=JJ;j++){
		pf[.,j+1] = rq_fnm(X,-Y,tau*l[j])
		tempy[(dd-1)*(j-1)+1..(dd-1)*j] = (pf[|2,j+1\.,j+1|] - pf[|2,1\.,1|]) * d0
		tempx1[j] = pf[1,j+1] - pf[1,1]
		pm[.,j+1] = rq_fnm(X,-Y,mm*tau*l[j])
	}
	pf[1,.] = pf[1,.] + J(1,JJ+1,b0)
	tempx = tempx1 # I(dd-1) 
	
	if (args()<=11){
		dtemp = invsym(tempx' * tempx) * tempx' * tempy
	}
	else{
		dtemp = delta_MSE
	}
	
	real matrix Tri, Xd, tempXd, idx, QH, omega_0, W, delta, mom, beta, invQH
	real scalar Nb
	
	for(i = 1; i <=2; ++i){
		Tri = (-dtemp:/d0, I(dd-1))
		Xd = X:/(J(1,dd,abs(X*(d0, dtemp')')):^.5)
		tempXd = Xd * (d0\dtemp)
	
		idx = abs(tempXd :- mm_median(tempXd)) :< 3*(mm_quantile(tempXd,1,.75) - mm_quantile(tempXd,1,.25))
		Xd = select(Xd,idx)
		QH = Xd' * Xd/rows(idx)
		invQH = invsym(QH)
		omega_0 = invQH * Qx * invQH
	
		W = (I(JJ) # Tri) * Gamma * (L # omega_0) * Gamma' * (I(JJ) # Tri')
		W = luinv(W)
		dtemp = luinv(tempx'*W*tempx) * tempx' * W * tempy	
	}
	
	delta = dtemp
	Var = luinv(G'*W*G)
	Out = delta
	
	mom = J(JJ*(dd-1), 1, 0)
	for(j = 1; j <= JJ; ++j){
		mom[(dd-1)*(j-1)+1..(dd-1)*j] = (pf[|2,j+1\.,j+1|] - pf[|2,1\.,1|]) :- (pf[1,j+1] - pf[1,1]) * delta/d0
	}
	
	Nb = log(mm) * sqrt(tau*T)/mm_median(abs(pm[|1,1\1,.|] - pf[|1,1\1,.|])')
	beta = - pf[|2,1\.,1|] + delta/d0 * pf[1,1]
	Out = (Out\beta)
	Nb1 = abs(sqrt(tau*T)/mm_median(pf[|1,1\1,.|]'))
	Dis = Nb:^2 * mom' * W * mom
	
	/*
	v.myout = Out
	v.myvar = Var
	v.mydis = Dis
	v.myNb1 = Nb1
	*/
}

//------------------------------------------------------------------------------
//myfun_hom
void myfun_hom(real scalar tau, real scalar mm, real matrix phi, real matrix X, ///
				real matrix Y, real matrix l, real matrix par, | real scalar dis_d, ///
				real scalar dis_b, real scalar Vd, real matrix Vb, ///
				real scalar Nb, real matrix delta_MSE)
{
	//struct myres scalar v
	real scalar T, dv, JJ,i, j, dbb, dd, kount
	real matrix psi, phitemp, psitemp 
	T = rows(X)
	dv = cols(X)
	JJ = cols(l)
	
	dbb = sum(phi)
	psi = J(1,dv-1,1) - phi
	dd = sum(psi)
	phitemp = J(dbb,dv-1,0)
	psitemp = J(dd,dv-1,0)
	kount = 1
	
	for(i=1;i<=dbb;i++){
		for(j=kount;j<=dv-1;j++){
			if (phi[j] == 1) {
				phitemp[i,j] = 1
				kount = j + 1
				break
			}
		}
	}
	
	kount = 1
	for(i=1;i<=dd;i++){
		for(j=kount;j<=dv-1;j++){
			if (psi[j] == 1) {
				psitemp[i,j] = 1
				kount = j + 1
				break
			}
		}
	}
	
	real matrix l1, L, Gamma, Gamma2, Gamma3
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
	
	real matrix tempy1, tempy2, Gd, Gb, Qx, pf, pm, tempx1, tempx, dtemp, dtemp2, /// 
				Tri, Xd, tempXd, idx, QH, omega_0, W1, W2, beta, mom1, mom2, GpG, ///
				invQH, psiTri, delta
	if (dd>0){
	
	tempy1 = J(dd*JJ,1,0)
	tempy2 = J(dbb*(JJ+1),1,0)
	Gd = log(l') # I(dd)
	Gb = - J(JJ+1,1,1) # I(dbb)
	Qx = X'*X/T
	pf = J(dv,JJ+1,0)
	pm = J(dv,JJ+1,0)
	pf[.,1] = rq_fnm(X,-Y,tau)
	pm[.,1] = rq_fnm(X,-Y,mm*tau)
	tempy2[|1,1\dbb,1|] = phitemp * pf[|2,1\.,1|]
	tempx1 = J(JJ,1,0)
	
	for(j=1;j<=JJ;j++){
		pf[.,j+1] = rq_fnm(X,-Y,tau*l[j])
		tempy1[dd*(j-1)+1..dd*j] = psitemp * (pf[|2,j+1\.,j+1|] - pf[|2,1\.,1|])
		tempy2[dbb*j+1..dbb*(j+1)] = phitemp * pf[|2,j+1\.,j+1|]
		tempx1[j] = pf[1,j+1] - pf[1,1]
		pm[.,j+1] = rq_fnm(X,-Y,tau*mm*l[j])
	}
	
	tempx = tempx1 # I(dd)	
	if (args() <= 12){
		dtemp = invsym(tempx' * tempx) * tempx' * tempy1
	}
	else{
		dtemp = delta_MSE	
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
	for(j=1; j<=JJ; j++){
		mom1[dd*(j-1)+1..dd*j] = psitemp * (pf[|2,j+1\.,j+1|] - pf[|2,1\.,1|]) :- (pf[1,j+1] - pf[1,1]) * delta
		mom2[dbb*(j-1)+1..dbb*j] = phitemp * pf[|2,j\.,j|] + beta
	}
	
	mom2 = mom2\(phitemp * pf[|2,JJ+1\.,JJ+1|] + beta)
	Nb = log(mm) * sqrt(tau*T)/mm_median(abs(pm[1,.] - pf[1,.])')
	
	dis_d = Nb^2 * mom1' * W1 * mom1
	dis_b = Nb^2 * mom2' * W2 * mom2
	par = delta
	par = par\beta	
	}
	
	else if (dd==0){
	
	Gb = - J(JJ+1,1,1) # I(dbb)
	Qx = X'*X/T
	pf = J(dv,JJ+1,0)
	pm = J(dv,JJ+1,0)
	pf[.,1] = rq_fnm(X,-Y,tau)
	pm[.,1] = rq_fnm(X,-Y,mm*tau)
	tempy2[|1,1\dbb,1|] = phitemp * pf[|2,1\.,1|]
	tempx1 = J(JJ,1,0)
	
	for(j=1;j<=JJ;j++){
		pf[.,j+1] = rq_fnm(X,-Y,tau*l[j])
		tempy2[dbb*j+1..dbb*(j+1)] = phitemp * pf[|2,j+1\.,j+1|]
		tempx1[j] = pf[1,j+1] - pf[1,1]
		pm[.,j+1] = rq_fnm(X,-Y,tau*mm*l[j])
	}
	
	GpG = Gamma3 # (phitemp * Gamma2)
	W2 = luinv(GpG * (L # invsym(Qx)) * GpG' ) 
	beta = - luinv(J(JJ+1,1,I(dbb))' * W2 * J(JJ+1,1,I(dbb))) * J(JJ+1,1,I(dbb))' * W2 * tempy2
	Vb = luinv(Gb' * W2 * Gb)
	
	mom2 = J(JJ*dbb, 1, 0)
	for(j = 1; j <= JJ; ++j){
		mom2[dbb*(j-1)+1..dbb*j,1] = phitemp * pf[|2,j\.,j|] + beta
	}
	mom2 = mom2\(phitemp * pf[|2,JJ+1\.,JJ+1|] + beta )
	Nb = log(mm) * sqrt(tau*T)/mm_median(abs(pm[1,.] - pf[1,.])')
	
	dis_d = 0
	dis_b = Nb^2 * mom2' * W2 * mom2
	
	par = beta
	Vd = 0
	}
	
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

mata: mata mosave myfun_combined(),dir(PERSONAL) replace
mata: mata mosave myfun_hetero(),dir(PERSONAL) replace
mata: mata mosave myfun_hom(),dir(PERSONAL) replace
mata: mata mosave rq_fnm(),dir(PERSONAL) replace
