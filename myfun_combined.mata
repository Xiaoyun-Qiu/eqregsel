*! program myfun_combined Xiaoyun Qiu 5Aug2015

/*
//------------------------------------------------------------------------------
//Stata program: creates a new command
capture program drop myfun_combined
program myfun_combined
	version 13
	syntax varlist[min=2 numeric] [if] [in] [, ]

end
*/

//------------------------------------------------------------------------------
//Mata part: main function
mata:
version 13
mata clear
mata set matastrict on

//void myfun_combined(string scalar varlist, string scalar touse)
void myfun_combined(string scalar varlist)
{
	real matrix Z, X, Y
	real scalar dv, N
	//st_view(Z=.,.,tokens(varlist),touse)
	st_view(Z=.,.,tokens(varlist))
	X = Z[|1,2\.,.|]
	N = rows(X)
	X = (J(N,1,1) , X)
	Y = Z[.,1]
	dv = cols(X)
	
	
	real scalar G, JJ, B, boots, mm, b0, d0, lower, upper, step, dd, dbb, i, j, b, gg
	real matrix l, ciboot1, phi, nb_mse, disd, disbd, disbb, disdbias, disbbbias, disbdbias
				
	G = 40
	l = (0.65,0.85,1.15,1.45)
	JJ = cols(l)
	B = 500
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
	
	pointer(real matrix) rowvector par_boot, par_bootb, V_homb, V_homd, V_mse
	par_boot = J(1,B,NULL)
	par_bootb = J(1,B,NULL)
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
	for(gg=1;gg<=G;gg++){
		tau = lower + step * gg
		myfun_hetero(tau,mm,b0,d0, X,Y,l,par_mse[.,gg],*V_mse[gg],chid[gg],nb_mse[gg])
		myfun_hom(tau,mm,phi,X,Y,l,par_hom[.,gg],*V_homd[gg],*V_homb[gg],Nb_hom[gg],tempdis_d,chibb[gg])
	}
	
	real matrix tempchi_bootd, tempchi_bootbd, tempchi_bootbb
	tempchi_bootd = J(B,G,0)
	tempchi_bootbd = J(B,G,0)
	tempchi_bootbb = J(B,G,0)
	
	real matrix Zz, idz, tempidz, seltempidz, tempboot, tempvmse, ///
				 tempvhomb
	real scalar  tempnbmse, tempvhomd,tempnbhom
	for(b=1;b<=B;b++){
		idz = 1::N
		tempidz = jumble(idz)
		seltempidz = tempidz[|1,1\boots,1|]
		Zz = Z[seltempidz,.]
		tempboot = J(3*(dv-1),G,0)
		for(gg=1;gg<=G;gg++){
			tau = lower + step * gg
			myfun_hetero(tau,mm,b0,d0,Zz[|1,2\.,.|],Zz[|1,1\.,1|],l,tempboot[|1,gg\2*(dv-1),gg|],tempvmse,tempchi_bootd[b,gg],tempnbmse)
		    myfun_hom(tau,mm,phi,Zz[|1,2\.,.|],Zz[|1,1\.,1|],l,tempboot[|2*(dv-1)+1,gg\.,gg|],tempvhomd,tempvhomb,tempnbhom,tempchi_bootbd[b,gg],tempchi_bootbb[b,gg])
		}
	par_boot[b] = &tempboot[|1,1\dv-1,.|]
	par_bootb[b] = &tempboot[|2*(dv-1)+1,1\.,.|]
	}
	
	//-----------------------------
	//the sizes of the matrices seem to be wrong here
	disdbias[|1,1\.,1|] = abs(mm_median(tempchi_bootd) :- median_d1)
	disbdbias[|1,1\.,1|] = abs(mm_median(tempchi_bootbd) :- median_d2)
	disbbbias[|1,1\.,1|] = abs(mm_median(tempchi_bootbb) :- median_b1)
	//-----------------------------
	
	real matrix disdvar, disbdvar, disbbvar
	disdvar = J(G,1,0)
	disbdvar = disdvar
	disbbvar = disdvar
	
	real matrix tempboot1, tempmse, tempboot2, tempbootd, tempbootb, temphom, ///
				temphomd, temphomb, w, beta_mse, ///
				delta_mse, delta_hom, beta_hom
	real scalar gstard, gstarbd, gstarbb
	w = (0,0)
	for(gg=1;gg<=G;gg++){
		tempboot1 = *par_boot[gg]
		tempmse = par_mse[|1,gg\dv-1,gg|]
		tempboot2 = *par_bootb[gg]
		tempbootd = tempboot2[|1,1\dd,.|]
		tempbootb = tempboot2[|dd+1,1\.,.|]
		temphom = par_hom[|1,gg\.,gg|]
		temphomd = temphom[|1,1\dd,.|]
		temphomb = temphom[|dd+1,1\.,.|]
		
		//-----------------------------
		//Here, kind of weird, too!!
		disdvar[gg] = sum((tempboot1 - J(1,B,mean(tempboot1')')):^2)
		disbdvar[gg] = sum((tempbootd - J(1,B,mean(tempbootd')')):^2)
		disbbvar[gg] = sum((tempbootb - J(1,B,mean(tempbootb')')):^2)
		//-----------------------------
		
		/*
		disdmse[gg] = sum((tempboot1 - J(1,B,tempmse)):^2)
		disbdmse[gg] = sum((tempbootd - J(1,B,temphomd)):^2)
		disbbmse[gg] = sum((tempbootb - J(1,B,temphomb)):^2)
		*/
		
		tau = lower + step * gg
		disd[gg] = disdvar[gg] * boots/N + disdbias[gg]/sqrt(tau*boots)
		disbd[gg] = disbdvar[gg] * boots/N + disbdbias[gg]/sqrt(tau*boots)
		disbb[gg] = disbbvar[gg] * boots/N + disbbbias[gg]/sqrt(tau*boots)
		
		/*
		disd[gg] = disdmse[gg] * boots/N + sqrt(disdbias[gg]/(tau*boots))
		disbd[gg] = disbdmse[gg] * boots/N + sqrt(disbdbias[gg]/(tau*boots))
		disbb[gg] = disbbmse[gg] * boots/N + sqrt(disbbbias[gg]/(tau*boots))
		*/
		
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
	par_bootstraphomb = J(2*(dv-1),B,0)
	par_bootstraphomd = J(2*(dv-1),B,0)
	tempboot = J(2*(dv-1),1,0)
	temp_bootstraphomb = (0)
	temp_bootstraphomd = (0)
	for(b=1;b<=B;b++){
		idz = 1::N
		tempidz = jumble(idz)
		Zz = Z[tempidz,.]
			tau = lower + step * gstard
			myfun_hetero(tau,mm,b0,d0,Zz[|1,2\.,.|],Zz[|1,1\.,1|],l,tempboot,tempvmse,tempdis,tempnbmse)
			tau = lower + step * gstarbb
			myfun_hom(tau,mm,phi,Zz[|1,2\.,.|],Zz[|1,1\.,1|],l,temp_bootstraphomb,tempvhomd,tempvhomb,tempnbhom,tempdis_d,tempdis_b)
			tau = lower + step * gstarbd
			myfun_hom(tau,mm,phi,Zz[|1,2\.,.|],Zz[|1,1\.,1|],l,temp_bootstraphomd,tempvhomd,tempvhomb,tempnbhom,tempdis_d,tempdis_b)
			
			//Some concern about this line!
			par_bootstrap[|1,b\.,b|] = tempboot
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
		beta_median_boot[|1,j\.,j|] = rq_fnm(Zz[|1,2\.,.|],Zz[|1,1\.,1|],0.5)
	}
	
	ci_med95 = (mm_quantile(beta_median_boot',1,0.025)',mm_quantile(beta_median_boot',1,0.975)')
	ci_med975 = (mm_quantile(beta_median_boot',1,0.0125)',mm_quantile(beta_median_boot',1,0.9875)')
	
		
	
}
end

mata: mata mosave myfun_combined(),dir(PERSONAL) replace
