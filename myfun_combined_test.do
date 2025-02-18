*! program myfun_combined Xiaoyun Qiu 5Aug2015


//------------------------------------------------------------------------------
//Stata program: creates a new command
capture program drop myfun_combined_test
program myfun_combined_test,eclass
	version 13
	syntax varlist(min=2 numeric) [if] [in] ,INDEPnum(real) Dnum(real) DFnum(real) Cnum(real) GRIDd(string)
	marksample touse
	tokenize `varlist'
	
	scalar quant = 0.5
	mata:myfun_combined("`varlist'","`touse'","`indep'","`dnum'","`dfnum'","`cnum'","quant","`gridd'")
	
	ereturn list

end


//------------------------------------------------------------------------------
//Mata part: main function
mata:
version 13
mata clear
mata set matastrict on

//void myfun_combined(string scalar varlist, string scalar touse)
void myfun_combined(string scalar varlist,string scalar touse, |string scalar indep, ///
					string scalar dnum, string scalar dfnum, string scalar cnum, ///
					string scalar qquant, string matrix gridd)
{
	real matrix Z, X, Y, Xdnp, Xcnp, gridnp, gridid, grid
	real scalar dv, N, xdnum, t_xdnum, xcnum, indepnum, quant
	st_view(Z=.,.,tokens(varlist),touse)
	
	indepnum = strtoreal(st_local("indep"))
	xdnum = strtoreal(st_local("dfnum"))
	t_xdnum = strtoreal(st_local("dnum"))
	xcnum = strtoreal(st_local("cnum"))
	quant = st_numscalar(qquant)
	grid = st_matrix(gridd)
	
	Y = Z[.,1]
	X = Z[|1,2\.,.|]

	N = rows(X)
	X = (J(N,1,1) , X)
	dv = cols(X)
	
	//gridid = ((grid[.,1] + grid[.,2]) :< 2)
	//grid = select(grid,gridid)
	grid = grid'
	gridnp = grid[|1,1\xdnum+xcnum,.|]
	
	rseed(13579)
	
	//--------------------------------------------------------------------------
	//defining parameters
	real scalar G, JJ, B, boots, mm, b0, d0, lower, upper, step, dd, dbb, i, j, b, gg, R
	real matrix l, ciboot1, phi, nb_mse, disd, disbd, disbb, disdbias, disbbbias, disbdbias
				
	G = 4
	l = (0.65,0.85,1.15,1.45)
	JJ = cols(l)
	B = 10
	boots = 550
	R = 200
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
	pointer(real matrix) rowvector  V_homb, V_mse
	//pointer(real matrix) rowvector  V_homd
	par_boot = J(B,1,NULL)
	par_bootb = J(B,1,NULL)
	V_homb = J(1,G,NULL)
	//V_homd = J(1,G,NULL)
	V_mse = J(1,G,NULL)
	
	for(gg=1;gg<=G;gg++){
		V_homb[gg] = &(0)
		//V_homd[gg] = &(0)
		V_mse[gg] = &(0)
	}
	
	real matrix par_mse, par_hom, Nb_hom, chid, tempV
	//real matrix  chibb
	real scalar median_d1, median_d2, median_b1, tempdis, tempdis_b, tempnb
	par_mse = J(2*(dv-1),G,0)
	par_hom = J(dv-1,G,0)
	Nb_hom = J(G,1,0)
	chid = J(G,1,0)
	//chibb = J(G,1,0)
	nb_mse = J(G,1,0)
	median_d1 = mm_median(rchi2(100000,1,(JJ-1)*(dv-1)))
	median_d2 = mm_median(rchi2(100000,1,(JJ-1)*dd))
	median_b1 = mm_median(rchi2(100000,1,JJ*dbb))
	
	real scalar nbtemp, tau, tempdis_d
	real matrix tempmsef, tempVd, tempVb, tempmse
	
	timer_clear()
	for(gg=1;gg<=G;gg++){
		timer_on(1)
		tempmsef = (0)
		tempmse = (0)
		tempV = (0)
		tempVd = (0)
		tempVb = (0)
	
		tau = lower + step * gg
		myfun_hetero(tau,mm,b0,d0, X,Y,l,tempmsef,tempdis,tempV,tempnb)
		par_mse[|1,gg\2*(dv-1),gg|] = tempmsef
		chid[gg] = tempdis
		V_mse[gg] = &tempV
		nb_mse[gg] = tempnb
		myfun_hom(tau,mm,phi,X,Y,l,tempmse,tempdis_d,tempdis_b,tempVd,tempVb,tempnb)
		par_hom[|1,gg\(dv-1),gg|] = tempmse
		Nb_hom[gg] = tempnb
		//chibb[gg] = tempdis_b
		//V_homd[gg] = &tempVd
		V_homb[gg] = &tempVb
		timer_off(1)
		
	}		
			
	//--------------------------------------------------------------------------
	//line 124
	real matrix tempchi_bootd, tempchi_bootbd, tempchi_bootbb
	tempchi_bootd = J(B,G,0)
	tempchi_bootbd = J(B,G,0)
	tempchi_bootbb = J(B,G,0)
	
	real matrix Zz, idz, tempidz, seltempidz, tempboot

	for(b=1;b<=B;b++){
		timer_on(2)
		idz = 1::N
		tempidz = jumble(idz)
		seltempidz = tempidz[|1,1\boots,1|]
		Zz = Z[seltempidz,.]
		X = Zz[|1,2\.,.|]
		X = (J(boots,1,1), X)
		Y = Zz[|1,1\.,1|]
		par_boot[b] = &(J(G,1,NULL))
		par_bootb[b] = &(J(G,1,NULL))
		
		tempmsef = (0)
		tempmse = (0)
		
		for(gg=1;gg<=G;gg++){	
			tau = lower + step * gg
			myfun_hetero(tau,mm,b0,d0,X,Y,l,tempmsef,tempdis)
			tempchi_bootd[b,gg] = tempdis
			(*par_boot[b])[gg] = &(tempmsef[|1,1\dv-1,1|])
		    myfun_hom(tau,mm,phi,X,Y,l,tempmse,tempdis_d,tempdis_b)
			(*par_bootb[b])[gg] = &(tempmse)
			tempchi_bootbd[b,gg] = tempdis_d
			tempchi_bootbb[b,gg] = tempdis_b
		}
		timer_off(2)

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
	tempboot1 = J(dv-1,B,0)
	tempboot2 = J(dv-1,B,0)
	for(gg=1;gg<=G;gg++){
		timer_on(3)
		for(b=1;b<=B;b++){
			tempboot1[|1,b\.,b|] = *(*par_boot[b])[gg]
			tempboot2[|1,b\.,b|] = *(*par_bootb[b])[gg]
		}
		
		//tempmse = par_mse[|1,gg\dv-1,gg|]
		tempbootd = tempboot2[|1,1\dd,.|]
		tempbootb = tempboot2[|dd+1,1\.,.|]
		temphom = par_hom[|1,gg\.,gg|]
		//temphomd = temphom[|1,1\dd,.|]
		//temphomb = temphom[|dd+1,1\.,.|]
		
		disdvar[gg] = mean(colsum((tempboot1 - J(1,B,mean(tempboot1')')):^2)')
		disbdvar[gg] = mean(colsum((tempbootd - J(1,B,mean(tempbootd')')):^2)')
		disbbvar[gg] = mean(colsum((tempbootb - J(1,B,mean(tempbootb')')):^2)')
		
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
		disdmse[gg] = mean(colsum((tempboot1 - J(1,B,tempmse)):^2)')
		disbdmse[gg] = mean(colsum((tempbootd - J(1,B,temphomd)):^2)')
		disbbmse[gg] = mean(colsum((tempbootb - J(1,B,temphomb)):^2)')
		
		tau = lower + step * gg
		disd[gg] = disdmse[gg] * boots/N + sqrt(disdbias[gg]/(tau*boots))
		disbd[gg] = disbdmse[gg] * boots/N + sqrt(disbdbias[gg]/(tau*boots))
		disbb[gg] = disbbmse[gg] * boots/N + sqrt(disbbbias[gg]/(tau*boots))
		*/
		//-----------------------------
		timer_off(3)
	}	
	
	minindex(disd,1,gstard,w)
	minindex(disbd,1,gstarbd,w)
	minindex(disbb,1,gstarbb,w)
	
	st_numscalar("e(gstard)",gstard)
	st_numscalar("e(gstarbd)",gstarbd)
	st_numscalar("e(gstarbb)",gstarbb)
	
	delta_mse = par_mse[|1,gstard\dv-1,gstard|]
	beta_mse = par_mse[|dv,gstard\2*(dv-1),gstard|]
	delta_hom = par_hom[|1,gstarbd\dd,gstarbd|]
	beta_hom = par_hom[|dd+1,gstarbb\dv-1,gstarbb|]
	real matrix V_mse_star
	real scalar nb_mse_star
	V_mse_star = *V_mse[gstard]
	nb_mse_star = nb_mse[gstard]
	
	st_matrix("e(delta_mse)",delta_mse)
	st_matrix("e(beta_mse)",beta_mse)
	st_matrix("e(V_mse_star)",V_mse_star)
	st_numscalar("e(nb_mse_star)",nb_mse_star)
	st_matrix("e(delta_hom)",delta_hom)
	st_matrix("e(beta_hom)",beta_hom)
	
	//--------------------------------------------------------------------------
	//bootstrap
	real matrix temp_bootstraphomb, temp_bootstraphomd, par_bootstrap, ///
				par_bootstraphomb, par_bootstraphomd
	
	par_bootstrap = J(2*(dv-1),B,0)
	par_bootstraphomb = J(dv-dd-1,B,0)
	par_bootstraphomd = J(dd,B,0)
	tempboot = (0)
	temp_bootstraphomb = (0)
	temp_bootstraphomd = (0)
	
	for(b=1;b<=B;b++){
		timer_on(4)
		//idz = 1::N
		//tempidz = jumble(idz)
		tempidz = rdiscrete(N,1,J(N,1,1/N))
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
		timer_off(4)		
	}
	
	//--------------------------------------------------------------------------
	real matrix cibootd95t, cibootb95t, cibootdlogt, cibootd975t, cibootb975t, ///
				cibootbb95t, cibootbb975t
	timer_on(5)			
	tau = lower + step * gstard		
	cibootd95t = (delta_mse :- 1.96*sqrt(diagonal((*V_mse[gstard]))/(N*tau)), delta_mse :+ 1.96*sqrt(diagonal((*V_mse[gstard]))/(N*tau)))
	cibootb95t = (beta_mse :- 1.96*sqrt(diagonal((*V_mse[gstard])):/nb_mse[gstard]),beta_mse :+ 1.96*sqrt(diagonal((*V_mse[gstard])):/nb_mse[gstard]))
	cibootdlogt = (delta_mse :- 7*sqrt(log(N))*sqrt(diagonal((*V_mse[gstard])):/(N*tau)),delta_mse :+ 7*sqrt(log(N))*sqrt(diagonal((*V_mse[gstard])):/(N*tau)))
	
	cibootbb95t = (beta_hom :- 1.96*sqrt(diagonal((*V_homb[gstarbb]))):/Nb_hom[gstarbb],beta_hom :+ 1.96*sqrt(diagonal((*V_homb[gstarbb]))):/Nb_hom[gstarbb])
	
	cibootd975t = (delta_mse :- 2.24*sqrt(diagonal((*V_mse[gstard])):/(N*tau)),delta_mse :+ 2.24*sqrt(diagonal((*V_mse[gstard])):/(N*tau)))
	cibootb975t = (beta_mse :- 2.24*sqrt(diagonal((*V_mse[gstard])):/nb_mse[gstard]),beta_mse :+ 2.24*sqrt(diagonal((*V_mse[gstard])):/nb_mse[gstard]))
	
	cibootbb975t = (beta_hom :- 2.24*sqrt(diagonal((*V_homb[gstarbb]))):/Nb_hom[gstarbb],beta_hom :+ 2.24*sqrt(diagonal((*V_homb[gstarbb]))):/Nb_hom[gstarbb])
	
	//--------------------------------------------------------------------------
	// bootstrap CI
	real matrix ciboot95p, ciboot975p,cibootbd95p,cibootbd975p, cibootbb95p, cibootbb975p
	
	ciboot95p = (mm_quantile(par_bootstrap',1,0.025)',mm_quantile(par_bootstrap',1,0.975)')
	ciboot975p = (mm_quantile(par_bootstrap',1,0.0125)',mm_quantile(par_bootstrap',1,0.9875)')
	
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
	timer_off(5)
	st_matrix("e(cibootd95t)",cibootd95t)
	st_matrix("e(cibootd975t)",cibootd975t)
	st_matrix("e(cibootb95t)",cibootb95t)
	st_matrix("e(cibootb975t)",cibootb975t)
	st_matrix("e(cibootbb95t)",cibootbb95t)
	st_matrix("e(cibootbb975t)",cibootbb975t)
	st_matrix("e(cibootdlogt)",cibootdlogt)
	
	st_matrix("e(ciboot95p)",ciboot95p)
	st_matrix("e(ciboot975p)",ciboot975p)
	st_matrix("e(cibootbd95p)",cibootbd95p)
	st_matrix("e(cibootbd975p)",cibootbd975p)
	st_matrix("e(cibootbb95p)",cibootbb95p)
	st_matrix("e(cibootbb975p)",cibootbb975p)
	
	st_matrix("e(cibootd95ps)",cibootd95ps)
	st_matrix("e(cibootd975ps)",cibootd975ps)
	st_matrix("e(cibootb95ps)",cibootb95ps)
	st_matrix("e(cibootb975ps)",cibootb975ps)
	st_matrix("e(cibootbd95ps)",cibootbd95ps)
	st_matrix("e(cibootbd975ps)",cibootbd975ps)
	st_matrix("e(cibootbb95ps)",cibootbb95ps)
	st_matrix("e(cibootbb975ps)",cibootbb975ps)
	
	//--------------------------------------------------------------------------
	
	real matrix beta_median, beta_median_boot, ci_med95, ci_med975
	
	beta_median = rq_fnm(X,Y,0.5)
	beta_median_boot = J(dv,500,0)
	
	for(j=1;j<=500;j++){
		timer_on(6)
		//idz = 1::N
		//tempidz = jumble(idz)
		tempidz = rdiscrete(N,1,J(N,1,1/N))
		Zz = Z[tempidz,.]
		X= Zz[|1,2\.,.|]
		X = (J(N,1,1) , X)
		Y = Zz[|1,1\.,1|]
		beta_median_boot[.,j] = rq_fnm(X,Y,0.5)
		timer_off(6)
	}
	
	ci_med95 = (mm_quantile(beta_median_boot',1,0.025)',mm_quantile(beta_median_boot',1,0.975)')
	ci_med975 = (mm_quantile(beta_median_boot',1,0.0125)',mm_quantile(beta_median_boot',1,0.9875)')
	
	st_matrix("e(ci_med95)",ci_med95)
	st_matrix("e(ci_med975)",ci_med975)
	
	
	//--------------------------------------------------------------------------
	//CLR bound for qth-QTE of the heteroskedastic variables
	
	real scalar M, r, rn, knv, knv2, temptheta1, temptheta2
	real matrix nor, y1s, y2s, row, y1, y2, Xgrid, Zu_ss, Zl_ss, thetau, thetal, ///
				Su, Sl, ll, theta0, e1, phiid, temp, ///
				temps, val1, val2, val3, val4, tempzss, row2, Zu_ss2, Zussmax, ///
				Zuss2max, row2id, Zlssmax, var1, var2, Zl_ss2, Zlss2max
	if(dbb<dv-1){
		
		M = cols(grid)
		nor = rnormal(dv-1,R,0,1)
		
		y1s = (0)
		y2s = (0)
		//quant is an input
		
		Y = Z[.,1]
		Xdnp = Z[|1,2\.,1+xdnum|]
		Xcnp = Z[|1,2+t_xdnum\.,1+t_xdnum+xcnum|]
		clr_bound(Y,Xdnp,Xcnp,quant,gridnp,y1s,y2s)
		st_matrix("e(y1s)",y1s)
		st_matrix("e(y2s)",y2s)
		
		real matrix rowid
		row = (y1s:!=-1000):*(y2s:!=1000):*(y1s:!=100):*(1..M)'
		rowid = (row:!=0)
		row = select(row,rowid)
		y1 = y1s[row]
		y2 = y2s[row]
		Xgrid = grid
		Xgrid = Xgrid[.,row]
		
		M = rows(row)
		Zu_ss = J(M,R,0)
		Zl_ss = J(M,R,0)
		thetau = J(M,1,0)
		thetal = J(M,1,0)
		Su = J(M,1,0)
		Sl = J(M,1,0)
		phiid = (phi :== 0) 
		ll = selectindex(phiid)
		theta0 = J(dd,2,0)
		
		real scalar tempscal
		for(j=1;j<=dd;j++){
			timer_on(7)
			e1 = J(dv-1,1,0)
			tempscal = ll[j]
			e1[tempscal] = 1
			
			for(mm=1;mm<=M;mm++){

				val1 = Xgrid[.,mm]'*delta_mse
				val2 = Xgrid[.,mm]'*beta_mse
				
				if(e1'*delta_mse>0 & 1+val1>0.05){ 
					temptheta1 = e1' * (beta_mse + delta_mse * (y2[mm] - val2 - b0)/(1 + val1))
					temptheta2 = -e1' * (beta_mse + delta_mse * (y1[mm] - val2 - b0)/(1 + val1))
				}
				else if(e1'*delta_mse<0 & 1+val1>0.05){
					temptheta1 = e1' * ( beta_mse  + delta_mse *(y1[mm] - val2 - b0)/(1 + val1))
					temptheta2 = -e1' * (beta_mse + delta_mse * (y2[mm] - val2 - b0)/(1 + val1))
				}
				else{
					temptheta1 = 1000
					temptheta2 = 1000
				}
				
				val3 = e1'*delta_mse
				temp = sign(val3) * (e1 - val3 * Xgrid[.,mm])/abs(1 + val1) 
				//real matrix lrec,squareroot
				//lrec = cholesky(V_mse_star)
				//squareroot = lrec * 2 - diag(lrec)
				val4 = temp'*matpowersym(V_mse_star,0.5)
				temps = sqrt(sum(val4:^2))/abs(nb_mse_star)
				
				for(r=1;r<=R;r++){
					tempzss = val4 * nor[.,r]/sqrt(sum(val4:^2))
					Zu_ss[mm,r] = tempzss
					Zl_ss[mm,r] = -tempzss
				}
				thetau[mm,1] = temptheta1
				Su[mm,1] = temps
				thetal[mm,1] = temptheta2
				Sl[mm,1] = temps
			}
			rn = 1 - 0.1/log(N)
			Zussmax = colmax(Zu_ss)
			knv = mm_quantile(Zussmax',1,rn)
			
			row2 = (thetau < (min(thetau + knv*Su) :+2*knv*Su))
			row2 = row2 :* (1 .. rows(row2))'
			row2id = (row2:!=0)
			row2 = select(row2,row2id)
			Zu_ss2 = Zu_ss[row2,.]
			Zuss2max = colmax(Zu_ss2)
			knv2 = mm_quantile(Zuss2max',1,0.975)			
			theta0[j,2] = min(thetau + knv2*Su)
			
			Zlssmax = colmax(Zl_ss)
			knv = mm_quantile(Zlssmax',1,rn)
			
			row2 = (thetal < min(thetal + knv*Sl):+2*knv*Sl)
			row2 = row2 :* (1 .. rows(row2))'
			row2id = (row2:!=0)
			row2 = select(row2,row2id)
			Zl_ss2 = Zl_ss[row2,.]
			Zlss2max = colmax(Zl_ss2)
			knv2 = mm_quantile(Zlss2max',1,0.975)
			theta0[j,1] = -min(thetal + knv2*Sl)
			timer_off(7)
		}
		
		real scalar qstard, qstarbd, qstarbb
		
		qstard = lower + step * gstard
		qstarbd = lower + step * gstarbd
		qstarbb = lower + step * gstarbb
		
		st_matrix("e(theta0)",theta0)
	}
	timer()
	
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
//rq_fnm.m

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
//---------------------------------------
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
//---------------------------------------
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

//------------------------------------------------------------------------------
//bound.m
void clr_bound(Y,Xdnp,Xcnp,tau,gridnp,y1,y2)
{
	real matrix D,tempx
	real scalar hh, nn, ll, dd, i, mm, quant
	D = (Y:>0)
	nn = rows(Y)
	ll = cols(Y)
	mm = cols(gridnp)
	dd = cols(Xdnp)
	hh = 1.06 * nn^(-1/(4 + cols(Xcnp)))
	y1 = J(mm,1,0)
	y2 = J(mm,1,0)
	quant = mm_quantile(Y,1,tau)
	
	rseed(13579)
	
	real matrix val1,val2,val3, xc, xd, Xc, Xd
	real scalar P1, temp1,temp2, rc
	for(i=1;i<=mm;i++){
		tempx = gridnp[|1,i\.,i|]
		val1 = (rowsum(Xdnp:==J(nn,1,tempx[1..dd]')):==dd)
		val2 = normalden(J(nn,1,tempx[|dd+1,1\.,1|]'):-Xcnp):/J(nn,1,sqrt(mm_colvar(Xcnp)))
		val2 = exp(rowsum(log(val2)))
		val3 = val1:*val2
		P1 = sum((D:==1):*val3)/sum(val3)
		val = D:*val3
			
		if((tau - 1 + P1)/P1>0){
			rc = mm_root(temp1=.,&fun1(),0,9,4.5,1000,Y,val,P1,tau)
		}
		else{
			temp1 = -1000
		}
		if(tau/P1<1){
			rc = mm_root(temp2=.,&fun2(),0,9,0,1000,Y,val,P1,tau)	
		}
		else{
			temp2 = 1000
		}
		y1[i,1] = temp1
		y2[i,1] = temp2
	}
}

function fun1(y,Y,val,P1,tau) return(sum((Y:<y):*val)/sum(val) - (tau-1+P1)/P1)
function fun2(y,Y,val,P1,tau) return(sum((Y:<y):*val)/sum(val) - tau/P1)
end

mata: mata mosave myfun_combined(),dir(PERSONAL) replace
mata: mata mosave myfun_hetero(),dir(PERSONAL) replace
mata: mata mosave myfun_hom(),dir(PERSONAL) replace
mata: mata mosave rq_fnm(),dir(PERSONAL) replace
mata: mata mosave clr_bound(),dir(PERSONAL) replace
