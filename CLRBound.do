*! CLRBound, Xiaoyun Qiu, 9Sep2015

//------------------------------------------------------------------------------
//Stata program: creates a new command
capture program drop CLRBound
program CLRBound,eclass
	version 13
	syntax varlist(min=2 numeric) [if] [in] ,INDEPnum(real) Dnum(real) DFnum(real) Cnum(real) GRIDd(string)
	marksample touse
	tokenize `varlist'
	
	scalar quant = 0.5
	mata:CLRBound("`varlist'","`touse'","`indep'","`dnum'","`dfnum'","`cnum'","quant","`gridd'")
	
	ereturn list

end

//------------------------------------------------------------------------------
//Mata part: main function
mata:
version 13
mata clear
mata set matastrict on

//void myfun_combined(string scalar varlist, string scalar touse)
void CLRBound(string scalar varlist,string scalar touse, |string scalar indep, ///
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
	
	real matrix phi
	real scalar dbb,dd,R,j,mm,b0,d0,lower,upper,boots,step,G
	
	G = 40
	boots = 550
	R = 200
	mm = 1.2
	b0 = 0
	d0 = 1
	lower = min((80/boots, 0.1))
	//lower = 0.05
	upper = 0.3
	step = (upper - lower)/G
	phi = (1,1,0,1,1)
	dbb = sum(phi)
	dd = dv - 1 - dbb
	
	rseed(13579)
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

		st_matrix("e(theta0)",theta0)
	}
	timer()
	
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
		val1 = (rowsum(Xdnp:==J(nn,1,tempx[1..dd]')):==3)
		val2 = normalden(J(nn,1,tempx[|dd+1,1\.,1|]'):-Xcnp):/J(nn,1,sqrt(mm_colvar(Xcnp)))
		val2 = exp(rowsum(log(val2)))
		val3 = val1:*val2
		P1 = sum((D:==1):*val3)/sum(val3)
		
		if((tau - 1 + P1)/P1>0){
			xc = tempx[|dd+1,1\.,1|]'
			xd = tempx[|1,1\dd,1|]'
			Xc = Xcnp
			Xd = Xdnp

			real scalar Xdrow, Xdcol
			real matrix xdtemp, xctemp, val
			Xdrow = rows(Xd)

			xdtemp = J(Xdrow,1,xd)
			xctemp = J(Xdrow,1,xc)
			val1 = (rowsum(Xd:==xdtemp):==3)
			val2 = exp(rowsum(log(normalden(xctemp-Xc/hh))))
			val3 = (D:==1)
			val = val3:*val1:*val2
			
			rc = mm_root(temp1=.,&fun1(),0,9,1e-5,1000,Y,val,P1,tau)

		}
		else{
			temp1 = -1000
		}
		if(tau/P1<1){
			xc = tempx[|dd+1,1\.,1|]'
			xd = tempx[|1,1\dd,1|]'
			Xc = Xcnp
			Xd = Xdnp
			Xdrow = rows(Xd)
	
			xdtemp = J(Xdrow,1,xd)
			xctemp = J(Xdrow,1,xc)
			val1 = (rowsum(Xd:==xdtemp):==3)
			val2 = exp(rowsum(log(normalden(xctemp-Xc/hh))))
			val3 = (D:==1)
			val = val3:*val1:*val2
			
			rc = mm_root(temp2=.,&fun2(),0,9,1e-5,1000,Y,val,P1,tau)
			
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

mata: mata mosave CLRBound(),dir(PERSONAL) replace
