*!Computing the CLR bound, Xiaoyun Qiu, 11Aug,2015
//Test what is wrong with the bound!

mata:
version 13
mata clear
mata set matastrict on


//CLR bound for qth-QTE of the heteroskedastic variables
void bound(Z,grid)	
{	real scalar M, r, rn, knv, knv2, temptheta1, temptheta2,dbb,dv,R,nb_mse_star, ///
				xdnum,xcnum,t_xdnum,quant, dd, j, mm, b0, N, gstard, gstardd, ///
				gstarbb,yl,yr
	real matrix nor, y1s, y2s, row, y1, y2, Xgrid, Zu_ss, Zl_ss, thetau, thetal, ///
				Su, Sl, ll, theta0, e1, phiid, temp, beta_mse, V_mse_star,delta_mse, ///
				temps, val1, val2, val3, val4, tempzss, row2, Zu_ss2, Zussmax, ///
				Zuss2max, row2id, Zlssmax, var1, var2, Zl_ss2, Zlss2max,Y,Xdnp,Xcnp, ///
				phi, gridnp, X,knv1, knv0, val44,up, lo

		beta_mse = ( -.21540589, .01351065, .21492722, .23769417, .03031442)'
		delta_mse = (.01938248,.00521279, -.02902246, -.0051029, -.00664483)'
		nb_mse_star = 3.989198859229658
		V_mse_star = ((.42270082, .16876999, .00591809, .06565046,-.01395695)\(.16876999, .45157217, .00027025, .04523325, .00157416)\( .00591809, .00027025,.00415519,-.00379549,-.00534896)\( .06565046, .04523325,-.00379549, .07046129,-.00751604)\( -.01395695,.00157416,-.00534896,-.00751604,.05711054))
		R = 500
		xdnum=3
		t_xdnum=3
		xcnum=1
		quant=0.5
		mm = 1.2
		b0 = 0
		grid = grid'
		gridnp = grid[|1,1\xdnum+xcnum,.|]
		
		Y = Z[.,1]
		X = Z[|1,2\.,.|]

		N = rows(X)
		X = (J(N,1,1) , X)
		dv = cols(X)
	
		phi = (1,1,0,1,1)
		dbb = sum(phi)
		dd = dv - 1 - dbb
		
		M = cols(grid)
		nor = rnormal(dv-1,R,0,1)
		
		yl = 0
		yr = 9
		y1s = (0)
		y2s = (0)
		//quant is an input
		
		Xdnp = Z[|1,2\.,1+xdnum|]
		Xcnp = Z[|1,2+t_xdnum\.,1+t_xdnum+xcnum|]
		clr_bound(Y,Xdnp,Xcnp,quant,gridnp,y1s,y2s,yl,yr)
		
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
		up = J(M,1,0)
		lo = J(M,1,0)
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
				val3 = e1'*delta_mse
				
				if(val3>0 & 1+val1>0.05){
					temptheta1 = e1' * (beta_mse + delta_mse * (y2[mm] - val2 - b0)/(1 + val1))
					temptheta2 = -e1' * (beta_mse + delta_mse * (y1[mm] - val2 - b0)/(1 + val1))
					up[mm] = 1
				}
				else if(val3<0 & 1+val1>0.05){
					temptheta1 = e1' * ( beta_mse  + delta_mse *(y1[mm] - val2 - b0)/(1 + val1))
					temptheta2 = -e1' * (beta_mse + delta_mse * (y2[mm] - val2 - b0)/(1 + val1))
					up[mm] = 2
				}
				else{
					temptheta1 = 1000
					temptheta2 = 1000
					up[mm] = 3
				}
				
				
				temp = sign(val3) * (e1 - val3 * Xgrid[.,mm]/abs(1 + val1)) 
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
				
				st_matrix("e(up)",up)
				st_matrix("e(e1)",e1)
				st_matrix("e(nor)",nor)
				st_matrix("e(Zl_ss)",Zl_ss)
				st_matrix("e(Zu_ss)",Zu_ss)
				thetau[mm,1] = temptheta1
				Su[mm,1] = temps
				thetal[mm,1] = temptheta2
				Sl[mm,1] = temps
			}
			rn = 1 - 0.1/log(N)
			Zussmax = colmax(Zu_ss)
			knv0 = mm_quantile(Zussmax',1,rn)
			st_numscalar("e(knv0)",knv0)
			
			row2 = (thetau < (min(thetau + knv0*Su) :+2*knv0*Su))
			row2 = row2 :* (1 .. rows(row2))'
			row2id = (row2:!=0)
			row2 = select(row2,row2id)
			Zu_ss2 = Zu_ss[row2,.]
			Zuss2max = colmax(Zu_ss2)
			knv1 = mm_quantile(Zuss2max',1,0.975)
			st_numscalar("e(knv1)",knv1)
			
			theta0[j,2] = min(thetau + knv1*Su)
			
			Zlssmax = colmax(Zl_ss)
			knv = mm_quantile(Zlssmax',1,rn)
			st_numscalar("e(knv)",knv)
			
			row2 = (thetal < (min(thetal + knv*Sl):+2*knv*Sl))
			row2 = row2 :* (1 .. rows(row2))'
			row2id = (row2:!=0)
			row2 = select(row2,row2id)
			Zl_ss2 = Zl_ss[row2,.]
			Zlss2max = colmax(Zl_ss2)
			
			knv2 = mm_quantile(Zlss2max',1,0.975)
			st_numscalar("e(knv2)",knv2)
			
			theta0[j,1] = -min(thetal + knv2*Sl)
			timer_off(7)
		}		
		
		st_matrix("e(Xgrid)",Xgrid)
		st_matrix("e(Zussmax)",Zussmax)
		st_matrix("e(Zlssmax)",Zlssmax)
		st_matrix("e(Zuss2max)",Zuss2max)
		st_matrix("e(Zlss2max)",Zlss2max)
		st_matrix("e(theta0)",theta0)
		st_matrix("e(thetau)",thetau)
		st_matrix("e(thetal)",thetal)
		st_matrix("e(Su)",Su)
		st_matrix("e(Sl)",Sl)
	timer()
	}
//------------------------------------------------------------------------------
//bound.m
void clr_bound(Y,Xdnp,Xcnp,tau,gridnp,y1,y2,yl,yr)
{
	real matrix D,tempx
	real scalar hh, nn, ll, dd, i, mm
	D = (Y:>0)
	nn = rows(Y)
	ll = cols(Y)
	mm = cols(gridnp)
	dd = cols(Xdnp)
	hh = 1.06 * nn^(-1/(4 + cols(Xcnp)))
	y1 = J(mm,1,0)
	y2 = J(mm,1,0)
	
	real matrix val1,val2,val3, xc, xd, Xc, Xd,val
	real scalar P1, temp1,temp2, rc
	for(i=1;i<=mm;i++){
		/*tempx = gridnp[|1,i\.,i|]
		val1 = (rowsum(Xdnp:==J(nn,1,tempx[1..dd]')):==3)
		rseed(13579)
		val2 = normal(J(nn,1,tempx[dd+1...]'):-Xcnp):/J(nn,1,sqrt(variance(Xcnp)))
		val2 = exp(rowsum(log(val2)))
		val3 = val1:*val2
		P1 = sum((D:==1):*val3)/sum(val3)
		
		if((tau - 1 + P1)/P1>0){
			xc = tempx[dd+1...]'
			xd = tempx[1..dd]'
			Xc = Xcnp
			Xd = Xdnp

			real scalar Xdrow, Xdcol
			real matrix xdtemp, xctemp, val
			Xdrow = rows(Xd)
			Xdcol = cols(Xd)
			xdtemp = J(Xdrow,1,xd)
			xctemp = J(Xdrow,1,xc)
			val1 = (rowsum(Xd:==xdtemp):==3)
			rseed(13579)
			val2 = exp(rowsum(log(normal(xctemp:-Xc/hh))))
			val3 = (D:==1)
			val = val3:*val1:*val2
			val = val/sum(val)
			rc = mm_root(temp1=.,&condQ1(),0,9,mm_quantile(Y,1,tau),1000,Y,val,P1,tau)

		}
		else{
			temp1 = -1000
		}
		if(tau/P1<1){
			xc = tempx[dd+1...]'
			xd = tempx[1..dd]'
			Xc = Xcnp
			Xd = Xdnp
			Xdrow = rows(Xd)
			Xdcol = cols(Xd)
			xdtemp = J(Xdrow,1,xd)
			xctemp = J(Xdrow,1,xc)
			val1 = (rowsum(Xd:==xdtemp):==3)
			rseed(13579)
			val2 = exp(rowsum(log(normal(xctemp:-Xc/hh))))
			val3 = (D:==1)
			val = val3:*val1:*val2
			val = val/sum(val)
			rc = mm_root(temp2=.,&condQ2(),0,9,mm_quantile(Y,1,tau),1000,Y,val,P1,tau)
			
		}
		else{
			temp2 = 1000
		}
		y1[i,1] = temp1
		y2[i,1] = temp2*/
		tempx = gridnp[|1,i\.,i|]
		val1 = (rowsum(Xdnp:==J(nn,1,tempx[1..dd]')):==dd)
		val2 = normalden(J(nn,1,tempx[|dd+1,1\.,1|]'):-Xcnp):/J(nn,1,sqrt(mm_colvar(Xcnp)))
		val2 = exp(rowsum(log(val2)))
		val3 = val1:*val2
		P1 = sum((D:==1):*val3)/sum(val3)
		val = (D:==1):*val3
			
		if((tau - 1 + P1)/P1>0){
			rc = mm_root(temp1=.,&fun1(),yl,yr,1e-5,1000,Y,val,P1,tau)
		}
		else{
			temp1 = -1000
		}
		if(tau/P1<1){
			rc = mm_root(temp2=.,&fun2(),yl,yr,1e-5,1000,Y,val,P1,tau)	
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

mata: mata mosave clr_bound(),dir(PERSONAL) replace

