*!Computing the CLR bound, Xiaoyun Qiu, 11Aug,2015

mata:
version 13
mata clear
mata set matastrict on

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

