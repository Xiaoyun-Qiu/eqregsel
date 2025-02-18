*! Pure Mata function myfun_hetero.mata, Xiaoyun Qiu, 8Aug2015

//------------------------------------------------------------------------------
//Mata version of myfun_hetero.m
//------------------------------------------------------------------------------
mata:
version 13
mata clear
mata set matastrict on

/*
struct myres {
	real matrix myout
	real matrix myvar
	real scalar mydis
	real scalar myNb1
}
*/

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

mata: mata mosave myfun_hetero(),dir(PERSONAL) replace
