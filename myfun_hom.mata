*! Mata version of myfun_hom, Xiaoyun Qiu,8Aug2015

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

mata: mata mosave myfun_hom(),dir(PERSONAL) replace
