version 12.0

mata:
real matrix myfun_MSE(real scalar tau,real scalar m,real scalar b0, real scalar d0, ///
						real matrix X, real matrix Y, real rowvector l,real rowvector delta_MSE)
{
	real scalar  T, d, J
	real matrix  pf, gamma, L, G, Qx, tempy, tempx1
	real vector  l1
	
	T = cols(X)
	d = rows(X)
	J = cols(l)
	pf = J(d,J+1,0)
	
	gamma = J(J,1,1)#(-I(d)), diag(1:/sqrt(l))#I(d)
	l1 = 1,l
	L = J(J+1,J+1,0)
	
	for (i =1; i <= J+1; i++){
		for(j=1; j<= J+1; j++){
			L[i,j] = min((l1[i],l1[j]))/sqrt(l1[i]*l1[j])
		}
	}
	
	G = log(l')#I(d-1)
	Qx = X'*X/T
	
	tempy = J((d-1)*J,1,0)
	tempx1 = J(J,1)
	
	
}
end






