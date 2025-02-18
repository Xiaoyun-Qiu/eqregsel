*! Test lq_fnm.mata XQ 6 Aug2015

version 13

cd "e:/Stata/duke/do"

import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

gen afqt_2 = afqt^2

gen y = - log_wage

local indepvar black hispanic age afqt afqt_2




/*end
mata
 st_view(X=.,.,"black hispanic age afqt afqt_2")
st_view(Y=.,.,"y")
X=(J(rows(X),1,1),X)

m = rows(X)
		
	u = J(m,1,1)
	q = 0.2
	a = (1 - q) :* u
	
	 k = -lq_fnm(X', -Y', X' * a, u, a)'
	
	A = X'
	c=-Y'
	b = X'*a
	x=a
	
	s =  u - x
	y = (invsym(A * A') * A * c')'
	r = c - y * A
	r = r + 0.001 * (r :== 0)
	z = r :* (r :> 0)
	w = z - r
	gap = c * x - y * b + w * u
	
	
	
	
*/
