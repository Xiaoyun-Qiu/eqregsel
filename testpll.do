*£¡ A test file for parallel

cd F:\onedrive\Stata\DMZ\version4_parallel_June2016
capture program drop qreg_simplified12
qui do qreg_simplied12.ado
qui do myfun_hom_new.mata
qui do stat.mata
qui do IC.mata
use cooked_12,clear
gen y = - log_wage
local G = 2
local B = 12
matrix l = (0.65, 0.85, 1.15, 1.45)
local J = colsof(l)
local d = 6
qui count
local N = r(N) 
local boots = 600
local lower = min(80/`boots',0.1)
local upper = 0.3
local step = (`upper' - `lower')/`G'
mat hom = (1,0,0,0,0)
local phi hom
local varlist "y black hispanic age AFQT0 AFQT0_2"

mat ss = J(1, rowsof(`phi'),1) * `phi'
local dbb = ss[1,1]
mat par_bootb = J(`dbb',`B',0)
mat chi_bootbb = J(`B',1,0)
scalar disbb = 100
mat thetatemp = J(`d',`B',0)	
di in gr "Begin big loop"
di "$S_TIME"

	local tau 0.2
	qui qreg_simplified12 `varlist',quantile(`tau') cformat(%10.0g)
	mat theta = e(b)'


		forvalues bb = 1/`B'{
			set seed `bb'
			preserve
			qui bsample 
			qui qreg_simplified12 `varlist',quantile(`tau') cformat(%10.0g)
			mat thetatemp[1,`bb'] = e(b)'
			restore
		}


    
		//mat thetatemp = e(thetatemp)
		* Adjust the order of coefficients of independent variables
		//mat thetahat = (theta[`d',1]\theta[1..`d'-1,1])
		//mat thetadagger = (thetatemp[`d',1..`B']\thetatemp[1..`d'-1,1..`B'])
		//mat Sigma=(thetadagger-thetahat*J(1,`B',1))*(thetadagger-thetahat*J(1,`B',1))'/`B'
		


	


