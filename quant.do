*!Quantile Regression, Xiaoyun Qiu 26july2015

capture program drop quant
program quant,eclass 
	version 13
	syntax varlist(min=2 numeric) [if] [in] [,  Tau(real 0.3)]
	marksample touse
	
	///parameters
	scalar m = 1.2
	scalar b0 = 0
	scalar d0 = 1
	matrix l = (0.65, 0.85, 1.15, 1.45)
	matrix l_col = colsof(l)
	local J = l_col[1,1]

	if "`tau'" == "" {
		display "Please specify a quantile index!"
	}
	else{
	
		local d: word count `varlist'

		tempvar trow
		local trow = (`d'-1)*`J'
		mat tempy = J(`trow',1,0)
		mat tempx1 = J(`J',1,0)
		
		qui qreg `varlist',quantile(`tau') cformat(%10.0g)
		mat beta1 = e(b)'
		mat pf = beta1
		
		tempvar quant tau1 tau2 
		local tau1 = `tau'	
		local quant = `tau'*m
		
		qui qreg `varlist',quantile(`quant') cformat(%10.0g)
		mat beta2 = e(b)'
		mat pm = beta2
		
		forvalues j = 1/`J'{
			local tau1 = `tau'*l[1,`j']
			qui qreg `varlist',quantile(`tau1') cformat(%10.0g)
			mat beta1`j' = e(b)'
			mat pf = (pf, beta1`j')			
			
			tempvar ind
			local ind = (`d'-1)*(`j'-1)+1
			mat tempy[`ind',1] = (pf[1..`d'-1,`j'+1] - pf[1..`d'-1,1])*d0
			mat tempx1[`j',1] = pf[`d',`j'+1] - pf[`d',1]
			
			local tau2 = `quant'* l[1,`j']
			qui qreg `varlist',quantile(`tau2') cformat(%10.0g)
			mat beta2`j' = e(b)'
			mat pm = (pm, beta2`j')
		} 
		
		mat pf[`d',1] = pf[`d',1...] + J(1,`J'+1,b0)
		
		
	}
end
