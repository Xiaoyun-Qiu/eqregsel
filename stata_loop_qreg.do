

	version 13

	cd "e:/Stata/duke/do"

	import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear
	gen afqt_2 = afqt^2
	qui count
	
	timer clear 1
	qui qreg log_wage black hispanic age afqt afqt_2
	mat beta_median = e(b)'
	mat beta_median_boot = J(6,500,0)
	
	
	
	forvalues j = 1/500{
		timer on 1
		bsample _N
		qui bsqreg log_wage black hispanic age afqt afqt_2
		mat beta_median_boot[1,`j'] = e(b)'
		timer off 1
	}
	
	timer list 1
	mat list beta_median_boot
	//mat ci_med95 = (mm_quantile(beta_median_boot',1,0.025)',mm_quantile(beta_median_boot',1,0.975)')
	//mat ci_med975 = (mm_quantile(beta_median_boot',1,0.0125)',mm_quantile(beta_median_boot',1,0.9875)')
	
