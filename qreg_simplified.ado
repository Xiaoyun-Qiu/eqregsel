*! version 4.1.5  23may2013
program qreg_simplified, eclass byable(recall) prop(sw mi)

	`vv' `BY' Estimate `0'
end

program Estimate, eclass byable(recall) sort prop(sw mi)
	//local vc = _caller()
	//local cmdline : copy local 0
	//version 13, missing

	//if `vc' <= 12 {
		syntax varlist(numeric fv) [fw aw] [if] [in] [,		///
				Quantile(real 0.5) WLSiter(integer 1)	///
				noLOg Level(cilevel) * ]
	//}
	//else {
		//syntax varlist(numeric fv) [fw aw iw pw] [if] [in] [,	///
		//		Quantile(real 0.5) WLSiter(integer 1)	///
		//		vce(string) noLOg Level(cilevel) * ]
	//}
	_get_diopts diopts options, `options'
	//local fvops = ("`s(fvops)'"=="true" | `vc'>12)

	tempname quant
	local quant = `quantile'
	//if (`quant'>=1) local quant = `quant'/100

	if `quant' <= 0 | `quant' >= 1 {
		di as err "{bf:quantile(`quantile')} is out of range"
		exit 198
	}

	//if ("`vce'"=="cluster" | "`vce'"=="robust") local vcetype Robust
	marksample touse
	qui count if `touse'
	local N = r(N)
	if (!`N') error 2000
	if (`N'<=2) error 2001

	gettoken dep indep : varlist

	//if (`fvops') local mse1 mse1
	//else local rmopt forcedrop

	_rmcoll `indep' `weights' if `touse', `rmopt'
	local indep `r(varlist)' 
	local varlist `dep' `indep'


	tempvar r i p
	qui gen long `i' = _n

	/* initial estimates via weighted least squares 		*/
	_qregwls `varlist' `weights' if `touse', r(`r') ///
		iterate(`wlsiter') quant(`quant') `log'

	if ("`log'"!="") local qui quietly

	/* stable sort							*/
	sort `r' `i'
	drop `r'
	drop `i'

	`qui' _qreg `varlist' if `touse' `weights', quant(`quant') ///
			`opopts'

	tempname b 

	mat `b' = e(b)


end

