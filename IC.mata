mata:
// version 13
mata clear
mata set matastrict on

void IC(string matrix Chi_bootbb, string matrix Par_bootb, real scalar J, ///
			real scalar dbb, real scalar B, real scalar boots, real scalar N, ///
			real scalar tau){
	real scalar median_b1, disbbbias, disbbvar, disbbtemp
	real matrix tempmean, chi_bootbb, par_bootb
	chi_bootbb = st_matrix(Chi_bootbb)
	par_bootb = st_matrix(Par_bootb)
	
	median_b1 = mm_median(rchi2(100000,1,J*dbb))
	disbbbias=abs(mm_median(chi_bootbb)*boots/N-median_b1)'
	tempmean = J(1,B,mean(par_bootb')')
	disbbvar = mean(colsum((par_bootb-tempmean):^2)')
	disbbtemp = disbbvar*boots/N + disbbbias/sqrt(tau*boots)
	st_numscalar("e(disbbtemp)",disbbtemp)
}
end
mata: mata mosave IC(),dir(PERSONAL) replace
