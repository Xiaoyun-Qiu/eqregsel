mata:
// version 13
mata clear
mata set matastrict on

void stat(string matrix Par_bootstraphomb, real scalar J, real scalar dbb, ///
		real scalar chibb	){
	real scalar  specificationtest
	real matrix std_b,par_bootstraphomb
	par_bootstraphomb = st_matrix(Par_bootstraphomb)
	std_b = mm_colvar(par_bootstraphomb')'                                      //Compute the standard deviation.
	std_b = sqrt(std_b)
	specificationtest = 1 - chi2(J*dbb,chibb)
	st_matrix("e(std_b)",std_b)
	st_numscalar("e(specificationtest)",specificationtest)
}

 
end
mata: mata mosave stat(),dir(PERSONAL) replace
