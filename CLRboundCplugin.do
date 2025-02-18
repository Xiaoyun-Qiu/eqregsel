*! CLR_bound with C_plugin solving the nonlinear equation, Xiaoyun Qiu,10Sep2015

capture program drop CLRboundCplugin
program CLRboundCplugin,eclass
	version 13
	syntax varlist(min=2,numeric), XDnum(real) XCnum(real) GRIDd(string)
	marksample touse
	tokenize `varlist'
	
	mat grid = (`gridd')'
	mat gridnp = grid[|1,1\`xdnum'+`xcnum',.|]
	
	
end
