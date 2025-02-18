*! Test myfun_combined.do Xiaoyun Qiu 22May2016

//version 13

//change the directory to the one containing the STATA codes and cooked.dta
cd F:\onedrive\Stata\DMZ\version3_June2016

use cooked_12,clear
	
qui do myfun_combined12.do
gen y = - log_wage

eqreg  y black hispanic age AFQT0 AFQT0_2


