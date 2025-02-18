*! Test myfun_combined.do Xiaoyun Qiu 22May2016

//version 13

//change the directory to the one containing the STATA codes and cooked.dta
cd D:\onedrive\Stata\DMZ\version3_June2016

use cooked_12,clear
	
qui do myfun_combined_v3_12.do
mat hom = (1,1,0,0,0)
gen y = - log_wage
eqreg  y black hispanic age AFQT0 AFQT0_2 , phi("hom") 


