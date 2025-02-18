*! Test myfun_hom.do Xiaoyun Qiu 31july2015

version 13

cd "e:/Stata/duke/do"

import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

do myfun_hom.do

gen afqt_2 = afqt^2

gen y = - log_wage

matrix l = (0.65, 0.85, 1.15, 1.45)

matrix phi = (1,1,0,1,1)

myfun_hom y black hispanic age afqt afqt_2, tau(0.4) m(1.2) l("l") phi("phi")

mat list e(par)
mat list e(Vb)



	
