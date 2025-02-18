*! Test quant.do Xiaoyun Qiu 26july2015

version 13

cd "e:/Stata/duke/do"

import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

do quant.do

gen afqt_2 = afqt^2

gen y = - log_wage

quant y black hispanic age afqt afqt_2, tau(0.1)


mat list pf, format(%10.0g)
mat list pm, format(%10.0g)
mat list tempy, format(%10.0g)
mat list tempx1, format(%10.0g)


