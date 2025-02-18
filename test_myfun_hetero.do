*! Test myfun_hetero.do Xiaoyun Qiu 26july2015

version 13

cd "e:/Stata/duke/do"

import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

do myfun_hetero.do

gen afqt_2 = afqt^2

gen y = - log_wage

myfun_hetero y black hispanic age afqt afqt_2, tau(0.2)




