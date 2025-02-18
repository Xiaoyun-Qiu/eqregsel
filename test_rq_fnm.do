version 13

cd "e:/Stata/duke/data"

import excel using test.xlsx, firstrow clear

gen afqt2=afqt^2

do "e:/Stata/duke/do/rq_fnm.do"

rq_fnm low_wage black hispanic age,quantile(0.25)

mat list e(b)
