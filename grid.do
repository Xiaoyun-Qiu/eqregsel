*! generating grid XQ 17Aug2015

version 13
cd "e:/Stata/duke/do"
import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

qui summarize afqt,d
local p5 = r(p5)
local p95 = r(p95)
local arange = `p95' - `p5'
local space = 200

clear

set obs 2400
egen tempc4 = fill(1(1)`space' 1(1)`space')
gen c4 =  `arange' * tempc4/(`space' + 1) + `p5'
drop tempc4
sort c4
egen c1 = fill(0 1 0 1) 
egen c2 = fill(0 0 1 1 0 0 1 1)
egen c3 = fill(26 26 26 26 27 27 27 27 28 28 28 28 26 26 26 26 27 27 27 27 28 28 28 28)
order c4,last


