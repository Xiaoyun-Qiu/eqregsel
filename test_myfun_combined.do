*! Test myfun_combined.do Xiaoyun Qiu 26july2015

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
egen c2 = fill(0 1 0 1) 
sort c4 c2
egen c3 = fill(26 27 28 26 27 28)
sort c4 c3 c2
egen c1 = fill(0 1 0 1)
order c4,last
order c1,first
gen c5 = c4 * c4
	
mata:st_view(gridd=.,.,"c1 c2 c3 c4 c5") 
mata:grid = gridd 
mata: gridid = ((grid[.,1] + grid[.,2]) :< 2)
mata: grid = select(grid,gridid)
mata:st_matrix("r(grid)",grid)
	
import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

//do myfun_combined.do

gen afqt_2 = afqt^2

qui do myfun_combined.do

//syntax
//myfun_combined log_wage black hispanic age afqt afqt_2 ,d(3) df(3) cf(1) lb(6) rb(9) quant(0.5) gridd("r(grid)") 

myfun_combined log_wage black hispanic age afqt afqt_2 , d(3) df(3) cf(1) gridd("r(grid)") 

