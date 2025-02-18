*! Test bound.do and clr_bound.do Xiaoyun Qiu 31july2015

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
gen c5 = c4 * c4
	
scalar quant = 0.5
	
mata:st_view(gridd=.,.,"c1 c2 c3 c4 c5") 
mata:grid = gridd 
mata:st_matrix("r(grid)",grid)

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
gen c5 = c4 * c4

mata:st_view(gridd=.,.,"c1 c2 c3 c4 c5") 
mata:grid = gridd 
mata: gridid = ((grid[.,1] + grid[.,2]) :< 2)
mata: grid = select(grid,gridid)
mata:st_matrix("r(grid)",grid)
	
scalar quant = 0.5
	
import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

gen afqt_2 = afqt^2

do bound.do

/*
mata
st_view(Xdnp=.,.,"black hispanic age")
st_view(Xcnp=.,.,"afqt")
st_view(Y=.,.,"log_wage")
tau = 0.5
xdnum = cols(Xdnp)
xcnum = cols(Xcnp)
grid = st_matrix("r(grid)")
gridid = ((grid[.,1] + grid[.,2]) :< 2)
grid = select(grid,gridid)
grid = grid'
gridnp = grid[|1,1\xdnum+xcnum,.|]
y1=(0)
y2=(0)
yl=0
yr=9
clr_bound(Y,Xdnp,Xcnp,tau,gridnp,y1,y2,yl,yr)
y1
y2
*/

/*
mata
st_view(Z=.,.,"log_wage black hispanic age afqt afqt_2")
grid = st_matrix("r(grid)")
gridid = ((grid[.,1] + grid[.,2]) :< 2)
grid = select(grid,gridid)
bound(Z,grid)
*/





	
