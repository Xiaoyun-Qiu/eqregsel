*! Test myfun_hcombined.do Xiaoyun Qiu 22May2016

//version 13

//change the directory to the one containing the STATA codes and cooked.dta
cd "F:\onedrive\Stata\duke\new code"

mata: phi = (1,0,0,0,0)
mata:st_matrix("r(phi)",phi)

use cooked,clear

//choose one of the versions 

/******************************************************************************\
the first version of STATA codes
\******************************************************************************/

do myfun_combined_v2.do

combined log_wage black hispanic age AFQT0 AFQT0_2 , phi("r(phi)") boots(600)

/******************************************************************************\
the second version of STATA codes
\******************************************************************************/

//do myfun_combined_v2.do
//combined log_wage black hispanic age AFQT0 AFQT0_2 , phi("r(phi)") boots(600)


/******************************************************************************\
the third version of STATA codes
\******************************************************************************/

//do myfun_combined_part1_v3.do
//estimator log_wage black hispanic age AFQT0 AFQT0_2 , phi("r(phi)") grid(3)
//do myfun_combined_part2_v3.do
//std_b log_wage black hispanic age AFQT0 AFQT0_2 , tau(0.3) phi("r(phi)") sigmahat("Sigmahat")
