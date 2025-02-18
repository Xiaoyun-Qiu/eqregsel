*! Test myfun_combined.do Xiaoyun Qiu 22May2016


//change the directory to the one containing cooked_12.dta
cd F:\onedrive\Stata\DMZ\version4_parallel_June2016

use cooked_12,clear

eqreg  log_wage black hispanic age AFQT0

