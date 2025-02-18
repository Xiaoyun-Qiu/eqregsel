Documentation of STATA code for “Extremal Quantile Regression for Selection Models and the Black-White Wage Gap”
Septempber 2015

Title
eqreg 

Syntax

eqreg varlist [if] [in] , [Grid(integer 20) BTP(integer 200) Boots(real 600)]

conditions	Description
Grid(#)	grid the interval of tau into # subintervals; default is grid(20)
BTP(#)	perform # bootstrap replications; default is btp(200)
Boots(#)	sample size when calculating criteria for selecting tau

Description

eqreg is an estimation method based on the paper Extremal Quantile Regression for Selection Models and the Black-White Wage Gap. Only the coefficient of homoskedastic variable can be estimated. This command estimates the coefficient of only one homoskedastic variable.

Output

Name	Description
e(tau0)	Optimal quantile index
e(specificationtest)	Value of specification test
e(beta_hom)	Coefficient of the 
e(std_b)	Standard error of the 

Remarks
1.	Remember to put the the homoskedastic varaible as the first independent variable. 
2.	The command generates a new variable y which stores the negative value of the dependent variable. After estimation, users can drop this variable.

Example

. eqreg log_wage black hispanic age AFQT0 AFQT0_2              
 (black is the homoskedatic variable to be estimated)

. eqreg log_wage black hispanic age AFQT0 AFQT0_2, grid(40) btp(300)

Version Requirements

This command requires Stata 12 or upper.

Methods and Formulas

See Xavier D'Haultfoeuille, Arnaud Maurel, Yichong Zhang (2013).

Supporting Package

This command requires -moremata- package. Before calling this command, please intall -moremata- by typing

 . ssc install moremata
