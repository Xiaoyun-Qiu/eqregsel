*! Test mata functions Xiaoyun Qiu 26july2015

//------------------------------------------------------------------------------
/*
Instructions:
This do file aims to test the mata functions, including myfun_hetero.mata, 
myfun_hom.mata and myfun_combined.mata. 

In order to test one or all of the above functions, please follow the steps.

Step1: Run this do-file.

Step2: Open the particular mata file that you want to test and run it.

Step3: Copy and paste the corresponding lines in this do-file into the windows.

Step4: After it stops(you would find a new blank line begins with :, if it stops,
or you can find the red cross turns dark), enter the names of matrices that you 
want to check, e.g. Out, or par, etc..

Step5: If you want to check the results of another mata function, enter "end"
to exist Mata before you repeat the above steps.

Note1: Step1 and Step2 can be exchanged.

Note2: If Stata returns error "lq_fnm() cannot find" or "rq_fnm() cannot find",
open rq_fnm() and run it (remember to exit Mata first by typing "end") before 
you do step3.

Note3: If other problems show up, please email me. Thanks!
*/
//------------------------------------------------------------------------------



version 13

cd "e:/Stata/duke/do"

import excel using "e:/Stata/duke/data/test1.xlsx", firstrow clear

do testfullpar10_7.mata

gen afqt_2 = afqt^2

gen y = - log_wage

//myfun_hetero y black hispanic age afqt afqt_2, tau(0.2)

//------------------------------------------------------------------------------
//Copy and paste the following lines into the windows to test myfun_hetero.mata
/*
mata
boots = 550
st_view(Z=.,.,"log_wage black hispanic age afqt afqt_2")
tempZ = jumble(Z)
Zz = tempZ[|1,1\boots,.|]
X = Zz[|1,2\.,.|]
X = (J(boots,1,1), X)
Y = Zz[|1,1\.,1|]
b0=0
d0=1
l = (0.65, 0.85, 1.15, 1.45)
tau = 0.2
Var=(0)
Out=(0)
Dis=0
Nb1=0
mm=1.2
myfun_hetero(tau, mm, b0, d0, X, Y, l, Out, Dis, Var, Nb1)
testfullpar( Z			)

*/

//------------------------------------------------------------------------------
//Copy and paste the following lines into the windows to test myfun_hom.mata
/*
mata
st_view(X=.,.,"black hispanic age afqt afqt_2")
st_view(Y=.,.,"y")
X=(J(rows(X),1,1),X)

l = (0.65, 0.85, 1.15, 1.45)
tau = 0.2
par=(0)
Vb=(0)
Vd=0
Nb=0
dis_b=0
dis_d=0
phi=(1,1,0,1,1)
mm=1.2
myfun_hom(tau, mm,phi,X,Y,l,par,dis_d,dis_b,Vd,Vb,Nb)
*/

//------------------------------------------------------------------------------
//Copy and paste the following lines into the windows to test myfun_combined
/*
mata
temp_myfun_combined("y black hispanic age afqt afqt_2")
*/


