{smcl}
{* *! version 1.2.1  26sep2015}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "examplehelpfile##syntax"}{...}
{viewerjumpto "Description" "examplehelpfile##description"}{...}
{viewerjumpto "Options" "examplehelpfile##options"}{...}
{viewerjumpto "Remarks" "examplehelpfile##remarks"}{...}
{viewerjumpto "Examples" "examplehelpfile##examples"}{...}
{title:Title}

{phang}
{bf:myfun_combined} {hline 2} Extremal Quantile Regression for Selection Models 


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:my:fun_combined}
{varlist}
{ifin}
{cmd:,} d(#) df(#) cf(#) grid(string) [tau(#)]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt d:num(#)}}total number of dependent variables{p_end}
{synopt:{opt df:num(#)}}number of free discrete independent variables{p_end}
{synopt:{opt cf:num(#)}}number of free continuous independent variables{p_end}
{synopt:{opth grid:d(matrix name)}}grid matrix of independent variables {it:newvar}{p_end}
{synopt:{opt tau(#)}}optional, the quantile of interest of the CLR bound (default=0.5){p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:myfun_combined} is an estimation method based on the paper {it:Extremal Quantile 
Regression for Selection Models and the Black-White Wage Gap}. In particular, the 
code conducts an automatic pretest for homoskedasiticity and yields the point 
estimates and its standard deviation of beta and delta corresponding to the homoskedastic 
variable(s) along with the CLR bound of the tau-th marginal quantile treatment effect 
for the heteroskedasitic variable(s) where the quantile index tau is user-specified.

{pstd}
In a special note, the grid generated in this function should satisfy the data restriction. 
For example, in our case of using NLSY79 data, black and hispanic cannot be both 1 
while AFQT^2 is fully determined on AFQT. For the detail, please check the program.

{pstd}
The main function it calls is the MATA function myfun_combined(). In myfun_combined(), 
the default spacing parameters we use are m = 1.2, l = (0.65,0.85,1.15,1.45). The number of 
subsamples is set to be 500. The subsample size is (150,300,500,600) for sample 
size (250,500,1000,2000) and the corresponding linear interpolation for sample 
size in between. For sample size N larger than 2000, the subsample size is set to 
be 600+0.2*(N-2000). Also, we use the rule of thumb tuning parameter and Gaussian 
kernel to estimate the nonparametric pieces when computing the CLR bound. 


{title:Output}

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt beta_mse}}the estimator of beta without imposing partial homoskedasticity{p_end}
{synopt:{opt delta_mse}}the estimator of delta without imposing partial homoskedasticity{p_end}
{synopt:{opt V_mse_star}}the asympotic variance covariance matrix for delta_mse(and also beta_mse).{p_end}
{synopt:{opt Nb_mse_star}}the convergence rate for beta_mse{p_end}
{synopt:{opt qstard}}the optimal quantile index used to estimate both delta_mse and delta_mse {p_end}
{synopt:{opt beta_hom}}the point estimates of beta under partial homoskedasticity for those variables who do not pass the interest{p_end}
{synopt:{opt delta_hom}}the point estimates of delta under partial homoskedasticity for those variables who do not pass the interest{p_end}
{synopt:{opt qstarbd}}the optimal quantile index used to estimate delta_hom {p_end}
{synopt:{opt qstarbb}}the optimal quantile index used to estimate beta_hom {p_end}
{synopt:{opt theta0}}the CLR bound for {it:tau}-th marginal quantile treatment effect for those variables who do not pass the pretest{p_end}
{synoptline}
{p2colreset}{...}




{marker remarks}{...}
{title:Remarks}

{pstd}

{cmd:grid}

{pstd}
There is no general rules for generating the grid of independent variables, since 
it varies with data sets. The best way is to remember what the grid matrix looks 
like for your data set and find a way to generate it in STATA. 

{pstd}
Grid is an m*n matrix, where m = space * # of val1 * # of val2 *бн and n is the 
number of free independent variables. # of val1 means how many different values 
the first independent variable can take. For binary variables, it should be two 
since binary variables can only take 0 and 1. For continuous variables, we use 
space to measure the number of values. In our example, we take space as 200. 
Users can specify a different value for space based on their data sets.

{pstd}
Basically, the grid matrix is a matrix containing all combination scenarios given 
different values of each independent variables. 

different cases

{pstd}
1.With continuous independent variables

{pstd}
Usually you should start with generating grid columns for continuous independent 
variables, since they usually have a large number of different values. The best 
way I can find is to use egen var = fill(1(1)space 1(1)space). And then do some 
transformation to this variable.
For categorical variables and binary variables, use the similar command 
egen var = fill(a b c a b c). 
Before you generate next columns, always remember to sort the previous column you 
have generated.

{pstd}
2.With more than one related binary variables

{pstd}
We sometimes use more than one binary variables to indicate different categories 
each observation belongs. In our case, we have black and Hispanic indicating the 
race of each observations. Then for the grid matrix, remember to drop those 
scenarios where black and Hispanic are both 1. 

{cmd:CLR_bound}

{pstd}
Users should specify the range of values the dependent variables can take. 

{marker examples}{...}
{title:Examples}

Generating the grid matrix.

{phang}{cmd:. set obs 2400}{p_end}
{phang}{cmd:. egen tempc4 = fill(1(1)`space' 1(1)`space')}{p_end}
{phang}{cmd:. gen c4 =  `arange' * tempc4/(`space' + 1) + `p5'}{p_end}
{phang}{cmd:. drop tempc4}{p_end}
{phang}{cmd:. sort c4}{p_end}
{phang}{cmd:. egen c2 = fill(0 1 0 1) }{p_end}
{phang}{cmd:. sort c4 c2}{p_end}
{phang}{cmd:. egen c3 = fill(26 27 28 26 27 28)}{p_end}
{phang}{cmd:. sort c4 c3 c2}{p_end}
{phang}{cmd:. egen c1 = fill(0 1 0 1)}{p_end}
 

{title:User Guide}

{pstd}
The following steps show users how to use this command. 

{phang}Step1: install moremata{p_end}

{phang}Step2: Put myfun_combined.ado into C:/ado/PERSONAL or put myfun_combined.do in current directory{p_end}

{phang}Step3: Start a new script calling myfun_combined.ado.{p_end}

{phang}Step3a: Load the data set.{p_end}

{phang}Step3b: Generate the grid for independent variables in MATA and return it as a matrix to STATA.{p_end}

{phang}Step3c: Load the data set again and call the command.{p_end}

{marker author}{...}
{title:Author}

{phang}Authors: Xavier D'Haultfoeuille, Arnaud Maurel, Yichong Zhang{p_end}

{phang}Coded by Xiaoyun Qiu, converted from the MATLAB code myfun_combined.m{break} {p_end}

{title:Also see}
{phang}test_myfun_combined.do{p_end}

{title:Other package}
{phang}moremata{p_end}
