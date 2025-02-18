{smcl}
{* *! version 1.2.1  18june2016}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "eqreg##syntax"}{...}
{viewerjumpto "Description" "eqreg##description"}{...}
{viewerjumpto "Options" "eqreg##options"}{...}
{viewerjumpto "Remarks" "eqreg##remarks"}{...}
{viewerjumpto "Examples" "eqreg##examples"}{...}
{title:Title}

{phang}
{bf:eqreg} {hline 5} Estimation Method Using Extremal Quantile Regression under Endogenous Selection


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:eqreg:}
[{varlist}]
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Options}
{synopt:{opt h:om(#)}}Estimate the coefficients of the first # independent variables; 
		default is hom(1){p_end}
{synopt:{opt b:oots(#)}}Subsampling size used for selecting optimal quantile index tau; 
		default is based on the size of the dataset{p_end}
{synopt:{opt g:rid(#)}}Discretize the interval of tau into # subintervals;  
		default is grid(40){p_end}
{synopt:{opt btp:(#)}}Perform # bootstrap replications; default is btp(150);{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:eqreg} estimates and provides inference for the coefficient of variables of interest under 
endogenous selection. The estimation method mainly assumes that under endogenous selection, the 
effect of the outcome on selection dominates those of covariates for sufficiently large values 
of the outcome, which is proposed by d'Haultfoeuille, X., Maurel, A., & Zhang, Y. (2016), based 
on the paper Extremal Quantile Regression for Selection Models and the Black-White Wage Gap. The 
variables of interest are assumed to have homogeneous effects on the outcome across its distribution. 
The number of variables of interest can be more than one.

{title:Saved results}

{phang}{cmd:eqreg} saves the following results in {cmd:e()}:

{phang}Scalars{p_end}
{col 10}{cmd:e(tau0)}{col 32}Optimal quantile index selected
{col 10}{cmd:e(specificationtest)}{col 32}Value of the specification test
{col 10}{cmd:e(boots)}{col 32}Subsampling size in the bootsratpping process
{col 10}{cmd:e(homvar)}{col 32}Number of variable(s) with homogeneous effect(s) on the outcome


{phang}Matrices{p_end}
{col 10}{cmd:e(beta_hom)}{col 32}Estimated coefficient(s) of interest
{col 10}{cmd:e(std_b)}{col 32}Standard error of the estimator(s)


{marker remarks}{...}
{title:Remarks}

{pstd}
Users should specify the number of variables of interest by using the option hom(). Apart from 
specifying number of variables whose coefficients are to be estimated, users should put names of 
those variables right after the dependent variable, followed by other independent variables. Users 
can specify more than one variables whose coefficients are to be estimated. During computation, the 
command will display the estimated time and progress bar. After estimation, results will be displayed.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. eqreg log_wage black hispanic age AFQT0 AFQT0_2} {space 13} (black is the variable of interest){p_end}

{phang}{cmd:. eqreg log_wage black hispanic age AFQT0 AFQT0_2, hom(2)} {space 5} (black and hispanic are variables of interest){p_end}

{phang}{cmd:. eqreg log_wage black hispanic age AFQT0 AFQT0_2, btp(300)}{p_end}

{title:Version requirements}

{p 4 4 2}This command requires Stata 12 or upper. 


{title:Methods and Formulas}

{p 4 6} See Xavier D'Haultfoeuille, Arnaud Maurel, Yichong Zhang (2016), Extremal Quantile Regression for Selection Models and the Black-White Wage Gap, Mimeo.

{title:Supporting Package}

{p 4 4 2}This command requires -moremata- package. Before calling this 
command, please intall -moremata- by typing

{phang}{cmd:. ssc install moremata}{p_end}
