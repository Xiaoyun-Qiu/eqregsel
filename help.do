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
{bf:eqreg} {hline 5} Conduct Extremal Quantile Regression for Selection Models


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:eqreg:}
[{varlist}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt g:rid(#)}}grid the interval of tau into # subintervals; 
		default is grid(20){p_end}
{synopt:{opt btp:(#)}}repeat # times during bootstraping; default is btp(200);{p_end}
{synopt:{opt b:oots(#)}}sample size when calculating criteria for selecting tau{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:eqreg} is an estimation method based on the paper "Extremal Quantile 
Regression for Selection Models and the Black-White Wage Gap". Only the 
coefficient of homoskedastic variable can be estimated. This command estimates
the coefficient of only one homoskedastic variable. 


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt grid(#)} ...

{phang}
{opt btp(#)} restricts the calculation to be based on only the means.  The
default is to use a trimmed mean.

{phang}
{opt boots(#)} ...

{title:Saved results}

{phang}{cmd:eqreg} saves the following results in {cmd:e()}:

{phang}Variables{p_end}
{col 10}{cmd:y}{col 25}Variable storing the negative value of dependent variable 

{phang}Scalars{p_end}
{col 10}{cmd:e(tau0)}{col 25}Optimal quantile index selected
{col 10}{cmd:e(specificationtest)}{col 25}Value of the specification test

{phang}Matrices{p_end}
{col 10}{cmd:e(beta_hom)}{col 25}Matrix containing the estimated coefficient of the homoskedastic varaible
{col 10}{cmd:e(std_b)}{col 25}Matrix containing the bootstapped standard error of the homoskedastic varaible


{marker remarks}{...}
{title:Remarks}

{pstd}
Remember to put the the homoskedastic varaible as the first independent variable

{marker examples}{...}
{title:Examples}

{phang}{cmd:. eqreg log_wage black hispanic age AFQT0 AFQT0_2}{p_end}

{phang}{cmd:. eqreg log_wage black hispanic age AFQT0 AFQT0_2, grid(40) btp(300)}{p_end}

{title:Version requirements}

{p 4 4 2}This command requires Stata 12 or upper. 


{title:Methods and Formulas}

{p 4 6} See Xavier D'Haultfoeuille, Arnaud Maurel, Yichong Zhang (2013).

{title:Supporting Package}

{p 4 4 2}This command requires -moremata- package. Before calling this 
command, please intall -moremata- by typing

{phang}{cmd:. ssc install moremata}{p_end}
