% script mytry.m is example to use the main function myfun_combined.m to
% produce column (4) of table 4 in the paper. The original dataset we use
% is the raw data for paper Neal & Johnson JPE 96. We want to thank Prof. Neal for his generosity of sharing the data with us.  
clear;
clc;
%% This part of the code loads the raw dataset and labels each column with their corresponding names. 
tic;
load('cooked79_2.mat');
female = cooked(:,10);
insamp1 = cooked(:,11);
insamp2 = cooked(:,12); %variable insamp2 is the indicator produced by Prof. Neal to indicate the sample used for their median regression. 
Mgrades = cooked(:,15);
Fgrades = cooked(:,16);
Rgrades = cooked(:,17);
momdad = max(Mgrades,Fgrades); % variable momdad stores the higher grade obained by respondent's parents. 
idx = (insamp2==1).*(female ==0);  % variable idx is a sample indicate. The sample we will choose is the one used to produce the second column of table 4 in Neal and Johnson JPE 96. 
samp2 = cooked(idx==1,:);

log_wage = samp2(:,2);
black = samp2(:,5);
hispanic = samp2(:,4);
age = samp2(:,3);
AFQTO = samp2(:,6); %variable AFQTO stores the original AFQT
AFQTC = samp2(:,8); %variable AFQTO stores the comparable AFQT merged from the dataset produced by Altonji et.al JOLE 2012
sch = samp2(:,13); % variable sch stores the year of schooling

Rgrades_high = (samp2(:,17)>11).*(samp2(:,17)<16); % respondent is high school graduates
Rgrades_col = (samp2(:,17)>15); % respondent is college graduates
Fgrades_high = (samp2(:,16)>11).*(samp2(:,16)<16); % respondent's father is high school graduates
Fgrades_col = (samp2(:,16)>15); % respondent's father is college graduates
Mgrades_high = (samp2(:,15)>11).*(samp2(:,15)<16); % respondent's mather is high school graduates
Mgrades_col = (samp2(:,15)>15); % respondent's mather is college graduates
momdadhigh = (momdad(idx==1,:)>11).*(momdad(idx==1,:)<16); % respondent's parents are at most high school graduates
momdadcol = (momdad(idx==1,:)>15); % at least one of respondent's parents is college graduates
parents = (samp2(:,14)==11);

%%

Y = log_wage;
Xc = [AFQTO,AFQTO.^2]; % the continuous independent variables  
Xd = [black,hispanic,age]; % the discrete independent variables


%% this section produces the grid which will be used to compute the CLR bound of QTE for variables whose correponding delta is nonzero. This part should be user-specific and are about 
% to change due to different data structrue.
tic;
space = 200; % space parameter is the number of grid use want to use for the continuous variable, which in this case is AFQT
AFQT = Xc(:,1); 
range = quantile(AFQT,0.95) - quantile(AFQT,0.05); % we use AFQT in the range of its 5% quantile to its 95% quantile
[black_gv,hispanic_gv,age_gv,AFQT_gv] = ndgrid(0:1,0:1,(26:28),quantile(AFQT,0.05)+(1:space)*range/(space+1)); % produce the grid 
grid = [black_gv(:),hispanic_gv(:),age_gv(:),AFQT_gv(:),AFQT_gv(:).^2];
grid = grid(grid(:,1)+grid(:,2)<2,:); % we remove the grid which has both variable black and hispanic equals 1 because it never happens in the data. 
grid = grid';
Xdnp = Xd; % variable Xdnp contains the discrete dependent variable used for nonparameter estimation. It can be potentially different from Xd because some of Xd can be fully determined by others.
Xcnp = Xc(:,1); % variable Xdnp contains the continuous dependent variable used for nonparameter estimation. In this case, Xcnp is not Xc because AFQT^2 is fully determined by AFQT
gridnp = grid(1:4,:); % varialbe gridnp constains the corresponding grid which is free of being determined. 
%%
quant = 0.5; % variable quant indicates at which quantile index we will compute the bound of QTE, by using CLR approach
[delta_hetero,beta_hetero,V_hetero_star,Nb_hetero_star,qstard,phi,delta_hom,beta_hom,Nb_hom_star,V_homd_star,V_homb_star,qstarbd,qstarbb,theta0] = myfun_combined(Xd,Xc,Y,Xdnp,Xcnp,quant,grid,gridnp);
%%
save(['mytry_AFQTO',date,'.mat']);
toc;