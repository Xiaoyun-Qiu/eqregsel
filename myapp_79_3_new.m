clear;
%clc;
%%
tic
load('cooked79_2.mat');
tic
female = cooked(:,10);
age = cooked(:,3);
insamp1 = cooked(:,11);
insamp2 = cooked(:,12);
Mgrades = cooked(:,15);
Fgrades = cooked(:,16);
parents = cooked(:,14);
Rgrades = cooked(:,17);
momdad = max(Mgrades,Fgrades);

idx = (insamp2==1).*(female ==0);
samp2 = cooked(idx==1,:);

log_wage = samp2(:,2);
black = samp2(:,5);
hispanic = samp2(:,4);
age = samp2(:,3);
AFQTO = samp2(:,6);
AFQTC = samp2(:,8);
sch = samp2(:,13);

Rgrades_high = (samp2(:,17)>11).*(samp2(:,17)<16);
Rgrades_col = (samp2(:,17)>15);
Fgrades_high = (samp2(:,16)>11).*(samp2(:,16)<16);
Fgrades_col = (samp2(:,16)>15);
Mgrades_high = (samp2(:,15)>11).*(samp2(:,15)<16);
Mgrades_col = (samp2(:,15)>15);
momdadhigh = (momdad(idx==1,:)>11).*(momdad(idx==1,:)<16);
momdadcol = (momdad(idx==1,:)>15);

Enroll90 = (samp2(:,18) == 2)+(samp2(:,18) == 3);
Enroll91 = (samp2(:,19) == 2)+(samp2(:,19) == 3);
parents = (samp2(:,14)==11);

%%



Y = log_wage;
X = [ones(length(Y),1),black,hispanic,age,AFQTO,AFQTO.^2];
boots = 600;
%%
phi = [1,0,0,0,0];
[beta_hom,std_b,specificationtest, tau0]=myfun_combined(X,Y,phi,boots);

cibootbb95ps = [beta_hom-1.96*std_b,beta_hom+1.96*std_b];
cibootbb975ps = [beta_hom-2.24*std_b,beta_hom+2.24*std_b];


toc;
%%
table = [beta_hom(1),std_b(1),specificationtest,tau0];

save(['app79_male_AFQTO_',date,'.mat']);






