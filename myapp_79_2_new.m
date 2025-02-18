clear;
%clc;
%%

load('cooked79_2.mat');
female = cooked(:,10);
age = cooked(:,3);
insamp1 = cooked(:,11);
insamp2 = cooked(:,12);
Mgrades = cooked(:,15);
Fgrades = cooked(:,16);
parents = cooked(:,14);
Rgrades = cooked(:,17);
momdad = max(Mgrades,Fgrades);
%idx = (insamp2==1).*(female ==0).*(momdad>-1).*(parents>-1);
%idx = (insamp2==1).*(female ==0).*(Rgrades<13).*(Rgrades>-1);
%idx = (insamp2==1).*(female ==0).*(age<28);
%idx = (insamp2==1).*(female ==0).*(age<28).*(momdad>-1).*(parents>-1);
%idx = (insamp2==1).*(female ==1);
%idx = (insamp2==1).*(female ==0).*(Rgrades>-1).*(momdad>-1).*(parents>-1);
%idx = (insamp2==1).*(female ==0).*(age<28);
idx = (insamp2==1).*(female ==0);
%idx = (insamp2==1).*(female ==0).*(age<28).*(momdad>-1).*(parents>-1);
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
%X = [ones(length(Y),1),black,hispanic,age,AFQTC,AFQTC.^2];
X = [ones(length(Y),1),black,hispanic,age,AFQTO,AFQTO.^2];
%X = [ones(length(Y),1),black,hispanic,age];
%X = [ones(length(Y),1),black,hispanic,age,AFQTC,AFQTC.^2,parents,momdadhigh,momdadcol];
%X = [ones(length(Y),1),black,hispanic,age,AFQTC,AFQTC.^2,parents,Rgrades_high,Rgrades_col,momdadhigh,momdadcol];

[N,d] = size(X);



%%
tic;
%G = 40; 
l = [0.65,0.85,1.15,1.45]; 
% l = [0.4, 0.7, 1.3, 2];
J = length(l);
%B = 500;
boots = 550;

B = 10;
G = 4;

m = 1.2;
b0 = 0;
d0 = 1;

ciboot1 = zeros(2*d-2,2);

%lower = 25/boots/mean(black); 
lower = min(80/boots, 0.1);
%lower = 0.05;
upper = 0.3;
step = (upper-lower)/G; 

phi = [1,1,0,1,1];
%phi = [1,1,0,1,1,1,1,1];
%phi = [1,1,0,1,1,1,1,1];
%phi = [1,1,0];
%phi = [1,1,1,1,1,1,1,1,1,1];
%phi = [1,1,0,1,1,1,1,1];
dd = d - 1 - sum(phi);
db = sum(phi);
V_mse = zeros(d-1,d-1,G);
nb_mse = zeros(G,1);
disd = zeros(G,1);disbd = disd;
disbb = disd;
disdbias = zeros(G,1);disbbbias = disdbias;disbdbias = disdbias;


par_boot = zeros(d-1,B,G);
par_mse = zeros(2*(d-1),G);
par_bootb = zeros(d-1,B,G);
par_hom = zeros(d-1,G);

options = statset('UseParallel','always') ;

display('Part 1: Setup')
toc;

median_d1 = median(random('chi2',(J-1)*(d-1),100000,1));
median_d2 = median(random('chi2',(J-1)*(d-1-sum(phi)),100000,1));
median_b1 = median(random('chi2',J*sum(phi),100000,1));
Nb_hom = zeros(G,1);
V_homb = zeros(db,db,G);
V_homd = zeros(dd,dd,G);
chid = zeros(G,1);
chibb = zeros(G,1);

tic;

for g = 1:G
warning off 

[tempmsef, V_mse(:,:,g), chid(g,1),nbtemp] = myfun_MSE_new(lower+step*g,m,b0,d0,X,Y,l);
tempmse = tempmsef(1:(d-1));
par_mse(:,g) = tempmsef;
nb_mse(g,1) = nbtemp;
[temphom, V_homd(:,:,g), V_homb(:,:,g),Nb_hom(g,1),~,chibb(g,1)] = myfun_hom(lower+step*g,m,phi,X,Y,l);
par_hom(:,g) = temphom;
end

display('Part 2: First loop (G)')
toc;
%%
Z = [Y,X];
tempchi_bootd = zeros(B,G);
tempchi_bootbd = zeros(B,G);
tempchi_bootbb = zeros(B,G);

tic;

parfor b = 1:B
Zz = Z(randperm(N,boots),:);      
tempboot = zeros(3*(d-1),G);
for g = 1:G
    tau = lower+step*g;
    [tempboot(1:2*(d-1),g), ~,tempchi_bootd(b,g)] = myfun_MSE_new(tau,m,0,1,Zz(:,2:end),Zz(:,1),l);
    [tempboot(2*(d-1)+1:end,g), ~, ~, ~,tempchi_bootbd(b,g),tempchi_bootbb(b,g)] = myfun_hom(tau,m,phi,Zz(:,2:end),Zz(:,1),l);
end
par_boot(:,b,:) = tempboot(1:(d-1),:);
par_bootb(:,b,:) = tempboot(2*(d-1)+1:end,:);
end

disdbias(:,1) = abs(median(tempchi_bootd) - median_d1);
disbdbias(:,1) = abs(median(tempchi_bootbd) - median_d2);
disbbbias(:,1) = abs(median(tempchi_bootbb) - median_b1);

display('Part 3: second loop(B*G)')
toc;






%%

tic;

disdvar = zeros(G,1);
disbdvar = disdvar;
disbbvar = disdvar;
for g = 1:G
        tempboot1 = par_boot(:,:,g);
        tempmse = par_mse(1:d-1,g);
        tempboot2 = par_bootb(:,:,g);
        tempbootd = tempboot2(1:dd,:);
        tempbootb = tempboot2(dd+1:end,:);
        temphom = par_hom(:,g);

        temphomd = temphom(1:dd,:);
        temphomb = temphom(dd+1:end,:);

disdvar(g,1) = mean(sum((tempboot1 - repmat(mean(tempboot1,2),1,B)).^2));
disbdvar(g,1) = mean(sum((tempbootd - repmat(mean(tempbootd,2),1,B)).^2));
disbbvar(g,1) = mean(sum((tempbootb - repmat(mean(tempbootb,2),1,B)).^2));
        
        
% disdmse(g,1) = mean(sum((tempboot1 - repmat(tempmse,1,B)).^2));
% disbdmse(g,1) = mean(sum((tempbootd - repmat(temphomd,1,B)).^2));
% disbbmse(g,1) = mean(sum((tempbootb - repmat(temphomb,1,B)).^2));
end
for g = 1:G
tau = lower+step*g;    
disd(g,1) = disdvar(g,1)*boots/N + disdbias(g,1)/sqrt((tau*boots));
disbd(g,1) = disbdvar(g,1)*boots/N+disbdbias(g,1)/sqrt((tau*boots));
disbb(g,1) = disbbvar(g,1)*boots/N+disbbbias(g,1)/sqrt((tau*boots));
% disd(g,1) = disdmse(g,1)*boots/N + sqrt(disdbias(g,1)/(tau*boots));
% disbd(g,1) = disbdmse(g,1)*boots/N+sqrt(disbdbias(g,1)/(tau*boots));
% disbb(g,1) = disbbmse(g,1)*boots/N+sqrt(disbbbias(g,1)/(tau*boots));


end

[A,gstard] = min(disd);
delta_mse = par_mse(1:d-1,gstard);
beta_mse = par_mse(d:2*(d-1),gstard);

[A,gstarbd] = min(disbd);
delta_hom = par_hom(1:dd,gstarbd);
[A,gstarbb] = min(disbb);
beta_hom = par_hom(dd+1:end,gstarbb);

display('Part4: two loops')
toc;

%% bootstrap

tic;

for b = 1:B    
    Zz = Z(randsample(N,N,true),:);          
     [tempboot, ~,~,~] = myfun_MSE_new(lower+step*gstard,m,0,1,Zz(:,2:end),Zz(:,1),l);
     [temp_bootstraphomb, ~, ~,~,~,~] = myfun_hom(lower+step*gstarbb,m,phi,Zz(:,2:end),Zz(:,1),l);
     [temp_bootstraphomd, ~, ~,~,~,~] = myfun_hom(lower+step*gstarbd,m,phi,Zz(:,2:end),Zz(:,1),l);
      %[tempboothom(b,g), ~, ~,~,tempchi_bootbd(b,g),tempchi_bootbb(b,g)] = myfun_hom(lower+step*g,m,phi,Zz(:,2:end),Zz(:,1),l);

par_bootstrap(:,b) = tempboot(:,1)';
par_bootstraphomb(:,b) = temp_bootstraphomb((dd+1):(d-1),1);
par_bootstraphomd(:,b) = temp_bootstraphomd(1:dd,1);

end



%%

cibootd95t = [delta_mse-1.96*sqrt(diag(V_mse(:,:,gstard))/(N*(lower+step*gstard))),delta_mse+1.96*sqrt(diag(V_mse(:,:,gstard))/(N*(lower+step*gstard)))];
cibootb95t = [beta_mse-1.96*sqrt(diag(V_mse(:,:,gstard))/(nb_mse(gstard))),beta_mse+1.96*sqrt(diag(V_mse(:,:,gstard))/(nb_mse(gstard)))];
cibootdlogt = [delta_mse-7*sqrt(log(N))*sqrt(diag(V_mse(:,:,gstard))/(N*(lower+step*gstard))),delta_mse+7*sqrt(log(N))*sqrt(diag(V_mse(:,:,gstard))/(N*(lower+step*gstard)))];

cibootbb95t = [beta_hom-1.96*sqrt(diag(V_homb(:,:,gstarbb)))/Nb_hom(gstarbb,1),beta_hom+1.96*sqrt(diag(V_homb(:,:,gstarbb)))/Nb_hom(gstarbb,1)];

cibootd975t = [delta_mse-2.24*sqrt(diag(V_mse(:,:,gstard))/(N*(lower+step*gstard))),delta_mse+2.24*sqrt(diag(V_mse(:,:,gstard))/(N*(lower+step*gstard)))];
cibootb975t = [beta_mse-2.24*sqrt(diag(V_mse(:,:,gstard))/(nb_mse(gstard))),beta_mse+2.24*sqrt(diag(V_mse(:,:,gstard))/(nb_mse(gstard)))];

cibootbb975t = [beta_hom-2.24*sqrt(diag(V_homb(:,:,gstarbb)))/Nb_hom(gstarbb,1),beta_hom+2.24*sqrt(diag(V_homb(:,:,gstarbb)))/Nb_hom(gstarbb,1)];
%% bootstrap CI
ciboot95p = [quantile(par_bootstrap,0.025,2),quantile(par_bootstrap,0.975,2)];
ciboot975p = [quantile(par_bootstrap,0.0125,2),quantile(par_bootstrap,0.9875,2)];


cibootbd95p = [quantile(par_bootstraphomd,0.025,2),quantile(par_bootstraphomd,0.975,2)];
cibootbd975p = [quantile(par_bootstraphomd,0.0125,2),quantile(par_bootstraphomd,0.9875,2)];


cibootbb95p = [quantile(par_bootstraphomb,0.025,2),quantile(par_bootstraphomb,0.975,2)];
cibootbb975p = [quantile(par_bootstraphomb,0.0125,2),quantile(par_bootstraphomb,0.9875,2)];


%% bootstrap std
%std_msedboot = (quantile(par_bootstrap(1:d-1,:),0.75,2) - quantile(par_bootstrap(1:d-1,:),0.25,2))/(norminv(0.75) - norminv(0.25));
std_msedboot = std(par_bootstrap(1:d-1,:),[],2);

cibootd95ps = [delta_mse-1.96*std_msedboot,delta_mse+1.96*std_msedboot];
cibootd975ps = [delta_mse-2.24*std_msedboot,delta_mse+2.24*std_msedboot];

std_msebboot = std(par_bootstrap(d:end,:),[],2);

cibootb95ps = [beta_mse-1.96*std_msebboot,beta_mse+1.96*std_msebboot];
cibootb975ps = [beta_mse-2.24*std_msebboot,beta_mse+2.24*std_msebboot];

std_homdboot = std(par_bootstraphomd,[],2);

cibootbd95ps = [delta_hom-1.96*std_homdboot,delta_hom+1.96*std_homdboot];
cibootbd975ps = [delta_hom-2.24*std_homdboot,delta_hom+2.24*std_homdboot];


std_hombboot = std(par_bootstraphomb,[],2);

cibootbb95ps = [beta_hom-1.96*std_hombboot,beta_hom+1.96*std_hombboot];
cibootbb975ps = [beta_hom-2.24*std_hombboot,beta_hom+2.24*std_hombboot];

display('part5: B*3')
toc;
%%

tic;

specificationtest = [cdf('chi2',chid(gstard),(J-1)*(d-1)),cdf('chi2',chibb(gstarbb),J*sum(phi))];

beta_median  = rq_fnm(X, Y, 0.5);
beta_median_boot = zeros(d,500);
for j = 1:500
    Zz = Z(randi(N,N,1),:);
    beta_median_boot(:,j) = rq_fnm(Zz(:,2:end), Zz(:,1), 0.5);
end
ci_med95 = [quantile(beta_median_boot,0.025,2), quantile(beta_median_boot,0.975,2)];
ci_med975 = [quantile(beta_median_boot,0.0125,2), quantile(beta_median_boot,0.9875,2)];

display('part6')
toc;
%%

save(['app79_male_AFQTO',num2str(boots/1000),'_',num2str(G),'_',num2str(upper),'_',date,'.mat']);

%%
% 
% D = (Y>0);
% White = 1 - black.*hispanic;
% Dblack = D(black==1,:);
% AFQTblack = AFQT(black==1,:);
% Nblack = length(Dblack);
% 
% hblack = 2.34*Nblack^(-1/5)*std(AFQTblack);
% 
% P(1,:) = ker(Dblack,AFQTblack,hblack);
% 
% Dhispanic = D(hispanic==1,:);
% AFQThispanic = AFQT(hispanic==1,:);
% Nhispanic = length(Dhispanic);
% 
% hhispanic = 2.34*Nhispanic^(-1/5)*std(AFQThispanic);
% 
% P(2,:) = ker(Dhispanic,AFQThispanic,hhispanic);
% 
% 
% Dwhite = D(white==1,:);
% AFQTwhite = AFQT(white==1,:);
% Nwhite = length(Dwhite);
% 
% hwhite = 2.34*Nwhite^(-1/5)*std(AFQTwhite);
% 
% P(3,:) = ker(Dwhite,AFQTwhite,hwhite);
% 
% AFQTmean = [mean(AFQTblack),mean(AFQThispanic),mean(AFQTwhite)] ;
% 
% %% 
% plot(-0.9+(1:180)*0.01,P(1,:),'black-',-0.9+(1:180)*0.01,P(2,:),'black--',-0.9+(1:180)*0.01,P(3,:),'black-.','lineWidth',2);
% xlabel('AFQT');
% % ylabel('Probability of Employment');
% legend(['Blacks (', num2str(AFQTmean(1)),')'],['Hispanics (', num2str(AFQTmean(2)),')'],['Whites (', num2str(AFQTmean(3)),')'],'Location','southeast')
% title('NLSY79: Conditional Probability of Employment on AFQT');
% saveas(gcf,['lew_conden_79',date,'.fig']);
% print(gcf,'-dpng','-r300')
% 
% 
% %%
% K = @(u) 3/4*(1-u.^2).*(abs(u)<1);
% pdfAFQT = zeros(3,180);
% for r = 1:180
%     
% pdfAFQT(1,r) = mean(K((AFQTblack - (-0.9)+r*0.01)/hblack))/hblack;
% pdfAFQT(2,r) = mean(K((AFQThispanic - (-0.9)+r*0.01)/hhispanic))/hhispanic;
% pdfAFQT(3,r) = mean(K((AFQTwhite - (-0.9)+r*0.01)/hwhite))/hwhite;
% end
% 
% plot(-0.9+(1:180)*0.01,pdfAFQT(1,:),'black-',-0.9+(1:180)*0.01,pdfAFQT(2,:),'black--',-0.9+(1:180)*0.01,pdfAFQT(3,:),'black-.','lineWidth',2);
% xlabel('AFQT');
% % ylabel('Probability of Employment');
% legend(['Blacks (', num2str(AFQTmean(1)),')'],['Hispanics (', num2str(AFQTmean(2)),')'],['Whites (', num2str(AFQTmean(3)),')'],'Location','southwest')
% title('NLSY79: Density of AFQT');
% saveas(gcf,['lew_79_den',date,'.fig']);
% print(gcf,'-dpng','-r300')

%% OLS 
DD = (log_wage>0);
Yols = log_wage(DD==1,:);
Xols = X(DD==1,:);
Beta_ols = (Xols'*Xols)\(Xols'*Yols);
beta_ols_boot = zeros(d,500);
Zols = [Yols,Xols];
beta_median_std = std(beta_median_boot,[],2);

for j = 1:500
    Zz = Zols(randi(sum(DD),sum(DD),1),:);
    beta_ols_boot(:,j) = (Zz(:,2:end)'*Zz(:,2:end))\(Zz(:,2:end)'*Zz(:,1));
end
beta_ols_std = std(beta_ols_boot,[],2);

ci_ols95 = [quantile(beta_ols_boot,0.025,2), quantile(beta_ols_boot,0.975,2)];
ci_ols975 = [quantile(beta_ols_boot,0.0125,2), quantile(beta_ols_boot,0.9875,2)];

ER(1,1) = sum(DD.*black)/sum(black);
ER(1,2) = sum(DD.*(1-black-hispanic))/sum(1-black-hispanic);







