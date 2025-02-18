function [delta_hetero,beta_hetero,V_hetero_star,Nb_hetero_star,qstard,phi,delta_hom,beta_hom,Nb_hom_star,V_homd_star,V_homb_star,qstarbd,qstarbb,theta0] = myfun_combined(Xd,Xc,Y,Xdnp,Xcnp,quant,grid,gridnp)
%MYFUN_COMBINED Summary of this function goes here
%%%  Input
%  Xd: a matrix of discrete dependent variables
%  Xc: a matrix of discrete dependent variables
%  Y: a vector of independent variable
%  Xdnp: a matrix of free discrete dependent variables (free means it is not determined by other components of X = [Xd,Xc]) 
%  Xcnp: a matrix of free continuous dependent variables
%  quant: the quantile of interest of the CLR bound
%  grid: the grid of X used to compute the CLR bound
%  gridnp: the grid of free X = [Xdnp,Xcnp] used in nonparametric estimation when computing CLR bound 
%%%  Output
%  delta_hetero: the estimator of delta without imposing partial homoskedasticity
%  beta_hetero: the estimator of beta without imposing partial homoskedasticity
%  V_hetero_star: the asymptotic variance covariance for delta_hetero (and also beta_hetero). The convergence rate for delta_hetero is sqrt(N tau) where N is the total sample size and tau is the quantile index used to compute delta_hetero 
%  Nb_hetero_star: the convergence for beta_hetero 
%  qstard: the optimal quantile index used to estimate both delta_hetero and beta_hetero 
%  phi: a 1 by d vector contains the results of pretest where d is the dimension of X (exclude intercept). If j-th entry is 1, then it indicates the j-th variable is homoskedastic. 
%  delta_hom: the point estimates of delta under partial homoskedasticity for those variables who do not pass the pretest 
%  beta_hom: the point estimates of beta under partial homoskedasticity for those variables who pass the pretest
%  Nb_hom_star: the convergence rate of beta_hom
%  V_homd_star: the variance covariance matrix of delta_hom
%  V_homb_star: the variance covariance matrix of beta_hom
%  qstarbd: the optimal quantile index used to estimate delta_hom
%  qstarbb: the optimal quantile index used to estimate beta_hom
%  theta0: the CLR bound for quant-th marginal quantile treatment effect for those variables who do not pass the pretest 

X = [ones(length(Y),1),Xd,Xc];
[N,d] = size(X);
b0 = 0;
d0 = 1;
G = 40; 
l = [0.65,0.85,1.15,1.45]; 
B = 500;
boots = floor(0.6*N*(N<500) + (300 + 0.4*(N-500))*(N<1000)+(500+1/10*(N-1000))*(N<1000)+(600+0.2*(N-2000)));
R = 200;
m = 1.2;

J = length(l);

lower = min(80/boots, 0.1);
upper = 0.3;
step = (upper-lower)/G; 
% notation-wise, everything ended with "_hetero" is computed without imputing zero and everything ended with "_hom" is computed when imputing zeros.  
V_hetero = zeros(d-1,d-1,G); % variance covariance matrix for delta_hetero  
nb_hetero = zeros(G,1); % convergence rate of beta_hetero for each quantile
disd = zeros(G,1); % the object function we use to select tau for the estimation of delta without imputing zero
disbd = disd; % the object function we use to select tau for the estimation of delta when zeros are imputed
disbb = disd; % the object function we use to select tau for the estimation of delta when zeros are imputed
% the object funtion takes the form of |bias| + variance, bias is stored in disdbias, disbdbias and disbbbias 
disdbias = zeros(G,1);disbbbias = disdbias;disbdbias = disdbias; % in terms of notation, ***d: something about delta without imputing zero, ***bd: delta when imputing zero, ***bb: beta when imputing zero
chid = zeros(G,1);
chibb = zeros(G,1);

par_boot = zeros(d-1,B,G); % it is a estimator of delta without imputing zero, for each subsample and each quantile,  
par_hetero = zeros(2*(d-1),G); % it stores the estimator of delta and beta computed from the whole sample and for each quantile, 
par_bootb = zeros(d-1,B,G); % it stores the estimator of nonzero delta and beta (with corresponding delta=0) for each subsample and each quantile
par_hom = zeros(d-1,G); % it stores the estimator of nonzero delta and beta (with corresponding delta=0) for the whole sample and each quantile

median_d1 = median(random('chi2',(J-1)*(d-1),100000,1)); % it is the median of chi square distribution with degree of freedom = (J-1)*(d-1)

%% compute the estimator of beta and delta without imputing zero 
for g = 1:G
warning off 
[tempmsef, V_hetero(:,:,g), chid(g,1),nbtemp] = myfun_hetero(lower+step*g,m,b0,d0,X,Y,l);
par_hetero(:,g) = tempmsef;
nb_hetero(g,1) = nbtemp;
end

%% compute the subsample bias
Z = [Y,X];
tempchi_bootd = zeros(B,G); % temporarily store minimal distance when delta is not imputed as zero
tempchi_bootbd = zeros(B,G); % temporarily store the minimal distance for delta when imputing zeros
tempchi_bootbb = zeros(B,G); % temporarily store the minimal distance for beta when imputing zeros

parfor b = 1:B
Zz = Z(randperm(N,boots),:);    % shuffle the data and draw a subsample of them   
tempboot = zeros(2*(d-1),G);    % temporarily store the estimates of delta and beta without imputing zero
for g = 1:G
    warning off 
    tau = lower+step*g;        % quantile 
    [tempboot(1:2*(d-1),g), ~,tempchi_bootd(b,g)] = myfun_hetero(tau,m,0,1,Zz(:,2:end),Zz(:,1),l);
end
par_boot(:,b,:) = tempboot(1:(d-1),:);
end
%% compute the optimal quantile index and its corresponding estimator of delta under heteroskedasticity
disdbias(:,1) = abs(median(tempchi_bootd) - median_d1);  %compute a proxy of bias, which is the difference between the median of subsample minimal distance and the median of chisquare((J-1)*(d-1))

disdvar = zeros(G,1); % disdvar stores the variance of the estimator for estimator of delta without imputing zero
disbdvar = disdvar;  % disbdvar stores the variance of the estimator of delta when imputing zero
disbbvar = disdvar; % disbbvar stores the variance of the estimator of beta when imputing zero
for g = 1:G
    tau = lower+step*g;    
    tempboot1 = par_boot(:,:,g); 
    disdvar(g,1) = mean(sum((tempboot1 - repmat(mean(tempboot1,2),1,B)).^2)); % compute variance 
    disd(g,1) = disdvar(g,1)*boots/N + disdbias(g,1)/sqrt((tau*boots)); % compute the value of the object function for each grid 
end

[A,gstard] = min(disd); % gstard is the optimal grid who minimizes the object function disd 
delta_hetero = par_hetero(1:d-1,gstard); 
beta_hetero = par_hetero(d:2*(d-1),gstard);
V_hetero_star = V_hetero(:,:,gstard);
Nb_hetero_star = nb_hetero(gstard,1);
%% pretest delta = 0 using critical value C_n = sqrt(log(n))
phi = (abs(delta_hetero)<log(N)^(1/2)*sqrt(diag(V_hetero_star)/(N*(lower+step*gstard))))';
%% estimate the parital homoskedasitic model
dd = d - 1 - sum(phi); % dd is the number of nonzero delta after the test
db = sum(phi);  % db is the number of zero deltas after the test (i.e. number of beta we can compute)
median_d2 = median(random('chi2',(J-1)*(d-1-sum(phi)),100000,1)); % median of chisquare with degree of freedom = (J-1)*(d-1-sum(phi))
median_b1 = median(random('chi2',J*sum(phi),100000,1)); % median of chisquare with degree of freedom = J*sum(phi)

Nb_hom = zeros(G,1); % store the convergence rate of beta estimated using full sample for each quantile when zero is imputed.  
V_homb = zeros(db,db,G); % store the convergence rate of beta estimated using full sample for each quantile when zero is imputed.
V_homd = zeros(dd,dd,G); % store the convergence rate of delta estimated using full sample for each quantile when zero is imputed.

for g = 1:G
warning off 
[temphom, V_homd(:,:,g), V_homb(:,:,g),Nb_hom(g,1),~,chibb(g,1)] = myfun_hom(lower+step*g,m,phi,X,Y,l);
par_hom(:,g) = temphom;
end

parfor b = 1:B
Zz = Z(randperm(N,boots),:);      
tempboot = zeros(d-1,G);
for g = 1:G
    warning off 
    tau = lower+step*g;
    [tempboot(:,g), ~, ~, ~,tempchi_bootbd(b,g),tempchi_bootbb(b,g)] = myfun_hom(tau,m,phi,Zz(:,2:end),Zz(:,1),l);
end
par_bootb(:,b,:) = tempboot;
end

disbdbias(:,1) = abs(median(tempchi_bootbd) - median_d2);  
disbbbias(:,1) = abs(median(tempchi_bootbb) - median_b1);
% compute the optimal quantile index for the partial homoskedastic model

for g = 1:G
    tau = lower+step*g;    
    tempboot2 = par_bootb(:,:,g);
    tempbootd = tempboot2(1:dd,:);
    tempbootb = tempboot2(dd+1:end,:); 
    disbdvar(g,1) = mean(sum((tempbootd - repmat(mean(tempbootd,2),1,B)).^2));
    disbbvar(g,1) = mean(sum((tempbootb - repmat(mean(tempbootb,2),1,B)).^2));
    disbd(g,1) = disbdvar(g,1)*boots/N+disbdbias(g,1)/sqrt((tau*boots));
    disbb(g,1) = disbbvar(g,1)*boots/N+disbbbias(g,1)/sqrt((tau*boots));
end
[A,gstarbd] = min(disbd);
delta_hom = par_hom(1:dd,gstarbd);
[A,gstarbb] = min(disbb);
beta_hom = par_hom(dd+1:end,gstarbb);
Nb_hom_star = Nb_hom(gstarbb,1);
V_homd_star = V_homd(:,:,gstarbd);
V_homb_star = V_homb(:,:,gstarbb);
%% CLR bound for qth-QTE of the heteroskedastic variables

if sum(phi)<d-1 % if there is nonzero delta, we do CLR bound to compute qth-QTE, 

M = size(grid,2);
nor = randn(d-1,R);


[y1s,y2s] = bound(Y,Xdnp,Xcnp,quant,gridnp); % compute the upper and lowre bound, see equation 2.3 in page 8 of the paper. corresponding to the two Q funcitons 

row = (y1s~=-1000).*(y2s~=1000).*(y1s~=100).*(1:M)'; % get rid of the bounds which we do not want to use
row = row(row~=0);
y1 = y1s(row,:);
y2 = y2s(row,:);
Xgrid = grid;
Xgrid = Xgrid(:,row);

M = length(row);
Zu_ss = zeros(M,R);   % Z_n^* for the upper bound where Z_n^* is defined in Appendix A.1
Zl_ss = zeros(M,R);   % Z_n^* for the lower bound where Z_n^* is defined in Appendix A.1
thetau = zeros(M,1);  % theta evaluated at each grid for the upper bound, see Appendix A.1
thetal = zeros(M,1);  % theta evaluated at each grid for the lower bound, see Appendix A.1
Su = zeros(M,1);   % s_n evaluated at each grid for the upper bound, see the first equation on page 38 for the definition of s_n(x)
Sl = zeros(M,1);   % s_n evaluated at each grid for the lower bound
ll = find(phi==0);
theta0 = zeros(d-1-sum(phi),2);
for j = 1:(d-1-sum(phi))
    e1 =zeros(d-1,1);
    e1(ll(j)) = 1;

parfor m = 1:M
   
    if e1'*delta_hetero > 0 && 1+Xgrid(:,m)'*delta_hetero>0.05
    
    temptheta1 = beta_hetero + delta_hetero*(y2(m) - Xgrid(:,m)'*beta_hetero -b0)/(1+Xgrid(:,m)'*delta_hetero);  % upper bound
    temptheta1 = e1'*temptheta1;
    temptheta2 = -(beta_hetero + delta_hetero*(y1(m) - Xgrid(:,m)'*beta_hetero-b0)/(1+Xgrid(:,m)'*delta_hetero)); % - lower bound
    temptheta2 = e1'*temptheta2;
    elseif  e1'*delta_hetero < 0 && 1+Xgrid(:,m)'*delta_hetero>0.05
    temptheta1 = beta_hetero + delta_hetero*(y1(m) - Xgrid(:,m)'*beta_hetero -b0)/(1+Xgrid(:,m)'*delta_hetero); 
    temptheta1 = e1'*temptheta1;
    temptheta2 = -(beta_hetero + delta_hetero*(y2(m) - Xgrid(:,m)'*beta_hetero-b0)/(1+Xgrid(:,m)'*delta_hetero)); 
    temptheta2 = e1'*temptheta2;    
    else
    temptheta1 = 1000;
    temptheta2 = 1000;
    end
    temp = sign(e1'*delta_hetero)*(e1-e1'*delta_hetero*Xgrid(:,m)/abs(1+Xgrid(:,m)'*delta_hetero));    
    temps = sqrt(sum((temp'*sqrtm(V_hetero_star)).^2))/abs(Nb_hetero_star);
        
 
    for r = 1:R
    tempzss = temp'*sqrtm(V_hetero_star)*nor(:,r)/(sqrt(sum((temp'*sqrtm(V_hetero_star)).^2)));
    Zu_ss(m,r) = tempzss;
    Zl_ss(m,r) = -tempzss;
    end
    thetau(m,1) = temptheta1;
    Su(m,1)=temps;
    thetal(m,1) = temptheta2;
    Sl(m,1)= temps;
end
rn = 1-0.1/log(N);
Zussmax = max(Zu_ss);
knv = quantile(Zussmax,rn);

row2 = (thetau < min(thetau + knv*Su)+2*knv*Su); %compute the contact set \hat{\mathcal{X}}_n in appendix A.1 for the upper bound

row2 = row2.*(1:length(row2))';
row2 = row2(row2~=0);
Zu_ss2 = Zu_ss(row2,:);
Zuss2max = max(Zu_ss2);
knv2 = quantile(Zuss2max,0.975); % K_{2n}(0.975) in the second last equation of appendix A.1 for the upper bound
theta0(j,2) = min(thetau + knv2*Su);

Zlssmax = max(Zl_ss);
knv = quantile(Zlssmax,rn);

row2 = (thetal < min(thetal + knv*Sl)+2*knv*Sl);  %compute the contact set \hat{\mathcal{X}}_n in appendix A.1 for the lower bound

row2 = row2.*(1:length(row2))';
row2 = row2(row2~=0);
Zl_ss2 = Zl_ss(row2,:);
Zlss2max = max(Zl_ss2);
knv2 = quantile(Zlss2max,0.975); % K_{2n}(0.975) in the second last equation of appendix A.1. 97.5% quantile of  -X = - 2.5% quantile of X 
theta0(j,1) = -min(thetal + knv2*Sl);
end
        
end
qstard = lower+step*gstard;
qstarbd = lower+step*gstarbd;
qstarbb = lower+step*gstarbb;

end

