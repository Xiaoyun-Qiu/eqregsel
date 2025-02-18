function [beta_hom,std_b,specificationtest,tau0] = myfun_combined(X,Y,boots,G,B)
%%% Input
% myfun_combined computes the homoskedastic beta and the J-test statistic.  
% X: Dependent variables
% Y: Independent variable
% phi: User supplied vector to indicate which variable (in our application, it is "Black") is assumed to be homoskedastic. (The code allows for more than one covariates to be homoskedastic.)
% boots: Number of subsample size. In simulation, we use (150,300,500,700) for sample size (300,500,1000,2000).
%%% Output
% beta_hom: The point estimator corresponding to the optimal quantile index tau_{n,2} selected based on the procedure described in the paper. 
% std_b: Standard deviation for the point estimator. 
% specificationtest: P-value for the J-test.
% tau0: Optimal quantile index selected. 
% B: The number of replications in our subsampling and bootstrap method.  
% G: The number of grids.
rng(14,'combRecursive')
[N,d] = size(X);
if nargin == 2        
    boots = 0.6*N*(N <=500) + (300 + 0.4*(N-500))*(N>500)*(N<=1000)+ (500 + 0.2*(N-1000))*(N>1000)*(N<=2000) + (700+0.1*(N-2000))*(N>2000);boots = floor(boots);
    G = 40; B = 150; 
elseif nargin == 3
    G = 40; B = 150; 
    elseif   nargin == 4
B = 150;
end
            
phi = zeros(d-1,1);phi(1)=1;
l = [0.9,1.1]; % l is corresponding to the moment we use. 
J = length(l); 
lower = min(80/boots, 0.1); % lower is the lower bound for tau_{n,0}'.
upper = 0.3; % upper is the upper bound for tau_{n,0}'.
step = (upper-lower)/G; 

db = sum(phi);

%%
disbb = zeros(G,1);disbbbias = zeros(G,1);


options = statset('UseParallel','always') ;

median_b1 = median(random('chi2',(J-1)*sum(phi),100000,1)); %Compute the median of the chi-square random variable with (J*sum(phi)) degree of freedom.

%%
tic
Z = [Y,X];
chi_bootbb = zeros(B,G);
par_bootb = zeros(db,B,G);
disbbvar = zeros(G,1);
for g=1:G
warning off
thetahat = rq_fnm(X, -Y, lower+step*g);
thetadagger = zeros(d,B);
tau = lower+step*g;
for s=1:B  % Bootstrap the first stage estimator, use it to compute omega_0.  
    Zz = Z(randsample(N,N,true),:);       
    thetadagger(:,s) = rq_fnm(Zz(:,2:end),-Zz(:,1), tau);
end
Sigma = (thetadagger - repmat(thetahat,1,B))*(thetadagger - repmat(thetahat,1,B))'/B; %Compute the optimal weighting matrix.

for s = 1:B %Use subsample to compute the point estimator and J-test statistic. 
Zz = Z(randperm(N,boots),:);      

pf = zeros(d,J);
for j = 1:J
       pf(:,j) = rq_fnm(Zz(:,2:end), -Zz(:,1), tau*l(j));       
end
[par_bootb(:,s,g),chi_bootbb(s,g)] = myfun_hom_new(d,phi,pf,l,Sigma);
end

end


%%
for g=1:G
tau = lower+step*g;
disbbbias(g,1) = abs(median(chi_bootbb(:,g))*boots/N - median_b1); %Compute an proxy of bias.

disbbvar(g,1) = mean(sum((par_bootb(:,:,g) - repmat(mean(par_bootb(:,:,g),2),1,B)).^2,1)); %Compute the variance.
disbb(g,1) = disbbvar(g,1)*boots/N+disbbbias(g,1)/sqrt(tau*boots); %Compute a proxy of MSE for the full sample with size N. 
end
[A,gstarbb] = min(disbb); %Select tau.


%% bootstrap the confidence interval
par_bootstraphomb = zeros(db,B);
tau = (lower+step*gstarbb);%Compute tau_{n,2}
thetahat = rq_fnm(X, -Y, tau);
thetadagger = zeros(d,B);

l1 = [1;0.2];J1 = length(l1);
parfor s=1:B  % Bootstrap the first stage estimator, use it to compute omega_0.   
    Zz = Z(randsample(N,N,true),:);       
    thetadagger(:,s) = rq_fnm(Zz(:,2:end),-Zz(:,1), tau);
end
    pf = zeros(d,J1);

for j = 1:J1
       pf(:,j) = rq_fnm(X, -Y, tau*l1(j));       
end

Sigmahat = (thetadagger - repmat(thetahat,1,B))*(thetadagger - repmat(thetahat,1,B))'/B; %Compute the optimal weighting matrix
[beta_hom,chibb] = myfun_hom_new(d,phi,pf,l1,Sigmahat); % Point estimator and the J-test statistic with tau_{n,2} 

parfor b = 1:B  %Bootstrap with tau_{n,2}, use the point estimator to compute the standard deviation.   
    Zz = Z(randsample(N,N,true),:); 
    pf = zeros(d,J1);
for j = 1:J1
       pf(:,j) = rq_fnm(Zz(:,2:end), -Zz(:,1), tau*l1(j));       
end

     [par_bootstraphomb(:,b),~] = myfun_hom_new(d,phi,pf,l1,Sigmahat);
end

%%
std_b = std(par_bootstraphomb,[],2); %Compute the standard deviation.

specificationtest = 1-cdf('chi2',chibb,(J1-1)*sum(phi));% Compute the p-value for the J-test. 
tau0 = tau; %Output tau_{n,2}
end
%%






