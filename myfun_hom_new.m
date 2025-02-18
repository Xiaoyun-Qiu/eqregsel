function [beta,dis_b] = myfun_hom_new(d,phi,pf,l,Sigma)
%%% Input
% myfun_hom_new computes the point estimator and the J-test statistic for a given quantile index tau.  
% tau: Quantile 
% phi: User supplied vector to indicate which variable (in our application, it is "Black") is assumed to be homoskedastic. (The code allows for more than one covariates to be homoskedastic.)
% X: Dependent variables
% Y: Independent variable
% l: Equations we will explore by minimum distance estimations are indexed by l. In general, we will use quantile level tau, tau*l to estimate beta and delta 
% Sigma: An estimator for \omega_0 in the paper
%%% Output
% beta: The point estimator
% dis_b: J-test statistic for this specification. 

% [T,d]= size(X);

J = length(l);
%% Convert vector phi into a matrix which picks out the covariates that are homoskedastic. 
db = sum(phi);
phitemp = zeros(db,d-1);
count = 1;
for i = 1:db
    counta = 1;
    for j = count:d-1
        if phi(j) == 1 && counta == 1
            phitemp(i,j) = 1;
            counta = counta+1;
            count = j +1;
        end
    end
end
%% Compute matrix L, Gamma_2, Gamma_3 in the paper.  
l1 = l; 
L = zeros(J,J);
for i = 1:J
    for j = 1:J
        L(i,j) = min(l1(i),l1(j))/sqrt(l1(i)*l1(j));
    end
end

gamma2 = [zeros(d-1,1),eye(d-1)];
gamma3 = diag(1./sqrt(l1));
%% The first step: extremal quantile regression
tempy2 = zeros(db*(J),1);
for j = 1:J
       tempy2(db*(j-1)+1:db*(j),1) = phitemp*pf(2:end,j);   
end
%% The second step: minimum distance estimation    
omega_0 = Sigma;
W2 = (kron(gamma3,phitemp*gamma2)*kron(L,omega_0)*kron(gamma3,phitemp*gamma2)')^(-1); % W2 is the optimal weighting matrix for homo beta
beta = -(repmat(eye(db),J,1))'*W2*(repmat(eye(db),J,1))\(repmat(eye(db),J,1))'*W2*tempy2;
for j = 1:J
mom2(db*(j-1)+1:db*j) = phitemp*pf(2:end,j) + beta; %mom2 is the value of moments for homo beta
end

dis_b = mom2*W2*mom2'; % here we compute the distance of homo beta evaluated at extremal quantile estimator of delta
 
    
