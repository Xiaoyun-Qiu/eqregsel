function [out,V,dis,Nb1] = myfun_MSE_new(tau,m,b0,d0,X,Y,l,delta_MSE)
% myfun_hetero computes beta and delta under the normalization that b_0 =
% 0 and d_0 = 1 without assuming partial homoskedasticity
%%%input
% tau: quantile 
% m: a tuning parameter that is used to normalize our estimator. should be different from but very close to 1. In simulation and applications, m is set to be 1.2 
% b0: the location normalizing factor
% d0: the scale normalizing factor
% X: dependent variables
% Y: independent variables
% l: equations we will explore by minimum distance estimations are indexed by l. In general, we will use quantile level tau, tau*l to estimate beta and delta 
% delta_hetero: the user supplied starting value of delta. If it is left unspecified, the program will use OLS estimator
%%% output
% out: 2*(d-1) by 1 vector collects (d-1) estimator of delta and (d-1) estimator of beta, where d is the dimension of x (including intercept).
% V: The asymptotic variance for delta
% dis: minimum distance of delta computed by plugging in the extremal estimator of delta
% Nb1: the convergence rate for beta
[T,d]= size(X);
J = length(l);
pf = zeros(d,J+1);

gamma = [-repmat(eye(d),J,1),kron(diag(1./sqrt(l)),eye(d))]; % Matrix Gamma is defined in section 3.2 of the paper
l1 = [1,l];
L = zeros(J+1,J+1); % Matrix L is defined in section 3.2 of the paper
for i = 1:J+1
    for j = 1:J+1
        L(i,j) = min(l1(i),l1(j))/sqrt(l1(i)*l1(j));
    end
end

G = kron(log(l'),eye(d-1)); % Matrix G is defined in section 3.2, above theorem 3.1
Qx = X'*X/T; 

tempy = zeros((d-1)*J,1);
tempx1 = zeros(J,1);

pf(:,1) = rq_fnm(X, -Y, tau);
pm = rq_fnm(X,-Y,m*tau);
pm(:,1) = rq_fnm(X,-Y,m*tau);

for j = 1:J
       pf(:,j+1) = rq_fnm(X, -Y, tau*l(j));       
       tempy((d-1)*(j-1)+1:(d-1)*j) = (pf(2:end,j+1)-pf(2:end,1))*d0; 
       tempx1(j) = pf(1,j+1)-pf(1,1);
       pm(:,j+1) = rq_fnm(X, -Y, m*tau*l(j));
end
pf(1,:) = pf(1,:) + repmat(b0,1,J+1);
tempx = kron(tempx1,eye(d-1));

if nargin == 7
    dtemp = (tempx'*tempx)\tempx'*tempy;
else
    dtemp = delta_MSE;
end
for r = 1:2 % we do a feasible two step estimator. But instead of stopping at one iteration, we do 3 iterations. So this can be also viewed as a simple version of CUE. 
    Tri = [-dtemp/d0,eye(d-1)];
    Xd = X./((repmat(abs(X*[d0,dtemp']'),1,d)).^(0.5));
    % X*delta is assumed to be positive, but it may not be the case when we
    % compute it with an estimator of delta, when it is close to zero, Xd
    % is very large and this will cause numerical problems. In order to
    % deal with this, we trim out the extremely large value of Xd. This
    % will provide robustness to the program.
    tempXd = Xd*[d0;dtemp];
    idx = (abs(tempXd-median(tempXd))<3*(quantile(tempXd,0.75) - quantile(tempXd,0.25)));
    % we trim out Xd whose distance to the median is 3 times larger than
    % the interquartile range.
    Xd = Xd(idx,:);    
    QH = Xd'*Xd/length(idx);
    omega_0 = QH^(-1)*Qx*QH^(-1);
    W = (kron(eye(J),Tri)*gamma*kron(L,omega_0)*gamma'*kron(eye(J),Tri'))^(-1);
    dtemp = (tempx'*W*tempx)\tempx'*W*tempy; 
end
   delta = dtemp;
    out(1:d-1,1) = delta;
    V = (G'*W*G)^(-1);
    % mom is the value of each moment evaluated at extremal quantile estimator of delta 
for j = 1:J
mom((d-1)*(j-1)+1:(d-1)*j) = (pf(2:end,j+1)-pf(2:end,1)) - (pf(1,j+1)-pf(1,1))*delta/d0;
end
Nb = log(m)*sqrt(tau*T)/median(abs(pm(1,:) - pf(1,:)));
dis = Nb^2*mom*W*mom'; % this is the value of distance 
beta = - pf(2:end,1) + delta/d0*pf(1,1);
out(d:2*d-2,1) = beta; 
Nb1 = abs(sqrt(tau*T)/median((pf(1,:))));
end
