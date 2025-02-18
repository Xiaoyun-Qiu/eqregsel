function [par, Vd, Vb, Nb, dis_d,dis_b] = myfun_hom(tau,m,phi,X,Y,l,delta_MSE)
%%% Input
% myfun_hom computes the homoskedastic beta and the heteroskedastic delta 
% tau: quantile 
% m: a tuning parameter that is used to normalize our estimator. should be different from but very close to 1. In simulation and applications, m is set to be 1.2 
% X: dependent variables
% Y: independent variables
% l: equations we will explore by minimum distance estimations are indexed by l. In general, we will use quantile level tau, tau*l to estimate beta and delta 
% phi: a 1 by d vector contains the results of pretest where d is the dimension of X (exclude intercept). If j-th entry is 1, then it indicates the j-th variable is homoskedastic. 
% delta_MSE: the user supplied starting value of delta. If it is left unspecified, the program will use OLS estimator
%%% Output
% par: a (d-1) by 1 vector collects estimators of heteroskedastic delta and estimators of homoskedastic beta, where d is the dimension of x (including intercept).
% Vd: The asymptotic variance for delta 
% Vb: The asymptotic variance for beta 
% dis_d: minimum distance for delta evaluated by plugging in the extremal estimator of delta under partial homoskedasticity 
% dis_b: minimum distance for beta evaluated by plugging in the extremal estimator of beta under partial homoskedasticity

[T,d]= size(X);
par = zeros(d-1,1);
J = length(l);
db = sum(phi);
psi = ones(1,d-1) - phi;
dd = sum(psi);
% based on the user supplied vector of phi, we divide our estimation
% procedure to two parts, one is the minimum distance estimation of homo beta and the
% other is the minimum distance estimation of hetero delta. phitemp and
% psitemp are created to separate variable to homo part and hetero part.
phitemp = zeros(db,d-1);
psitemp = zeros(dd,d-1);
count = 1;
for i = 1:db
    for j = count:d-1
        if phi(j) == 1
            phitemp(i,j) = 1;
            count = j +1;
            break;
        end
    end
end

count = 1;
for i = 1:dd
    for j = count:d-1
        if psi(j) == 1
            psitemp(i,j) = 1;
            count = j +1;
            break;
        end
    end
end

l1 = [1,l];
L = zeros(J+1,J+1);
for i = 1:J+1
    for j = 1:J+1
        L(i,j) = min(l1(i),l1(j))/sqrt(l1(i)*l1(j));
    end
end


Gamma = [-repmat(eye(d),J,1),kron(diag(1./sqrt(l)),eye(d))]; % matrix Gamma is defined above theorem 3.1 in section 3.2
Gamma2 = [zeros(d-1,1),eye(d-1)]; % matrix Gamma2 is defined above theorem 3.2 in section 3.2 
Gamma3 = diag(1./sqrt(l1)); % matrix Gamma3 is defined above theorem 3.2 in section 3.2
if dd > 0
tempy1 = zeros(dd*J,1);
tempy2 = zeros(db*(J+1),1);
Gd = kron(log(l'),eye(dd));
Gb = -kron(ones(J+1,1),eye(db));
Qx = X'*X/T;
pf = zeros(d,J+1);
pm = zeros(d,J+1);
pf(:,1) = rq_fnm(X, -Y, tau);
pm(:,1) = rq_fnm(X,-Y,m*tau);
tempy2(1:db,1) = phitemp*pf(2:end,1);  
tempx1 = zeros(J,1);
for j = 1:J
       pf(:,j+1) = rq_fnm(X, -Y, tau*l(j));       
       tempy1(dd*(j-1)+1:dd*j,1) = psitemp*(pf(2:end,j+1)-pf(2:end,1)); 
       tempy2(db*j+1:db*(j+1),1) = phitemp*pf(2:end,j+1);   
       tempx1(j) = pf(1,j+1)-pf(1,1);
       pm(:,j+1) = rq_fnm(X, -Y, tau*m*l(j));
end
    
tempx = kron(tempx1,eye(dd));

if nargin == 6
     dtemp = (tempx'*tempx)\tempx'*tempy1;
     dtemp2 = dtemp'*psitemp;
else
    dtemp = delta_MSE;
    dtemp2 = dtemp'*psitemp;
end
for r = 1:2
    Tri = [-dtemp2',eye(d-1)];
    Xd = X./((repmat(abs(X*[1,dtemp2]'),1,d)).^(0.5));
    tempXd = Xd(:,end);
    % Trimming is applied to guard against very large Xd. For detail,
    % pleasee see the description of myfun_mse_new
    idx = (abs(tempXd-median(tempXd))<3*(quantile(tempXd,0.75) - quantile(tempXd,0.25)));
    Xd = Xd(idx,:);    
    QH = Xd'*Xd/length(idx);
    omega_0 = QH^(-1)*Qx*QH^(-1);
    W1 = (kron(eye(J),psitemp*Tri)*Gamma*kron(L,omega_0)*Gamma'*kron(eye(J),Tri'*psitemp'))^(-1); % W1 is the optimal weighting matrix for hetero delta
    delta = (tempx'*W1*tempx)\tempx'*W1*tempy1;
    dtemp = delta;
    dtemp2 = dtemp'*psitemp;
end

Vd = (Gd'*W1*Gd)^(-1);
W2 = (kron(Gamma3,phitemp*Gamma2)*kron(L,omega_0)*kron(Gamma3,phitemp*Gamma2)')^(-1); % W2 is the optimal weighting matrix for homo beta
beta = -(repmat(eye(db),J+1,1))'*W2*(repmat(eye(db),J+1,1))\(repmat(eye(db),J+1,1))'*W2*tempy2;
Vb = (Gb'*W2*Gb)^(-1);
for j = 1:J
mom1(dd*(j-1)+1:dd*j) = psitemp*(pf(2:end,j+1)-pf(2:end,1)) - (pf(1,j+1)-pf(1,1))*delta; %mom1 is the value of moments for hetero delta
mom2(db*(j-1)+1:db*j) = phitemp*pf(2:end,j) + beta; %mom2 is the value of moments for homo beta
end
mom2(db*J+1:db*(J+1)) = phitemp*pf(2:end,J+1) + beta;
Nb = log(m)*sqrt(tau*T)/median(abs(pm(1,:) - pf(1,:))); %Nb is the normalizing factor for homo beta, the normalizing factor for hetero delta is \sqrt{N*tau}

dis_d = Nb^2*mom1*W1*mom1'; % here we compute the distance of hetero delta evaluated at extremal quantile estimator of delta
dis_b = Nb^2*mom2*W2*mom2'; % here we compute the distance of homo beta evaluated at extremal quantile estimator of delta
par(1:dd,1) = delta;
par((dd+1):(d-1),1) = beta;
elseif dd == 0
   
    Gb = -kron(ones(J+1,1),eye(db));
    Qx = X'*X/T;
    pf = zeros(d,J+1);
    pm = zeros(d,J+1);
    pf(:,1) = rq_fnm(X, -Y, tau);
    pm(:,1) = rq_fnm(X,-Y,m*tau);
tempy2(1:db,1) = phitemp*pf(2:end,1);  
tempx1 = zeros(J,1);
for j = 1:J
       pf(:,j+1) = rq_fnm(X, -Y, tau*l(j));       
       tempy2(db*j+1:db*(j+1),1) = phitemp*pf(2:end,j+1);   
       tempx1(j) = pf(1,j+1)-pf(1,1);
       pm(:,j+1) = rq_fnm(X, -Y, tau*m*l(j));
end
    

    W2 = (kron(Gamma3,phitemp*Gamma2)*kron(L,Qx^(-1))*kron(Gamma3,phitemp*Gamma2)')^(-1); 
    beta = -(repmat(eye(db),J+1,1))'*W2*(repmat(eye(db),J+1,1))\(repmat(eye(db),J+1,1))'*W2*tempy2;
    Vb = (Gb'*W2*Gb)^(-1);
    for j = 1:J
    mom2(db*(j-1)+1:db*j) = phitemp*pf(2:end,j) + beta; %mom2 is the value of moments for homo beta
    end
mom2(db*J+1:db*(J+1)) = phitemp*pf(2:end,J+1) + beta;
Nb = log(m)*sqrt(tau*T)/median(abs(pm(1,:) - pf(1,:))); %Nb is the normalizing factor for homo beta, the normalizing factor for hetero delta is \sqrt{N*tau}

dis_d = 0;
dis_b = Nb^2*mom2*W2*mom2'; % here we compute the distance of homo beta evaluated at extremal quantile estimator of delta

par((dd+1):(d-1),1) = beta;
Vd = 0;
end    
    
