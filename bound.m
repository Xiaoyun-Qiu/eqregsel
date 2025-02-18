function [y1,y2] = bound(Y,Xdnp,Xcnp,tau,gridnp)
%bound compute the upper and lower bounds of the QTE 
%%% Input
%  Y: the independent variable
%  Xdnp: a matrix of free discrete dependent variables (free means it is not determined by other components of X = [Xd,Xc]) 
%  Xcnp: a matrix of free continuous dependent variables
%  tau: the quantile of interest of the CLR bound
%  gridnp: the grid of free X = [Xdnp,Xcnp] used in nonparametric estimation when computing CLR bound 
%%% Output
%  y1: the lower bound 
%  y2: the upper bound
D=(Y>0); 
h = 1.06*length(Y)^(-1/(4+size(Xcnp,2)));
m = size(gridnp,2);
y1 = zeros(m,1);
y2 = zeros(m,1);
[n,l] = size(Y);
dd = size(Xdnp,2);
parfor i = 1:m;
       tempx = gridnp(:,i);
       P1 = sum((D==1).*(sum(Xdnp==repmat(tempx(1:dd)',n,1),2)==3).*prod(normpdf(repmat(tempx(dd+1:end)',n,1)-Xcnp)./repmat(std(Xcnp),n,1),2))...
           /sum((sum(Xdnp==repmat(tempx(1:dd)',n,1),2)==3).*prod(normpdf(repmat(tempx(dd+1:end)',n,1)-Xcnp)./repmat(std(Xcnp),n,1),2));
       if (tau-1+P1)/P1 > 0
       temp1 = fzero(@(y) condQ(y,tempx(dd+1:end)',tempx(1:dd)',Y,D,Xcnp,Xdnp,h)-(tau-1+P1)/P1,quantile(Y,tau));
       else 
       temp1 = - 1000;
       end
       if tau/P1 < 1
       temp2 = fzero(@(y) condQ(y,tempx(dd+1:end)',tempx(1:dd)',Y,D,Xcnp,Xdnp,h)-tau/P1,quantile(Y,tau));
       else 
       temp2 = 1000;    
       end    
       y1(i,1) = temp1;
       y2(i,1) = temp2;        
end
       
end    




