function q = condQ(y,xc,xd,Y,D,Xc,Xd,h)
%q3 is the kernel density estimator of the conditional probability P(Y<y|D=1,X=x)
%%% Input
% y: the value at which the conditional probability is evaluated
% xc: the continuous variable value at which the conditional probability is conditional on
% xd: the discrete variable value at which the conditional probability is conditional on
% Y: the independent variable
%Xc: the continuous elements of X
%Xd: the discrete elements of X
%h: the tunning parameter
%%% Output
% q: estimator of P(Y<y|D=1,X=x)

[n,m]=size(Xd);
xdtemp = repmat(xd,n,1);
xctemp = repmat(xc,n,1);
q = sum((Y<y).*(D==1).*(sum((Xd==xdtemp),2)==3).*prod(normpdf(xctemp-Xc/h),2))/sum((sum((Xd==xdtemp),2)==3).*prod(normpdf(xctemp-Xc/h),2).*(D==1));



end

