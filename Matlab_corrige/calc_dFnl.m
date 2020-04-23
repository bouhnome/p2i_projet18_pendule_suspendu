function [dFX dFdX] = calc_dFnl(X,dX,ddX,t)

global l a OMEGA OMEGA2 beta

lambdapp=-a*OMEGA^2*sin(OMEGA*t);
dFX=zeros(length(X));
dFdX=zeros(length(X));
dFddX = zeros(length(X));

dFX(1,1) = 0;
dFX(1,2) = -l*alpha*ddX(2)*cos(X(2)) + l*alpha*sin(X(2))*(dX(2).^2);
dFX(2,1) = 0;
dFX(2,2) = -beta*cos(X(2))*ddX(1) + (OMEGA2^2)*l*cos(x(2))  ;  %%% plus terme en lambda seconde

dFdX(1,1) = 0;
dFdX(1,2) = -2*l*alpha*cos(X(2))*dX(2);
dFdX(2,1) =0;
dFdX(2,2) =0;


dFddX(1,1) = 0;
dFddX(1,2) = -l*sin(X(2))*alpha;
dFddX(2,1) = -beta*sin(X(2));
dFddX(2,2) = 0;

end

