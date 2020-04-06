function [dFX dFdX] = calc_dFnl(X,dX,t)

global l0 a omega1 omega2 epsilon1 epsilon2 omega;

dFX=zeros(length(X));
dFdX=zeros(length(X));
lambdapp=-a*omega^2*sin(omega*t);

dFX(1,1)=-dX(2)^2+omega1^2;
dFX(1,2)=omega2^2*sin(X(2))-sin(X(2))*lambdapp/l0;
dFX(2,1)=-(2*dX(1)*dX(2)+2*epsilon2*omega2*dX(2)+omega2^2*sin(X(2))-sin(X(2))*lambdapp/l0)/(X(1)^2);
dFX(2,2)=(omega2^2*cos(X(2))-cos(X(2))*lambdapp/l0)/X(1);

dFdX(1,1)=2*epsilon1*omega1;
dFdX(1,2)=-2*X(1)*dX(2);
dFdX(2,1)=2*dX(2)/X(1);
dFdX(2,2)=(2*dX(1)+2*epsilon2*omega2)/X(1);

end

