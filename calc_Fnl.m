function Fnl = calc_Fnl(X,dX,t)

global l0 a omega1 omega2 eps1 eps2 omega;
Fnl=zeros(size(X));
lambdapp=-a*omega^2*sin(omega*t);
%
Fnl(1)=-X(1)*dX(2)^2+2*eps1*omega1*dX(1)+omega1^2*(X(1)-1)-omega2^2*cos(X(2))+cos(X(2))*lambdapp/l0;
Fnl(2)=(2*dX(1)*dX(2)+2*eps2*omega2*dX(2)+omega2^2*sin(X(2))-sin(X(2))*lambdapp/l0)/X(1);

end

