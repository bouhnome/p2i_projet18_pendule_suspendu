function dxdt=Pendule_elastique(t,x)
%
global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K
%
lambdapp=-a*omega0^2*sin(omega0*t);
gammma=1-alpha*BETA*sin(x(3))^2;

%
dxdt=zeros(4,1);
%
dxdt(1)=x(2);
dxdt(2)=-lambdapp+l*alpha*cos(x(3))*x(4)^2/gammma-2*l*alpha*eps2*omega2*sin(x(3))*x(4)/gammma-2*eps1*omega1*x(2)/gammma-omega1^2*x(1)/gammma-(1-gammma)*l*omega2^2/(gammma*BETA);
dxdt(3)=x(4);
dxdt(4)=(alpha*BETA*sin(x(3))*cos(x(3))*x(4)^2-2*eps2*omega2*x(4)-2*BETA*eps1*omega1*sin(x(3))*x(2)/l-BETA*omega1^2*sin(x(3))*x(1)/l-omega2^2*sin(x(3)))/gammma;
%
end
