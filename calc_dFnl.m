
function [dFX dFdX dFddX] = calc_dFnl(X,dX,ddX,lambdapp)
%%les constantes du problème
global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K

%%
%calcul de la jacobienne 
dFX(1,1) = 0;
dFX(1,2) = -l*alpha*ddX(2)*cos(X(2)) + l*alpha*sin(X(2))*(dX(2).^2);
dFX(2,1) = 0;
dFX(2,2) = -BETA*cos(X(2))*(ddX(1)+lambdapp) + (omega2^2)*l*cos(X(2))  ;



dFdX(1,1) = 0;
dFdX(1,2) = -2*l*alpha*cos(X(2))*dX(2);
dFdX(2,1) =0;
dFdX(2,2) =0;



dFddX(1,1) = 0;
dFddX(1,2) = -l*sin(X(2))*alpha;
dFddX(2,1) = -BETA*sin(X(2));
dFddX(2,2) = 0;


end
