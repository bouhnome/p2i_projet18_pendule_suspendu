function Fnl = calc_Fnl(X,dX,ddX,lambdapp)
%%les constantes du problème
global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K
%%

Fnl = [lambdapp-l*alpha*sin(X(2))*ddX(2)-l*alpha*cos(X(2))*dX(2)*dX(2);-(BETA*sin(X(2)))*(ddX(1)+lambdapp)+omega2^2*l*sin(X(2))];


end 
