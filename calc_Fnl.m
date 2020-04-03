function Fnl = calc_Fnl(X,dX,ddX,lambdapp)
global k m Ma g l Io omega1 omega2 omega eps1 eps2 alpha beta lamda0
global M C K %matrices-coefficients du formalisme de newmark 

Fnl = zeros(size(X));

Fnl = [lambdapp-l*alpha*sin(X(2))*ddX(2)-l*alpha*cos(X(2))*dX(2)*dX(2)  ;     -(beta*sin(X(2)))*(ddX(1)+lambdapp)+omega2^2*l*sin(X(2))];


end 
