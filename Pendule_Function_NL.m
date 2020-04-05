function dX = Pendule_Function_NL(t,X)
global k m Ma g l Io omega1 omega2 omega eps1 eps2 alpha beta a lambdapp gamma

dX = zeros (4,1);
lambdapp = -a*(omega^2)*sin(omega*t);
gamma = 1 - alpha*beta*sin(X(3));

dX(1) = X(2);
dX(2) = lambdapp + ( (1*alpha*cos(X(3))*(X(4)^2)) - (2*l*alpha*eps2*omega2*sin(X(3))*X(4)) - (2*eps1*omega1*X(2)) - ((omega1^2)*X(1)) - ((1 - gamma)*l*(omega2^2)/beta) )/gamma;
dX(3) = X(4);
dX(4) = ( (alpha*beta*sin(X(3))*cos(X(3))*(X(4)^2)) - (2*eps2*omega2*X(4)) - ( ( (2*beta*eps1*omega1*sin(X(3))*X(2)) + (beta*(omega1^2)*sin(X(3))*X(1)) )/l ) - ((omega2^2)*sin(X(3))) )/gamma;
end