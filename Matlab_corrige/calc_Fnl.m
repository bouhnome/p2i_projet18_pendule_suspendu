function Fnl = calc_Fnl(X,dX,ddX,t)

%%%add global variables

lambdapp = -a*(OMEGA^2)*sin(OMEGA*t)

Fnl=zeros(size(X));
Fnl(1)= lambdapp - l*alpha*sin(X(2))*ddX(2) - (l*alpha*cos(X(2))*ddX(2);
Fnl(2) = -(beta*sin(X(2))*(ddX(1)+ lambdapp) + OMEGA2^2*l*sin(X(2));

end

