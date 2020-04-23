% Forced Duffing Oscillator

function xdot=duffing(t,x)

global knl p0 OMEGA
global M C K

xdot(1)=x(2);
xdot(2)=-C*x(2)-K*x(1)-knl*x(1)^3+p0*cos(OMEGA*t);

xdot=xdot';

end
