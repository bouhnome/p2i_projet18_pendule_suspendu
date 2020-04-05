function sol = Pendule_Function_L(t)
%Resolution pour la seconde equa diff
eps1 = 0.005;
eps2 = 0.005;
k = 0.2;
m = 0.5;
Ma = 5;
g =9.81 ;
l = 5;
Io = 4*m*(l^2);
w1 = (k/(m+Ma))^0.5;
w2 = sqrt(m*g*l/Io);
a = 0.02;
thetap_0 = 0.1;


a2 = 1;
b2 = 2*eps2*w2;
c2 = w2^2;

delta2 = b2^2 - 4*a2*c2;

X2_1 = (-b2 - 1i*sqrt(-delta2))/(2*a2);
X2_2 = conj(X2_1);

A2 = 0;
B2 = thetap_0/imag(X2_1);

% sol 2 = theta et sol 4 = thetapoint
sol2 = exp(-real(X2_1)*t)*(A2*cos(imag(X2_1)*t)+B2*sin(imag(X2_2)*t));
sol4 = -real(X2_1)*imag(X2_1)*exp(-real(X2_1)*t)*B2*cos(imag(X2_1)*t);

%Resolution pour la premiere equa-diff
a1 = 1;
b1 = 2*eps1*w1;
c1 = w1^2;

delta1 = b1^2 - 4*a1*c1 ;

X1_1 = (-b1 - 1i*sqrt(-delta1))/(2*a1);
X1_2 = conj(X1_1);

w = (w1 + w2)/2;
zsp = ((a*w^2)/((w1^2 - w^2)^2 + 4*(eps1*w1*w)^2)) * ((w1^2 - w2^2)*sin(w*t) + 2*eps1*w1*w*cos(w*t));

A1 = -zsp;
B1 = -A1/imag(X1_1);
% sol1 = z et sol3 = zpoint
sol1 = exp(-real(X1_1)*t)*(A1*cos(imag(X1_1)*t)+B1*sin(imag(X1_2)*t));
sol3 = -real(X1_1)*exp(-real(X1_1)*t)*(A1*cos(imag(X1_1)*t)+B1*sin(imag(X1_2)*t)) + exp(-real(X1_1)*t)*(-A1*imag(X1_1)*sin(imag(X1_1)*t)+B1*imag(X1_2)*cos(imag(X1_2)*t));

sol = [sol1 sol2 sol3 sol4];

end