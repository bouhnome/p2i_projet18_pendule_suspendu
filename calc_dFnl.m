
function [dFX dFdX dFddX] = calc_dFnl(X,dX,ddX,lambdapp)
%%les constantes du problème
k = 0.2;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
m = 0.5;%masse de la tige(solide S2)
Ma = 5;%masse du solide S1
g =9.81 ;%valeur du champs de gravité
l = 5;%la moitie la longueure de la tige(solide S2)
Io = 4*m*(l^2)/3;%moment d'inertie de la tige par rapport à son axe de rotation
eps1 = 0.005;%facteur d'amortissement visqueux lié à omega1
eps2 = 0.005;%facteur d'amortissement visqueux li2 à omega2
a = 0;%amplitude de l'exitation d'entree 
omega1 = sqrt(k /(m+Ma));%pulsation du système masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de résonnance du pendule seul
omega0 = (omega1 + omega2)/2;%pulsation de l'exitation  
alpha = m/(m+Ma) ;%coefficient adimensionnel 
BETA = m*l^2/Io;%coefficient adimensionnel 


M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];

%%
%calcul de la jacobienne 
dFX(1,1) = omega1^2;
dFX(1,2) = -l*alpha*ddX(2)*cos(X(2)) + l*alpha*sin(X(2))*(dX(2).^2);
dFX(2,1) = 0;
dFX(2,2) = -BETA*cos(X(2))*(ddX(1)+lambdapp) + (omega2^2)*l*cos(X(2))  ;



dFdX(1,1) = 2*eps1*omega1;
dFdX(1,2) = -dX(2)*2*l*alpha*cos(X(2));
dFdX(2,1) =0;
dFdX(2,2) =2*eps2*omega2*l;



dFddX(1,1) = 1;
dFddX(1,2) = -l*sin(X(2))*alpha;
dFddX(2,1) = -BETA*sin(X(2));
dFddX(2,2) = l;


end
