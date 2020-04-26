function [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot)

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
tt=[t_init:dt:t_tot];%vecteur temps à retourner
precNR=1.e-9;
lambdapp=-a*omega0^2*sin(tt*omega0);
lambdapp=lambdapp';

%% sol init
t=t_init;n=1;
X=X0; dX=dX0;
%on détermine la valeur de l'accélération initiale
ddX=-inv([1,-l*alpha*sin(X(2,1));-BETA*sin(X(2,1)),l])*[-lambdapp(1,1)+(l*alpha*cos(X(2,1)))*dX(2,1)^2-2*eps1*omega1*dX(1,1)-omega1^2*X(1,1);BETA*sin(X(2,1))*lambdapp(1,1)-2*eps2*omega2*l*dX(2,1)-omega2^2*l*sin(X(2,1))];



%% integration temporelle
for t=t_init+dt:dt:t_tot;
    n=n+1;
    
    % prediction
    X=X+dt*dX+(dt^2/2)*ddX;
    dX=dX+dt*ddX;
    ddX=-inv([1,-l*alpha*sin(X(2,1));-BETA*sin(X(2,1)),l])*[-lambdapp(n,1)+(l*alpha*cos(X(2,1)))*dX(2,1)^2-2*eps1*omega1*dX(1,1)-omega1^2*X(1,1);BETA*sin(X(2,1))*lambdapp(n,1)-2*eps2*omega2*l*dX(2,1)-omega2^2*l*sin(X(2,1))];
    % Calcul du residu
      res=-calc_Fnl(X,dX,ddX,lambdapp(n,1));
    
    while (norm(res)>precNR);    %Newton Raphson
    
        % Calcul du residu
         res=-calc_Fnl(X,dX,ddX,lambdapp(n,1));
        % Calcul de la matrice effective
        [dFX dFdX dFddX]=calc_dFnl(X,dX,ddX,lambdapp(n,1));
        J=(4/dt^2)*(dFddX)+(2/dt)*(dFdX)+dFX;
        % Calcul de la correction
        deltaX=J\res;
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
       
        
    end
    Xt(:,n)=X;
    dXt(:,n)=dX;
end
