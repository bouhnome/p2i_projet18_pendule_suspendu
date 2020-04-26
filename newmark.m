function [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot)

%%les constantes du problème
global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K
%%
tt=[t_init:dt:t_tot];%vecteur temps à retourner
precNR=1.e-9;
lambdapp=-a*omega0^2*sin(tt*omega0);
lambdapp=lambdapp';

%% sol init
t=t_init;n=1;
X=X0; dX=dX0;
%on détermine la valeur de l'accélération initiale
ddX=inv([1,-l*alpha*sin(X(2,1));-BETA*sin(X(2,1)),l])*[-lambdapp(1,1)+(l*alpha*cos(X(2,1)))*dX(2,1)^2-2*eps1*omega1*dX(1,1)-omega1^2*X(1,1);BETA*sin(X(2,1))*lambdapp(1,1)-2*eps2*omega2*l*dX(2,1)-omega2^2*l*sin(X(2,1))];



%% integration temporelle
for t=t_init+dt:dt:t_tot;
    n=n+1;
    
    % prediction
    X=X+dt*dX+(dt^2/2)*ddX;
    dX=dX+dt*ddX;
    % Calcul du residu
      Fnl=calc_Fnl(X,dX,ddX,lambdapp(n,1));
      res=-M*ddX-C*dX-K*X-Fnl;

    
    
    while (norm(res)>precNR);    %Newton Raphson
    
        % Calcul du residu
        Fnl=calc_Fnl(X,dX,ddX,lambdapp(n,1));
        res=-M*ddX-C*dX-K*X-Fnl;
        % Calcul de la matrice effective
        [dFX dFdX dFddX]=calc_dFnl(X,dX,ddX,lambdapp(n,1));
        J=(4/dt^2)*(M+dFddX)+(2/dt)*(C+dFdX)+K+dFX;
        % Calcul de la correction
        deltaX=J\res;
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
       
        
    end
    Xt(:,n)=X;
    dXt(:,n)=dX;
end
