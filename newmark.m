function [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot)

%%les constantes du problème

%%
tt=[t_init:dt:t_tot];%vecteur temps à retourner

precNR=1.e-9;

%% sol init
lambdapp0=0;%valeur initialle de lambdapp
t=t_init;n=1;
X=X0; dX=dX0;
%on détermine la valeur de l'accélération initiale
ddX=inv([1,-l*alpha*sin(X(2,1));-BETA*sin(X(2,1)),l])*[-lambdapp0+(l*alpha*cos(X(2,1)))*dX(2,1)^2-2*eps1*omega1*dX(1,1)-omega1^2*X(1,1);BETA*sin(X(2,1))*lambdapp0-2*eps2*omega2*l*dX(2,1)-omega2^2*l*sin(X(2,1))];

Fnl=calc_Fnl(X,dX,ddX);

%% integration temporelle
for t=t_init+dt:dt:t_tot;
    n=n+1;
    
    % prediction
    X=X+dt*dX+(dt^2/2)*ddX;
    dX=dX+dt*ddX;
    ddX=M\(P-C*dX-K*X-Fnl);
    % Calcul du residu
    Fnl=calc_Fnl(X,dX,ddX);
    res=P-M*ddX-C*dX-K*X-Fnl;
    
    
    while (norm(deltaX)>precNR);    %Newton Raphson
    
        % Calcul de la matrice effective
        [dFX dFdX dFddX]=calc_dFnl(X,dX,ddX);
        J=(4/dt^2)*(M+dFddX)+(2/dt)*(C+dFdX)+K+dFX;
        % Calcul de la correction
        deltaX=J\res;
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
        % Calcul du residu
        Fnl=calc_Fnl(X,dX,ddX);
        res=-M*ddX-C*dX-K*X-Fnl;
        
    end
    Xt(:,n)=X;
    dXt(:,n)=dX;
end
