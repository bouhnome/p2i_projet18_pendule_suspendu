function [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot)

global OMEGA OMEGA1 OMEGA2
global M C K
%%%%%STARTED HERE



precNR=1.e-9;

% sol init
t=t_init;n=1;
X=X0; dX=dX0; ddX = ddX0;

Fnl=zeros(size(X));
dFX=zeros(length(X));
dFdX=zeros(length(X));
dFddX = zeros(length(X));



Fnl=calc_Fnl(X,dX,ddX);
ddX=M\(P-C*dX-K*X-Fnl);
tt(n)=t;
Xt(:,n)=X;
dXt(:,n)=dX;
% integration temporelle
for t=t_init+dt:dt:t_tot;   %Boucle sur les pas de temps
    n=n+1;
    % prediction
    iter=0;
    X=X+dt*dX+(dt^2/2)*ddX;
    dX=dX+dt*ddX;
    ddX=ddX;
    % Calcul du residu
    P=calc_P(t);
    Fnl=calc_Fnl(X,dX);
    res=P-M*ddX-C*dX-K*X-Fnl;
    while (norm(res)>precNR);    %Newton Raphson
        iter=iter+1;
        % Calcul de la Jacobienne
        [dFX dFdX]=calc_dFnl(X,dX);
        J=(4/dt^2)*(M+dFddX)+(2/dt)*(C+dFdX)+K+dFX;
        % Calcul de la correction
        deltaX=J\res;
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
        % Calcul du residu
        Fnl=calc_Fnl(X,dX);
        res=P-M*ddX-C*dX-K*X-Fnl;
%         res=deltaX;
    end
    tt(n)=t;
    Xt(:,n)=X;
    dXt(:,n)=dX;
end
