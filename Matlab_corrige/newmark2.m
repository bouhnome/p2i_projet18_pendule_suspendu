function [tt,Xt,dXt]=newmark2(X0,dX0,t_init,dt,t_tot)

global knl p0 OMEGA
global M C K

precNR=1.e-9;

% sol init
t=t_init;n=1;
X=X0;dX=dX0;

Fnl=zeros(size(X));
dFX=zeros(length(X));
dFdX=zeros(length(X));
P=zeros(size(X));

P=p0*cos(OMEGA*t);
Fnl=knl*X.^3;
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
    P=p0*cos(OMEGA*t);
    Fnl=knl*X.^3;
    res=P-M*ddX-C*dX-K*X-Fnl;
    normres=norm(res)/norm(P);
%     pause
    while (normres>precNR);    %Newton Raphson
        iter=iter+1;
        % Calcul de la Jacobienne
        dFX=3*knl*X.^2; dFdX=0;
        J=(4/dt^2)*M+(2/dt)*(C+dFdX)+K+dFX;
        % Calcul de la correction
        deltaX=J\res;
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
        % Calcul du residu
        Fnl=knl*X.^3;
        res=P-M*ddX-C*dX-K*X-Fnl;
%         normres=norm(res)/norm(P);
        normres=norm(deltaX)/norm(X);
    end
    tt(n)=t;
    Xt(:,n)=X;
    dXt(:,n)=dX;
end
