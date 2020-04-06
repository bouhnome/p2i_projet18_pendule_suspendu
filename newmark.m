function [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot)

global l0 a omega1 omega2 epsilon1 epsilon2 omega;
% global knl f0 OMEG
global M C K

precNR=1.e-11;
% sol init
t=t_init;n=1;
X=X0;dX=dX0;

Fnl=zeros(size(X));
dFX=zeros(length(X));
dFdX=zeros(length(X));
P=zeros(size(X));
% P=f0*cos(OMEG*t);
Fnl=calc_Fnl(X,dX,t);
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
    P=0;
%     P=f0*cos(OMEG*t);
    Fnl=calc_Fnl(X,dX,t);
    res=P-M*ddX-C*dX-K*X-Fnl;
    normres=norm(res)/norm(P);
    while (normres>precNR);    %Newton Raphson
        iter=iter+1;
        % Calcul de la Jacobienne
        [dFX dFdX]=calc_dFnl(X,dX,t);
        J=(4/dt^2)*M+(2/dt)*(C+dFdX)+K+dFX;
        % Calcul de la correction
        deltaX=J\res;
        X=X+deltaX;
        dX=dX+(2/dt)*deltaX;
        ddX=ddX+(4/dt^2)*deltaX;
        % Calcul du residu
        Fnl=calc_Fnl(X,dX,t);
        res=P-M*ddX-C*dX-K*X-Fnl;
%         normres=norm(res)/norm(P);
        normres=norm(deltaX)/norm(X);
    end
    tt(n)=t;
    Xt(:,n)=X;
    dXt(:,n)=dX;
end
