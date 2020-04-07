clear all; close all;

global k m Ma g l Io omega1 omega2 omega0 eps1 eps2 alpha Beta a



periode=2*pi/omega0;      % periode de l'excitation et de la reponse
nb_pts_per=30;          % nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per;  % taille du pas de temps
nb_per=70;              % nb de periodes pour le calcul temporel
t_tot=nb_per*periode;   % temps final
t_init=0;               % temps initial

%Le vecteur X est le vecteur des inconnues il contient z et theta 
z0 = 1;%valeur de z initiale 
zp0 = 2;%valeur de z_point initiale
theta0 = pi;%valeur de theta initiale 
thetap0 = 2;% valeur de theta_point initiale
X0=[z0;theta0];dX0=[zp0;thetap0];             % conditions initiales

[tt,Xt,dXt]=newmarklin(X0,dX0,t_init,dt,t_tot);   % Integration par newmark

%% Comparaison avec solution analytique

%SHABAAZ COMPLETE ICI 
thetaAnaly = ;
zAnaly = ;
thetapAnaly = ;
zpAnaly = ;



% Plot de z
figure(1)
plot(tt,zAnaly(1,:))
hold on
plot(tt,Xt(1,:))
legend('Z- analytique','Z- newmark')
%Plot de theta
figure(2)
plot(tt,thetaAnaly(1,:));
hold on
plot(tt,Xt(2,:));
legend('theta- analytique','theta- newmark')
%Plot de Zpoint
figure(3)
plot(tt,zpAnaly(1,:));
hold on
plot(tt,dXt(1,:));
legend('Zpoint- analytique','Zpoint- newmark')
%Plot de thetapoint
figure(4)
plot(tt,thetapAnaly(1,:));
hold on
plot(tt,dXt(2,:));
legend('thetapoint- analytique','thetapoint- newmark')
