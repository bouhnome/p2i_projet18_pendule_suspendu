clear all; close all;

global k m Ma g l Io omega1 omega2 omega eps1 eps2 w0 a

periode=2*pi/omega;      % periode de l'excitation et de la reponse
nb_pts_per=30;          % nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per;  % taille du pas de temps
nb_per=30;              % nb de periodes pour le calcul temporel
t_tot=nb_per*periode;   % temps final
t_init=0;               % temps initial

%Le vecteur X est le vecteur des inconnues il contient z et theta 
z0 = 0;%valeur de z initiale 
zp0 = 0.1;%valeur de z_point initiale
theta0 = 0;%valeur de theta initiale 
thetap0 = 0.1;% valeur de theta_point initiale
X0=[z0;theta0];dX0=[zp0;thetap0];             % conditions initiales

[t,Xt,dXt]=newmarklin(X0,dX0,t_init,dt,t_tot);   % Integration par newmark

%% Comparaison avec solution analytique
%N = nombre d'iteration et taille de la matrice
N =1000;
Z = zeros(N,1);
Th = zeros(N,1);
Zp = zeros(N,1);
Thp = zeros(N,1);
for i = 1:N
    Sol = Pendule_Function_L(i);
    Z(i,1) = Sol(1);
    Th(i,1) = Sol(2);
    Zp(i,1) = sol(3);
    Thp(i,1) = Sol(4);
end
t = 1:N;
% Plot de z
figure(1)
plot(t,Z(:,1))
hold on
plot(t,Xt(1,:))
legend('Z- analytique','Z- newmark')
%Plot de theta
figure(2)
plot(t,Th(:,1));
hold on
plot(t,Xt(2,:));
legend('theta- analytique','theta- newmark')
%Plot de Zpoint
figure(3)
plot(t,Zp(:,1));
hold on
plot(t,dXt(1,:));
legend('Zpoint- analytique','Zpoint- newmark')
%Plot de thetapoint
figure(4)
plot(t,Thp(:,1));
hold on
plot(t,dXt(2,:));
legend('thetapoint- analytique','thetapoint- newmark')


 

