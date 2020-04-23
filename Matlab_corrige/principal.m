clear all; close all;

global knl p0 OMEGA OMEGA1 OMEGA2 EPS1 EPS2
global M C K m



knl=0.25;
p0=.5;
OMEGA=1.5;
% OMEGA=1.5;

OMEGA1 = sqrt(k/(m+M));
%%%%OMEGA2 = sqrt(m*g*l)/

%%%STARTED HERE
M = [1 0;0 l];
C = [2*OMEGA1*EPS1, 0;0 , 2*EPS2*OMEGA2*l];
K = [OMEGA1^2, 0;0, 0];



periode=2*pi/OMEGA;      % periode de l'excitation et de la reponse
nb_pts_per=30;          % nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per;  % taille du pas de temps
nb_per=30;              % nb de periodes pour le calcul temporel
t_tot=nb_per*periode;   % temps final
t_init=0;               % temps initial

X0=1;dX0=0;             % conditions initiales
% X0=1.5;dX0=0;           % conditions initiales
% X0=2;dX0=0;             % conditions initiales
[tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot);   % Integration par Newmark
plot(tt,Xt,'b+-')        % On trace le deplacement au cours du temps

% Pour OCTAVE uniquement, il faut charger le package odepkg pour utiliser ode45
% pkg load odepkg 

%% Comparaison avec ode45 a decommenter
% hold on;
% options = odeset('RelTol',1e-8);
% z0=[X0 dX0];
% [t45 z45]=ode45(@duffing,0:dt:t_tot,z0,options);
% plot(t45,z45(:,1),'r-')
