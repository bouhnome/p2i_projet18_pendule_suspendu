clear all; close all;

global knl p0 OMEGA
global M C K

omega0=1.;
M=1;
C=0.;
K=M*omega0^2;
knl=0.;
p0=0.;
OMEGA=1.;

periode=2*pi/OMEGA;      % periode de l'excitation et de la reponse
nb_pts_per=30;          % nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per;  % taille du pas de temps
nb_per=10;              % nb de periodes pour le calcul temporel
t_tot=nb_per*periode;   % temps final
t_init=0;               % temps initial

X0=1;dX0=0;             % conditions initiales
[tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot);   % Integration par Newmark
plot(tt,Xt,'b+-')        % On trace le deplacement au cours du temps

%% Comparaison avec sol analytique a decommenter
% hold on;
% t=t_init:dt:t_tot;
% plot(t,cos(omega0*t),'r-')
% title('Cas lineaire non force')
