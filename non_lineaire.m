clear all; close all;

global k m Ma g l Io omega1 omega2 omega eps1 eps2 alpha beta lamda0

k = ;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
m = ;%masse de la tige(solide S2)
Ma = ;%masse du solide S1
g =9.8 ;%valeur du champs de gravité
l = ;%2 fois la longeure de la tige(solide S2)
Io = ;%moment d'inertie de la tige par rapport à son axe de rotation
eps1 = ;%facteur d'amortissement visqueux lié à omega1
eps2 = ;%facteur d'amortissement visqueux li2 à omega2
omega = ;%pulsation de l'exitation  
lambda0 = ;%amplitude de l'exitation 
omega1 = sqrt(k /(m+Ma));%pulsation du système masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de résonnance du pendule seul
alpha = m/(m+Ma) ;%coefficient adimensionnel 
beta = m*l^2/Io;coefficient adimensionnel 

global M C K %matrices-coefficients du formalisme de newmark 

M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];


periode=2*pi/omega;      % periode de l'excitation et de la reponse
nb_pts_per=30;          % nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per;  % taille du pas de temps
nb_per=30;              % nb de periodes pour le calcul temporel
t_tot=nb_per*periode;   % temps final
t_init=0;               % temps initial

%Le vecteur X est le vecteur des inconnues il contient z et theta 
z0 = ;%valeur de z initiale 
zp0 = ;%valeur de z_point initiale
zpp0 = ;%valeur de z_point_point initiale
theta0 = ;%valeur de theta initiale 
thetap0 = ;% valeur de theta_point initiale
thetapp0 = ;%valeur de theta_point_point initiale
X0=[z0;theta0];dX0=[zp0;thetap0];ddX0 =[zpp0;thetapp0];             % conditions initiales

[tt,Xt,dXt]=newmark(X0,dX0,ddX0,t_init,dt,t_tot);   % Integration par Newmark

plot(tt,Xt,'b+-')        % On trace le deplacement au cours du temps

%% Comparaison avec ode45 a decommenter à effectuer

