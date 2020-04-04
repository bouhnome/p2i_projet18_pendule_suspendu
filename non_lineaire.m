clear all; close all;

global k m Ma g l Io omega1 omega2 omega eps1 eps2 alpha beta a

k = 0.2;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
m = 0.5;%masse de la tige(solide S2)
Ma = 5;%masse du solide S1
g =9.81 ;%valeur du champs de gravit�
l = 5;%la moitie la longueure de la tige(solide S2)
Io = 4*m*(l^2);%moment d'inertie de la tige par rapport � son axe de rotation
eps1 = 0.005;%facteur d'amortissement visqueux li� � omega1
eps2 = 0.005;%facteur d'amortissement visqueux li2 � omega2
a = 0.02;%amplitude de l'exitation d'entree 
omega1 = sqrt(k /(m+Ma));%pulsation du syst�me masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de r�sonnance du pendule seul
omega = (omega1 + omega2)/2;%pulsation de l'exitation  
alpha = m/(m+Ma) ;%coefficient adimensionnel 
beta = m*l^2/Io;%coefficient adimensionnel 

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
timespan = t_init:dt:t_tot;

%Le vecteur X est le vecteur des inconnues, il contient z et theta 
z0 = 0;%valeur de z initiale 
zp0 = 0.1;%valeur de z_point initiale
theta0 = 0;%valeur de theta initiale 
thetap0 = 0.1;% valeur de theta_point initiale
X0=[z0;theta0];dX0=[zp0;thetap0];           % conditions initiales

[tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot);   % Integration par Newmark

%% Comparaison avec ode45
Y0 = [z0,zp0, theta0, thetap0]; 
[tt, Y] = ode45(@Pendule_Function_NL, timespan, Y0);


