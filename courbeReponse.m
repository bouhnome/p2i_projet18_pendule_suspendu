close all
clear vars

global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K

%%les constantes du problème
k = 300;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
m = 0.1;%masse de la tige(solide S2)
Ma = 0.5;%masse du solide S1
g =9.81 ;%valeur du champs de gravité
l = 0.1;%la moitie la longueure de la tige(solide S2)
Io = 4*m*(l^2)/3;%moment d'inertie de la tige par rapport à son axe de rotation
eps1 = 0.005;%facteur d'amortissement visqueux lié à omega1
eps2 = 0.005;%facteur d'amortissement visqueux li2 à omega2
a = 0.05;%amplitude de l'exitation d'entree 
omega1 = sqrt(k /(m+Ma));%pulsation du système masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de résonnance du pendule seul
omega0 = (omega1 + omega2)/2;%pulsation de l'exitation  
alpha = m/(m+Ma) ;%coefficient adimensionnel 
BETA = m*l^2/Io;%coefficient adimensionnel 

M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];

nb_pts_per=30;          % nb de points par periode pour l integration temporelle
nb_per=30;              % nb de periodes pour le calcul temporel
t_init=0;               % temps initial  

omegaf=5;domega=0.05;

% conditions initiales
z0 = 0.5;%valeur de z initiale 
zp0 = 0;%valeur de z_point initiale
theta0 = 30*pi/180;%valeur de theta initiale 
thetap0 = 0;% valeur de theta_point initiale
X0=[z0;theta0];
dX0=[zp0;thetap0]; 
omega=omega0;

%% Montee en frequence
% boucle sur Omega 
j=0;
while (omega<omegaf)
  OME(j+1)=omega;
  periode=2*pi/omega;      % periode de l'excitation et de la reponse
  dt=periode/nb_pts_per;  % taille du pas de temps
  t_tot=nb_per*periode;   % temps final
  [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot);   % Integration par Newmark
  W(j+1)=max(Xt(1,end-3*nb_pts_per:end));
  X0=Xt(:,end);
  dX0=dXt(:,end);   % nouvelles CI
  omega=omega+domega;
  j=j+1;
end
plot(OME,W,'r-o')    % courbe de reponse
title('Courbe de reponse')
xlabel('Omega');ylabel('Amplitude max');

