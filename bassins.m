close all
clear vars

global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K

k = 0.2;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
m = 0.5;%masse de la tige(solide S2)
Ma = 5;%masse du solide S1
g =9.81 ;%valeur du champs de gravité
l = 5;%la moitie la longueure de la tige(solide S2)
Io = 4*m*(l^2)/3;%moment d'inertie de la tige par rapport à son axe de rotation
eps1 = 0.005;%facteur d'amortissement visqueux lié à omega1
eps2 = 0.005;%facteur d'amortissement visqueux li2 à omega2
a = 0;%amplitude de l'exitation d'entree 
omega1 = sqrt(k /(m+Ma));%pulsation du système masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de résonnance du pendule seul
omega0 = (omega1 + omega2)/2;%pulsation de l'exitation  
alpha = m/(m+Ma) ;%coefficient adimensionnel 
BETA = m*l^2/Io;%coefficient adimensionnel 

M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];

periode=2*pi/omega0;      % periode de l'excitation et de la reponse
nb_pts_per=30;          % nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per;  % taille du pas de temps
nb_per=30;              % nb de periodes pour le calcul temporel
t_tot=nb_per*periode;   % temps final
t_init=0;               % temps initial  

% Conditions initiales
dz0=0.25; dv0=0.5;
x0=-3:dz0:4;   
dx0=-7:dv0:7;

W=zeros(length(dx0),length(x0));

for i=1:length(x0)      % boucle sur les CI
    for j=1:length(dx0)
          [tt,Xt,dXt]=newmarklin(x0,dx0,t_init,dt,t_tot);
          W(j,i)=max(Xt(1,end-2*nb_pts_per:end));
    end
    figure(1)
    pcolor(x0,dx0,W)
    caxis([0 max(W(:))])
    drawnow
%   shading('flat')
end

