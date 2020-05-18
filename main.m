clear all; close all;
%% Constantes du problème
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
a = 0.001;%amplitude de l'exitation d'entree 
omega1 = sqrt(k /(m+Ma));%pulsation du système masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de résonnance du pendule seul
omega0 = (omega1 + omega2)/2;%pulsation de l'exitation  
alpha = m/(m+Ma) ;%coefficient adimensionnel 
BETA = m*l^2/Io;%coefficient adimensionnel 


M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];



%% Données de simulation
periode=2*pi/omega0;      % periode de l'excitation et de la reponse
nb_pts_per=40;          % nb de points par periode pour l integration temporelle
%il faut noter que la valeur du pas de temps change dramatiquement les
%solutions calculés en régime chaotique(effet papillon?)
%aussi augmenter le pas de temps donne des résultats plus précis mais
%ralentit l'animation
dt=periode/nb_pts_per;  % taille du pas de temps
nb_per=70;              % nb de periodes pour le calcul temporel
t_tot=nb_per*periode;   % temps final
t_init=0;               % temps initial     

%Le vecteur X est le vecteur des inconnues il contient z et theta 
z0 = 0.003;%valeur de z initiale 
zp0 = 0;%valeur de z_point initiale
theta0 = 2*pi/180;%valeur de theta initiale 
thetap0 = 0;% valeur de theta_point initiale
X0=[z0;theta0];dX0=[zp0;thetap0];             % conditions initiales
         

[tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot);   % Integration par Newmark

%% Comparaison avec ode45 à effectuer(non linéaire)
[t,x]=ode45(@Pendule_suspendu,tt,[z0;zp0;theta0;thetap0]);
%% Comparaison avec la solution Analytique à effectuer(linéaire)
%Constantes relatives à la solution analytique
v=-eps2*omega2;
mu=sqrt(omega2^2-(eps2*omega2)^2);
A=theta0;
B=(thetap0-A*v)/mu;


 thetaAnaly = A*exp(v*t).*cos(mu*t)+B*exp(v*t).*sin(mu*t);
 thetapAnaly = A*v*exp(v*t).*cos(mu*t)-mu*A*exp(v*t).*sin(mu*t)+B*v*exp(v*t).*sin(mu*t)+B*mu*exp(v*t).*cos(mu*t);
 
v=-eps1*omega1;
mu=sqrt(omega1^2-(eps1*omega1)^2);
coeffsVW=inv([omega1^2-omega0^2,2*eps1*omega1*omega0;-2*eps1*omega1*omega0,omega1^2-omega0^2])*[0;a*omega0^2];
V=coeffsVW(1);
W=coeffsVW(2);
A=z0-V;
B=(zp0-A*v-W*omega0)/mu;
 

 zAnaly =A*exp(v*t).*cos(mu*t)+B*exp(v*t).*sin(mu*t)+V*cos(omega0*t)+W*sin(omega0*t);
 zpAnaly =A*v*exp(v*t).*cos(mu*t)-mu*A*exp(v*t).*sin(mu*t)+B*v*exp(v*t).*sin(mu*t)+B*mu*exp(v*t).*cos(mu*t)-V*omega0*sin(omega0*t)+W*omega0*cos(omega0*t);
%% Plot et animation

if(theta0<=2*pi/180 && z0<=0.003 && a<=0.001)
    choix=2;
else
choix = 1; %variable choix prennant des valeurs de 1 à 3 pour choisir d'afficher la solution non-linéaire et la solution linéaire
end

if choix ==1 %Affichage de la solution non-linéaire et comparaison avec ODE45
%Plot de z
figure(1)
plot(tt,Xt(1,:),tt,x(:,1))
legend('Newmark','ODE45')
title('Evolution de z')

%Plot de theta
figure(2)
plot(tt,Xt(2,:),tt,x(:,3));
legend('Newmark','ODE45')
title('Evolution de theta')

%Plot de Zpoint
figure(3)
plot(tt,dXt(1,:),tt,x(:,2));
legend('Newmark','ODE45')
title('Evolution de z point')

%Plot de thetapoint
figure(4)
plot(tt,dXt(2,:),tt,x(:,4));
legend('Newmark','ODE45')
title('Evolution de z theta point')


%Animation
figure(5)
animation(tt,Xt,a*sin(omega0*tt)) 
end

if choix ==2 % Affichage de la solution linéaire et comparaison avec la solution analytique
%A NOTER que l'on observe la superposition de toutes les figures pour des
%valeurs initiales inférieures ou égales à 3mm et 2°
%on peut prendre une valeur d'amplitude de l'exitation inférieur ou égale à
%1mm

%Plot de z
figure(1)
plot(tt,Xt(1,:),tt,zAnaly)
legend('Newmark','Analytique')
title('Evolution de z')

%Plot de theta
figure(2)
plot(tt,Xt(2,:),tt,thetaAnaly);
legend('Newmark','Analytique')
title('Evolution de theta')

%Plot de Zpoint
figure(3)
plot(tt,dXt(1,:),tt,zpAnaly);
legend('Newmark','Analytique')
title('Evolution de z point')

%Plot de thetapoint
figure(4)
plot(tt,dXt(2,:),tt,thetapAnaly);
legend('Newmark','Analytique')
title('Evolution de theta point')


%Animation
figure(5)
animation(tt,Xt,a*sin(omega0*tt)) 
end

%%