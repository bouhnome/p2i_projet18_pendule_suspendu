clear vars; close all;

global Io l a omega1 omega2 eps1 eps2 omega alpha beta;
global Ma m g k;
global M C K;

M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];
k ; m ; Ma ; g=9.8 ; l ; Io ;
eps1 ; eps2 ;
omega ; lambda0 ;
omega1 = sqrt(k /(m+Ma));
omega2 = sqrt(m*g*l/Io);
alpha = m/(m+Ma) ;

periode=2*pi/omega;
nb_pts_per=30;
dt=periode/nb_pts_per;
ttot=40*periode;

OMEGA_debut=0.5;OMEGA_fin=2.5;dOMEGA=0.05;

nb_pts_per=30;          % nb de points par periode pour l integration temporelle
nb_per=50;              % nb de periodes pour le calcul temporel
t_init=0;               % temps initial

% conditions initiales
z0=0.5;dz0=0;
k=0; omega=OMEGA_debut;

%% Montee en frequence
% boucle sur Omega 
while (omega<OMEGA_fin)
  OME(k+1)=omega;
  periode=2*pi/omega;      % periode de l'excitation et de la reponse
  dt=periode/nb_pts_per;  % taille du pas de temps
  t_tot=nb_per*periode;   % temps final
  % calcul de z(t)
  [tt,Xt,dXt]=newmark(z0,dzoui0,t_init,dt,t_tot);   % Integration par Newmark
  % recherche de l'amplitude max
  W(k+1)=max(Xt(1,end-3*nb_pts_per:end));
  txt=sprintf('ome=%7.5f x=%0.5g',OME(k+1),W(k+1));   % affichage resultat
  disp(txt);
  z0=Xt(1,end);dz0=dXt(1,end);   % nouvelles CI
  omega=omega+dOMEGA;
  k=k+1;
end
plot(OME,W,'r-o')    % courbe de reponse
title('Courbe de reponse')
xlabel('Omega');ylabel('Amplitude max');

