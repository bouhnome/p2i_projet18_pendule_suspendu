clear all; close all;

global knl p0 OMEGA
global M C K

M=1;
C=0.1;
K=1;
knl=0.25;
p0=.5;

OMEGA_debut=.5;OMEGA_fin=2.5;dOMEGA=0.05;

nb_pts_per=50;          % nb de points par periode pour l integration temporelle
nb_per=50;              % nb de periodes pour le calcul temporel
t_init=0;               % temps initial

% conditions initiales
X0=0.5;dX0=0;
k=0; OMEGA=OMEGA_debut;

%% Montee en frequence
% boucle sur Omega 
while (OMEGA<OMEGA_fin)
  OME(k+1)=OMEGA;
  periode=2*pi/OMEGA;      % periode de l'excitation et de la reponse
  dt=periode/nb_pts_per;  % taille du pas de temps
  t_tot=nb_per*periode;   % temps final
  % calcul de z(t)
  [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot);   % Integration par Newmark
  % recherche de l'amplitude max
  W(k+1)=max(Xt(1,end-3*nb_pts_per:end));
  txt=sprintf('ome=%7.5f x=%0.5g',OME(k+1),W(k+1));   % affichage resultat
  disp(txt);
  X0=Xt(1,end);dX0=dXt(1,end);   % nouvelles CI
  OMEGA=OMEGA+dOMEGA;
  k=k+1;
end
plot(OME,W,'r-o')    % courbe de reponse
title('Courbe de reponse')
xlabel('Omega');ylabel('Amplitude max');

%% Descente en frequence  a decommenter
% % conditions initiales
% X0=0.5;dX0=0;
% k=0; OMEGA=OMEGA_fin;
% % boucle sur Omega 
% while (OMEGA>OMEGA_debut)
%   OME2(k+1)=OMEGA;
%   periode=2*pi/OMEGA;      % periode de l'excitation et de la reponse
%   dt=periode/nb_pts_per;  % taille du pas de temps
%   t_tot=nb_per*periode;   % temps final
%   % calcul de z(t)
%   [tt,Xt,dXt]=newmark(X0,dX0,t_init,dt,t_tot);   % Integration par Newmark
%   % recherche de l'amplitude max
%   W2(k+1)=max(Xt(1,end-3*nb_pts_per:end));
%   txt=sprintf('ome=%7.5f x=%0.5g',OME(k+1),W(k+1));
%   disp(txt);
%   X0=Xt(1,end);dX0=dXt(1,end);
%   OMEGA=OMEGA-dOMEGA;
%   k=k+1;
% end
% hold on;
% plot(OME2,W2,'b-^')
