close all
clear vars
clf
%clc

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

% Conditions initiales
% pas de la grille
dz0=0.5; dv0=1;    % grille grossiere
% dx0=0.2; dv0=.4;   % grille moyenne
% bornes de la grille
z0=-6:dz0:7;   
v0=-15:dv0:15;

W=zeros(length(v0),length(z0));

for i=1:length(x0)      % boucle sur les CI
    for j=1:length(v0)
%         [tt,Xt,dXt]=newmark(x0(i),v0(j),0,dt,ttot);
          [tt,Xt,dXt]=newmarklin(z0(i),v0(j),0,dt,ttot);
          W(j,i)=max(Xt(1,end-2*nb_pts_per:end));
    end
    figure(1)
    pcolor(z0,v0,W)
    caxis([0 max(W(:))])
    drawnow
%     shading('flat')
end
