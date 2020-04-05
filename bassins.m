close all
clear vars
clf
%clc

global knl p0 OMEGA
global M C K

M=1;
C=0.1;
K=1;
knl=0.25;
p0=.5;

% OMEGA=1.35;
OMEGA=1.5;
% OMEGA=1.62;
periode=2*pi/OMEGA;
nb_pts_per=31;
dt=periode/nb_pts_per;
ttot=40*periode;

% Conditions initiales
% pas de la grille
dx0=0.5; dv0=1;;    % grille grossiere
% dx0=0.2; dv0=.4;   % grille moyenne
% bornes de la grille
x0=-6:dx0:7;   
v0=-15:dv0:15;

W=zeros(length(v0),length(x0));

for i=1:length(x0)      % boucle sur les CI
    for j=1:length(v0)
%         [tt,Xt,dXt]=newmark(x0(i),v0(j),0,dt,ttot);
          [tt,Xt,dXt]=newmark2(x0(i),v0(j),0,dt,ttot);
          W(j,i)=max(Xt(1,end-2*nb_pts_per:end));
    end
    figure(1)
    pcolor(x0,v0,W)
    caxis([0 max(W(:))])
    drawnow
%     shading('flat')
end
