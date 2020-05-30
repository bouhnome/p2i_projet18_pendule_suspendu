close all
clear vars

global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K

k = 300;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
m = 0.1;%masse de la tige(solide S2)
Ma = 0.5;%masse du solide S1
g =9.81 ;%valeur du champs de gravité
l = 0.1;%la moitie la longueure de la tige(solide S2)
Io = 4*m*(l^2)/3;%moment d'inertie de la tige par rapport à son axe de rotation
eps1 = 0.005;%facteur d'amortissement visqueux lié à omega1
eps2 = 0.005;%facteur d'amortissement visqueux li2 à omega2 
omega1 = sqrt(k /(m+Ma));%pulsation du système masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de résonnance du pendule seul
omega0 = (omega1 + omega2)/2;%pulsation de l'exitation  
alpha = m/(m+Ma) ;%coefficient adimensionnel 
BETA = m*l^2/Io;%coefficient adimensionnel 

M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];

%
scz=get(0,'screensize');
%
% ne nombre de point
% amplitude initiale, finale et increment d'amplitude pour les diagrammes des bifuraction
% ne nombre de point d'incrementation de l'amplitude
ai=0.05; af=0.1; ne=101; da=(af-ai)/(ne-1); 
%
a=ai;
% 
hfg=zeros(1,10);
for ii=1:6
   hfg(ii)=figure(ii);
end
%
% Position initiale
z0 = 0.5;%valeur de z initiale 
zp0 = 0;%valeur de z_point initiale
theta0 = 45*pi/180;%valeur de theta initiale 
thetap0 = 0;% valeur de theta_point initiale
X0=[z0;theta0];
dX0=[zp0;thetap0]; 
%
% pulsation et periode d'excitation 
T=2*pi/omega0;

% np nombre de periode d'integration pour chaque solution temporelle
% ns nombre de periode d'integration a partir de laquelle la reponse est consideree comme stable
% ni nombre de point d'integration pour la plus petite periode
% nth nombre de periode de trace de la reponse temporelle en theta
np=500; ns=400; ni=40; nth=20;
% balayage en temps
t0=0; tf=np*T; dt=T/ni;
%
%

for i=1:ne
    % incrementaion de a de 1 mm
    a=a+da;
    %
    % integration des equations de mouvement
    [tt,Xt,dXt]=newmark(X0,dX0,t0,dt,tf);
    %
    % z et position angulaire du pendule
    z=Xt(1,:);
    theta=Xt(2,:);
    %
    % nt nombre de points d'integration 
    % na ordre du point e partir duquel on considere que la solution est stable
    nt=length(z);
    na=nt-(np-ns)*ni;
    %
    %
%   position
    figure(hfg(1));
    set(hfg(1),'position',[10 50 scz(3)/4 scz(4)/3]);
    plot(tt(nt-nth*ni:nt),z(nt-nth*ni:nt));
    xlabel('t (s)')
    ylabel('z (m)')
    title('z');
    drawnow
    %
%   theta sur les nth dernieres periodes
    figure(hfg(2));
    set(hfg(2),'position',[10 150+scz(4)/3 scz(3)/4 scz(4)/3]);
    plot(tt(nt-nth*ni:nt),theta(nt-nth*ni:nt))
    xlabel('t (s)')
    ylabel('\theta (°)')
    title('position angulaire')
    drawnow
%
    thetap=dXt(2,:);
    zp=dXt(1,:);
%
%   construction des sections de Poincare par echantillonnage tous les ni
%   points des ns derniere periodes
    xp=theta(na:ni:nt);
    xpp=thetap(na:ni:nt);
    yp=z(na:ni:nt);
    ypp=zp(na:ni:nt);
% 
% trace des sections de Poincare pour theta et z
    figure(hfg(3));
    set(hfg(3),'position',[3*scz(3)/4-10 50 scz(3)/4 scz(4)/3]);
    plot(yp(:),ypp(:),'.b');
    xlabel('z (m)')
    ylabel('dz/dt (m/s)')
    title('section de Poincare "en z"')
    drawnow
%  
    figure(hfg(4));
    set(hfg(4),'position',[3*scz(3)/4-10 150+scz(4)/3 scz(3)/4 scz(4)/3]);
    plot(xp(:),xpp(:),'.b');
    xlabel('\theta (°)')
    ylabel('d\theta/dt (°/s)')
    title('section de Poincare "en theta"')
    drawnow
    
    amp(1:length(ypp))=a;
    % trace du diagramme de bifurcation
    figure(hfg(5));
    set(hfg(5),'position',[2*scz(3)/4-100 50 scz(3)/4 scz(4)/3]);
    plot(amp,ypp,'.b');
    xlabel('a (m)')
    ylabel('z (m)')
    title('diagramme de bifurcation "en z"')
    drawnow
    hold on
%
    figure(hfg(6));
    set(hfg(6),'position',[2*scz(3)/4-100 150+scz(4)/3 scz(3)/4 scz(4)/3]);
    plot(amp,xpp,'.b');
    xlabel('a (m)')
    ylabel('\theta (rd)')
    title('diagramme de bifurcation "en theta"')
    drawnow
    hold on
%
end
