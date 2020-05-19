close all
clear vars

global k m Ma g l Io eps1 eps2 a omega1 omega2 omega0 alpha BETA M C K

k = 0.2;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
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
ai=0.05; af=0.1; ne=200; da=(af-ai)/(ne-1); 
%
a=ai;
f1=omega1/(2*pi);
T1=1/f1;
%

f2=omega2/(2*pi);
T2=1/f2;
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
omega=omega2*5/4;
% omega=omega2*5.1/4;
T=2*pi/omega0;
%
% np nombre de periode d'integration pour chaque solution temporelle
% ns nombre de periode d'integration a partir de laquelle la reponse est consideree comme stable
% ni nombre de point d'integration pour la plus petite periode
% nth nombre de periode de trace de la reponse temporelle en theta
np=1500; ns=800; ni=100; nth=100;
% balayage en temps
t0=0; tf=np*T; dt=T/ni; t_balayage=t0:dt:tf;
%
%
for k=1:ne
    % incrementaion de a de 1 mm
    a=a+da;
    %
    % integration des equations de mouvement
    [tt,Xt,dXt]=newmarklin(X0,dX0,t0,dt,tf);
    %
    % longueur et position angulaire du pendule
    z=Xt(1,:);
    theta=Xt(2,:);
    %
    % nt nombre de points d'integration 
    % na ordre du point e partir duquel on considere que la solution est stable
    nt=length(z);
    na=nt-(np-ns)*ni;
    %
    xh=z.*tan(theta);
    %
%     % trajectoire
    figure(hfg(1));
    set(hfg(1),'position',[10 50 scz(3)/4 scz(4)/3]);
    plot(xh(na:nt),-z(na:nt),'Ob');
    xlabel('x (m)')
    ylabel('z (m)')
    axis([-0.1 0.1 -0.5 0.5]);
    title('positions successives de la masse du pendule');
    drawnow
    %
%     % theta sur les nth dernieres periodes
    figure(hfg(2));
    set(hfg(2),'position',[10 150+scz(4)/3 scz(3)/4 scz(4)/3]);
    plot(tt(nt-nth*ni:nt),theta(nt-nth*ni:nt))
    xlabel('t (s)')
    ylabel('\theta (rd)')
    axis([(nt-nth*ni)*dt tf -0.00025 0.00025]);
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
% trace des sections de Poincare pour theta et l
    figure(hfg(3));
    set(hfg(3),'position',[3*scz(3)/4-10 50 scz(3)/4 scz(4)/3]);
    plot(yp(:),ypp(:),'.b');
    xlabel('z (m)')
    ylabel('dz/dt (m/s)')
    % cacul de l'amplitude de la reponse e partir du modele lineaire
    % ls longuer statique et la amplitude dynamique
    % amplitude de dl/dt en multipliant par -omega
    axis([-0.1 1.6 -4.5 0]);
    title('section de Poincare "en z"')
    drawnow
%  
    figure(hfg(4));
    set(hfg(4),'position',[3*scz(3)/4-10 150+scz(4)/3 scz(3)/4 scz(4)/3]);
    plot(xp(:),xpp(:),'.b');
    xlabel('\theta (rd)')
    ylabel('d\theta/dt (rd/s)')
    axis([-1.2 1.2 -6.0 6.0]);
    title('section de Poincare "en theta"')
    drawnow
    
    amp(1:length(yp))=a;
    % trace du diagramme de bifurcation
    figure(hfg(5));
    set(hfg(5),'position',[2*scz(3)/4-100 50 scz(3)/4 scz(4)/3]);
    plot(amp,yp,'.b');
    xlabel('a (m)')
    ylabel('z (m)')
    axis([ai ai+ne*da -0.6 0.6]);
    title('diagramme de bifurcation "en z"')
    drawnow
    hold on
%
    figure(hfg(6));
    set(hfg(6),'position',[2*scz(3)/4-100 150+scz(4)/3 scz(3)/4 scz(4)/3]);
    plot(amp,xp,'.b');
    xlabel('a (m)')
    ylabel('\theta (rd)')
    axis([ai ai+ne*da -0.04 0.04]);
    title('diagramme de bifurcation "en theta"')
    drawnow
    hold on
%
end


