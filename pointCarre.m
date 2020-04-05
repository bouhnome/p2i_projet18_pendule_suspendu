clear vars
close all
% parametres
% acceleration pesanteur g, longueur a vide du ressort, omega1 pulsation de
% resonance de systeme masse M ressort k, omega2 pulsation de resonance du
% pendule non elastique (longueur l0), epsilon1 et 2 facteurs d'amortissements
% associe a omega1 et omega2, omega pulsation d'excitation et a amplitude
% d'excitation
%
global l0 a omega1 omega2 eps1 eps2 omega;
global M C K m g k;
M=eye(2);C=zeros(2);K=zeros(2);
%
g=9.81; m=0.02; l0=0.6;
%
scz=get(0,'screensize');
%
% ne nombre de point
ne=200;
a=0.01;
da=0.0004522613;
%
% amortissement
eps1=0.03;
eps2=0.02;
%
% raideur des ressots 
k=0.5;
%
% pulsation et frequence de resonance du systeme masse ressort
omega1=sqrt(k/m);
f1=omega1/(2*pi);
T1=1/f1;
%
% pulsation de resonance du pendule seul
omega2=sqrt(g/l0);
f2=omega2/(2*pi);
T2=1/f2;
% 
hfg=zeros(1,10);
for ii=1:6
   hfg(ii)=figure(ii);
end
%
% Position initiale
theta0=0.0;
thetap0=0.01;
rho0=1.0;
rhop0=0;
x0=[rho0,rhop0,theta0,thetap0];
X0=[rho0;theta0];
dX0=[rhop0;thetap0];
%
% pulsation et periode d'excitation 
omega=omega2*5/4;
% omega=omega2*5.1/4;
T=2*pi/omega;
%
% np nombre de periode d'integration pour chaque solution temporelle
% ns nombre de periode d'integration a partir de laquelle la reponse est consideree comme stable
% ni nombre de point d'integration pour la plus petite periode
% nth nombre de periode de trace de la reponse temporelle en theta
np=1500; ns=300; ni=100; nth=100;
% balayage en temps
t0=0; tf=np*T; dt=T/ni; t_balayage=t0:dt:tf;
%
%
for k=1:ne
    % incrementaion de a de 1 mm
    a=a+da;
    %
    % integration des equations de mouvement
%     [t,x]=ode45(@Pendule_elastique,t_balayage,x0);
    [t,x,dx]=newmarklin(X0,dX0,t0,dt,tf);
    %
    % longueur et position angulaire du pendule
%     l=x(:,1)*l0;
%     theta=x(:,3);
    l=x(1,:)*l0;
    theta=x(2,:);
    %
    % nt nombre de points d'integration 
    % na ordre du point e partir duquel on considere que la solution est stable
    nt=length(l);
    na=nt-(np-ns)*ni;
    %
    xh=l.*sin(theta);
    zv=l.*cos(theta);
    %
    % trajectoire
    figure(hfg(1));
    set(hfg(1),'position',[10 50 scz(3)/4 scz(4)/3]);
    plot(xh(na:nt),-zv(na:nt),'Ob');
    xlabel('x (m)')
    ylabel('z (m)')
    axis([-1.2 1.2 -2. 0.2]);
    title('positions successives de la masse du pendule');
    drawnow
    %
    % theta sur les nth dernieres periodes
    figure(hfg(2));
    set(hfg(2),'position',[10 150+scz(4)/3 scz(3)/4 scz(4)/3]);
    plot(t(nt-nth*ni:nt),theta(nt-nth*ni:nt))
    xlabel('t (s)')
    ylabel('\theta (rd)')
    axis([(nt-nth*ni)*dt tf -1.5 1.5]);
    title('position angulaire')
    drawnow
%
%     thetap=x(:,4);
%     lp=x(:,2);
    thetap=dx(2,:);
    lp=dx(1,:);
%
%   construction des sections de Poincare par echantillonnage tous les ni
%   points des ns derniere periodes
    xp=theta(na:ni:nt);
    xpp=thetap(na:ni:nt);
    yp=l(na:ni:nt);
    ypp=lp(na:ni:nt);
% 
% trace des sections de Poincare pour theta et l
    figure(hfg(3));
    set(hfg(3),'position',[3*scz(3)/4-10 50 scz(3)/4 scz(4)/3]);
    plot(yp(:),ypp(:),'.b');
    xlabel('l (rd)')
    ylabel('dl/dt (rd/s)')
    % cacul de l'amplitude de la reponse e partir du modele lineaire
    % ls longuer statique et la amplitude dynamique
    % amplitude de dl/dt en multipliant par -omega
    axis([0.2 1.6 -4.5 -0.5]);
    title('section de Poincare "en l"')
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
end
