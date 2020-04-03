

t_tot=25;
dt=1e-3;
X0=[1;1];
dX0 = [2;2];
tt=0:dt:t_tot;


k = 1;%2 fois la constante de raideur. (chaque ressort a une raideur de k/2)
m = 1;%masse de la tige(solide S2)
Ma = 1;%masse du solide S1
g =9.8 ;%valeur du champs de gravité
l = 1;%2 fois la longeure de la tige(solide S2)
Io = 1;%moment d'inertie de la tige par rapport à son axe de rotation
eps1 = 1;%facteur d'amortissement visqueux lié à omega1
eps2 = 1;%facteur d'amortissement visqueux li2 à omega2
omega = 1;%pulsation de l'exitation  
lamda0 = 1;%amplitude de l'exitation 
omega1 = sqrt(k /(m+Ma));%pulsation du système masse ressort (pendule immobile)
omega2 = sqrt(m*g*l/Io);%pulsation de résonnance du pendule seul
alpha = m/(m+Ma) ;%coefficient adimensionnel 
beta = m*l^2/Io;%coefficient adimensionnel 

M = [1,0;0,l];
C = [2*eps1*omega1,0;0,2*eps2*omega2*l];
K = [omega1^2,0;0,0];
  
  tt=[t_init:dt:t_tot];%vecteur tt à retourner 
  Y = zeros(2*size(X0,1),size(tt,2));%création de la matrice utilisé dans les calculs cf commentaire de l'issue #8
  ddlambda = -omega^2*lamda0*sin(omega*tt);%terme d'exitation 

  %initialisation
  Y(1,1)=X0(1);
  Y(2,1)=X0(2);
  Y(3,1)=dX0(1);
  Y(4,1)=dX0(2);

  %on  pose les variables qui permettent de simplifer les expressions cf commentaire de l'issue #8
  E = [-2*eps1*omega1,0;0,-2*eps2*omega2];
  F = [-omega1^2,0;0,-omega2^2];
  G = [-ddlambda;zeros(1,size(tt,2))];
  
  %on inverse la matrice en dehors de la boucle pour ne pas gaspiller du temps de calcul 
  A = [dt^2*F/4-eye(2) , dt^2*E/4   ;   dt*F/2 , dt*E/2-eye(2)];
  A = inv(A);


  for i=1:size(tt)-1%boucle d'intégration temporelle
  
    Y(:,i+1)=A*[-Y(1:2,i)-dt*Y(3:4,i)-dt^2*E/4*Y(3:4,i)-dt^2*F/4*Y(1:2,i)-dt^2/4*G(:,i)-dt^2/4*G(:,i+1);-Y(3:4,i)-dt*E/2*Y(3:4,i)-dt*F/2*Y(1:2,i)-dt*G(:,i)/2-dt/2*G(:,i+1)];%cf dernière équation dans le commentaire de l'issue#8

  end
  %on met les résultats dans les variables à renvoyer
  Xt=[Y(1,:);Y(2,:)];
  dXt=[Y(3,:);Y(4,:)];
  plot(tt,Xt)

