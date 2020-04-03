[tt,Xt,dXt]=newmarklin(X0,dX0,t_init,dt,t_tot);  
  global k m Ma g l Io omega1 omega2 omega eps1 eps2 alpha beta lamda0
  
  tt=[t_init:dt:t_tot];%vecteur tt à retourner 
  Y = zeros(2*size(X0),size(tt));%création de la matrice utilisé dans les calculs cf commentaire de l'issue #8
  ddlambda = -omega^2*lamda0*sin(omega*tt);%terme d'exitation 

  %initialisation
  Y(1,1)=X0(1);
  Y(2,1)=X0(2);
  Y(3,1)=dX0(1);
  Y(4,1)=dX0(2);

  %on  pose les variables qui permettent de simplifer les expressions cf commentaire de l'issue #8
  E = [-2*eps1*omega1,0;0,-2*eps2*omega2];
  F = [-omega1^2,0;0,-omega2^2];
  G = [-ddlambda;0];
  
  %on inverse la matrice en dehors de la boucle pour ne pas gaspiller du temps de calcul 
  A = [dt^2*F/4 -eyes(4) , dt^2*E/4   ;   dt*F/2 , dt*E/2 - eyes(4)];
  A = inv(A);


  for i=1:size(tt)-1%boucle d'intégration temporelle
  
    Y(:,i+1)=A*[  -Y(1:2,i)-dt*Y(3:4,i)-dt^2*E/4*Y(3:4,i)-dt^2*F/4*Y(1:2,i)-dt^2/4*G(:,i)-dt^2/4*G(:,i+1) ;  -Y(3:4,i) -dt*E/2*Y(3:4,i) - dt*F/2*Y(1:2,i) -dt*G(:,i)/2 -dt/2*G(:,i+1) ];%cf dernière équation dans le commentaire de l'issue#8

  end
  %on met les résultats dans les variables à renvoyer
  Xt=[Y(1,:);Y(2,:)];
  dXt=[Y(3,:);Y(4,:)];

end
