[tt,Xt,dXt]=verlet(X0,dX0,t_init,dt,t_tot);  
  global k m Ma g l Io omega1 omega2 omega eps1 eps2 alpha beta lamda0
  
  tt=[t_init:dt:t_tot];%vecteur tt à retourner 
 
  Xt = zeros(size(X0),size(tt))%déclaration du vecteur contenant les valeur de X au cours du temps
  dXt = zeros(size(X0),size(tt))%idem pour dXt
  ddXt = zeros(size(X0),size(tt))%idem pour ddXt

  %initialisation
  Xt(:,1)= X0;
  dXt(:,1)=dX0;




  for i=2:size(tt)%boucle d'intégration temporelle
  
  %partie à remplir avec le schéma de verlet

  end

end
