function [tt,Xt,dXt]=newmark(X0,dX0,ddX0,t_init,dt,t_tot)

global k m Ma g l Io omega1 omega2 omega eps1 eps2 alpha beta lamda0
global M C K %matrices-coefficients du formalisme de newmark 
