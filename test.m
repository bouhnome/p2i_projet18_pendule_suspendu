%taille de la figure 
size = 1;
%position générale de la figure
offsetx = 0;
offsety = 0;



%%%%=================================tracé du solide 1==============%
%position de S1 relative à la figure
offsetS1x=0;
offsetS1y=0; 
%nous dessinons seulement la partie haute et droite de S1 et nous
%utilisons les symmetries
XS1= [540 740 740 800 800 840 840];
YS1= [220 220 160 160 220 220 160];


%%%%=================================tracé du solide 0*==============%
%position de S0* relative à la figure
offsetS0stx=0;
offsetS0sty=0; 
%nous dessinons seulement la partie droite de S0* et nous
%utilisons les symmetries
XS0st= [544 824 824 824 904 904 824 824 824 544];
YS0st= [648 648 343 648 648 683 683 725 683 683];



%=======================================Le plot====================% 
for(j=1:2000)  
    
offsetS1x=0;
offsetS1y=10*sin(j);
offsetS0stx=0;
offsetS0sty=0; 

a= 32*(10*size-10)+320;
xS1 = (XS1-offsetx-offsetS1x)*size;
yS1 = (YS1-offsety-offsetS1y)*size; 
xS1oppo= 2*(XS1(1)-offsetx-offsetS1x)*size-(XS1-offsetx-offsetS1x)*size;
yS1oppo= a-(YS1+offsety+offsetS1y)*size;




xS0st=(XS0st-offsetx-offsetS0stx)*size;
yS0st= -YS0st*size+343*size-offsety*size-offsetS0sty*size;
yS0st(3)=yS1oppo(5);
xS0stoppo = 2*(XS0st(1)-offsetx-offsetS0stx)*size-size*(XS0st-offsetx-offsetS0stx);

plot(xS1,yS1,'-b',  xS1oppo,yS1,'-b',   xS1,yS1oppo,'-b',     xS1oppo,yS1oppo,'-b',xS0st,yS0st,xS0stoppo,yS0st)

axis([0 1000 -400 400])
pause(0.1)
end
