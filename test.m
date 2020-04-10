%taille de la figure 
size = 1.5;
%position générale de la figure
offsetx = 200;
offsety = -400;


%%%%=================================tracé du solide 1==============%
%position de S1 relative à la figure
offsetS1x=0;
offsetS1y=0; 
%nous dessinons seulement la partie haute et droite de S1 et nous
%utilisons les symmetries
xS1= [540 740 740 800 800 840 840];
yS1= [220 220 160 160 220 220 160];

a= 32*(10*size-10)+320;

S1 = plot((xS1-offsetx-offsetS1x)*size,(yS1-offsety-offsetS1y)*size,'-b',  2*(xS1(1)-offsetx-offsetS1x)*size-(xS1-offsetx-offsetS1x)*size,(yS1-offsety-offsetS1y)*size,'-b',   (xS1-offsetx-offsetS1x)*size,a-(yS1+offsety+offsetS1y)*size,'-b',     2*(xS1(1)-offsetx-offsetS1x)*size-(xS1-offsetx-offsetS1x)*size,a-(yS1+offsety+offsetS1y)*size,'-b')
for j = 1:2000
    
    %%%%on dessine S1%%%%
    offsetS1x=0;
    offsetS1y=10*sin(j);%donner ici le mouvement en z
    S1 =plot((xS1-offsetx-offsetS1x)*size,(yS1-offsety-offsetS1y)*size,'-b',  2*(xS1(1)-offsetx-offsetS1x)*size-(xS1-offsetx-offsetS1x)*size,(yS1-offsety-offsetS1y)*size,'-b',   (xS1-offsetx-offsetS1x)*size,a-(yS1+offsety+offsetS1y)*size,'-b',     2*(xS1(1)-offsetx-offsetS1x)*size-(xS1-offsetx-offsetS1x)*size,a-(yS1+offsety+offsetS1y)*size,'-b')
    
    
    axis([0 1000 0 1000])
    pause(0.01)
end

