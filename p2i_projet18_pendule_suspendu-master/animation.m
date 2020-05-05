
function animation(tt,Xt,lambda) 

Size = 0.01;
%position générale de la figure
offsetx = 0;
offsety = -400;



%%%%=================================tracé du solide 1==============%

%nous dessinons seulement la partie haute et droite de S1 et nous
%utilisons les symmetries
XS1= [540 740 740 800 800 840 840];
YS1= [220 220 160 160 220 220 160];


%%%%=================================tracé du solide 0*==============%

%nous dessinons seulement la partie droite de S0* et nous
%utilisons les symmetries
XS0st= [544 824 824 824 904 904 824 824 824 544];
YS0st= [648 648 343 648 648 683 683 725 683 683];


%%%%===============================tracé du solide 0===================%s

%nous dessinons seulement la partie droite de S0 et nous
%utilisons les symmetries
XS0=[545 788 788 810 828  828 788 788]+14;
XS0(1)=XS0(1)-14;
YS0=[739 739 702 702 702  780 780 739];



%=======================================Le plot====================% 
for j=1:size(tt,2)-1
    

    
%position de S1 relative à la figure
offsetS1x=0;
offsetS1y=-Xt(1,j)/Size; %mettre ici le mouvement en z
theta = Xt(2,j);%mettre le mouvement en theta ici
%position de S0* relative à la figure
offsetS0stx=0;
offsetS0sty=-lambda(j)/Size; %mettre ici le mouvement en lambda

%position de S0 relative à la figure
offsetS0x=0;
offsetS0y=0; 

a= 32*(10*Size-10)+320;
xS1 = (XS1-offsetx-offsetS1x)*Size;
yS1 = (YS1-offsety-offsetS1y)*Size; 
xS1oppo= 2*(XS1(1)-offsetx-offsetS1x)*Size-(XS1-offsetx-offsetS1x)*Size;
yS1oppo= a-(YS1+offsety+offsetS1y)*Size;




xS0st=(XS0st-offsetx-offsetS0stx)*Size;
yS0st= -YS0st*Size+343*Size-offsety*Size-offsetS0sty*Size;
yS0st(3)=yS1oppo(5);
xS0stoppo = 2*(XS0st(1)-offsetx-offsetS0stx)*Size-Size*(XS0st-offsetx-offsetS0stx);

xS0=(XS0-offsetx-offsetS0x)*Size;
yS0= -YS0*Size-offsety*Size-offsetS0y*Size+320*Size;
yS0st(8)=yS0(3);
xS0oppo = 2*(XS0(1)-offsetx-offsetS0x)*Size-Size*(XS0-offsetx-offsetS0x);



figure(5)
plot(xS1,yS1,'-b',  xS1oppo,yS1,'-b',   xS1,yS1oppo,'-b',     xS1oppo,yS1oppo,'-b',xS0st,yS0st,'-g',xS0stoppo,yS0st,'-g', xS0,yS0,'-bla',xS0oppo,yS0,'-bla' , [xS0st(2),xS0st(2)],[2*yS1(5)-yS1(5),yS0st(2)+600*Size],'-g',[xS0stoppo(2),xS0stoppo(2)],[2*yS1(5)-yS1(5),yS0st(2)+600*Size],'-g',[xS0st(2),xS0st(2)],[yS0(6),yS0st(2)-200*Size],'-g',[xS0stoppo(2),xS0stoppo(2)],[yS0(6),yS0st(2)-200*Size],'-g',[xS1(1) xS1(1)+Size*324*sin(theta)],[yS1(3) yS1(3)-324*Size*cos(theta)],'-r')
axis([100-offsetx 1000-offsetx -600-offsety 400-offsety]*Size)

pause(0.01)
end
end
