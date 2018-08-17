%Caso de estudio 3: generación de movimiento delimitado por un conjunto de pares de puntos - CE3
function [FO,SVR]=funcionMecanismoCuatroBarras3(p)
format long;
r1=p(1);r2=p(2);r3=p(3);r4=p(4);rcx=p(5);rcy=p(6);teta0=p(7);x0=p(8);y0=p(9);
teta12=p(10);teta22=p(11);teta32=p(12);teta42=p(13);teta52=p(14);teta62=p(15);
teta72=p(16);teta82=p(17);teta92=p(18);teta102=p(19);
puntos1=[1.768 2.3311;1.947 2.6271;1.595 2.7951;1.019 2.7241;0.479 2.4281;0.126 2.0521;-0.001 1.720;0.103 1.514;0.442 1.549;1.055 1.905];
puntos2=[1.9592 2.44973;2.168 2.675;1.821 2.804;1.244 2.720;0.705 2.437;0.346 2.104;0.195 1.833;0.356 1.680;0.558 1.742;1.186 2.088];

SVR=0;

g1=r1+r2-r3-r4;
if(g1>0)
    SVR=SVR+g1;
end
g2=r2-r3;
if(g2>0)
    SVR=SVR+g2;
end
g3=r3-r4;
if(g3>0)
    SVR=SVR+g3;
end
g4=r4-r1;
if(g4>0)
    SVR=SVR+g4;
end
g5=teta12-teta22;
if(g5>0)
    SVR=SVR+g5;
end
g6=teta22-teta32;
if(g6>0)
    SVR=SVR+g6;
end
g7=teta32-teta42;
if(g7>0)
    SVR=SVR+g7;
end
g8=teta42-teta52;
if(g8>0)
    SVR=SVR+g8;
end
g9=teta52-teta62;
if(g9>0)
    SVR=SVR+g9;
end
g10=teta62-teta72;
if(g10>0)
    SVR=SVR+g10;
end
g11=teta72-teta82;
if(g11>0)
    SVR=SVR+g11;
end
g12=teta82-teta92;
if(g12>0)
    SVR=SVR+g12;
end
g13=teta92-teta102;
if(g13>0)
    SVR=SVR+g13;
end
%Los angulos de los puntos se encuentran en las posiciones 10-15
FO=0;
j=1;
for i=10:size(p,2)
    A1=2*r3*((r2*cos(p(i)))-(r1*cos(0)));
    B1=2*r3*((r2*sin(p(i)))-(r1*sin(0)));
    C1=(r1^2)+(r2^2)+(r3^2)-(r4^2)-(2*r1*r2*cos(0-p(i)));
    teta3=2*atan((-B1+sqrt((B1^2)+(A1^2)-(C1^2)))/(C1-A1));
    if(isreal(teta3)==1)
        Cxr=(r2*cos(p(i)))+(rcx*cos(teta3))-(rcy*sin(teta3));
        Cyr=(r2*sin(p(i)))+(rcx*sin(teta3))+(rcy*cos(teta3));
        
        C=[cos(teta0) -sin(teta0);sin(teta0) cos(teta0)]*[Cxr; Cyr]+[x0; y0];
        C=C';
        diferencia1=(puntos1(j,1)-C(1,1))^2+(puntos1(j,2)-C(1,2))^2;
        diferencia2=(puntos2(j,1)-C(1,1))^2+(puntos2(j,2)-C(1,2))^2;
        FO=FO+(diferencia1+diferencia2);
        j=j+1;
    else
        FO=1000;
        %SVR=SVR+1000;
        break;
    end
end
