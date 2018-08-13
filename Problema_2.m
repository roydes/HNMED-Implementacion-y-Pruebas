%Caso de estudio 2: seguimiento de una trayectoria no alineada, con sincronización previa - CE2
function [FO,SVR]=funcionMecanismoCuatroBarras2(p)
format long;
r1=p(1);r2=p(2);r3=p(3);r4=p(4);rcx=p(5);rcy=p(6);teta0=0;x0=0;y0=0;
angulos=[((2*pi)/12) ((3*pi)/12) ((4*pi)/12) ((5*pi)/12) ((6*pi)/12)];
%teta12=((2*pi)/12);teta22=((3*pi)/12);teta32=((4*pi)/12);teta42=((5*pi)/12);teta52=((6*pi)/12);
puntos=[3 3;2.759 3.363;2.372 3.663;1.890 3.862;1.355 3.943];

SVR=0;

g1=r1+r2-r3-r4;
if(g1>0)
    SVR=SVR+g1;
    %display(g1);
end
g2=r2-r3;
if(g2>0)
    SVR=SVR+g2;
    %display(g2);
end
g3=r3-r4;
if(g3>0)
    SVR=SVR+g3;
    %display(g3);
end
g4=r4-r1;
if(g4>0)
    SVR=SVR+g4;
    %display(g4);
end

%%if(SVR==0)
FO=0;
for i=1:size(puntos,1)
    A1=2*r3*((r2*cos(angulos(i)))-(r1*cos(0)));
    B1=2*r3*((r2*sin(angulos(i)))-(r1*sin(0)));
    C1=(r1^2)+(r2^2)+(r3^2)-(r4^2)-(2*r1*r2*cos(0-angulos(i)));
    teta3=2*atan((-B1+sqrt((B1^2)+(A1^2)-(C1^2)))/(C1-A1));
    
    if(isreal(teta3)==1)
        Cxr=(r2*cos(angulos(i)))+(rcx*cos(teta3))-(rcy*sin(teta3));
        Cyr=(r2*sin(angulos(i)))+(rcx*sin(teta3))+(rcy*cos(teta3));
        
        C=[cos(teta0) -sin(teta0);sin(teta0) cos(teta0)]*[Cxr; Cyr]+[x0; y0];
        C=C';
        diferencia=(puntos(i,1)-C(1,1))^2+(puntos(i,2)-C(1,2))^2;
        FO=FO+diferencia;        
    else
        FO=1000;
        %SVR=SVR+1000;
        break;
    end
end


