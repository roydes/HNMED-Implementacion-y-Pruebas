%Caso de estudio 1: seguimiento de una trayectoria lineal vertical, sin sincronización previa - (CE1)
function [FO,SVR]=Problema_1(p)
format long;
r1=p(1);r2=p(2);r3=p(3);r4=p(4);rcx=p(5);rcy=p(6);teta0=p(7);x0=p(8);y0=p(9);
teta12=p(10);teta22=p(11);teta32=p(12);teta42=p(13);teta52=p(14);teta62=p(15);
puntos=[20 20;20 25;20 30;20 35;20 40;20 45];

% if((r2<r3 && r3<r4 && r4<r1) && (teta12<teta22 && teta22<teta32 && teta32<teta42 && teta42<teta52 && teta52<teta62))
SVR=0;
%evaluar restricciones condicion de Grashof

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
%evaluacion de las restricciones de la secuencia de angulos
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

% if(SVR==0)
%Los angulos de los puntos se encuentran en las posiciones 10-15
sumatoria=0;
j=1;
for i=10:size(p,2)
    A1=2*r3*((r2*cos(p(i)))-(r1*cos(0)));
    B1=2*r3*((r2*sin(p(i)))-(r1*sin(0)));
    C1=(r1^2)+(r2^2)+(r3^2)-(r4^2)-(2*r1*r2*cos(p(i)-0));
    teta3=2*atan((-B1+sqrt((B1^2)+(A1^2)-(C1^2)))/(C1-A1));
    if(isreal(teta3)==1)
        Cxr=(r2*cos(p(i)))+(rcx*cos(teta3))-(rcy*sin(teta3));
        Cyr=(r2*sin(p(i)))+(rcx*sin(teta3))+(rcy*cos(teta3));
        
        C=[cos(teta0) -sin(teta0);sin(teta0) cos(teta0)]*[Cxr; Cyr]+[x0; y0];
        C=C';
        diferencia=(puntos(j,1)-C(1,1))^2+(puntos(j,2)-C(1,2))^2;
        sumatoria=sumatoria+diferencia;
        j=j+1;
    else
        sumatoria=10000;
        %SVR=SVR+1000;
        break;
    end
end
FO=sumatoria;

