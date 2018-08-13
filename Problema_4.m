function [fit, sumr]=funcionGriperTresDedos(p)
format long;

%Esta función evalua la funcion objetivo del problema de optimización
%del mecanismo de seis barras. Dicha funcion es el error cuadratico
%acumulado mas la normalizacion de las longitudes de los eslabones y es el valor que regresa esta rutina
punpres=[-10 160;10 170;40 165]; 
peso1=1;
peso2=0;

ctr=[];
ctr(1)=p(10)+p(9)-p(8)-p(7);
ctr(2)=p(7)-p(10);
ctr(3)=p(8)-p(10);
ctr(4)=p(9)-p(10);
ctr(5)=secuencia(p(4:6));
ctr(6)=p(11)-p(9);
ctr(7)=p(2)-(p(3)/3);
ctr(8)=(p(3)/5)-p(2);
%La entrada de la función es el conjunto de las 14 variables de diseño, las
%5 restricciones y los puntos de precision.
% y tiene el siguiente orden: p=[r1, r2, r3, r41,r42,r43,r0,r2p,r5,r6,rex, rey,alfa,teta0]
% ctr=[4-grashof, 1-secuencia], punpres=[(c1x, c1y), (c2x, c2y),(c3x, c3y)]

% Edgar Alfredo Portilla Flores
%7 de octubre de 2014

%Esta es la variable donde se almacenan las violaciones de las restricciones
contx=0;

%Se forma un vector de las dimensiones de las barras
barras=p(1:3);
barras(4:9)=p(7:12);

%Se evaluan las violaciones de las restricciones
if ctr(1)>0
    contx=contx+ctr(1); 
end
if ctr(2)>0
    contx=contx+ctr(2);
end
if ctr(3)>0
    contx=contx+ctr(3);
end
if ctr(4)>0
    contx=contx+ctr(4);
end
if ctr(5)>0
    contx=contx+ctr(5);
end
if ctr(6)>0
    contx=contx+ctr(6);
end
if ctr(7)>0
    contx=contx+ctr(7);
end
if ctr(8)>0
    contx=contx+ctr(8);
end

    error=0;
    for i=1:3
        if p(3+i)==0
           teta2=acos((p(1)^2+p(2)^2-p(3)^2)/(2*p(1)*p(2)));
           phi=acos((p(2)^2+p(3)^2-p(1)^2)/(2*p(2)*p(3)));
           teta3=pi+teta2+phi;
        else
            L=sqrt(p(1)^2+p(3+i)^2);
            gamma=acos((p(2)^2+L^2-p(3)^2)/(2*p(2)*L));
            beta=atan(p(3+i)/p(1));
            teta2=beta+gamma;
            phi=acos((p(2)^2+p(3)^2-L^2)/(2*p(2)*p(3)));
            teta3=pi+teta2+phi;
        end
        teta2p=teta2+p(13);
        A=2*p(9)*(p(8)*cos(teta2p)-p(7)*cos(p(14)));
        B=2*p(9)*(p(8)*sin(teta2p)-p(7)*sin(p(14)));
        C=p(7)^2+p(8)^2+p(9)^2-p(10)^2-2*p(7)*p(8)*cos(p(14)-teta2p);
        teta5=2*atan((-B-sqrt(B^2+A^2-C^2))/(C-A));
        
        D=2*p(10)*(p(7)*cos(p(14))-p(8)*cos(teta2p));
        E=2*p(10)*(p(7)*sin(p(14))-p(8)*sin(teta2p));
        F=p(7)^2+p(8)^2+p(10)^2-p(9)^2-2*p(7)*p(8)*cos(p(14)-teta2p);
        teta6=2*atan((-E+sqrt(E^2+D^2-F^2))/(F-D));
        
%         r7=sqrt(p(11)^2+p(12)^2); %ASÍ LO HIZO PEDRO
%         delta=atan(p(12)/p(11));
%         teta7=delta+teta5-pi;
%         
%         Ex=p(7)*cos(p(14))+p(10)*cos(teta6)+r7*cos(teta7);
%         Ey=p(7)*sin(p(14))+p(10)*sin(teta6)+r7*sin(teta7);

        tetaey = teta5-pi/2; %ASÍ LO HIZO EDITH
        teta2p = teta2+p(13);
        Ex=p(8)*cos(teta2p)+p(11)*cos(teta5)+p(12)*cos(tetaey);
        Ey=p(8)*sin(teta2p)+p(11)*sin(teta5)+p(12)*sin(tetaey);
        
        error=error+(punpres(i,1)-Ex)^2+(punpres(i,2)-Ey)^2;
    end
    verifica=isreal(error); %verifica si produce imaginarios la funcion objetivo.
    if verifica==0
        fit=1000;
        contx=contx+1000;
    else
		if(peso2==0)%No tiene por qué hacerce el cálculo de normalización de barras, cuando el peso dos es 0 ya que daría 0
			fit=error;
		else
			longs=0;
			difs=0;
			for i=1:9
				for j=1:9
				   difs=difs+(barras(i)-barras(j))^2; 
				end
			end
			for i=1:9
				longs=longs+barras(i)^2;
			end
			
			fit=peso1*error+peso2*sqrt(difs/longs);
		end
    end
    sumr=contx;
function [fit]=secuencia(KC)
%Esta función evalua la secuencia de la posicion de la corredera del
%mecanismo de entrada.Verifica que dicha secuencia sea ascendente o descendente.
%El valor que regresa la funcion es cero si dicha secuencia es valida y
%diferente de cero si no lo es. En este ultimo caso el valor representa la
%suma de violacion de restricciones de la secuencia. Un valor mayor indica
%que la secuencia es menos ordenada en ascendencia.

%La entrada de la función es el conjunto de valores de la posicion de la
%corredera y tiene el siguiente orden: [r41,r42,r43]


%Esta es la variable donde se almacenan las violaciones de las restricciones
contx=0;


if (KC(1)<KC(2) && KC(2)<KC(3))
    fit=0;
else
    if (KC(1)>KC(2) && KC(2)>KC(3))
        fit=0;
    else
       contx=abs(KC(1)-KC(2))+abs(KC(1)-KC(3))+abs(KC(2)-KC(3));
       fit=contx;
    end
end

    
 
    



    
    



   