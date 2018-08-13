%% Problema de optimización para un mecanismo de 6 barras
% * Esta función calula el valor de la funcion objetivo definida en el
% problema y la suma de violación de restricciones
% * La función tiene un único parámetro, el cual es el vector de variables
% de diseño definido para el problema

function [Fdp,sVR]=PO_SIXBAR_NoG(p)

Fdp=funcionobjetivo(0,p);
r=restricciones(0,0,0,p);
sVR=SVR(0,0,28,r,1);

%% Evalua la funcion objetivo

function f=funcionobjetivo(nvariables,x)
Dx=[10.32 10.59 11.82 14.69 19.06 24.35 30.61 40.73 44.00 45.63];
Dy=[92.98 91.17 86.30 79.96 74.04 69.38 65.85 63.24 63.10 63.09];

NDP=size(Dx,2);
R=size(x,1);
wp=1;
for r=1:R
    sum=0;
    for n=1:NDP
        [OUTvar]=modelSIXBARKNEE(x(r,1:13),x(r,13+n));
        sum=sum+wp*((Dx(n)-OUTvar(5))^2+(Dy(n)-OUTvar(6))^2);
    end
    f=sum; 
end


%% Evalua las restricciones

function c=restricciones(nvariables,fsob,epsilon,x)

%Orden
c(:,nvariables+fsob+1)=x(:,14)-x(:,15);
c(:,nvariables+fsob+2)=x(:,15)-x(:,16);
c(:,nvariables+fsob+3)=x(:,16)-x(:,17);
c(:,nvariables+fsob+4)=x(:,17)-x(:,18);
c(:,nvariables+fsob+5)=x(:,18)-x(:,19);
c(:,nvariables+fsob+6)=x(:,19)-x(:,20);
c(:,nvariables+fsob+7)=x(:,20)-x(:,21);
c(:,nvariables+fsob+8)=x(:,21)-x(:,22);
c(:,nvariables+fsob+9)=x(:,22)-x(:,23);

R=size(x,1);
for r=1:R
    r1=x(r,1);
    r2=x(r,2);
    r3=x(r,3);
    r4=x(r,4);
    r5=x(r,5);
    r6=x(r,6);
    r7=x(r,7);
    r8=x(r,8);
    r0=x(r,9);
    theta0=x(r,10);
    theta1=x(r,11);
    theta5=x(r,12);
    theta8=x(r,13);

    for v1=1:10
        theta2=x(r,13+v1);

        A=2*r4*(r1*cos(theta1)-r2*cos(theta2));
        B=2*r4*(r1*sin(theta1)-r2*sin(theta2));
        C=r1^(2)+r2^(2)-r3^(2)+r4^(2)-2*r1*r2*cos(theta1-theta2);
        D=2*r3*(r2*cos(theta2)-r1*cos(theta1));
        E=2*r3*(r2*sin(theta2)-r1*sin(theta1));
        F=r1^(2)+r2^(2)+r3^(2)-r4^(2)-2*r1*r2*cos(theta2-theta1);
   
        %x(r,nvariables+fsob+7+2*v1)=-(D.^(2)+E.^(2)-F.^(2));

        theta3=2*atan((-E+sqrt(1*(D.^(2)+E.^(2)-F.^(2))))./(F-D));
        theta4=2*atan((-B-sqrt(1*(A.^(2)+B.^(2)-C.^(2))))./(C-A));
        DD=2*r6*(r5*cos(theta3+theta5)-r8*cos(theta4-theta8));
        EE=2*r6*(r5*sin(theta3+theta5)-r8*sin(theta4-theta8));
        FF=r8^(2)+r5^(2)+r6^(2)-r7^(2)-2*r8*r5*cos((theta3+theta5)-(theta4-theta8));
        
        if v1 == 10
            if (D^(2)+E^(2)-F^(2)) <= 0
                c(r,nvariables+fsob+28)=(D^(2)+E^(2)-F^(2));
            elseif (DD^(2)+EE^(2)-FF^(2)) <=0
                c(r,nvariables+fsob+28)=(DD^(2)+EE^(2)-FF^(2));
            else
                c(r,nvariables+fsob+28)=(D^(2)+E^(2)-F^(2))+(DD^(2)+EE^(2)-FF^(2));
            end
        else
            c(r,nvariables+fsob+8+2*v1)=-(D.^(2)+E.^(2)-F.^(2));
            c(r,nvariables+fsob+9+2*v1)=-(DD.^(2)+EE.^(2)-FF.^(2));
        end
    end
end


%% Evalua las restricciones y calcula SVR para cada individuo

function s=SVR(nv,nf,ng,x,ni)

for nia=1:ni
    sumaSVR=0;
    for v=(nv+nf+1):(nv+nf+ng)
        %% Criterio de suma
        % * Si g_{nia} es positivo, es decir que esta fuera de la
        % restriccion actual, entonces se incrementa SVR con el valor de
        % g_{nia}
        if x(nia,v)>0
            sumaSVR=sumaSVR+x(nia,v);
        end
    end
    %% Asigna el resultado a cada individuo en la posicion nv+nf+ng+1
%     x(nia,nv+nf+ng+1)=sumaSVR;
    s=sumaSVR;
end

function [OUTvar]=modelSIXBARKNEE(DESIGNvar,INvar)

r1=DESIGNvar(1);
r2=DESIGNvar(2);
r3=DESIGNvar(3);
r4=DESIGNvar(4);
r5=DESIGNvar(5);
r6=DESIGNvar(6);
r7=DESIGNvar(7);
r8=DESIGNvar(8);
r0=DESIGNvar(9);
theta0=DESIGNvar(10);
theta1=DESIGNvar(11);
theta5=DESIGNvar(12);
theta8=DESIGNvar(13);


theta2=INvar(:,1);
%% ========================================================================
%           ANALISIS DE POSICION
%==========================================================================
A=2*r4*(r1*cos(theta1)-r2*cos(theta2));
B=2*r4*(r1*sin(theta1)-r2*sin(theta2));
C=r1^(2)+r2^(2)-r3^(2)+r4^(2)-2*r1*r2*cos(theta1-theta2);
D=2*r3*(r2*cos(theta2)-r1*cos(theta1));
E=2*r3*(r2*sin(theta2)-r1*sin(theta1));
F=r1^(2)+r2^(2)+r3^(2)-r4^(2)-2*r1*r2*cos(theta2-theta1);
% theta3=2*atan((-E+sqrt(1*(D.^(2)+E.^(2)-F.^(2))))./(F-D));
% theta4=2*atan((-B-sqrt(1*(A.^(2)+B.^(2)-C.^(2))))./(C-A));
theta3=real(2*atan((-E+sqrt(1*(D.^(2)+E.^(2)-F.^(2))))./(F-D)));
theta4=real(2*atan((-B-sqrt(1*(A.^(2)+B.^(2)-C.^(2))))./(C-A)));
%% ========================================================================
%           ANALISIS DE POSICION
%==========================================================================
% % theta1-->theta8
% % theta2-->theta5
% % theta3-->theta6
% % theta4-->theta7

AA=2*r7*(r8*cos(theta4-theta8)-r5*cos(theta3+theta5));
BB=2*r7*(r8*sin(theta4-theta8)-r5*sin(theta3+theta5));
CC=r8^(2)+r5^(2)-r6^(2)+r7^(2)-2*r8*r5*cos((theta4-theta8)-(theta3+theta5));
DD=2*r6*(r5*cos(theta3+theta5)-r8*cos(theta4-theta8));
EE=2*r6*(r5*sin(theta3+theta5)-r8*sin(theta4-theta8));
FF=r8^(2)+r5^(2)+r6^(2)-r7^(2)-2*r8*r5*cos((theta3+theta5)-(theta4-theta8));
% theta6=2*atan((-EE+sqrt(1*(DD.^(2)+EE.^(2)-FF.^(2))))./(FF-DD));
% theta7=2*atan((-BB-sqrt(1*(AA.^(2)+BB.^(2)-CC.^(2))))./(CC-AA));
theta6=real(2*atan((-EE+sqrt(1*(DD.^(2)+EE.^(2)-FF.^(2))))./(FF-DD)));
theta7=real(2*atan((-BB+sqrt(1*(AA.^(2)+BB.^(2)-CC.^(2))))./(CC-AA)));
%==========================================================================+
x=r0*cos(theta0)+r2*cos(theta2)+r3*cos(theta3)+r8*cos(theta4-theta8)+r7*cos(theta7);
y=r0*sin(theta0)+r2*sin(theta2)+r3*sin(theta3)+r8*sin(theta4-theta8)+r7*sin(theta7);


OUTvar=[theta3,theta4,theta6,theta7,x,y];



