function [Ahorro_Total, E_F]=MR_DE_rand_1_bin (E,NP,e)

format long
%sFO=zeros(ejecuciones,2); %FO y Número de corrida
    E_F=[];
    rng('shuffle', 'twister');
    %% INICIALIZACIÓN DE VARIABLES A TRAVÉS DE LECTURA EN EXCEL
    P1=zeros(24,1);
    P2=zeros(24,1);
    P3=zeros(24,1);
    P4=zeros(24,1);
    FO=zeros(24,1);
	SVR=zeros(24,1);
    Pcarga=zeros(24,1);
    %%P1--> Diesel P2-->Battery, P3-->PV, P4--> wind; Potencia[W], SOC
    celdas_variables24=xlsread('Potencias_redes.xls','C5:G28');
    PL_VEC=celdas_variables24(:,1);                 %X de potencias de la carga a suplir
    SOC_VEC=celdas_variables24(:,3);              %X de SOC de la bateria
    P1_VEC=celdas_variables24(:,2);                 %X potencias Diesel
    P3_VEC=celdas_variables24(:,4);                 %X potencias de PV
    P4_VEC=celdas_variables24(:,5);                 %X potencias viento
    
    %%parámetros bateria para restricción SOC
    n_B=0.85; % Round trip efficiency
    n_c=0.85;
    n_d=1;
    
    % celdas_variables3=xlsread('Potencias_redes.xls','K5:M5');
    % Bc_max=celdas_variables3(3); %Capacidad máxima
    % alpha_bateria=n_B/Bc_max;
    % alpha_c=n_c/Bc_max;
    % alpha_d=n_d/Bc_max;
    % SOC_max=celdas_variables3(2);    %Estado de carga mínimo de la bateria
    % SOC_min=celdas_variables3(1);    %Estado de carga máximo de la bateria
    
    Bc_max=2000; %Capacidad máxima
    alpha_bateria=n_B/Bc_max;
    alpha_c=n_c/Bc_max;
    alpha_d=n_d/Bc_max;
    SOC_max=0.95;    %Estado de carga mínimo de la bateria
    SOC_min=0.4;    %Estado de carga máximo de la bateria
    
    %%parámetros de inicialización para cumplir con la sumatoria de la restricción de SOC para las 24 horas
    sumP2=0;
    sumP3=0;
    sumP4=0;
    HORAS=24;
    %%Parámetros generales
    var=4; %Cantidad de variables de diseño.
    %% INICIA LA EJECUCIÓN DE LOS 24 CASOS DE OPTIMIZACIÓN
    %%%%%%%%%%%% CREAR DIRECTORIOS NECESARIOS PARA GUARDAR ARCHIVOS%%%%%%%%%%%%%%%%%%%%%
    [stat,struc] = fileattrib;
    PathCurrent = struc.Name;
    GlobalFolderName='/Ejecuciones_Micro_Red';
    mkdir(PathCurrent, GlobalFolderName);
    FolderName='/ED';
    PathFolder = [PathCurrent GlobalFolderName FolderName];
    mkdir([PathCurrent GlobalFolderName], FolderName);  
    
    for hora=1:HORAS    %Cada hora representa una problema de optimización diferente      

        
%         display(['HORA ',num2str(hora)]);
        %%Inicialización de variables y cotas necesarias para cada hora
        PL=PL_VEC(hora);                    %Potencia de la carga en la hora
        SOC_0=SOC_VEC(hora);                %SOC actual en la hora
        P1_max=P1_VEC(hora);                %Potencia Diesel
        P2_max=SOC_0*Bc_max-SOC_min*Bc_max; %Potencia de la bateria que depende de capacidad máxima de la misma y el SOC(0)
        P3_max=P3_VEC(hora);                %Potencia de PV
        P4_max=P4_VEC(hora);                %Potencia wind
        
        P1_min=0;
        P2_min=0;
        %P2_min=SOC_min*Bc_max;
        P3_min=0;
        P4_min=0;
        
        cotas=[P1_min P1_max;P2_min P2_max;P3_min P3_max;P4_min P4_max];
        
        %% INICIA EL ALGORITMO ED  
        
        %%Generando la población inicial
        X=zeros(NP,size(cotas,1));
        for j=1:size(cotas,1)
            X(:,j)=(cotas(j,2)-cotas(j,1)).*rand(NP,1) + cotas(j,1);
        end
        
        %Evaluando la función objetivo y las restricciones de la población inicial
        
        %Restricción de desigualdades distintas a máximos y mínimos--> ej:SOC
        %A_NNO=[0 alpha_d -alpha_c -alpha_c];
        b_NNO=-SOC_min+SOC_0-alpha_d*sumP2+alpha_c*sumP3+alpha_c*sumP4;
        %A_DO=[0 -alpha_d alpha_c alpha_c];
        b_DO=SOC_max-SOC_0+alpha_d*sumP2-alpha_c*sumP3-alpha_c*sumP4;
        
        for i=1:NP
            %Función objetivo
            [X(i,var+1)]=funcionObjetivo(X(i,1:var)');
            %Calculando la restricción de igualdad
            %P1_k+P2_k+P3_k+P4_k-PL_k=0
            %utilizaremos un epsilon "cero gordo" de 1E-6
            X(i,var+2)=sum(X(i,1:var))-PL;
            %Fin de Calculando la restricción de igualdad.
            
            %Evaluación de las restricciones de desigualdad
            Lizq1=alpha_d*X(i,2)-alpha_c*X(i,3)-alpha_c*X(i,4);
            X(i,var+3)=Lizq1-b_NNO;
            Lizq2=-alpha_d*X(i,2)+alpha_c*X(i,3)+alpha_c*X(i,4);
            X(i,var+4)=Lizq2-b_DO;
            %Fin de evaluación de las restricciones de desigualdad
            
            %Calculando la SVR (Suma de violación de restricciones) de cada
            %individuo, Si la casilla tiene un cero, el individuo es factible.
            X(i,var+5)=0;
            if abs(X(i,var+2))>0.000001
                X(i,var+5)=abs(X(i,var+2));
            end
            for j=var+3:var+4
                if X(i,j)>0
                    X(i,var+5)=X(i,var+5)+X(i,j);
                end
            end
            %Fin de Calculando la SVR
        end
        %Fin de generando la primera población
        
        %Se inicializa el mejor individuo con el primer vector X, será
        %reemplazado posteriormente en la primera generación
        best=X(1,:);
        
        NFE=0;
        mejoresGen=[];
        algo=0;
        g=0;
        while NFE<E
            g=g+1; % Número de la generación
            % Cálculo del factor de cruza. Este factor es por ejecución o corrida
            a=0.8;
            b=1.0;
            
            CR=a+(b-a)*rand(1);
            
            %Cálculo del factor de escalamiento o mutación
            a=0.3;
            b=0.9;
            F=a+(b-a)*rand(1);
            
            %Generación de vectores para cruza
            for i=1:NP
                r0=randi(NP,1);
                while(r0==i)
                    r0=randi(NP,1);
                end
                r1=randi(NP,1);
                while(r1==i || r1==r0)
                    r1=randi(NP,1);
                end
                r2=randi(NP,1);
                while(r2==i || r2==r0 || r2==r1)
                    r2=randi(NP,1);
                end
                
                %Se genera el vector ruido
                ruido=X(r0,1:var)+(F*(X(r1,1:var)-X(r2,1:var)));
                
                
                %Se verifican los límites (cotas) para cada variable del individuo
                for n=1:var
                    while(ruido(n)<cotas(n,1) || ruido(n)>cotas(n,2))
                        if(ruido(n)<cotas(n,1))
                            ruido(n)=2.0*cotas(n,1)-ruido(n);
                        else
                            ruido(n)=2.0*cotas(n,2)-ruido(n);
                        end
                    end
                end
                
                %Generando el vector trial (u)
                u=zeros();
                jrand=randi(var,1);
                for j=1:var
                    randj=rand();
                    if(randj<CR || jrand==j)
                        u(j)=ruido(j);
                    else
                        u(j)=X(i,j);
                    end
                end
                %Fin de generando el vetor trial (u)
                
                %Función objetivo
                [u(var+1)]=funcionObjetivo(u');
                NFE=NFE+1;
                %Calculando la restricción de igualdad
                %P1_k+P2_k+P3_k+P4_k-PL_k=0
                %utilizaremos un epsilon "cero gordo" de 1E-6
                u(var+2)=sum(u(1:var))-PL;
                %Fin de Calculando la restricción de igualdad.
                
                %Evaluación de las restricciones de desigualdad
                Lizq1=alpha_d*u(2)-alpha_c*u(3)-alpha_c*u(4);
                u(var+3)=Lizq1-b_NNO;
                Lizq2=-alpha_d*u(2)+alpha_c*u(3)+alpha_c*u(4);
                u(var+4)=Lizq2-b_DO;
                %Fin de evaluación de las restricciones de desigualdad
                
                %Calculando la SVR (Suma de violación de restricciones) de cada
                %individuo, Si la casilla tiene un cero, el individuo es factible.
                u(var+5)=0;
                if abs(u(var+2))>0.000001
                    u(var+5)=abs(u(var+2));
                end
                for j=var+3:var+4
                    if u(j)>0
                        u(var+5)=u(var+5)+u(j);
                    end
                end
                %Fin de Calculando la SVR
                
                %Reemplazo del padre por el trial (u) siguiendo reglas de Deb
                if (X(i,var+5)>0) && (u(var+5)>0) %Si ambos violan restricciones
                    if (X(i,var+5)>=u(var+5)) %Se toma al que tiene menor SVR.
                        X(i,:)=u;
                    end
                elseif ((X(i,var+5)==0) && (u(var+5)==0)) %Si ambos no violan restricciones
                    if X(i,var+1)>=u(var+1)% Se toma al que tiene mejor F.O.
                        X(i,:)=u;
                    end
                elseif ((X(i,var+5)>0) && (u(var+5)==0)) %El padre viola y el trial(u) no viola
                    X(i,:)=u;
                end
                
                %Se reemplaza al mejor individuo con criterios de Deb
                if(best(var+5)==0 && X(i,var+5)==0) %Si ambos son factibles se selecciona el individuo con mejor valor de aptitud (función objetivo)
                    if(X(i,var+1)<=best(var+1))
                        best=X(i,:);
                    end
                elseif(best(var+5)~=0 && X(i,var+5)==0)  %Si sólo un individuo es factible, éste es seleccionado
                    best=X(i,:);
                elseif(best(var+5)~=0 && X(i,var+5)~=0) %Si ambos no son factibles, selecciona aquel con la menor suma de violaciones
                    if(X(i,var+5)<= best(var+5))
                        best=X(i,:);
                    end
                end                
                
            end
           %Si el individuo es factible se almacena en la matriz
            %mejoresGen para que pase a formar parte de la grafica de
            %convergencia
            if hora==11
               [E_F]=uptdate_E_F(NFE,E_F,best);
            end
            
        end
        x_k=best; 
        %%Termina el algoritmo DE
        
        %% Generación vectores de Potencia para gráfica
        P1(hora)=x_k(1);
        P2(hora)=x_k(2);
        P3(hora)=x_k(3);
        P4(hora)=x_k(4);
        FO(hora)=x_k(5);
		SVR(hora)=x_k(9);
        Pcarga(hora)=sum(P1(hora)+P2(hora)+P3(hora)+P4(hora));
        %Time(hora)=hora;
        %Recálculo del SOC de la bateria a partir de las potencias P2,P3,P4 entregadas a la carga
        P_descargar=x_k(2)/Bc_max;                          %valor a descargar de la batería a partir de P2 entregado
        P_cargar=((P3_max-x_k(3))+(P4_max-x_k(4)))/Bc_max;   %valor a cargar de la batería a partir de P3, P4 excedentes
        SOC_VEC(hora+1)=SOC_0-P_descargar+P_cargar;    %actualización de SOC a partir de P2 entregada y valor a cargar de P3,P4
        if SOC_VEC(hora+1)>=0.95
            SOC_VEC(hora+1)=0.95;
        end
        sumP2=sumP2+x_k(2);
        sumP3=sumP3+x_k(3);
        sumP4=sumP4+x_k(4);
    end
    
    
%     display(['Evaluaciones=',num2str(NFE)]);
% %     display(['Suma FO=',num2str(sum(FO(:)))])
% %     display(size(SOC_VEC));
    
   %%%%%%%%%%%% GENERAR EXCEL CON DATOS DE LAS MEJORES SOLUCIONES POR GENERACIÓN%%%%%%%%%%%%%%%%%%%%%
    Name = ['/Potencias_' num2str(e) '.xls'];
    Nameexcel = [PathFolder Name];
    encabezado = {'Hora' 'P1'  'P2'  'P3'  'P4'  'Carga'  'FO' 'SVR' 'SOC de bateria (%)'};
    vector_horas=[0:23]';
    xlswrite(Nameexcel, encabezado, 1, 'A1');
    xlswrite(Nameexcel, vector_horas, 1, 'A2:A25')
    xlswrite(Nameexcel,sum(FO(:)),1,'G26');
    
    %se agrupan los vectores calculados en una matriz. Esto con el fin de
    %escribir la matriz completa y no hacer una escritura por vector, lo
    %cual alenta la ejecución del programa
     P_matriz=[P1,P2,P3,P4,Pcarga,FO,SVR,SOC_VEC(1:24,:)];
     xlswrite(Nameexcel,P_matriz,'B2:I25');
    
    %Gráfica de los valores de potencia en cada hora
    Time=0:1:23;
    figure(e);
    p=plot(Time,Pcarga,Time,P1,Time,P2,Time,P3,Time,P4);
    set(p,'linewidth',2);
    title('Gráfica P vs t')
    xlabel('Tiempo [h]');
    ylabel('Potencia [W]');
    legend('Carga','P1','P2','P3','P4');
    saveas(e,fullfile(PathFolder, ['Grafica_' num2str(e)]), 'png');
    sFO(e,1)=sum(FO(:));
    sFO(e,2)=e;
    Ahorro_Total=sum(FO(:));
    close all
    beep;

%encabezado={'FO' 'Corrida'};
%xlswrite('1-Resultados_DERand1Bin_12k.xlsx',encabezado,1,'A1');
%xlswrite('1-Resultados_DERand1Bin_12k.xlsx',sFO,1,'A2');
%display(mat2str(sFO));
%close all
end
%% Función Objetivo
function y=funcionObjetivo(x)
%Gestión de generación para microred en modo isla
cf=1;
w1=0.25;
w2=0.25;
w3=0.25;
w4=0.25;
y=[0.000435*w1*cf 0 0 0]*x.^2+[0.3*w1*cf 119*w2 -545.016*w3 -152.616*w4]*x+[14.88*w1*cf 0 0 0]*ones(size(x,1),1);
end
function[E_F]=uptdate_E_F(E,E_F,best)
EF=[];
n=get_problem_n();
if(best(n+5)==0)
EF(1)=best(n+1); 
EF(2)=E;
s= size(E_F);
if(isempty(E_F))
s=[0 0];
end

E_F(s(1)+1,:,:)=EF;
end
end
