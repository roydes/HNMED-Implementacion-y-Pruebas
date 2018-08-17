%L_SHADE
function [Ahorro_Total,E_F]= MR_CLSHADE(e)
format long

    rng('shuffle', 'twister');
    %% INICIALIZACI�N DE VARIABLES A TRAV�S DE LECTURA EN EXCEL
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
    
    %%par�metros bateria para restricci�n SOC
    n_B=0.85; % Round trip efficiency
    n_c=0.85;
    n_d=1;
    
    % celdas_variables3=xlsread('Potencias_redes.xls','K5:M5');
    % Bc_max=celdas_variables3(3); %Capacidad m�xima
    % alpha_bateria=n_B/Bc_max;
    % alpha_c=n_c/Bc_max;
    % alpha_d=n_d/Bc_max;
    % SOC_max=celdas_variables3(2);    %Estado de carga m�nimo de la bateria
    % SOC_min=celdas_variables3(1);    %Estado de carga m�ximo de la bateria
    
    Bc_max=2000; %Capacidad m�xima
    alpha_bateria=n_B/Bc_max;
    alpha_c=n_c/Bc_max;
    alpha_d=n_d/Bc_max;
    SOC_max=0.95;    %Estado de carga m�nimo de la bateria
    SOC_min=0.4;    %Estado de carga m�ximo de la bateria
    
    %%par�metros de inicializaci�n para cumplir con la sumatoria de la restricci�n de SOC para las 24 horas
    sumP2=0;
    sumP3=0;
    sumP4=0;
    HORAS=24;
    %%Par�metros del algoritmo de evoluci�n diferencial
    var=4; %Cantidad de variables de dise�o.
    %% INICIA LA EJECUCI�N DE LOS 24 CASOS DE OPTIMIZACI�N
    
    
    %%%%%%%%%%%% CREAR DIRECTORIOS NECESARIOS PARA GUARDAR ARCHIVOS%%%%%%%%%%%%%%%%%%%%%
    [stat,struc] = fileattrib;
    PathCurrent = struc.Name;
    GlobalFolderName='/Ejecuciones_LSHADE';
    mkdir(PathCurrent, GlobalFolderName);
    FolderName=['/Ejecucion' num2str(e)];
    PathFolder = [PathCurrent GlobalFolderName FolderName];
    mkdir([PathCurrent GlobalFolderName], FolderName);
    E_F=[];
    
    for hora=1:HORAS    %Cada hora representa una problema de optimizaci�n diferente     
        
       %display(['HORA ',num2str(hora)]);
        %%Inicializaci�n de variables y cotas necesarias para cada hora
        PL=PL_VEC(hora);                    %Potencia de la carga en la hora
        SOC_0=SOC_VEC(hora);                %SOC actual en la hora
        P1_max=P1_VEC(hora);                %Potencia Diesel
        P2_max=SOC_0*Bc_max-SOC_min*Bc_max; %Potencia de la bateria que depende de capacidad m�xima de la misma y el SOC(0)
        P3_max=P3_VEC(hora);                %Potencia de PV
        P4_max=P4_VEC(hora);                %Potencia wind
        
        P1_min=0;
        P2_min=0;
        %P2_min=SOC_min*Bc_max;
        P3_min=0;
        P4_min=0;
        
        cotas=[P1_min P1_max;P2_min P2_max;P3_min P3_max;P4_min P4_max];
        
        %% INICIA EL ALGORITMO ED        
        
        %SHADE
        H=6; %L-SHADE
        mFCR=zeros(H,2); mFCR(:,:)=0.5;
        A=[]; %%Es el archivo manejado en JADE
        ks=1;
        p=0.11;

        %L-SHADE
        TV=0.001;
        rNinit=18;
        Ninit=round(var*rNinit);
        Nmin=4;
        NP=Ninit;
        MAX_NFE=10000;
        NFE=0;
        %
        
        
        %%Generando la poblaci�n inicial
        X=zeros(NP,size(cotas,1));
        for j=1:size(cotas,1)
            X(:,j)=(cotas(j,2)-cotas(j,1)).*rand(NP,1) + cotas(j,1);
        end
        
        %Evaluando la funci�n objetivo y las restricciones de la poblaci�n inicial
        
        %Restricci�n de desigualdades distintas a m�ximos y m�nimos--> ej:SOC
        %A_NNO=[0 alpha_d -alpha_c -alpha_c];
        b_NNO=-SOC_min+SOC_0-alpha_d*sumP2+alpha_c*sumP3+alpha_c*sumP4;
        %A_DO=[0 -alpha_d alpha_c alpha_c];
        b_DO=SOC_max-SOC_0+alpha_d*sumP2-alpha_c*sumP3-alpha_c*sumP4;
        
        for i=1:NP
            %Funci�n objetivo
            [X(i,var+1)]=funcionObjetivo(X(i,1:var)');
            %Calculando la restricci�n de igualdad
            %P1_k+P2_k+P3_k+P4_k-PL_k=0
            %utilizaremos un epsilon "cero gordo" de 1E-6
            X(i,var+2)=sum(X(i,1:var))-PL;
            %Fin de Calculando la restricci�n de igualdad.
            
            %Evaluaci�n de las restricciones de desigualdad
            Lizq1=alpha_d*X(i,2)-alpha_c*X(i,3)-alpha_c*X(i,4);
            X(i,var+3)=Lizq1-b_NNO;
            Lizq2=-alpha_d*X(i,2)+alpha_c*X(i,3)+alpha_c*X(i,4);
            X(i,var+4)=Lizq2-b_DO;
            %Fin de evaluaci�n de las restricciones de desigualdad
            
            %Calculando la SVR (Suma de violaci�n de restricciones) de cada
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
      
        %Fin de generando la primera poblaci�n

        %Se busca al mejor individuo de la poblaci�n inicial usando
        %criterios de Deb
        best=X(1,:);%Para iniciar la busqueda del mejor, el best se inicializa con el primer individuo de la poblaci�n
        for i=2:NP
            if(best(var+5)==0 && X(i,var+5)==0) %Si ambos son factibles se selecciona el individuo con mejor valor de aptitud (funci�n objetivo)
                if(X(i,var+1)<=best(var+1))
                    best=X(i,:);
                end
            elseif(best(var+5)~=0 && X(i,var+5)==0)  %Si s�lo un individuo es factible, �ste es seleccionado
                best=X(i,:);
            elseif(best(var+5)~=0 && X(i,var+5)~=0) %Si ambos no son factibles, selecciona aquel con la menor suma de violaciones
                if(X(i,var+5)<= best(var+5))
                    best=X(i,:);
                end
            end
        end        
        
        %TRIALS, MUTACI�N Y CRUZA
        mejoresGen=[];
        g=0;
        algo=0;
        while (NP>=Nmin && NFE<MAX_NFE)
            g=g+1;
            %%SHADE
            sF=[];
            sCR=[];
            DeltaF=[]; 

            top=round(NP*p);
            if(top==0)
                top=1;
            end
            X_Order=sortrows(X,[var+5,var+1]);
            
            %Generaci�n de vectores para cruza
            for i=1:NP
                %%SHADE
                ri=randi(H,1);
                F=cauchyrnd(mFCR(ri,1),0.1);
                while(F<=0 || F>1)
                    if(F>1)
                        F=1;
                    end
                    if(F<=0)
                        F=cauchyrnd(mFCR(ri,1),0.1);
                    end
                end

                if(mFCR(ri,2)==TV)  %L-SHADE
                    CR=0;           %L-SHADE 
                else
                    CR=normrnd(mFCR(ri,2),0.1);
                    if(CR<0)
                        CR=0;
                    elseif(CR>1)
                        CR=1;
                    end
                end

                %Se crea el vector mutante de acuerdo a current-to-pbest
                rbest=randi(top,1);
                x_pbest=X_Order(rbest,:);
    %             display(rbest);

                r1=randi(NP,1);
                while(r1==i)
                    r1=randi(NP,1);
                end
                r2=randi(NP,1);
                while(r2==i || r2==r1)
                    r2=randi(NP,1);
                end

                %Se genera el vector ruido
                ruido=X(i,1:var)+F*(x_pbest(1:var)-X(i,1:var))+F*(X(r1,1:var)-X(r2,1:var));
                
                
                %Se verifican los l�mites (cotas) para cada variable del individuo
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
                
                %Funci�n objetivo
                [u(var+1)]=funcionObjetivo(u');
                NFE=NFE+1;
                %Calculando la restricci�n de igualdad
                %P1_k+P2_k+P3_k+P4_k-PL_k=0
                %utilizaremos un epsilon "cero gordo" de 1E-6
                u(var+2)=sum(u(1:var))-PL;
                %Fin de Calculando la restricci�n de igualdad.
                
                %Evaluaci�n de las restricciones de desigualdad
                Lizq1=alpha_d*u(2)-alpha_c*u(3)-alpha_c*u(4);
                u(var+3)=Lizq1-b_NNO;
                Lizq2=-alpha_d*u(2)+alpha_c*u(3)+alpha_c*u(4);
                u(var+4)=Lizq2-b_DO;
                %Fin de evaluaci�n de las restricciones de desigualdad
                
                %Calculando la SVR (Suma de violaci�n de restricciones) de cada
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
                if (X(i,var+5)~=0) && (u(var+5)~=0) %Si ambos violan restricciones
                    if (X(i,var+5)>=u(var+5)) %Se toma al que tiene menor SVR.
                        if(u(var+1)<X(i,var+1))
                            %                                 if(size(A,1)>=NP)
                            %                                     A(randi(NP,1),:)=[];
                            %                                 end
                            %                                 A=vertcat(A,X(i,:));
                            sF=[sF,F];
                            sCR=[sCR,CR];
                            DeltaF=[DeltaF,abs(u(var+1)-X(i,var+1))];
                        end
                        X(i,:)=u;
                    end
                elseif ((X(i,var+5)==0) && (u(var+5)==0)) %Si ambos no violan restricciones
                    if X(i,var+1)>=u(var+1)% Se toma al que tiene mejor F.O.
                        if(u(var+1)<X(i,var+1))
                            %if(size(A,1)>=NP)
                            %    A(randi(NP,1),:)=[];
                            %end
                            %A=vertcat(A,X(i,:));
                            sF=[sF,F];
                            sCR=[sCR,CR];
                            DeltaF=[DeltaF,abs(u(var+1)-X(i,var+1))];
                        end
                        X(i,:)=u;
                    end
                elseif ((X(i,var+5)>0) && (u(var+5)==0)) %El padre viola y el trial(u) no viola
                    sF=[sF,F];
                    sCR=[sCR,CR];
                    DeltaF=[DeltaF,abs(u(var+1)-X(i,var+1))];
                    X(i,:)=u;
                end
                
                %Se reemplaza al mejor individuo con criterios de Deb
                if(best(var+5)==0 && X(i,var+5)==0) %Si ambos son factibles se selecciona el individuo con mejor valor de aptitud (funci�n objetivo)
                    if(X(i,var+1)<=best(var+1))
                        best=X(i,:);
                    end
                elseif(best(var+5)~=0 && X(i,var+5)==0)  %Si s�lo un individuo es factible, �ste es seleccionado
                    best=X(i,:);
                elseif(best(var+5)~=0 && X(i,var+5)~=0) %Si ambos no son factibles, selecciona aquel con la menor suma de violaciones
                    if(X(i,var+5)<= best(var+5))
                        best=X(i,:);
                    end
                end                
                
            end
            
            %%SHADE
            if(~isempty(sF) && ~isempty(sCR))
                %%Se actualiza la memoria
                %display('SHADE');
                sumF1=0;
                sumF2=0;
                sumCR1=0;
                sumCR2=0;
                for n=1:length(sCR)
                    wk=DeltaF(n)/sum(DeltaF);
                    sumCR1=sumCR1+(wk*(sCR(n)^2));
                    sumCR2=sumCR2+(wk*sCR(n));
                    sumF1=sumF1+(wk*(sF(n)^2));
                    sumF2=sumF2+(wk*sF(n));
                end
                
                
                meanWF=sumF1/sumF2;
                meanWCR=sumCR1/sumCR2;
                
                mFCR(ks,1)=meanWF;
                
                if((mFCR(ks,2)==TV) || (max(sCR)==0))
                    mFCR(ks,2)=TV;
                else
                    mFCR(ks,2)=meanWCR;
                end
                
                %display(mat2str(mFCR(ks,:)));
                
                ks=ks+1;
                if(ks>H)
                    ks=1;
                end
            end
            NG=NP;
            NP=round((((Nmin-Ninit)/MAX_NFE)*NFE)+Ninit);
            
            if(NP<NG)
                X=sortrows(X,[var+2,var+1]);
                borrar=NG-NP;
                for n=1:borrar
                    X(size(X,1),:)=[];
                    
                    %A(randi(NG,1),:)=[];
                end
            end
            
            %         display(size(X,1));
            %         display(NFE);
            
            %Si el individuo es factible se almacena en la matriz
            %mejoresGen para que pase a formar parte de la grafica de
            %convergencia
            if(best(var+5)==0)
                if(algo==0)
                    algo=g;
                end
                if hora==13
                   E_F= uptdate_E_F(NFE,E_F,best);
                end
            end
            
        end
%         display(['Mejor Soluci�n: ',mat2str(best)]);
%         display(['Generaciones: ',num2str(g)]);
%         display(['Evaluaciones: ',num2str(NFE)]);
%         display(['Primer_Factible_En_Generacion: ',num2str(algo)]);
        
        
        %Gr�fica de convergencia
%         figure(hora);
%         q=plot(mejoresGen(:,var+6),mejoresGen(:,var+1));
%         %set(p,'linewidth',2);
% %         set(gca,'YScale','log');
% %         set(gca,'XScale','log');
%         
%         set(gca,'YScale','linear');
%         set(gca,'XScale','linear');
%         
%         
%         title('Gr�fica f(x) vs NFE')
%         xlabel('NFE');
%         ylabel('f(x)');
%         saveas(hora,fullfile(PathFolder, ['Grafica_Convergencia_Hora' num2str(hora)]), 'png');
        x_k=best;
        %%Termina el algoritmo DE
        
        %% Generaci�n vectores de Potencia para gr�fica
        P1(hora)=x_k(1);
        P2(hora)=x_k(2);
        P3(hora)=x_k(3);
        P4(hora)=x_k(4);
        FO(hora)=x_k(5);
		SVR(hora)=x_k(9);
        Pcarga(hora)=sum(P1(hora)+P2(hora)+P3(hora)+P4(hora));
        %Time(hora)=hora;
        %Rec�lculo del SOC de la bateria a partir de las potencias P2,P3,P4 entregadas a la carga
        P_descargar=x_k(2)/Bc_max;                          %valor a descargar de la bater�a a partir de P2 entregado
        P_cargar=((P3_max-x_k(3))+(P4_max-x_k(4)))/Bc_max;   %valor a cargar de la bater�a a partir de P3, P4 excedentes
        SOC_VEC(hora+1)=SOC_0-P_descargar+P_cargar;    %actualizaci�n de SOC a partir de P2 entregada y valor a cargar de P3,P4
        if SOC_VEC(hora+1)>=0.95
            SOC_VEC(hora+1)=0.95;
        end
        sumP2=sumP2+x_k(2);
        sumP3=sumP3+x_k(3);
        sumP4=sumP4+x_k(4);
    end
%     display(['NFE=',num2str(NFE)]);
%     display(['Suma FO=',num2str(sum(FO(:)))])
%     display(size(SOC_VEC));
    
    %%%%%%%%%%%% GENERAR EXCEL CON DATOS DE LAS MEJORES SOLUCIONES POR GENERACI�N%%%%%%%%%%%%%%%%%%%%%
%     Name = ['/Potencias_' num2str(e) '.xls'];
%     Nameexcel = [PathFolder Name];
%     encabezado = {'Hora' 'P1'  'P2'  'P3'  'P4'  'Carga'  'FO' 'SVR' 'SOC de bateria (%)'};
%     vector_horas=[0:23]';
%     xlswrite(Nameexcel, encabezado, 1, 'A1');
%     xlswrite(Nameexcel, vector_horas, 1, 'A2:A25')
%     xlswrite(Nameexcel,sum(FO(:)),1,'G26');
%     
%     %se agrupan los vectores calculados en una matriz. Esto con el fin de
%     %escribir la matriz completa y no hacer una escritura por vector, lo
%     %cual alenta la ejecuci�n del programa
%     P_matriz=[P1,P2,P3,P4,Pcarga,FO,SVR,SOC_VEC(1:24,:)];
%     xlswrite(Nameexcel,P_matriz,'B2:I25');
    
    %Gr�fica de los valores de potencia en cada hora
%     Time=0:1:23;
%     figure(e);
%  %   p=plot(Time,Pcarga,Time,P1,Time,P2,Time,P3,Time,P4);
%     set(p,'linewidth',2);
%     title('Gr�fica P vs t')
%     xlabel('Tiempo [h]');
%     ylabel('Potencia [W]');
%     legend('Carga','P1','P2','P3','P4');
%     saveas(e,fullfile(PathFolder, ['Grafica_Final' num2str(e)]), 'png')   
     Ahorro_Total=sum(FO(:));
    close all
%xlswrite('1-Resultados_L-SHADE.xlsx',sFO,1,'A1');
%display(mat2str(sFO));
end
%% Funci�n Objetivo
function y=funcionObjetivo(x)
%Gesti�n de generaci�n para microred en modo isla
cf=1;
w1=0.25;
w2=0.25;
w3=0.25;
w4=0.25;
y=[0.000435*w1*cf 0 0 0]*x.^2+[0.3*w1*cf 119*w2 -545.016*w3 -152.616*w4]*x+[14.88*w1*cf 0 0 0]*ones(size(x,1),1);
end
function[E_F]=uptdate_E_F(E,E_F,best)
EF=[];
n=4;
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
function r = cauchyrnd(mu,c)
    r = c.*tan(pi*(rand(1)-0.5))+mu;
end
