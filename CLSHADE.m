function [best,E_F] = CLSHADE(D,cotas,MAX_NFE,c,H,p,rNinit,Nmin)
close all; format long;
if nargin < 4
  error('argumentos_entrada :  D,cotas,MAX_NFE y corridas son argumentos obligatorios')
end

if(D<1)
    error('El argumento D debe ser mayor a 1');
end

if(size(cotas,1)~=D || size(cotas,2)~=2)
    error('La matriz cotas debe tener D filas por 2 columnas');
end

if(MAX_NFE<1)
    error('El argumento MAX_NFE debe ser mayor a 1');
end


if (~exist('H','var'))
        H = 6;
else
    if H<1
        error('El parámetro de H debe ser mayor a 1');
    end

end

if (~exist('p','var'))
        p = 0.11;
else
    if p<=0 || p>1
        error('El parámetro p debe ser mayor a 0 y menor que 1');
    end
end

if (~exist('rNinit','var'))
        rNinit = 18;
else
    if rNinit<1
        error('El parámetro de rNinit debe ser mayor a 1');
    end
end

if (~exist('Nmin','var'))
        Nmin = 4;
else
    if Nmin<4 || Nmin>round(D*rNinit)
        error('El parámetro Nmin debe ser mayor o igual a 4 y menor que D*rNinit');
    end        
end

%fprintf('====================================== C-LSHADE ======================================\n');
%fprintf('                      -------------Valores ingresados-------------   \n\n');
%display(['D=',num2str(D)]);
%display(['cotas=',mat2str(cotas)]);
%display(['MAX_NFE=',num2str(MAX_NFE)]);
%fprintf('\n');
%display(['H=',num2str(H)]);
%display(['p=',num2str(p)]);
%display(['rNinit=',num2str(rNinit)]);
%display(['Nmin=',num2str(Nmin)]);




%Parámetros constantes de CLSHADE
Ninit=round(D*rNinit);
tv=-1;
mejores=zeros(30,D+3);
E_F=[];
%fprintf('                  ----------Inicia ejecución del algoritmo----------   \n\n');
%INICIAN LAS EJECUCIONES DEL ALGORITMO

    rng('shuffle', 'twister');
 %   display(['CORRIDA ',num2str(c)]);
    k=1; %contador de posición de memoria
    g=0; %contador de generaciones
    NFE=0; %contador de evaluaciones
    NP=Ninit; %población
    mFCR=zeros(H,2); %Se crea una matriz de H filas por 2 columnas
    mFCR(:,:)=0.5; %inicialización de memoria MF y MCR. 
                   %Los valores de F se almacenarán en la columna 1 de mFCR y los de CR en la columna 2
    
    %INICIA POBLACIÓN INICIAL
    %Creación de la población inicial
    X=zeros(NP,size(cotas,1));
    for j=1:size(cotas,1)
        X(:,j)=(cotas(j,2)-cotas(j,1)).*rand(NP,1) + cotas(j,1);
    end
    
    %Evalación de la población inicial
    for i=1:NP
        [X(i,D+1),X(i,D+2)]=funcionObjetivo(X(i,1:D)); %La función regresa el valor de FO y SVR
    end
    NFE=NP;
    %Se inicializa al mejor individuo
    best=X(1,:);
    %TERMINA POBLACIÓN INICIAL
    
    %INICIA CICLO EVOLUTIVO
    while (NP>=Nmin && NFE<MAX_NFE)
        g=g+1; %Se incrementa el contador de generaciones
        
        %Se declaran como vacías las estructuras temporales sF, sCR y sDF
        sF=[]; 
        sCR=[];
        sDF=[]; 
        
        %Cálculos necesarios para la mutación curren-to-pbest
        top=round(NP*p); %Tamaño del subconjunto de población para current-to-pbest
        
        %Se evita que haya un error si el subconjunto es de 0 individuos
        if(top==0)
            top=1;
        end
        
        %Se ordena la población dando prioridad a la SVR (D+2) y luego a la
        %%FO(D+1). La función regresa sortrows regresa tanto los vectores como los índices ordenados 
        [vals,indices] = sortrows(X,[D+2,D+1]);
    
        %INICIA MUTACIÓN, CRUZA Y SELECCIÓN POR INDIVIDUO
        for i=1:NP
            %Inicia generación de los valores de F y CR que se utilizarán por individuo
            ri=randi(H,1);
            F=cauchyrnd(mFCR(ri,1),0.1); %mFCR(ri,1) es el valor de F en la memoria
            while(F<=0 || F>1)
                if(F>1)
                    F=1;
                end
                if(F<=0)
                    F=cauchyrnd(mFCR(ri,1),0.1); %distribución de Cauchy
                end
            end

            if(mFCR(ri,2)==tv)  
                CR=0;           
            else
                CR=normrnd(mFCR(ri,2),0.1); %distribución normal. mFCR(ri,2) es el valor de CR en la memoria
                if(CR<0)
                    CR=0;
                elseif(CR>1)
                    CR=1;
                end
            end
            %Termina generación de F y CR
            
            %Inicia creación de vector mutante de acuerdo a current-to-pbest
            randomNum=randi(top,1); %índice de uno de los mejores individuos
            x_pbest=X(indices(randomNum),:); %individuos seleccionado
            
            r1=randi(NP,1); %índice para vector r1
            while(r1==i)
                r1=randi(NP,1);
            end
            r2=randi(NP,1); %índice para vector r2
            while(r2==i || r2==r1)
                r2=randi(NP,1);
            end
            
            v=X(i,1:D)+F*(x_pbest(1:D)-X(i,1:D))+F*(X(r1,1:D)-X(r2,1:D)); %vector mutante
            %Termina creación del vector mutante
            
            %Inicia verificación de los límites (cotas) para cada variable del individuo
            for n=1:D
                    if(v(n)<cotas(n,1))
                        %v(n)=2.0*cotas(n,1)-v(n);
                        v(n)=(cotas(n,1)+X(i,n))/2;
                    elseif (v(n)>cotas(n,2))
                        %v(n)=2.0*cotas(n,2)-v(n);
                        v(n)=(cotas(n,2)+X(i,n))/2;
                    end
            end
            %Termina verificación de límites
            
            %Inicia creación de vector trial (u)
            u=zeros();
            jrand=randi(D,1);
            for j=1:D
                randj=rand();
                if(randj<CR || jrand==j)
                    u(j)=v(j);
                else
                    u(j)=X(i,j);
                end
            end
            
            [u(D+1),u(D+2)]=funcionObjetivo(u); %Se evalúa el trial con la función objetivo
            NFE=NFE+1; %se aumenta el contador de evaluaciones
            %Termina creación del vector trial
            
            %Inicia selección entre padre e hijo
            if(u(D+2)==0 && X(i,D+2)==0) %si ambos son factibles se selecciona el individuo con mejor valor de aptitud (función objetivo)
                if(u(D+1)<=X(i,D+1))
                    if(u(D+1)<X(i,D+1))
                        sF=[sF,F]; %se guarda el valor de F en sF
                        sCR=[sCR,CR]; %se guarda el valor de CR en sCR
                        sDF=[sDF,abs(u(D+1)-X(i,D+1))]; %se guarda diferencia de aptitud en sDF
                    end
                    X(i,:)=u; %se reemplaza el target por el trial
                end
            elseif(u(D+2)==0 && X(i,D+2)~=0)  %si sólo el trial es factible, éste es seleccionado
                sF=[sF,F]; %se guarda el valor de F en sF
                sCR=[sCR,CR]; %se guarda el valor de CR en sCR
                sDF=[sDF,abs(u(D+1)-X(i,D+1))]; %se guarda diferencia de aptitud en sDF
                X(i,:)=u;
            elseif(u(D+2)~=0 && X(i,D+2)~=0) %si ambos no son factibles, selecciona aquel con la menor SVR
                if(u(D+2)<=X(i,D+2))
                    if(u(D+2)<X(i,D+2))
                        sF=[sF,F]; %se guarda el valor de F en sF
                        sCR=[sCR,CR]; %se guarda el valor de CR en sCR
                        sDF=[sDF,abs(u(D+1)-X(i,D+1))]; %se guarda diferencia de aptitud en sDF
                    end
                    X(i,:)=u;
                end
            end
            %Termina selección entre padre e hijo
            
            %Se actualiza el mejor individuo segun reglas de Deb 
            if(best(D+2)==0 && X(i,D+2)==0) %Si ambos son factibles se selecciona el individuo con mejor valor de aptitud (función objetivo)
                if(X(i,D+1)<=best(D+1))
                    best=X(i,:);
                end
            elseif(best(D+2)~=0 && X(i,D+2)==0)  %Si sólo un individuo es factible, éste es seleccionado
                best=X(i,:);
            elseif(best(D+2)~=0 && X(i,D+2)~=0) %Si ambos no son factibles, selecciona aquel con la menor SVR
                if(X(i,D+2)<= best(D+2))
                    best=X(i,:);
                end
            end
        end
        %TERMINA MUTACIÓN, CRUZA Y SELECCIÓN POR INDIVIDUO
        
        %INICIA ACTUALIZACIÓN DE MEMORIA 
        if(~isempty(sDF)) %si la estructura no es vacía, se actualiza la memoria
            %Se inicializan las variables que se ocuparán para calcular el promedio de Lehmer
            numeradorF=0;
            denominadorF=0;
            numeradorCR=0;
            denominadorCR=0;
            
            %Se calculan los promedios de Lehmer utilizando las estructuras temporales sCR y sF
            for n=1:length(sCR)
                wk=sDF(n)/sum(sDF); %peso ponderado
                numeradorCR=numeradorCR+(wk*(sCR(n)^2));
                denominadorCR=denominadorCR+(wk*sCR(n));
                numeradorF=numeradorF+(wk*(sF(n)^2));
                denominadorF=denominadorF+(wk*sF(n));
            end
            meanWF=numeradorF/denominadorF; %promedio de Lehmer para F
            meanWCR=numeradorCR/denominadorCR; %promedio de Lehmer para CR
            
            %Se actualiza la posición k de la memoria mFCR
            mFCR(k,1)=meanWF; %se actualiza el valor de F en la posición k
            
            %Se actualiza el valor de CR en la posición k
            if((mFCR(k,2)==tv) || (max(sCR)==0))
                mFCR(k,2)=tv;
            else
                mFCR(k,2)=meanWCR;
            end
            
            k=k+1; %se actualiza el contador de posición de la memoria
            if(k>H)
                k=1; %se reinicia el contador de posición si excede el tamaño H
            end
        end
        %TERMINA ACTUALIZACIÓN DE LA MEMORIA
        
        %INICIA AJUSTE DE POBLACIÓN (NP) DE INDIVIDUOS
        NG=NP;
        NP=round((((Nmin-Ninit)/MAX_NFE)*NFE)+Ninit); %función lineal para reducción de población
        
        %si el nuevo tamaño de la población disminuye, se eliminan los peores individuos 
        if(NP<NG) 
            peores=NG-NP;
            for n=1:peores
                %Se ordena la población dando prioridad a la SVR (D+2) y luego a la
                %%FO(D+1). La función sortrows regresa tanto los vectores como los índices ordenados
                [vals, indices]=sortrows(X,[D+2,D+1]);
                X(indices(end),:)=[];
            end
        end
        %TERMINA EL AJUSTE DE POBLACIÓN DE INDIVIDUOS
    [E_F]=uptdate_E_F(NFE,E_F,best);

    end

    %TERMINA CICLO EVOLUTIVO
    
    %se muestra en pantalla mejor solución encontrada
   % display(['Mejor Solución: ',mat2str(best(1:D))]);
   % display(['FO: ',mat2str(best(D+1))]);
   % display(['SVR: ',num2str(best(D+2))]);
   % display(['Evaluaciones realizadas: ',num2str(NFE)]);
   % fprintf('\n');
    
    %Se escriben los mejores individuos de cada corrida en una matriz
    %mejores(c,1:D+2)=best;
    %mejores(c,D+3)=c;
    

%TERMINAN LAS EJECUCIONES DEL ALGORITMO

%Se escriben los mejores individuos de cada generación en un archivo
%encabezado={};
%for i=1:D
 %   encabezado=[encabezado,['x' num2str(i)]];
%end
%encabezado=[encabezado,'FO'];
%encabezado=[encabezado,'SVR'];
%encabezado=[encabezado,'Corrida'];
%nombreArchivo=['/CLSHADE_' num2str(30) 'Ejecuciones_' num2str(NFE) 'Eval.xlsx'];
%xlswrite(nombreArchivo,encabezado,1,'A1');
%xlswrite(nombreArchivo,mejores,1,'A2'); 


end

function r = cauchyrnd(mu,c)
    r = c.*tan(pi*(rand(1)-0.5))+mu;
end

function [f_value, constrain_violation_sum] =  funcionObjetivo(x)
 
[f_value, constrain_violation_sum]=f(x);
 
end