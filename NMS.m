function [bestsol]=NMS(NE)
%Esta función aplica el metodo de Nelder Mead  con operador de encogimiento
% para busqueda del mínimo de la función
%del probema 2 de síntesis de un mecanismo de cuatro barras
%La entrada de la función es una matriz con las N+1 soluciones candidatas
% N es la cantidad de variables de diseño.
alfa=1;
beta=0.5;
gamma=2;

%Generacion aleatoria de los miembros del simplex.
rand('twister',sum(100*clock));  %semilla del generador de aleatorios
N= get_problem_n();
rango=get_problem_bound();
S=[ ]; %Vector que almacena el simplex del método.
for i=1:30  %generar población inicial aleatoria
       for j=1:N
            S(i,j)=rango(j,1)+(rango(j,2)-rango(j,1))*rand(1);
       end
end
%Generación aleatoria de los miembros del simplex

%epsilon=1e-50;
vector=[];
xbar=[];
for i=1:30
   [z(i,1) z(i,2)]=f(S(i,:));
   vector(i,:)=[S(i,:) z(i,:)];
end
vector;
vector=clasificagral(vector,N);
S=vector(1:N+1,:);
vector=[];
paro=1; 
E=N+1;
%paro=abs(S(1,N+1)-S(N+1,N+1))
%while (paro>epsilon)
while (E<NE)
    %Computing the centroid
    xbar(1,1:N)=0;
    for i=1:N
       for j=1:N 
            xbar(i)=xbar(i)+S(j,i)/N;
       end
    end
    %Computing the centroid
    %Computing the reflected vertex
    for i=1:N
        sr(i)=(1+alfa)*xbar(i)-alfa*S(N+1,i);
    end
    sr=limita(sr);
    [sr(N+1) sr(N+2)]=f(sr(1:N));
    E=E+1;
    %Computing the reflected vertex
    sk=sr;
    band=compara(sr,S(N,:));
     if band==1     %sr(N+1)<S(N,N+1)
         band=compara(sr,S(1,:));
         if band==1 %sr(N+1)<S(1,N+1)
            %Computing the expanded vertex 
            for i=1:N
               se(i)=gamma*sr(i)+(1-gamma)*xbar(i);
            end
            se=limita(se);
            [se(N+1) se(N+2)]=f(se(1:N));
            %Computing the expanded vertex
            band=compara(se,S(1,:));
            if band==1   %se(N+1)<S(1,N+1)
                sk=se;
            end
         end
     else
         st=S(N+1,:);
         st=limita(st);
         [st(N+1) st(N+2)]=f(st(1:N));
         band=compara(sr,st);
         if band==1   %sr(N+1)<st(N+1)
             st=sr;
         end
         %Computing the contraction vertex
         for i=1:N
            sc(i)=beta*st(i)+(1-beta)*xbar(i);
         end
         sc=limita(sc);
         [sc(N+1) sc(N+2)]=f(sc(1:N));
         E=E+1;
         %Computing the contraction vertex
         band=compara(sc,S(N,:));
         if band==1   %sc(N+1)<S(N,N+1)
             sk=sc;
         else
             %Computing the shrinkage operation
             for j=2:N
                 for i=1:N
                     S(j,i)=(S(1,i)+S(j,i))/2;
                 end
                 S(j,1:N)=limita(S(j,1:N));
                 [S(j,N+1) S(j,N+2)]=f(S(j,1:N));
             end             
             for i=1:N
                 sk(i)=(S(1,i)+S(N+1,i))/2;
             end
             sk=limita(sk);
             [sk(N+1) sk(N+2)]=f(sk(1:N));
             E=E+1;
             %Computing the shrinkage operation
         end        
     end
     for i=1:N
         vector(i,:)=S(i,:);
     end
     vector(N+1,:)=sk; 
     paro=paro+1;
     S=clasificagral(vector,N);
     %paro=abs(S(1,N+1)-S(N+1,N+1))
end %fin del mientras
bestsol=S(1,:);
end
 
function [HN]=clasificagral(H,Num)

%Función para clasificar lo elementos del simplex
%La entrada es el arreglo del vertex con el siguiente orden
%[variables_diseño, F.O., SVR]
%Num es el numero de variables de diseño.

[ren col]=size(H);
j=1;
k=1;
band1=0;
band2=0;
for i=1:ren
    if H(i,Num+2)==0
        temp1(j,:)=H(i,:);
        j=j+1;
        band1=1;
    else
        temp2(k,:)=H(i,:);
        k=k+1;
        band2=1;
    end
end
%temp1
%temp2
if band1==1 && band2==1 %Caso donde hubo factibles y no factibles
   HN1=sortrows(temp1,Num+1);
   HN2=sortrows(temp2,Num+2);
   HN=[HN1;HN2];
else
    if band1==1 && band2==0 %Caso donde solo hubo factibles
        HN=sortrows(temp1,Num+1);
    else
        HN=sortrows(temp2,Num+2); %Cso donde solo hubo no factibles
    end
end
end
function [bandera] = compara(sol1,sol2)

% Esta función evalua si la sol1 es mejor que la sol2
%utilizando las reglas de DEB
%Si sol1 es mejor que sol2 la bandera regresa el valor 1
%en caso contrario regresa el valor 0.

 var=get_problem_n();
 
 if sol1(1,var+2)>0 && sol2(1,var+2)>0 %ambos violan restricciones
     if sol1(1,var+2)>=sol2(1,var+2)
         bandera=0;  %sol2 es mejor
     else
         bandera=1;  %sol1 es mejor
     end
 else
     if sol1(1,var+2)<=0 && sol2(1,var+2)<=0 %ambos NO violan restricciones
         if sol1(1,var+1)>=sol2(1,var+1) % se evalua la funcion objetivo
            bandera=0;  %sol2 es mejor
         else
            bandera=1;  %sol1 es mejor
         end
     else
         if sol1(1,var+2)>0 && sol2(1,var+2)==0 %sol1 viola restricciones y sol2 NO viola
             bandera=0;
         else
             bandera=1;
         end
     end
 end
end
function [acotado]=limita(soln)

% Función que verifica los limites de las variables
% del vector candidato a solución
% La entrada es el vector candidato y la salida es
% el vector acotado.

rango=get_problem_bound();

var=get_problem_n;
vector=soln;
for k=1:var
  while vector(1,k)<rango(k,1) || vector(1,k)>rango(k,2)
    if vector(1,k) < rango(k,1)
           vector(1,k)=(rango(k,1)*2)-vector(1,k);
    else
         if vector(1,k) > rango(k,2)
            vector(1,k)=(rango(k,2)*2)-vector(1,k);
         end
    end  
  end
end

acotado=vector;
    
end


 