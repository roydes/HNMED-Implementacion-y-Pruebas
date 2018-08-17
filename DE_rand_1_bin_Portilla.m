% Variante del Algoritmo Evolución Diferencial DE/rand/1/bin
function [best,E_F]=DE_rand_1_bin_Portilla(NE,ind) %semilla del generador de aleatorios
rand('twister',sum(100*clock));
%rng(0,'twister'); 
amt0=clock;
n= get_problem_n();
rango=get_problem_bound();
E=ind;
E_F=[];

 
w = warning ('off','all');

hijo=zeros(1,n+2);
vector=zeros(ind,n+2); %Vector qu almacena la poblacion de la ED. [6 variables, 1 F.O., 1 SVR]

var=n; %num de variables


restric=4; %num de restricciones
amt0=clock;
  for i=1:ind  %generar población inicial aleatoria
       for j=1:var
            vector(i,j)=rango(j,1)+(rango(j,2)-rango(j,1))*rand(1);
       end
       %evaluar a los ind en la funcion objetivo y las restricciones
        vector(i,:)=evaluate_x_i( vector(i,:),n);
  end
best=vector(1,:);
while E<NE%Ciclo del algoritmo para el numero maximo de generaciones
     %disp('generación')
     %g
     % Calculo del factor de mutación o escalamiento. Este factor es
     % calculado por generación.
     a=0.3;
     b=0.9;
     F=a+(b-a)*rand(1);
     a=0.8;
     b=1;
     CR=a+(b-a)*rand(1);
     %Ciclo del algoritmo de ED.
     for i=1:ind
         %Seleccionar tres individuos diferentes
           r1=randint(1,1,[1,ind]);
          while(r1==i)
                r1=randint(1,1,[1,ind]);   
          end
          r2=randint(1,1,[1,ind]); 
          while(r2==r1 || (r2==i))
                r2=randint(1,1,[1,ind]);   
          end
          r3=randint(1,1,[1,ind]);
          while(r3==r1 || (r3==r2) || r3==i)
               r3=randint(1,1,[1,ind]);   
          end   
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
                                       
         %Se genera el vector de ruido r3+F(r1-r2)
         for k=1:var
            vectRuido(1,k)=vector(r1,k)+F*(vector(r2,k)-vector(r3,k));
             %Verificando limites de las variables en el hijo
             if vectRuido(1,k) < rango(k,1)
                vectRuido(1,k)=(rango(k,1)*2)-vectRuido(1,k);
             else
                 if vectRuido(1,k) > rango(k,2)
                     vectRuido(1,k)=(rango(k,2)*2)-vectRuido(1,k);
                 end
             end
         end   
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         jrand=randint(1,1,[1,var]); %garantiza que al menos un elemento
                                       %del hijo sera diferente q el padre
         
        %Generar hijo a partir del target y del vector ruido
        for k=1:var
            randj=rand(1);  %para la cruza ver q % va a tener de variación
           if (randj<CR || jrand==k)
               hijo(1,k)=vectRuido(1,k);        
           else
               hijo(1,k)=vector(i,k);
           end
        end     
             
       %evaluar a los ind en la funcion objetivo y las restricciones
       hijo(1,:)=evaluate_x_i( hijo(1,:),n);
       E=E+1;
       %Reemplazo del padre por el hijo si es mejor el hijo de acuerdo a
       %las reglas de Deb.
       
       if (vector(i,var+2)>0) && (hijo(1,var+2)>0) %si ambos violan restricciones
             if (vector(i,var+2) >= hijo(1,var+2) ) %tomo el q tenga menor suma de SVR
                 vector(i,:)=hijo(1,:);
             end
             
        else 
              if (vector(i,var+2)<=0) && (hijo(1,var+2)<=0) %si ninguno viola restricciones
                  if (vector(i,var+1)>hijo(1,var+1)) % tomo al que tiene mejor valor fx
                       vector(i,:)=hijo(1,:);
                  end
              end
        end
          
        if vector(i,var+2)>0 && hijo(1,var+2)==0 % si el padre viola restricciones y el hijo no
              % tomo al hijo
              vector(i,:)=hijo(1,:);
        end
        % Actualizando el mejor
        if (best(1,var+2)>0) && (hijo(1,var+2)>0) %si ambos violan restricciones
             if (best(1,var+2) >= hijo(1,var+2) ) %tomo el q tenga menor suma de SVR
                 best(1,:)=hijo(1,:);
             end
             
        else 
              if (best(1,var+2)<=0) && (hijo(1,var+2)<=0) %si ninguno viola restricciones
                  if (best(1,var+1)>hijo(1,var+1)) % tomo al que tiene mejor valor fx
                       best(1,:)=hijo(1,:);
                  end
              end
        end
          
        if best(1,var+2)>0 && hijo(1,var+2)==0 % si el padre viola restricciones y el hijo no
              % tomo al hijo
              best(1,:)=hijo(1,:);
        end
        
    end %del for de ind
    [E_F]=uptdate_E_F(E,E_F,best);

end %del for de genmax

end