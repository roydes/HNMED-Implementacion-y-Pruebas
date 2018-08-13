%**************** HOOKE-JEEVES METHOD ***********************
function [best]= HHJED1(NE,NP)
%--------------STEP-1------------------------
xk=[];
delta_i=0.5;

alpha=2;
epsilon= 0.0000001;
n= get_problem_n();
delta(1:NP,1:n)=delta_i;
rango=get_problem_bound();
X=zeros(NP,n+2);
X=generate_random_pop(NP,n,rango);
X=evaluate_pop(X,NP,n);
E=NP;
i=n;
best=X(1,:);
while E<NE
  %--------------STEP-2------------------------
  % Exploratory Move
  %disp('Xk inicial')
 
      [F,CR]=generate_F_and_CR();

  for l=1:NP 
       if i>n
      i=1;
  end
   xk=X(l,:);   
   xk= evaluate_x_i(xk,n);
   x= Exploratory_Move(xk,delta(l,:),rango,i);
   E=E+3;
  %--------------STEP-3------------------------
    
    if are_vectors_equals(x,xk)==false
        xk1=x;
        %disp('Exploration sussces')
        EP=1;
        while compare_by_rules(xk1,xk) && EP<3
  %--------------STEP-4------------------------
  %Pattern Move
           [xr1,xr2,xr3]=generate_xr1_xr2_xr3(l,X);
            xp=xk1(1:n)+F*(best(1:n)-xr1(1:n))+(1-F)*(xk(1:n)-xr3(1:n));
            xp=limit(xp(1:n),rango);
            xp=evaluate_x_i(xp,n);
 %--------------STEP-5------------------------
  % Exploratory Move
             %disp('Secon EM')
            x=Exploratory_Move(xp,delta(l,:),rango,i);
            E=E+3;
            xk=xk1;
            xk1=x;
            EP=EP+1;
        end
    else
        delta(l,:)=delta(l,:)/alpha;
    end
    X(l,:)=xk;  
   if compare_by_rules(xk,best)
            best=xk;
   end
    % DIFERENTIAL EVOLUTION FOR SMALL POPULATION (XL) THE BESTS POINTS
        %_________________________________________________________________
        
        [xr1,xr2,xr3]=generate_xr1_xr2_xr3(l,X);
        parent=X(l,:);
        v=xr1(1:n)+F*(xr2(1:n)-xr3(1:n));        
        v=limit(v,rango);
        u=zeros();
        jrand=randi(n,1);
        
        for j=1:n
            randj=rand();
            if(randj<CR || jrand==j)
                u(j)=v(j);
            else
                u(j)=parent(j);
            end
        end
        u= evaluate_x_i(u,n);
        E=E+1;
        
        if compare_by_rules(u,parent)
            X(l,:)=u;
        end
        if compare_by_rules(u,best)
            best=u;
        end
    i=i+1;
  end
E
best(n+1)
end
end
function [x]= Exploratory_Move(x,delta,rango,i)
n= get_problem_n();
   xp=x;
   xn=x;
  xp(i)=xp(i)+delta(i);
  xn(i)=xn(i)-delta(i);
   
   xp= limit(xp,rango);
   xn= limit(xn,rango);
    fx= evaluate_x_i(x,n);%vector with the evaluation of f(x) and the cvs
    fxp= evaluate_x_i(xp,n);%vector with the evaluation of +f(x) and the cvs
    fxn= evaluate_x_i(xn,n);%vector with the evaluation of -f(x) and the cvs
    
     x=find_minimun([fx;fxp;fxn]);
%}
end



function [population]=generate_population(size,n,rango)
 population=[];
 for i=1:size
  for j=1:n
    population(i,j)=rango(j,1)+(rango(j,2)-rango(j,1))*rand(1);
  end
 end
end


function [best, population]=operador_ED(population,population_size,gen,var)
     a=0.3;
     b=0.9;
     F=a+(b-a)*rand(1);
     a=0.8;
     b=1;
     CR=0.5;
rango=get_problem_bound();
for g=1:gen
 for i=1:population_size
  r1=randi([1,population_size],1,1);
          while(r1==i)
                %r1=randint(1,1,[1,ind]);  
                r1=randi([1,population_size],1,1);
          end
          %  r2=randint(1,1,[1,ind]); 
          r2=randi([1,population_size],1,1);
          while(r2==r1 || (r2==i))
                %r2=randint(1,1,[1,ind]);   
                r2=randi([1,population_size],1,1);
          end
         % r3=randint(1,1,[1,ind]);
         r3=randi([1,population_size],1,1);
          while(r3==r1 || (r3==r2) || r3==i)
              % r3=randint(1,1,[1,ind]);   
              r3=randi([1,population_size],1,1);
          end   
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
                                       
         %Se genera el vector de ruido r3+F(r1-r2)
         for k=1:var
             vectRuido(1,k)=population(r1,k)+F*(population(r2,k)-population(r3,k));
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
         
         %jrand=randint(1,1,[1,var]); %garantiza que al menos un elemento
                                       %del hijo sera diferente q el padre
         jrand=randi([1,var],1,1);
         
        %Generar hijo a partir del target y del vector ruido
        for k=1:var
            randj=rand(1);  %para la cruza ver q % va a tener de variación
           if (randj<CR || jrand==k)
               hijo(1,k)=vectRuido(1,k);        
           else
               hijo(1,k)=population(i,k);
           end
        end 
        hijo=evaluate_x_i(hijo,var);
        
        if compare_by_rules(hijo,population(i,:))
            population(i,:)=hijo;
        end
 end
end
 population=sort_by_contrain_first(population,population_size);
 best=population(1,:);
 end


function [x_evaluated]=evaluate(x,size,n)
        for i = 1:size
            [error_i,cvs_i] =  f(x(i,1:n));
             x(i,n+1)=error_i;
             x(i,n+2)=cvs_i;
        end
        x_evaluated=x;
end 
function [x_i]=evaluate_x_i(x_i,n)
            [error_i,cvs_i] =  f(x_i(1:n));
             x_i(n+1)=error_i;
             x_i(n+2)=cvs_i;
end 
function [x_i]=limit(x_i,rango)
        for j=1:length(rango)
            while x_i(j)<rango(j,1) || x_i(j)>rango(j,2)
                if  x_i(j)<rango(j,1)
                     x_i(j)=rango(j,1)*2- x_i(j);
                end
                if  x_i(j)>rango(j,2)
                     x_i(j)=rango(j,2)*2- x_i(j);
                end
            end
    
        end
 
 end
% compare x_j is selected over x_j1 according the Deb's rules
function [selected]=compare_by_rules(x_j,x_j1)
         n=length(x_j);   
         selected=false;

        if (x_j(n)>0) && (x_j1(n)>0) %si ambos violan restricciones
             if (x_j(n) < x_j1(n) ) %tomo el q tenga menor suma de SVR
                 selected=true;
             end
             
        else 
              if (x_j(n)<=0) && (x_j1(n)<=0) %si ninguno viola restricciones
                  if x_j(n-1)< x_j1(n-1) % tomo al que tiene mejor valor fx
                      selected=true;
                  end
              end
        end
          
        if x_j1(n)>0 && x_j(n)==0 % si el padre viola restricciones y el hijo no
            selected=true; % tomo al hijo
        end
  
end

function [sorted_x] = sort_by_contrain_first(x,n)
s=size(x);
    for i = 1:s(1)
        for j = 1:s(1)-i 
            if compare_by_rules(x(j+1,:),x(j,:))
                z= x(j,:);
                x(j,:) = x(j+1,:);
                x(j+1,:) = z;
            end
        end
    end
    sorted_x=x;
end    


function [equals] =are_vectors_equals(x1,x2)
n=get_problem_n();
equals=true;
  for i=1:n
      if x1(i)~=x2(i)
          equals=false;
          break;
      end
  end
end 

 function [minimun]=find_minimun(X)
s=size(X);
minimun=X(1,:);
for i=2:s(1)
    if compare_by_rules(X(i,:),minimun)
      minimun=X(i,:);
    end
end
  
 end

