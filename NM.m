%**************** NELDER-MEAD METHOD ***********************
 function [result]= NM(E)

rand('twister',sum(100*clock));  %semilla del generador de aleatorios

%STEP-1
%Choosing gamma>1, 0<beta<1
gamma=1.05;
beta=0.5;
alpha=2;

n= get_problem_n()+1;%(N+1)
rango=get_problem_bound();

% Generating a random initial simplex
for i=1 :n
    for j=1:n-1
       x(i,j)=rango(j,1)+(rango(j,2)-rango(j,1))*rand(1);
    end
end

%STEP-2
x=evaluate(x,n);

Evaluaciones=0;
while Evaluaciones<E
    %Descending sorting points by the objective function value
    x=sort_by_contrain_first(x,n);
    %Highest function evaluation-Worst point
    xh=x(1,:);
    %Next to the worst point
    xg=x(2,:);
    %Lowest-function evaluation-best point
    xl=x(n,:);
    %Calculating the centroid
    xc=sum(x(2:n,:))/(n-1);
    %STEP 3
    %Reflection
    %disp('Reflection');
    xr=alpha*xc(1:n-1)-xh(1:n-1);
    % Evaluate xr
    [error_r,cvs_r] =  f(xr);
    xr(n)=error_r;
    xr(n+1)=cvs_r;
    Evaluaciones=Evaluaciones+1;
    xnew=xr;
    if compare_by_rules(xr,xl)
        %Expansion
        %disp('Expansion');
        xnew= (1+gamma)*xc(1:n-1)-gamma*xh(1:n-1);
    elseif ~compare_by_rules(xr,xh)
        %Contraction
        %disp('Contraction outside');
        xnew= (1-beta)*xc(1:n-1)+beta*xh(1:n-1);
    elseif compare_by_rules(xg,xr) && compare_by_rules(xr,xh)
        %Contraction
        %disp('Contraction inside');
        xnew= (1+beta)*xc(1:n-1)-beta*xh(1:n-1);
    end
    % evaluate xnew
      xnew(1:n-1)=limit(xnew(1:n-1),rango);
      xnew=evaluate_x_i(xnew(1:n-1),n);
      Evaluaciones=Evaluaciones+1;
      x(1,:)=xnew;
     
end
    x=sort_by_contrain_first(x,n);
    result=x(n,:);
 end

 
function [x_evaluated]=evaluate(x,n)
        for i = 1:n
            [error_i,cvs_i] =  f(x(i,1:n-1));
             x(i,n)=error_i;
             x(i,n+1)=cvs_i;
        end
        x_evaluated=x;
end 
function [x_i]=evaluate_x_i(x_i,n)
            [error_i,cvs_i] =  f(x_i(1:n-1));
             x_i(n)=error_i;
             x_i(n+1)=cvs_i;
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
             if (x_j(n) <= x_j1(n) ) %tomo el q tenga menor suma de SVR
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
    for i = 1:n
        for j = 1:n-i 
            if compare_by_rules(x(j,:),x(j+1,:))
                z= x(j,:);
                x(j,:) = x(j+1,:);
                x(j+1,:) = z;
            end
        end
    end
    sorted_x=x;
end    





