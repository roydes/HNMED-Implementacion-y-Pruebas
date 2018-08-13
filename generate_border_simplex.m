function [X]=generate_border_simplex(NS)
n=get_problem_n();
bound=get_problem_bound();
% Generador aleatorio con distribucón uniforme
rng(0,'v5uniform');
% Límites para q_k
a=0.05;
b=3*n/100;
if b>0.95
   b=0.95;
end
% Límites para p_j
c=0.05;
d=5*n/100;
if d>0.95
   d=0.95;
end
v=1;
k=1;
S=zeros(n+1,n+2);
X=zeros(NS*(n+1),n+2);
while k<=NS
for i=1:n+1
% Generar q_k para cada simplex.
       q_k=a+(b-a)*rand(1);
% Situar para cada vertice i en la vecindad 
% del límite supieror de la dimensión i
for j=1:n
    if i==j
       p_j=c+(d-c)*rand(1); 
       S(i,j)=((bound(j,2)-bound(j,1))/2+((bound(j,2)-bound(j,1))/2)*p_j);
    else
       if rand()<0.5  
          r=(NS*(n+1)/100*rand());
          S(i,j)=bound(j,1)+(((bound(j,2)-bound(j,1))/2)*q_k) +r;
       else
          r=(NS*(n+1)/100*rand());   
          S(i,j)=bound(j,1)+(((bound(j,2)-bound(j,1))/2)*q_k) -r;
       end
    end
end
% Limitar y evaluar el vértice. 
S(i,:)=limit(S(i,:),bound); 
S(i,:)=evaluate_x_i(S(i,:),n);
end
% Agregar el simplex a la población
X(v:v+n,:)=S;
v=v+n+1;
k=k+1;
end
end
