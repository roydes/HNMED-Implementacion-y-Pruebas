function [X]=generate_regular_simplex(NS)
n=get_problem_n();
rango=get_problem_bound();
rng(0,'v5uniform');
%rng('shuffle','twister');
v=1;
k=1;
w=(NS/2)*(n+1)+1;
S=zeros(n+1,n+2);
SO=zeros(n+1,n+2);
X=zeros(NS*(n+1),n+2);


while k<=NS
    c=generate_random_vector(rango,n);
    
    for i=1:n+1
        for j=1:n
            d=rango(j,2)-rango(j,1)/n;
            p=d*((sqrt(n+1)+n-1)/n*sqrt(2));
            q=d*((sqrt(n+1)-1)/n*sqrt(2));
            if i==1
                x1=c-(p+(n-1)*q)/(n+1);
                x1=limit(x1,rango);
                x1=evaluate_x_i(x1,n);
                S(1,:)=x1;
            else
                if i==j
                    S(i,j)=x1(j)+p+d*0.7*rand();
                else
                    S(i,j)=x1(j)+q-d*0.7*rand();
                end
            end
            %SO(i,j)=rango(j,2)+rango(j,1)-S(i,j);
        end
        S(i,:)=limit(S(i,:),rango);
        S(i,:)=evaluate_x_i(S(i,:),n);
        %SO(i,:)=evaluate_x_i(SO(i,:),n);
    end
    X(v:v+n,:)=S;
    %X(w:w+n,:)=SO;
    %w=w+n+1;
    v=v+n+1;
    k=k+1;
end
%}
rng('shuffle','twister');
%S(n+1,1)=rango(j,2)-((rango(j,2)-rango(j,1))*A)/n;