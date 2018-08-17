% Variante del Algoritmo Evolución Diferencial DE/rand/1/bin
function [best,E_F]=DE_rand_1_bin(NE,NP)
rand('twister',sum(100*clock));
amt0=clock;
n= get_problem_n();
rango=get_problem_bound();
X=generate_random_pop(NP,n,rango);
X=evaluate_pop(X,NP,n);
best=X(1,:);
E=NP;
E_F=[];

while E<NE
    %por generacion generar F y CR
    [F,CR]=generate_F_and_CR();
    for i=1:NP
        [xr1,xr2,xr3]=generate_xr1_xr2_xr3(i,X);
        v=xr1(1:n)+F*(xr2(1:n)-xr3(1:n));
        v=limit(v,rango);
        u=zeros();
        jrand=randi(n,1);
        for j=1:n
            randj=rand();
            if(randj<CR || jrand==j)
                u(j)=v(j);
            else
                u(j)=X(i,j);
            end
        end
        u=evaluate_x_i(u,n);
        E=E+1;
        
        if compare_by_rules(u,X(i,:))
            X(i,:)=u;
        end
        if compare_by_rules(u,best)
            best=u;
        end
        
    end
    E_F=uptdate_E_F(E,E_F,best);

end

end
