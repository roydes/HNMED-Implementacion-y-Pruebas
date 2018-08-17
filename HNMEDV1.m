function [best,E_F]= HNMEDV1(NE,NS)

rand('twister',sum(100*clock));
gamma=1.05;
beta=0.5;
alpha=2;
n= get_problem_n();
rango=get_problem_bound();
NP=(n+1)*NS;
X=generate_random_pop(NP,n,rango);
X=evaluate_pop(X,NP,n);
best=X(1,:);
E=NP;
XL=zeros(NS,n+2);
E_F=[];
while E<NE
    k=1;
    w=1;
    [F,CR]=generate_F_and_CR();
    while w<NP
        % NELDER MEAD FOR SIMPLEX K
        %__________________________________________________________________
        %Select Simplex Sk
        Sk=X(w:(w+n),:);
        %Sort
        Sk=sort_by_rules(Sk);
        %Lowest-function evaluation-best point
        xl=Sk(1,:);
        xl(n+1);
        %highest-function value is the worst point
        xh=Sk(n+1,:);
        %Next to the worst point
        xg=Sk(n,:);
        % Calculate centroid
        xc= sum(Sk(1:n,:))/(n);
        
        %________Reflection____________________
        xr=alpha*xc(1:n)-xh(1:n);
        % Evaluate xr
        xr=evaluate_x_i(xr,n);
        xnew=xr;
        E=E+1;
        
        
        if compare_by_rules(xr,xl)
        %________Expansion_____________________
            xnew= (1+gamma)*xc(1:n)-gamma*xh(1:n);
        elseif ~compare_by_rules(xr,xh)
        %________Contraction outside___________
            xnew= (1-beta)*xc(1:n)+beta*xh(1:n);
        elseif compare_by_rules(xg,xr) && compare_by_rules(xr,xh)
        %________Contraction inside______________
            xnew= (1+beta)*xc(1:n)-beta*xh(1:n);
        end
        % evaluate xnew
        xnew(1:n)=limit(xnew(1:n),rango);
        xnew=evaluate_x_i(xnew(1:n),n);
        E=E+1;
        
        
        if compare_by_rules(xnew,xh)
            Sk(n+1,:)=xnew;
        else
        %******* xnew no better than worst then apply V1 operator********
            xnew= (1+rand(1)*n)*xc(1:n)-rand(1)*xl(1:n);
            xnew=limit(xnew(1:n),rango);
            xnew=evaluate_x_i(xnew,n);
            E=E+1;
            Sk(n+1,:)=xnew;
        end
        
        if compare_by_rules(xnew,xl)    xl=xnew; end 
                XL(k,:)=xl; 
        if compare_by_rules(xnew,best)   best=xnew;end

               
        % DIFERENTIAL EVOLUTION FOR SMALL POPULATION (XL) THE BESTS POINTS
        %_________________________________________________________________
        
        [xr1,xr2,xr3]=generate_xr1_xr2_xr3(w,X);
        parent=XL(k,:);
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
            Sk(1,:)=u;
            XL(k,:)=u;
        end
        if compare_by_rules(u,best)
            best=u;
        end
        X(w:w+n,:)=Sk;
        k=k+1;
        w=w+n+1;
    end
    [E_F]=uptdate_E_F(E,E_F,best);
end
end


