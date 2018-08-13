function [best,E_F]= ED(NE,NP)

rng('shuffle','twister');
n= get_problem_n();
bounds=get_problem_bound();
X=zeros(NP,n+2);
X=generate_random_pop(NP,n,bounds);
X=evaluate_pop(X,NP,n);
best=find_minimun(X);
E=NP;
E_F=[];
while E<NE
    
    for i=1:NP
        % DIFERENTIAL EVOLUTION FOR SMALL POPULATION (XL) THE BESTS POINTS
        %_________________________________________________________________
        [F,CR]=generate_F_and_CR();
        [xr1,xr2,xr3]=generate_xr1_xr2_xr3(i,X);
        parent=X(i,:);
        v=xr1(1:n)+F*(best(1:n)-xr1(1:n))+(1-F)*(xr2(1:n)-xr3(1:n));
       %         v=xr1(1:n)+F*(xr2(1:n)-xr3(1:n));        

        v=limit(v,bounds);
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
            X(i,:)=u;
        end
        if compare_by_rules(u,best)
            best=u;
        end
    end
 X
    [E_F]=uptdate_E_F(E,E_F,best);
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