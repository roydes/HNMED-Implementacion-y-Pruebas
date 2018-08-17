function [best,E_F]= HNMEDV5(NE,NS)
% Inciando generador de números aleatorios.
rng('shuffle','twister');
% Seleccionando parámetros de expansión,
% contracción y reflexión.
gamma=1.05;
beta=0.5;
alpha=2;
% Asignando la n para el problema actual.
n= get_problem_n();
% Asignando los límites de las variables.
bounds=get_problem_bound();
% Calculando el tamaño de la población.
NP=(n+1)*NS;
X=zeros(NP,n+2);
% Generando y evaluando la población.  
X=generate_random_pop(NP,n,bounds);
X=evaluate_pop(X,NP,n);
best=X(1,:);
E=NP;
XL=zeros(NS,n+2);
E_F=[];
% Comienza el ciclo hasta llegar al número 
% máximo de evaluaciones
while E<NE
    k=1;
    w=1;
    % Se generan los factores 
    % de escala y crizamiento
    [F,CR]=generate_F_and_CR();
    while w<NP
        % NELDER MEAD SONBRE SIMPLEX SK
        % Seleccionar Sk de X
        Sk=X(w:(w+n),:);
        % Ordenea ascendentemente de acuerdo 
        % a las reglas de Deb
        Sk=sort_by_rules(Sk);
        % Seleccionar el mejor punto
        xl=Sk(1,:);
        % Seleccionar el segundo peor punto
        xg=Sk(n,:);
        % Seleccionar el peor punto
        xh=Sk(n+1,:);
        % Calcular el centroide
        xc= sum(Sk(1:n,:))/(n);
        %________Reflexión____________________
        xr=alpha*xc(1:n)-xh(1:n);
        xr=limit(xr(1:n),bounds);
        xr=evaluate_x_i(xr,n);
        xnew=xr;
        E=E+1;
        if compare_by_rules(xr,xl)
        %_______Expansión_____________________
            xnew= (1+gamma)*xc(1:n)-gamma*xh(1:n);
        elseif ~compare_by_rules(xr,xh)
        %_______Contracción  hacia fuera______
            xnew= (1-beta)*xc(1:n)+beta*xh(1:n);
        elseif compare_by_rules(xg,xr) && compare_by_rules(xr,xh)
        %_______Contracción hacia dentro______
            xnew= (1+beta)*xc(1:n)-beta*xh(1:n);
        end
        xnew(1:n)=limit(xnew(1:n),bounds);
        xnew=evaluate_x_i(xnew(1:n),n);
        E=E+1;
        if compare_by_rules(xnew,xh)
           Sk(n+1,:)=xnew;
        else
        % Si ninguno de los operadores originales 
        % mejora a xnew se aplica el operador de V5
           [xr1,xr2,xr3]=generate_xr1_xr2_xr3(w,X);
           xnew=xr1(1:n)+F*(best(1:n)-xr1(1:n))+(1-F)*(xr2(1:n)-xr3(1:n));
           xnew=limit(xnew(1:n),bounds);
           xnew=evaluate_x_i(xnew,n);
           E=E+1;
           Sk(n+1,:)=xnew;
        end
        if compare_by_rules(xnew,xl)
            xl=xnew; 
        end 
        XL(k,:)=xl; 
        if compare_by_rules(xnew,best)
            best=xnew;
        end
        % EVOLUCIÓN DIFERENCIAL SE APLICA  A XL(K)
        [xr1,xr2,xr3]=generate_xr1_xr2_xr3(w,X);
        parent=XL(k,:);
        v=xr1(1:n)+F*(xr2(1:n)-xr3(1:n));        
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
            Sk(1,:)=u;
            XL(k,:)=u;
        end
        if compare_by_rules(u,best)
            best=u;
        end
        % SE ACTUALIZA POBLACIÓN con Sk
        X(w:w+n,:)=Sk;
        k=k+1;
        w=w+n+1;
    end
    % Se almacena el valor de la F(x) de best y E para
    % Gráficas de convergencia
    [E_F]=uptdate_E_F(E,E_F,best);
end

end


