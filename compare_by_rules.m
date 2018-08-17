function [selected]=compare_by_rules(x1,x2)
n=get_problem_n();
selected=false;
% Si ambos violan restricciones 
% seleccionar el que tenga menor suma de SVR
if x1(n+2)~=0 && x2(n+2)~=0 
    if x1(n+2) < x2(n+2)   %  
        selected=true;
    end
    
else
% Si ninguno viola restricciones  
% seleccionar al que tiene menor valor f(x)
    if x1(n+2)==0 && x2(n+2)==0 
        if x1(n+1) < x2(n+1)  
            selected=true;
        end
    end
end
% si un inviduo viola restricciones
% y el otro no seleccionar este último
if x2(n+2)>0 && x1(n+2)==0 
    selected=true; 
end
end