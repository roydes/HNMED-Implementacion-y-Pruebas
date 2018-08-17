function [f_value, constrain_violation_sum] =  f(x)

switch get_problem_i()
    case 1
         [f_value,constrain_violation_sum]= Problema_1(x);
    case 2
         [f_value,constrain_violation_sum]= Problema_2(x);
    case 3
         [f_value,constrain_violation_sum]= Problema_3(x);
    case 4
         [f_value,constrain_violation_sum]= Problema_4(x);
    case 5
         [f_value,constrain_violation_sum]= Problema_5(x);  
    case 6
         [f_value,constrain_violation_sum]= Problema_6(x);  
    case 7
         [f_value,constrain_violation_sum]= Problema_7(x);  
    case 8
         [f_value,constrain_violation_sum]= Problema_8(x);  
    otherwise
        disp('Finish')
end
