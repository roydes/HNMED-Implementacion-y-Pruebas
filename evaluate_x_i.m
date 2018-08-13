function [x_i]=evaluate_x_i(x_i,n)
[f_value,cvs_i] =  f(x_i(1:n));
x_i(n+1)=f_value;
x_i(n+2)=cvs_i;
end