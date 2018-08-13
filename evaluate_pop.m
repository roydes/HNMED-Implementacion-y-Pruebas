function [X]=evaluate_pop(X,size,n)
for i = 1:size
    [error_i,cvs_i] =  f(X(i,1:n));
    X(i,n+1)=error_i;
    X(i,n+2)=cvs_i;
end
end
