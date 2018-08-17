function [X]=generate_random_pop(NP,n,bounds)
X=zeros(NP,n+2);
for i=1:NP
    X(i,1:n)=generate_random_vector(bounds,n);
end
end