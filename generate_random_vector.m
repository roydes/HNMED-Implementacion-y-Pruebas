function [x]= generate_random_vector(bounds,n)
for j=1:n
    x(j)=bounds(j,1)+(bounds(j,2)-bounds(j,1))*rand(1);
end
end