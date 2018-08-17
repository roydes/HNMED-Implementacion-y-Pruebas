function [x_i]=limit(x_i,bounds)
for j=1:length(bounds)
    while x_i(j)<bounds(j,1) || x_i(j)>bounds(j,2)
        if  x_i(j)<bounds(j,1)
            x_i(j)=2.0*bounds(j,1)-x_i(j);
        else    
            x_i(j)=2.0*bounds(j,2)-x_i(j);
        end
    end
    
end

end