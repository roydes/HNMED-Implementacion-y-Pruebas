
function [X]= quicksort(X, low, high)
    if low < high 
      [pi,X] = partition(X, low, high);
      X= quicksort(X, low, pi - 1 );
      X= quicksort(X, pi + 1, high);
    end
end
function [p,X]= partition(X, low, high) 
    pivot = X(high,:);
    i = low-1;  
    for j = low: high
        if compare_by_rules(X(j,:),pivot)
            i = i + 1;
            X=swap(X,i,j);
        end
    end    
    X=swap(X,i+1,high);
    p=i+1;
end

function X = swap(X,i,j)
    val = X(i,:);
    X(i,:) = X(j,:);
    X(j,:) = val;
end