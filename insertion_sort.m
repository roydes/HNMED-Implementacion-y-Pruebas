% Insertion sort en orden ascendente con criterio de comparación. Recibe
% una matriz mxn y la devuelve ordenada
function sorted = insertion_sort(X,compare)
if(isempty(X))
    error('¡¡Arreglo vacio!!');
end
s= size(X);
m=s(1);
for i = 2:m
    d = i;    
    while((d > 1) && compare(X(d,:),X(d-1,:)))
        temp = X(d,:);
        X(d,:) = X(d-1,:);
        X(d-1,:) = temp;
        d = d-1;
    end
end
sorted = X;
end