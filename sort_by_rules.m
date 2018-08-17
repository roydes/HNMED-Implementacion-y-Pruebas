function [X] = sort_by_rules(X)
X=insertion_sort(X,@compare_by_rules);
%s=size(X);
%X=quicksort(X, 1, s(1));
end