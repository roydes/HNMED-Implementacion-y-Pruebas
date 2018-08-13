%**************** NELDER-MEAD METHOD ***********************
 function [result]= NM2ELA(NE)

rng('shuffle', 'twister');  

%STEP-1
%Choosing gamma>1, 0<beta<1
gamma=1.05;
%input('Input a value for gamma')
beta=0.5;
alpha=2;
%Creating an initial simplex
n= get_problem_n();%(N+1)
bounds=get_problem_bound();
x=generate_random_pop(n+1,n,bounds);
x=evaluate_pop(x,n+1,n);
E=0;
while E<NE
    %Sorting points by according the Debs Rules
    x=sort_by_rules(x);
    %Highest function evaluation-Worst point
    xh=x(n+1,:);
    %Next to the worst point
    xg=x(n,:);
    %Lowest-function evaluation-best point
    xl=x(1,:);

    %Calculating the centroid
    xc=sum(x(1:n,:))/(n);
    %STEP 3
    %Reflection
    %disp('Reflection');
    xr=alpha*xc(1:n)-xh(1:n);
    xr=limit(xr(1:n),bounds);
    % Evaluate xr
    xr=evaluate_x_i(xr,n);
    E=E+1;

    xnew=xr;
    if compare_by_rules(xr,xl)
        %Expansion
        %disp('Expansion');
        xnew= (1+gamma)*xc(1:n)-gamma*xh(1:n);
    elseif ~compare_by_rules(xr,xh)
        %Contraction
        %disp('Contraction outside');
        xnew= (1-beta)*xc(1:n)+beta*xh(1:n);
    elseif compare_by_rules(xg,xr) && compare_by_rules(xr,xh)
        %Contraction
        %disp('Contraction inside');
        xnew= (1+beta)*xc(1:n)-beta*xh(1:n);
    end
    % evaluate xnew
      xnew(1:n)=limit(xnew(1:n),bounds);
      xnew=evaluate_x_i(xnew(1:n),n);
      E=E+1;

    if compare_by_rules(xnew,xh)
     x(n+1,:)=xnew;
   else 
      %no better then large random size expansion 
      xnew= (1+rand(1)*n)*xc(1:n)-rand(1)*xh(1:n);
      xnew=evaluate_x_i(xnew(1:n),n);
      E=E+1;
      x(n+1,:)=xnew;
    end
   if ~compare_by_rules(xnew,xh)
     % Make a random large expasion to ashive better points. If a bad point
     % is found a reflexion or expantion will be made in the next iteration
     % take the best as the worst to jump out of the local minima
      xnew= (1+n*rand)*xc(1:n)-rand*xl(1:n);
      xnew(1:n)=limit(xnew(1:n),bounds);
      xnew=evaluate_x_i(xnew(1:n),n);
      x(n+1,:)=xnew;
      E=E+1;
     %
   end
    
end
    result=x(1,:);
 end



