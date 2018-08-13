function[E_F]=uptdate_E_F(E,E_F,best)
EF=[];
n=get_problem_n();
if(best(n+2)==0)
EF(1)=best(n+1); 
EF(2)=E;
s= size(E_F);
if(isempty(E_F))
s=[0 0];
end

E_F(s(1)+1,:,:)=EF;
end

end