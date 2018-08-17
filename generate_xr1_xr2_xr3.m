function[xr1,xr2,xr3]=generate_xr1_xr2_xr3(i,X)
S=size(X);
NP=S(1);
r1=randi(NP,1);
while(r1==i)
    r1=randi(NP,1);
end
r2=randi(NP,1);
while(r2==i || r2==r1)
    r2=randi(NP,1);
end
r3=randi(NP,1);
while(r3==i || r3==r2 || r3==r1)
    r3=randi(NP,1);
end
xr1=X(r1,:);
xr2=X(r2,:);
xr3=X(r3,:);
end