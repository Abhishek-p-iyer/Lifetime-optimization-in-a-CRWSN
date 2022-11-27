function FindProbabilities(i)
global n N LuciferinLevel pb

LL_sum = 0;
for j=1:n
    LL_sum=LL_sum+N(i,j)*(LuciferinLevel(j)-LuciferinLevel(i));
end 

if(LL_sum==0)
    pb(i,:)=zeros(1,n);
else 
    for j=1:n
        pb(i,j)=(N(i,j)*(LuciferinLevel(j)-LuciferinLevel(i)))/LL_sum;
    end 
end 
