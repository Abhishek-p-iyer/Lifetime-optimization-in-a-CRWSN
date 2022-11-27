function FindNeighbors(i)
global n m A N rd Na LuciferinLevel 
n_sum=0;
for j=1:n
    if(j~=i)
        square_sum=0;
        for k=1:m
            square_sum=square_sum+(A(i,k)-A(j,k))^2;
        end 
        d=sqrt(square_sum);
        if(d<=rd(i) && (LuciferinLevel(i)<LuciferinLevel(j)))
            N(i,j)=1;
            n_sum=n_sum+1;
        end
    end 
    Na(i)=n_sum;
end
