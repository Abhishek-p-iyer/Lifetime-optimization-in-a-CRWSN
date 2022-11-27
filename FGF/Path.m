function Del=Path(i,j)
global A m
square_sum=0;
for k=1:m
    square_sum=square_sum + (A(i,k)-A(j,k))^2;
end 
hyp=sqrt(square_sum);
for k=1:m
    Del(:,k)=(A(j,k)-A(i,k))/hyp;
end 
