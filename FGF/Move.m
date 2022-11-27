function Move(i,j)
global A m step1 LuciferinLevel bound 

if((j~=0) && (LuciferinLevel(i)<LuciferinLevel(j)))
    temp(i,:)=A(i,:)+step1*Path(i,j);
    flag=0;
    for k=1:m
        if(temp(i,k)<-bound)||(temp(i,k)>bound)
            flag=1;
            break
        end 
    end 
    if(flag==0)
        A(i,:)=temp(i,:);
    end 
end 
