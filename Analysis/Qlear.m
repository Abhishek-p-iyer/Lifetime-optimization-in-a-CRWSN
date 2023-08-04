function Qlear
i=zeros(1,1000);
count1=0;
count2=0;
count3=0;
count4=0;
count0=0;
for c=1:1000
    if(c<132)
        i(c)=randi(5);
    elseif(c>132 && c<207)
        i(c)=randi(6);
    elseif(c>207 && c<281)
        i(c)=randi(5);
    elseif(c>281 && c<532)
        i(c)=randi(4);
    elseif(c>532 && c<561)
        i(c)=randi(3);
    elseif(c>561 && c<621)
        i(c)=randi(2);
    elseif(c>621 && c<660)
        i(c)=randi(1);
    elseif(c>660 && c<711)
        i(c)=randi(3);
    elseif(c>711 && c<722)
        i(c)=randi(2);
    elseif(c>722 && c<882)
        i(c)=randi(4);
    elseif(c>882 && c<861)
        i(c)=randi(1);
    elseif(c>861 && c<921)
        i(c)=0;
    elseif(c>921 && c<936)
        i(c)=randi(2);
    elseif(c>936 && c<1000)
        i(c)=0;
    end
end

for c=1:1000
    if(i(c)>0)
        i(c)=i(c)-1;
    end
end
for c=1:1000
    if i(c)==1
        count1=count1+1;
    elseif i(c)==2
        count2=count2+1;
    elseif i(c)==3
        count3=count3+1;
    elseif i(c)==0
        count0=count0+1;
    else 
        count4=count4+1;
    end
end
formatSpec = '%1i channels is being interfered %1i times\n';
fprintf('Zero channels are being interfered %1i times\n',count0);
fprintf(formatSpec,1,count1);
fprintf(formatSpec,2,count2);
fprintf(formatSpec,3,count3);
fprintf(formatSpec,4,count4);
figure('Name','Number of interfering channls')   
plot(1:1000,i,'r');
axis([0,1000,0,8]);
xlabel('Time');
ylabel('Number of interfering channels');
title('Number of Interfering channels Vs Time');


        

