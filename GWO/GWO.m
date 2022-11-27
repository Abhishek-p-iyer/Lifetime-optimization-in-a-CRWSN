clc
clear all 
close all 

global m n rs rd rmin gamma ro beta bound nt LuciferinLevel step1 A A_init 
m=2;  %Number of dimensions 
n=200; %Number of agents
rs=3;  %Sensor range
rd=rs*ones(n,100);  %local decision range 
rmin=0; %threshold decision range 
gamma=0.6; %luciferin enhancement constant 
ro=0.4; %luciferin decay constant 
step1=0.03; %Distance moved by each glowworm after decision taken
beta=0.08; %decision range gain 
nt=5; %desired number of neighbours 
 
bound=3; %workspace range 
DeployAgents; %deploy GW randomly 
LuciferinLevel=5*ones(n,1); %Initialization of luciferin levels 
j=1; 
iter=250;
Ave_d=zeros(iter,1); %Average distance 

while(j<=iter)
    UpdateLuciferin; %To update the luciferin level at the current position
    Act; %Select a direction to move
    for k=1:n  
        agent_x(k,j,:)=A(k,1); %store the values 
        agent_y(k,j,:)=A(k,2);
    end 
    j=j+1;
end 

plot(A_init(:,1),A_init(:,2),'x');
xlabel('X');
ylabel('Y');
hold on;
DefineAxis;
for k=1:n
    plot(agent_x(k,:,:),agent_y(k,:,:));
end
DefineAxis;
grid on;
hold on;
plot([-0.009982;-0.4792;1.285],[1.562;-0.609;0.02765],'ok');

figure(2);
plot(A(:,1),A(:,2),'.');
DefineAxis;
grid on;
hold on;
plot([-0.009982;-0.4792;1.285],[1.562;-0.609;0.02765],'.');

function DeployAgents
global n m A_init A bound
B=-bound*ones(n,m);
A_init = B+2*bound*rand(n,m);
A=A_init;
end

function UpdateLuciferin
global n A J LuciferinLevel gamma ro 
for i=1:n
    x=A(i,1);
    y=A(i,2);
    J(i,:)=3*(1-x)^2*exp(-(x^2)-(y+1)^2)-10*(x/5-x^3-y^5)*exp(-x^2-y^2)-1/3*exp(-(x+1)^2-y^2);
    LuciferinLevel(i,:)=(1-ro)*LuciferinLevel(i,:)+gamma*J(i,:);
end 
end

function Act 
global n rs rd N Na beta nt 
N(:,:)=zeros(n,n);
Na(:,:)=zeros(n,1);

for i=1:n
    FindNeighbors(i);
    FindProbabilities(i);
    Leader(i)=SelectAgents(i);
end 

for i=1:n
    j=Leader(i);
    Move(i,j);
    rd(i)=max(0,min(rs,rd(i)+beta*(nt-Na(i))));
end 
end

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
end 

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
end 

function j=SelectAgents(i)
global n pb
bound_lower=0;
bound_upper=0;
toss=rand;
j=0;
for k=1:n
    bound_lower=bound_upper;
    bound_upper=bound_upper+pb(i,k);
    if((toss>bound_lower) && (toss<bound_upper))
        j=k;
        break;
    end 
end 
end

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
end

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
end 


function DefineAxis
global bound
axis([-bound bound -bound bound]);
grid on;
end 

