
global m n rs rd rmin gamma ro beta bound nt LuciferinLevel step1 
m=2;  %Number of dimensions 
%dim = 30;
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
DeployAgents;
LuciferinLevel=5*ones(n,1); %Initialization of luciferin levels 
max_iter = 250;
t=1;

for i = 1:n
    X(i,:)=A(i,1);
    Y(i,:)=A(i,2);
    D(i,:)=(X(i,:).^2+Y(i,:).^2).^0.5;
    Sol(i,:)=1./D(i,:);
    Fitness(i)=fun(Sol(i,:));
end

[bestSmell,index] = min(Fitness); % Get the min fitness and its index

if(index>5)
    new_X = X(index,:); % the X axis of min fitness
    new_Y = Y(index,:); % the Y axis of min fitness
    Smellbest = bestSmell;
    best = Sol(index,:);
else 
    while(t<=max_iter)
        for i=1:n
            x=A(i,1);
            y=A(i,2);
            J(i,:)=3*(1-x)^2*exp(-(x^2)-(y+1)^2)-10*(x/5-x^3-y^5)*exp(-x^2-y^2)-1/3*exp(-(x+1)^2-y^2);
            LuciferinLevel(i,:)=(1-ro)*LuciferinLevel(i,:)+gamma*J(i,:);
        end 
        N(:,:)=zeros(n,n);
        Na(:,:)=zeros(n,1);
        for i=1:n
            FindNeighbors(i);
            FindProbabilites(i);
            Leader(i)=SelectAgents(i); 
        end 
        for i=1:n
            j=Leader(i);
            Move(i,j);
            rd(i)=max(0,min(rs,rd(i)+beta*(nt-Na(i))));
        end 
        t=t+1;
    end
end
function z = fun(u)
z = sum(u.^2);
end 
