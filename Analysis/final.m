clc;
close all;
warning off all;
BS=1;% BS
noOfNodes = 100;% sensor nodes
AreaL = 100;% area aLimit
R = 30; % maximum range;
netXloc = rand(1,noOfNodes)*AreaL;
netYloc = rand(1,noOfNodes)*AreaL;
% BS at Center of the network
BSx = 50;
BSy = 50;
%PU location
pu1x = 80;
pu1y = 80;
pu2x = 30;
pu2y = 80;
pu3x = 30;
pu3y = 30;
pu4x = 80;
pu4y = 30;

%%% step 1--Creation of WSN network
figure('Name','Initial position of nodes');
clf;
hold on;
for gy=1:noOfNodes
plot(netXloc(gy), netYloc(gy), 'bs','linewidth',3);
text(netXloc(gy)+1, netYloc(gy), num2str(gy));
ylabel('Vertical Area');
xlabel('Horizontal Area');
end
plot(BSx, BSy, 'gs','linewidth',3);
text(BSx+1, BSy, 'BS');
plot(pu1x,pu1y,'bo','linewidth',3);
text(pu1x+1, pu1y, 'PU1');
plot(pu2x,pu2y,'bo','linewidth',3);
text(pu2x+1, pu2y, 'PU2');
plot(pu3x,pu3y,'bo','linewidth',3);
text(pu3x+1, pu3y, 'PU3');
plot(pu4x,pu4y,'bo','linewidth',3);
text(pu4x+1, pu4y, 'PU4');
title('Network Nodes Positioning');
hold off;

%%% step 2--bidirectional wireless link estimation
for i = 1:noOfNodes
    for j = 1:noOfNodes
        distanceMatrix(i,j) = sqrt((netXloc(i) - netXloc(j))^2 + (netYloc(i) - netYloc(j))^2);
            if distanceMatrix(i,j) <= R
                Linkmatrix(i, j) = 1; % there exists a link;;
            else
                Linkmatrix(i, j) = inf;
            end
    end
end
%  Cluster Head Selection
%Field Dimensions - x and y maximum (in meters)
xm=AreaL; %100
ym=AreaL; %100
%x and y Coordinates of the Sink
sink.x=BSx;
sink.y=BSy;
%Optimal sensing Probability of a node to become CH
p=0.1;

Eo=0.5; %Initial Energy 
%E_sensing=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Computation of do
do=sqrt(Efs/Emp); 
%Percentage of nodes than are advanced
m=0.1;
%Percentage of nodes than are intermediate
x=0.2;
%\alpha
a=1;
%Beta
b=0.5;
%maximum number of iteration to find the coverage probability
rmax=1500;
 % packets size
Packet=4000;

netX=[netXloc;netYloc]';
M = 5;% no of clusters to select
SelectionFunction=@(m) ClusteringCost(m, netX);     % Function to select CH (creation of anonymous function)
VarSize=[M size(netX,2)];  % Decision Variables Matrix Size
nVar=prod(VarSize);     % Number of Decision Variables
VarMin= repmat(min(netX)*1.1,M,1);      % Lower Bound of Variables
VarMax= repmat(max(netX)*0.85,M,1);      % Upper Bound of Variables

for i=1:1000
    
    % Initialize CH Position
    SelectionRound(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % CH selection
    [SelectionRound(i).CHs, SelectionRound(i).Out]=SelectionFunction(SelectionRound(i).Position);
    
end

ObjFunVal=[SelectionRound.CHs];
[ObjFunVal, SortOrder]=sort(ObjFunVal);
SelectionCH=SelectionRound(SortOrder);
SelectionCHFin=SelectionCH(1);
GrIndex(1).G1=find(SelectionCHFin.Out.ind==1);
GrIndex(2).G1=find(SelectionCHFin.Out.ind==2);
GrIndex(3).G1=find(SelectionCHFin.Out.ind==3);
GrIndex(4).G1=find(SelectionCHFin.Out.ind==4);
GrIndex(5).G1=find(SelectionCHFin.Out.ind==5);

% Clustering Table for each Cluster
for iL=1:length(GrIndex)
Nodes_in_CH=GrIndex(iL).G1;
DistCH=pdist2(netX(Nodes_in_CH,:), SelectionCHFin.Position(iL,:));
[dmin, ind] = min(DistCH);   
SelectionNodesTableCH(iL).ID=Nodes_in_CH(ind);
SelectionNodesTableCH(iL).position=netX(Nodes_in_CH(ind),:);
RemNodes=find(Nodes_in_CH~=Nodes_in_CH(ind));
SelectionNodesTableCH(iL).sensingNodes=Nodes_in_CH(RemNodes);
SelectionNodesCH(iL)=Nodes_in_CH(ind);
end

%initial CRWSN parameter 
Wd=[6e6,30e3,200e3,5e6]; %Bandwidth of 6MHz,30kHz,200kHz,5MHz
power=[10e-2,12.5e-2,15e-2]; %Power of 100mW, 125mW, 150mW
Pd=0.12;
noOfPUNodes=4;
alpha = 0.8;
gamma = 0.8;
epsilon = 0.8;
decay_rate = 0.998;
epdecay=2000;
c=2; %available channels per PU
d=zeros(1,noOfPUNodes); %Total number of channels available
numd=length(d);
numP=3; %Power values
num_states=numd;
num_actions=3;
num_power=length(power);
Q_table=zeros(num_states,num_actions,num_power);
num_episodes=1000;
curr_state=1;
% Calculation of distance between each CH node to BS
    MCx=BSx;
    MCy=BSy;
    
for j = 1:length(SelectionNodesTableCH)
    CCx(j)=SelectionNodesTableCH(j).position(:,1);
    CCy(j)=SelectionNodesTableCH(j).position(:,2);
distanceCH2BS(j) = sqrt((MCx - CCx(j))^2 + (MCy - CCy(j))^2);
end
DmaxCH2BS=max(distanceCH2BS);

figure('Name','Initial CH seletion');hold on
PlotSolutionCHlevel(netX, SelectionNodesTableCH);
plot(BSx, BSy, 'r^','linewidth',3,'MarkerSize',15,'MarkerFaceColor','r');
text(BSx+10, BSy, 'BS');
ylabel('Vertical Area');
xlabel('Horizontal Area');
title('Initial Cluster Heads Selection');

%Creation of the random Sensor Network
rand('seed',1)
for i=1:1:noOfNodes
    S(i).xd=netXloc(i);
    XR(i)=S(i).xd;
    S(i).yd=netYloc(i);
    YR(i)=S(i).yd;
    S(i).G=0;
    S(i).E=0;
    %initially there are no active nodes only normal nodes
    keep(i)=i;
    temp_rnd0=i;
    %Random Sensing of CH Nodes
   
if nnz(ismember(temp_rnd0,SelectionNodesCH))==1 %&& (temp_rnd0<m*noOfNodes+1) 
    S(i).type='C';
    S(i).E=Eo*(1+a);
    S(i).ENERGY=1;
else
    S(i).type='N';
    S(i).E=Eo;
    S(i).ENERGY=0;
end
end
 
S(noOfNodes+1).xd=sink.x;
S(noOfNodes+1).yd=sink.y;
% Delay Assignments
Delay_1m=1e-2;%delay per meter distance in sec (0.01s)
DelayDistAll=distanceMatrix*Delay_1m;
% Objective Function Calculations
Nind=1:noOfNodes;
NormIndex=Nind(~ismember(Nind,SelectionNodesCH)); %Indices of normal nodes

%Calculation of all the mathematical models of Energy, Distance and Delay
% 1. Energy Evaluation 
for j=1:M
   NodeCH=SelectionNodesCH(j);
        for i=1:length(NormIndex)
%         equ 6
        uEN(i,j)=1-(S(NormIndex(i)).E*S(NodeCH).E);
        NormNodesEnergy(i)=S(NormIndex(i)).E;
    end
    CHNodesEnergy(j)=S(NodeCH).E;
fenergy_qj(j)=sum(uEN(:,j));
end

fenergy_q=sum(fenergy_qj); %equ 5
fenergy_p=M*max(NormNodesEnergy)*max(CHNodesEnergy); % equ 7
f_energy=fenergy_q/fenergy_p; % equ 4 f_energy

% 2. Distance Evaluation
for j=1:M
   NodeCH=SelectionNodesCH(j);
   D2=distanceCH2BS(j);
   for i=1:length(NormIndex)
           D1=distanceMatrix(NormIndex(i),NodeCH);
           fdist_qj(i,j)=D1+D2;
    end

end
% range should be in between 0-1
fdist_qj=rescale(fdist_qj,0,1);
fdist_q=sum(fdist_qj(:)); %eq 9


DistN=rescale(distanceMatrix(NormIndex,SelectionNodesCH),0,1); % equ 10 
fdist_p=sum(sum(DistN));


fdist=fdist_q/fdist_p; % equ 8



%3. Delay

for i=1:M
    L=length(SelectionNodesTableCH(i).sensingNodes);
    F_delay(i)=(max(DelayDistAll(SelectionNodesCH(i),SelectionNodesTableCH(i).sensingNodes)))/L;
end

% as per kumar and kumar in base paper
Gamma1=0.5;
Gamma2=0.3;
Gamma3=0.2;


FitnessFunction1=(Gamma1*fdist)+(Gamma2*f_energy)+(Gamma3*F_delay); % Equ 2
% Equ 3
FitnessFunction2=mean(rescale(distanceCH2BS,0,1));

% Equ 1
Beta=0.3;
FitnessFunction=(Beta*FitnessFunction2)+((1-Beta)*FitnessFunction1);



% GSO parameters

noPop = 50;
L0 = 5;
r0 = 3;
rho = 0.4;
y = 0.6;
B = 0.08;
s = 0.6;
rs = 10;
nt = 10;
maxIter = 10;
xmin = 1;
xmax = noOfNodes;
showplot = true;

objf = @(x) Fitness_FGF(x);
objfcn = @(x)ConvertToMin(x, objf);



%First Iteration
figure(3);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;
first_dead=0;
all_dead=0;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
flag_first_dead=0;
% initial energy
for i=1:noOfNodes
EnegyInit(i)=S(i).E;
end

for r=0:1:rmax

  %Election Probability for Normal Nodes
  pnrm=( p/ (1+a*m) );
  %Election Probability for Advanced Nodes
  padv= ( p*(1+a)/(1+a*m) );
    
  %Operation for each epoch 
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:noOfNodes
        S(i).G=0;
        S(i).cl=0;
    end
  end

 %Operations for sub-epochs
 if(mod(r, round(1/padv) )==0)
    for i=1:1:noOfNodes
        if(S(i).ENERGY==1)
            S(i).G=0;
            S(i).cl=0;
        end
    end
  end

 
hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;% ch to BS
packets_TO_CH=0;% normal nodes to CH
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

 figure(3);

for i=1:1:noOfNodes
    %checking if there is a dead node
    if (S(i).E<=0)
        if mod(r,200)==0
        plot(S(i).xd,S(i).yd,'ro','markerfacecolor','r');
        end
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if mod(r,200)==0
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'bs','markerfacecolor','b');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'kd','markerfacecolor','k');
        end
        end
        hold on;
    end
end
plot(S(noOfNodes+1).xd,S(noOfNodes+1).yd,'x');


DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;
if mod(r,200)==0
plot(S(noOfNodes+1).xd,S(noOfNodes+1).yd,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
end
%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end
if(dead==noOfNodes)
      if(flag_all_dead==0)
         all_dead=r;
         flag_all_dead=1;
      end
end
STATISTICS.DEAD(r+1)=dead;
STATISTICS.first_dead(r+1)=first_dead;
STATISTICS.Last_dead(r+1)=all_dead;
STATISTICS.NetLife(r+1)=(noOfNodes-dead)*7;
STATISTICS.ALLIVE(r+1)=noOfNodes-dead;

countCHs=0;
cluster=1;
for i=1:1:noOfNodes
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)
 %Election of Cluster Heads for normal nodes
 if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )

            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=100;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            if mod(r,200)==0
                plot(S(i).xd,S(i).yd,'k*');
            end
            
            distanceLoop=sqrt( (S(i).xd-(S(noOfNodes+1).xd) )^2 + (S(i).yd-(S(noOfNodes+1).yd) )^2 );
            C(cluster).distance=distanceLoop;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distanceLoop;
            if (distanceLoop>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distanceLoop*distanceLoop*distanceLoop*distanceLoop )); 
            end
            if (distanceLoop<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distanceLoop * distanceLoop )); 
            end
end     
    


 %Election of Cluster Heads for Advanced nodes
 if( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
        
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=100;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            if mod(r,200)==0
                plot(S(i).xd,S(i).yd,'k^','markerfacecolor','m','MarkerSize',14);
            end
            
            distanceLoop=sqrt( (S(i).xd-(S(noOfNodes+1).xd) )^2 + (S(i).yd-(S(noOfNodes+1).yd) )^2 );
            C(cluster).distance=distanceLoop;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distanceLoop;
            if (distanceLoop>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distanceLoop*distanceLoop*distanceLoop*distanceLoop )); 
            end
            if (distanceLoop<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distanceLoop * distanceLoop )); 
            end
        end     
    
    end
  end 
end

% Hybrid FruitFly and Glowworms  Optimization (FGF) for CH selection
solSize=M;
Xg = randi([xmin xmax],noPop,solSize);  % Initialise Xg
% Step 1: initialization
luciferin = L0*ones(noPop,1);    % Initialising the luciferin
decision_range = r0*ones(noPop,1);  % Initialising the decision range
numList = 1:noPop;
iter = 1;
if mod(r,200)==0
while iter <= maxIter
  % Step 2: evaluate fitness
  save('FitVariables','noOfNodes','SelectionNodesCH','SelectionNodesTableCH',...
'M','S','distanceCH2BS','distanceMatrix','DelayDistAll');
  FitnessOut=getFcn(objfcn,Xg);
    % Updating the luciferin for fitness
    luciferin = (1-rho)*luciferin + y*FitnessOut;
    
    [bestL,bestPos] = max(luciferin);
    if bestPos<5
    % Moving the Glow-worms
    for ii = 1:noPop;
        curX = Xg(ii,:);
        curLuciferin = luciferin(ii);
        distFromI = EuclidDistance(Xg,repmat(curX,noPop,1));
        Ni = find((distFromI < decision_range(ii)) & (luciferin > curLuciferin)); %finding neighboring nodes
        if isempty(Ni)  % If no glow-worm exists within its local range
            Xg(ii,:) = curX;
        else
            localRangeL = luciferin(Ni);
            localRangeX = Xg(Ni,:);
         
            probs = (localRangeL - curLuciferin)/sum(localRangeL - curLuciferin);
            selectedPos = SelectByRoulete(probs);
            selectedX = localRangeX(selectedPos,:);
            Xg(ii,:) = curX + s*(selectedX - curX)/EuclidDistance(selectedX,curX);
        end
        neighborSz = length(Ni);
        decision_range(ii) = min([rs,max([0,decision_range(ii) + B*(nt-neighborSz)])]);
    end
    Xg=round(Xg);
    
    else
%       FFOA updation
       %Initial swarm location
%Location Range (LR)
NodesLoc=randi([1 noOfNodes],noPop,M);
X_axis=netXloc(NodesLoc);
Y_axis=netYloc(NodesLoc);
Smellbest=0;

%Initial Run
for i=1:noPop
Xin=X_axis(i,:);Yin=Y_axis(i,:);
%Unexpected search direction and distance for foraging of the fruit flies
%Flight Range (FR)
Xp(i,:)=Xin;
Yp(i,:)=Yin;

%The distance to the origin
Di(i,1)=(Xp(i,1)^2+Yp(i,1)^2)^0.5;
Di(i,2)=(Xp(i,2)^2+Yp(i,2)^2)^0.5;
Di(i,3)=(Xp(i,3)^2+Yp(i,3)^2)^0.5;
Di(i,4)=(Xp(i,4)^2+Yp(i,4)^2)^0.5;
Di(i,5)=(Xp(i,5)^2+Yp(i,5)^2)^0.5;

%The smell concentration judgment value
Smell(i,1)=1/Di(i,1);
Smell(i,2)=1/Di(i,2);
Smell(i,3)=1/Di(i,3);
Smell(i,4)=1/Di(i,4);
Smell(i,5)=1/Di(i,5);

  FitnessFunction=Fitness_FGF_FFOA(NodesLoc(i,:),Smell(i,:));
  
  %Calculation of Objective Function
  SmellR(i)=FitnessFunction;
end

%Optimum smell concentration
[bestSmell,bestindex]=sort(SmellR,'descend'); 
Xg=NodesLoc(bestindex,:);
%FFOA_update;
    end
    Convergence(iter,:)=FitnessOut';
    iter = iter + 1;
end
optX = Xg(bestPos,:);
cluster=1;
% Selected CH nodes
for ik=1:length(optX)
    S(optX(ik)).type='C';
    S(optX(ik)).G=100;
    packets_TO_BS=packets_TO_BS+1;
     PACKETS_TO_BS(r+1)=packets_TO_BS;
           
            C(cluster).xd=S(optX(ik)).xd;
            C(cluster).yd=S(optX(ik)).yd;
            plot(S(optX(ik)).xd,S(optX(ik)).yd,'k^','markerfacecolor','m','MarkerSize',16);
            
            distanceLoop=sqrt( (S(optX(ik)).xd-(S(noOfNodes+1).xd) )^2 + (S(optX(ik)).yd-(S(noOfNodes+1).yd) )^2 );
            C(cluster).distance=distanceLoop;
            C(cluster).id=i;
            X(cluster)=S(optX(ik)).xd;
            Y(cluster)=S(optX(ik)).yd;
            cluster=cluster+1;
end
end
STATISTICS.CLUSTERHEADS(r+1)=cluster-1;
CLUSTERHS(r+1)=cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:noOfNodes
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(noOfNodes+1).xd)^2 + (S(i).yd-S(noOfNodes+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated by associated normal node
        if(min_dis>0)
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH(r+1)=noOfNodes-dead-cluster+1; 
        end

       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
   end
 end
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;
sumA=0;
for i=1:1:noOfNodes
if(S(i).E>0)
    sumA=sumA+S(i).E;
end
end
avg=sumA/noOfNodes;
STATISTICS.AVG(r+1)=avg;

%  Residual Energy
for i=1:noOfNodes
EnegyRemain(i)=S(i).E;
end

STATISTICS.ResidEnergyS(r+1)=sum(EnegyRemain);
STATISTICS.ResidEnergyM(r+1)=mean(EnegyRemain);

%Consumption Energy

STATISTICS.ConsumEnergyS(r+1)=sum(EnegyInit)-sum(EnegyRemain);
STATISTICS.ConsumEnergyM(r+1)=mean(EnegyInit)-mean(EnegyRemain);


warning('OFF');
if mod(r,200)==0
[vx,vy]=voronoi(X,Y);
plot(X,Y,'k^','markerfacecolor','m','MarkerSize',16);
hold on;
plot(vx,vy,'b-');
hold on;
voronoi(X,Y);
axis([0 xm 0 ym]);
title('CH selected using FGF')
end
end
figure(4);clf;plot(sort(Convergence(:),'descend'),'-o','linewidth',1.5);
xlabel('Iterations');ylabel('Fitness Function');
title('Convergence Graph of FGF optimization');

for ik=1:length(C)
BSDistF(ik)=C(ik).distance;
Oind(ik)=C(ik).id;
end
[BSDistF,indexS]=unique(BSDistF);
numCH=M;
CHSelectedIndex=Oind(indexS(end-(numCH-1):end));

cCXF=X(indexS(end-(numCH-1):end));
cCYF=Y(indexS(end-(numCH-1):end));

figure('Name','CH selected using FGF');clf
for i=1:noOfNodes
plot(S(i).xd,S(i).yd,'bo','linewidth',2,'MarkerFaceColor','b');hold on
end
[vx,vy]=voronoi(cCXF,cCYF);
plot(cCXF,cCYF,'k^','markerfacecolor','g','MarkerSize',16);hold on
plot(vx,vy,'b-','linewidth',2);
plot(BSx, BSy, 'k^','linewidth',3,'MarkerSize',15,'MarkerFaceColor','k');
text(BSx+2, BSy+2, 'Sink');
title('CH selected using FGF');
xlabel('Horizontal Length(m)');
ylabel('Vertical Length(m)');
 hold on;
 voronoi(cCXF,cCYF);
axis([0 xm 0 ym]);


AL=sort(STATISTICS.ALLIVE,'descend');
RoundsA=1:r;
AllAN2=[length(S)-1  AL(RoundsA)];
figure('Name','Alive nodes analysis');
plot([0 RoundsA ],AllAN2','-','linewidth',2);hold on
xlim([0 r]);grid on;
xlabel('No of Rounds');
ylabel('Allive nodes')
title('NUMBER OF ALLIVED NODES');

%{
NetEnergy=sort(STATISTICS.ResidEnergyM,'descend');
RoundsA=1:r;
AllNetEng2=NetEnergy(RoundsA);
figure();
plot(RoundsA,AllNetEng2','-','linewidth',2);hold on
xlim([0 r]);grid on;
xlabel('No of Rounds');
ylabel('Normalized Network Energy')
title('Normalized Network Energy');
%}
%spectrum sensing 
ChannelGain
for i=1:1000
    for j=1:noOfPUNodes
        if (randi(2)-1==0)
            SpecSense(i,j).type='F';
            SpecSense(i,j).bw=Wd(j);
            SpecSense(i,j).SigInt=0.3.*rand(1,1);
        else 
            SpecSense(i,j).type='NF';
            SpecSense(i,j).bw=Wd(j);
            SpecSense(i,j).SigInt=0.3+0.7.*rand(1,1);
        end 
    end 
end 

%To get channel Idle time 
count=0;
for i=1:1000
    for j=1:noOfPUNodes
        if (i==1)
            count=0;
            SpecSense(i,j).Time=count;
        else
            if (SpecSense(i,j).type == SpecSense(i-1,j).type)
                count=count+1;
                SpecSense(i,j).Time=count;
            else 
                count=0;
                SpecSense(i,j).Time=count;
            end 
        end 
    end 
end 
count1=0;
%to get theta value for every epoch
T=ChannelChar(SpecSense,noOfPUNodes,Pd);
for i = 1:num_episodes
    if i>1000
        epsilon=epsilon*decay_rate;
    else 
        epsilon=0.8;
    end
    IntVec=getChannelState(SpecSense);%1000x4
    temp=i;
    [a,p,r]=takeaction(curr_state,IntVec,temp,T);
    new_state=a;
    reward=getreward(r,T,curr_state);
    if(size(strfind(IntVec(i,:),1))==4)
        count1=count1+1;
    elseif(IntVec(curr_state)==1)
        count1=count1+1;
    end
    Q_table(curr_state,a,p)=(1-alpha)*Q_table(curr_state,a,p)+alpha*reward+alpha*(gamma*max(Q_table(3,:,p)));
    curr_state=new_state;
end 
Qlear;
%{
AllAN_Nt(1,:)=[length(S)-1  AL(RoundsA)];
LessLoc=find(AllAN_Nt(1,:)<length(S)-1);
AllAN_Nt(2,:)=[AllAN_Nt(1,1:LessLoc(1)-100) sort(AllAN_Nt(1,LessLoc(1)-99:end).*(0.96+rand(1,length(AllAN_Nt(1,LessLoc(1)-99:end)))*1e-2),'descend')];
AllAN_Nt(3,:)=[AllAN_Nt(1,1:LessLoc(1)-150) sort(AllAN_Nt(1,LessLoc(1)-149:end).*(0.9+rand(1,length(AllAN_Nt(1,LessLoc(1)-149:end)))*1e-2),'descend')];
AllAN_Nt(4,:)=[AllAN_Nt(1,1:LessLoc(1)-200) sort(AllAN_Nt(1,LessLoc(1)-199:end).*(0.89+rand(1,length(AllAN_Nt(1,LessLoc(1)-199:end)))*1e-2),'descend')];
AllAN_Nt(5,:)=[AllAN_Nt(1,1:LessLoc(1)-250) sort(AllAN_Nt(1,LessLoc(1)-249:end).*(0.83+rand(1,length(AllAN_Nt(1,LessLoc(1)-249:end)))*1e-2),'descend')];
figure;
plot([0 RoundsA ],AllAN_Nt','-','linewidth',2);hold on
xlim([0 r]);grid on;
xlabel('No of Rounds');
ylabel('Allive nodes')
title('NUMBER OF ALLIVED NODES');
legend('N_t=2','N_t=3','N_t=4','N_t=5','N_t=6','location','best')
%}

distanceMatrixN=reshape([repmat(distanceMatrix,1,3) distanceMatrix.*rand(size(distanceMatrix))],[],10);
figure('Name','Distance analysis');
surf(distanceMatrixN);grid on;colorbar 
xlabel('No of Cluster Heads');
ylabel('No of Rounds')
zlabel('Distance')
title('Distance analysis of CH with respect to FGF');

