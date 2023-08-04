function FitnessFunctionN=Fitness_FGF_FFOA(Xg,Smell)

% Objective Function Calculations
load('FitVariables')
SelectionNodesCH=Xg;
  MCx=S(noOfNodes+1).xd;
  MCy=S(noOfNodes+1).yd;
    
for j = 1:length(Xg)
    CCx(j)=S(Xg(j)).xd;
    CCy(j)=S(Xg(j)).yd;
distanceCH2BS(j) = sqrt((MCx - CCx(j))^2 + (MCy - CCy(j))^2);
end

% Normal nodes index
Nind=1:noOfNodes;
NormIndex=Nind(~ismember(Nind,SelectionNodesCH));

% 1. Energy Evaluation

% equ 5 & 6
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
%         equ 5
fenergy_q=sum(fenergy_qj);

% equ 7
fenergy_p=M*max(NormNodesEnergy)*max(CHNodesEnergy);

% equ 4 f_energy

f_energy=fenergy_q/fenergy_p;


% 2. Distance Evaluation

% equ 9
for j=1:M
   NodeCH=SelectionNodesCH(j);
%    CH-BS distance
   D2=distanceCH2BS(j);

    for i=1:length(NormIndex)
% normal node to CH distance
           D1=distanceMatrix(NormIndex(i),NodeCH);
%         equ 9
           fdist_qj(i,j)=D1+D2;
    end

end
% range sould be 0-1
fdist_qj=rescale(fdist_qj,0,1);
%         equ 9 sum
fdist_q=sum(fdist_qj(:));

% equ 10
DistN=rescale(distanceMatrix(NormIndex,SelectionNodesCH),0,1);
% DistN=distanceMatrix(NormIndex,SelectionNodesCH);
fdist_p=sum(sum(DistN));

% equ 8
fdist=fdist_q/fdist_p;



%  delay

for i=1:M
    F_delay(i)=(max(DelayDistAll(SelectionNodesCH(i),NormIndex)))/length(NormIndex);

end
F_delayM=min(F_delay);
% as per (kumar and kumar 2015 reference)
Gamma1=0.5;
Gamma2=0.3;
Gamma3=0.2;

% Equ 2
FitnessFunction1=(Gamma1*fdist)+(Gamma2*f_energy)+(Gamma3*F_delayM);

% Equ 3
FitnessFunction2=mean(rescale(distanceCH2BS,0,1));

% Equ 1
Beta=0.3;
FitnessFunction=(Beta*FitnessFunction2)+((1-Beta)*FitnessFunction1);
FitnessFunctionN=FitnessFunction*mean(Smell(:));
