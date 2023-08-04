function theta = ChannelChar(SpecSense,noOfPUNodes,Pd)
epsilon1 = 0.2;
epsilon2 = 0.5;
epsilon3 = 0.3;

for i=1:1000
    for j=1:noOfPUNodes
        theta(i,j) = (SpecSense(i,j).bw.*epsilon1+SpecSense(i,j).SigInt.*epsilon2+(1-Pd).*epsilon3.*SpecSense(i,j).Time);   
    end 
end
