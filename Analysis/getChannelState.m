function [b,acs]=getChannelState(SpecSense)

for i=1:1000
    for j=1:4
        if (SpecSense(i,j).SigInt > 0.3)
            b(i,j)=1;
        else 
            b(i,j)=0;
        end
    end
end