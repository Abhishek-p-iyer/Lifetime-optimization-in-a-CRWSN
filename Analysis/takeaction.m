function [a,p,r]=takeaction(curr_state,IntVec,temp,T)
    if (curr_state==0 || IntVec(temp,curr_state)==1) %Switching channels
        available_channel=strfind(IntVec(temp,:),0); 
        tf=isempty(available_channel);
        if(tf)
            curr_state=randi(4);
            a= 1;
            r=-1;
            p=randi(3);
        else
            curr_state=available_channel(1); 
            a=2;
            r=-0.1;
            p=randi(3);
        end
    else
        a=curr_state; %Normal communication
        r=T(temp,curr_state);
        p=randi(3);
    end
end