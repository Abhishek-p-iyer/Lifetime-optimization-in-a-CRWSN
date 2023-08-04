function reward = getreward(r,T,curr_state)

if (r==1)
    reward=-1;
elseif (r==2)
    reward=-0.1;
else
    reward=T(curr_state);
end
