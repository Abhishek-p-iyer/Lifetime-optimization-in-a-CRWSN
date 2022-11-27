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
