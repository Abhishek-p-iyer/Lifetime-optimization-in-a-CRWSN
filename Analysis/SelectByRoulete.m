function ret = SelectByRoulete(allProb)
cumProb = cumsum(allProb);
rn = rand;
hd = find(cumProb >= rn);
ret = hd(1);
