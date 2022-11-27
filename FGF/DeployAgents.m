function DeployAgents
global n m A_init A bound
B=-bound*ones(n,m);
A_init = B+2*bound*rand(n,m);
A=A_init;
