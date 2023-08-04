
function Fx = getFcn(objfcn,Xs)
n = size(Xs,1);
Fx = ones(n,1);
for k = 1:n
    Fx(k) = objfcn(Xs(k,:));
end