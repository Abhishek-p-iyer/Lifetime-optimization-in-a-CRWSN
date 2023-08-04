function ret = EuclidDistance(pos1,pos2)
ret = sqrt(sum((pos1-pos2).^2,2));
