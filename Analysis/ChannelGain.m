FGFy_throughput=pv(:,1);
xaxis=(1:1:1500);
figure
plot(xaxis,FGFy_throughput,'-r','linewidth',2, 'DisplayName','FGF'); 
xlabel('Number of rounds');
ylabel('Packets sent to BS');
title('Throughput vs Number of rounds');




