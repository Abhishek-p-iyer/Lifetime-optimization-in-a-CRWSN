
function GroupsAll=PlotSolutionCHlevel(X, sol)

    % Cluster Centers
    k = size(sol,2);
    
    Colors2=hsv(k);
      
    for j=1:k
        m = sol(j).position;
    % Cluster Indices
    ind = sol(j).sensingNodes;
    
        Xj = X(ind,:);
        GroupsAll(j).Coords=Xj;
         plot(Xj(:,1),Xj(:,2),'ko','LineWidth',2,'MarkerFaceColor',Colors2(j,:));
         text(Xj(:,1)+1,Xj(:,2)+1, num2str(ind));

        for xx=1:size(Xj,1)
            plot([m(1) Xj(xx,1)],[m(2) Xj(xx,2)],'--','LineWidth',1.2,'Color',Colors2(j,:));
        end
        hold on;
         plot(m(1),m(2),'s','LineWidth',2,'Color','k','MarkerSize',10,'MarkerFaceColor','m');
  text(m(1)+5, m(2), ['CH' num2str(j) '(' num2str(sol(j).ID) ')']);

    end
    
     
%     hold off;
    grid on;
    
end