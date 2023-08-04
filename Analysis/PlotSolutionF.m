
function GroupsAll=PlotSolutionF(X, sol)

    % Cluster Centers
    m = sol.Position;
    k = size(m,1);
    
    % Cluster Indices
    ind = sol.Out.ind;
    
    Colors = hsv(k);
    
    for j=1:k
        Xj = X(ind==j,:);
        GroupsAll(j).Coords=Xj;
        plot(Xj(:,1),Xj(:,2),'o','LineWidth',2,'Color',Colors(j,:));
        hold on;

    end
    
     
%     hold off;
    grid on;
    
end