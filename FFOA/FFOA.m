function [Smellbest,X,Y] = FOA(n,maxt,lb,ub,dim)   

% Parameters setting
if nargin < 1
    n = 20; % Population size
    maxt = 5e2; % Max iterations
    dim = 30; % Dimension of test function
    lb = -100 * ones(1,dim); % Lower bound of test function
    ub = 100 * ones(1,dim); % Upper bound of test function
end
% X = zeros(1 * dim);
% Y = zeros(1 * dim);
% new_X = zeros(1 * dim);
% new_Y = zeros(1 * dim);
% D = zeros(1 * dim);
% Sol = zeros(1 * dim);
% Fitness = zeros(n * 1);

% Initialize the original position
for i = 1:n
    X(i,:) = lb+(ub-lb).*rand(1,dim); % the position of X axis
    Y(i,:) = lb+(ub-lb).*rand(1,dim); % the position of Y axis
    D(i,:) = (X(i,:).^2 + Y(i,:).^2).^0.5; % Caculate the distance
    Sol(i,:) = 1./D(i,:); % the solution set
    Fitness(i) = fun(Sol(i,:)); % Caculate the fitness
end


[bestSmell,index] = min(Fitness); % Get the min fitness and its index
new_X = X(index,:); % the X axis of min fitness
new_Y = Y(index,:); % the Y axis of min fitness
Smellbest = bestSmell;
best = Sol(index,:);

% Start main loop
for t = 1:maxt
    for i = 1:n
        % Refer to the process of initializing
        X(i,:) = new_X + (ub - lb).*rand(1,dim);
        Y(i,:) = new_Y + (ub - lb).*rand(1,dim);
        D(i,:) = (X(i,:).^2 + Y(i,:).^2).^0.5;
        Sol(i,:) = 1./D(i,:);
        Fitness(i) = fun(Sol(i,:));
    end
    [bestSmell,index] = min(Fitness);
    % If the new value is smaller than the best value,update the best value
    if (bestSmell < Smellbest)
        X(i,:) = X(index,:);
        Y(i,:) = Y(index,:);
        Smellbest = bestSmell;
    end
    
    % Out put result each 100 iterations
    if round(t/100) == (t/100)
        Smellbest;
    end
    
    cg_curve(t) = Smellbest;
end

% Output/display
disp(['Number of evaluations: ',num2str(maxt)]);
disp(['Best solution=',num2str(best),'   fmin=',num2str(Smellbest)]);

% Draw the picture
semilogy((1:25:maxt),cg_curve(1:25:maxt),'k-o','markersize',5);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fruit fly (score) obtained so far');

hold on
axis tight
grid off
box on
legend('FOA')

% This is a classcial test function,namely Sphere function,which range is
% from -100 to 100.The dimension can be defined as you want.
function z = fun(u)
z = sum(u.^2);
