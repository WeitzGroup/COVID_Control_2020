%% ga algo to find 2 params - y intercept and slope; y -mx -c = 0

% c must be less than 1
% 
Aeq = [];
beq = [];

A = [];
b = [];

lb = 10;
ub = 500;
%% no intercept
parsM.cB = 5;
isolation_efficiency = 0.25:0.25:0.75;
%isolation_efficiency = 0.65;
shield_efficiency = 2:5;
x_f = zeros(length(isolation_efficiency),length(shield_efficiency),1);
%load('ga_line_params_constrained.mat')
% x_f(2,1,:) = [0;0];
% x_f(5,1,:)=[0;0];
% x_f(6,1,:)=[0;0];
vector_test = [1, 2];
for i = vector_test
    c_min = (1-isolation_efficiency(i))*parsM.cB;
    i
    for j = 1:3:4
        j
%         if (x_f(i,j,1)^2 + x_f(i,j,2)^2)
%             continue
%         end
        c_max = shield_efficiency(j)*parsM.cB;
        fun = @(line_params) heuristic_cost_new(line_params, c_min, c_max);
        % set max generation
        max_iter = 30;
        options = optimoptions('ga','MaxGenerations',max_iter,'PlotFcn', @gaplotbestf);
        rng default % For reproducibility
        x_f(i,j,:) = ga(fun,1,A,b,Aeq,beq,lb,ub,[],options);
        %x_f(i,j,:) = ga(fun,1,A,b,Aeq,beq,[],[],[],options); % no intercept, one variable
        save('ga_line_params_constrained_no_intercept_360_days.mat','x_f')
    end   
end