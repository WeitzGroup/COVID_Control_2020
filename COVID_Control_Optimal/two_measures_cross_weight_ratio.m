%% Parameters
% model parameters
parsM.Tinc = 4; % length of incubation period
parsM.Tinf = 6; % duration patient is infectious
parsM.etaI = 0.1; % transmission effectivenesss, note transmission rate beta ~ etaI*contactRate
parsM.mu = 1e-2; % case fatality ratio
parsM.cB = 5; % baseline contact rate
parsM.Ntot = 1e6; % total number of population (neglect death)
parsM.numSVar = 5; % number of state variables
parsM.numCVar = 3; % number of control variables

% control regulator parameters
parsC.W1 = 1; % weight - social connections
parsC.Wd = 1e5/parsM.Ntot; % scaled weight; % scaled weight
parsC.WI = 1e4/parsM.Ntot; % scaled weight; % scaled weight
parsC.WIf = 0; % final infected cases weight
parsC.WEf = 0; % final exposed cases weight


parsC.W2 = (1/parsM.cB^2).*eye(4);  % weight matrix
parsC.Ws = parsC.W2(1,1);
parsC.We = parsC.W2(2,2);
parsC.Wi = parsC.W2(3,3);
parsC.Wr = parsC.W2(4,4);


parsC.K = (10/(parsM.cB*parsM.Ntot)); % scale of exponential cost

% contacting rate range (policy efficiency)
multiple = 5;
parsC.cmax = parsM.cB*multiple; % constrain;
parsC.cmin = parsM.cB/multiple; % constrain; if cmin is very low (0.1), decreasing I is enough

% optimization algorithm parameters (backtracking step size)
parsC.alpha = 0.1; % alpha is *very* important 
parsC.beta = 0.5;
% numerical solver parameters;
parsT.dt = 5e-2;



log_weight_ratio_list = -4:0.5:4;
total_deaths_list = zeros(length(log_weight_ratio_list), 1);
work_frac_list = zeros(length(log_weight_ratio_list), 1);

for i = 1:length(log_weight_ratio_list)
    
% weight balancing
parsC.wr = 10^log_weight_ratio_list(i); % weight ratio
parsC.Wd = parsC.wr*(1e5/parsM.Ntot);
parsC.WI = parsC.wr*(1e4/parsM.Ntot);
parsC.WIf = 0;


% piecewise constant strategy or continouse open loop strategy
parsC.discrete = 1; % 0 - continous, 1 - discrete

%% Simulation of outbreak
% outbreak simulation
parsT.t0 = 0;
parsT.tf = 60;
initial_state = [parsM.Ntot - 1, 0, 1, 0, 0];
control_sol_outbreak = parsM.cB.*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3); % baseline
state_sol_outbreak = state_solver(parsM, parsC, parsT, initial_state, control_sol_outbreak);


%% No intervention SEIRD dynamics
% reset time
parsT.t0 = 60;
parsT.tf = 240; % after 6 months intervention
parsT.tseq = parsT.t0:parsT.dt:parsT.tf; % sequence of times
% set partition (policy update timepoints)
parsT.uniformDelta_t = 30;
parsT.tpartition = parsT.t0:parsT.uniformDelta_t:parsT.tf;
if parsT.tpartition(end) < parsT.tf
    parsT.tpartition = [parsT.tpartition, parsT.tf];
end

% initialize control
control_sol = parsM.cB.*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
initial_state = state_sol_outbreak(end, 2:end);

state_sol_no_intervention = state_solver(parsM, parsC, parsT, initial_state, control_sol);
state_sol_no_intervention_stack = [state_sol_outbreak; state_sol_no_intervention];

%% Projected GD implementation
% terminal condition of costate
terminal_costate = [0, parsC.WEf, parsC.WIf, 0, parsC.Wd];
u_init_1 = (parsC.cmin + 1).*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_init_2 = (parsC.cmax - 1).*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_random = parsC.cmin + (parsC.cmax - parsC.cmin).*rand(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_cB = control_sol;

u_old = u_cB.*1.2; % initialize control sol
% projected GD iterations
num_iter = 100;
J_saver = zeros(num_iter, 1);
kk_saver = zeros(num_iter, 1);

for iter = 1:num_iter
    state_sol_old = state_solver(parsM, parsC, parsT, initial_state, u_old);
    costate_sol_old = costate_solver(parsM, parsC, parsT, terminal_costate, state_sol_old, u_old);
    J_old = cost(parsM, parsC, parsT, state_sol_old, u_old);
    J_saver(iter) = J_old;
    
    gradJ = dJdu(parsM, parsC, state_sol_old, costate_sol_old, u_old);
    proj_u = proj_operator(parsC, parsT, u_old - gradJ); % projection, step size 1, kk = 0
    L2rhs = L2innerProduct(parsT, gradJ, u_old - proj_u); % kk = 0 rhs
   
    state_sol_curr = state_solver(parsM, parsC, parsT, initial_state, proj_u);
    J_curr = cost(parsM, parsC, parsT, state_sol_curr, proj_u);
    
    % steepest descent algorithm with Armijo step size
    kk = 0; J_curr_kk = J_curr;
    while J_curr_kk - J_old >= -parsC.alpha*L2rhs
        kk = kk + 1;
        proj_u_kk = proj_operator(parsC, parsT, u_old - (parsC.beta^kk).*gradJ); % projection
        state_sol_kk = state_solver(parsM, parsC, parsT, initial_state, proj_u_kk);
        J_curr_kk = cost(parsM, parsC, parsT, state_sol_kk, proj_u_kk);
        if kk > 30
            break;
        end
        L2rhs = L2innerProduct(parsT, gradJ, ...
            u_old - proj_operator(parsC, parsT, u_old - parsC.beta^kk.*gradJ));
    end
    
    % update
    if J_curr_kk - J_old < 0 % descent
        u_next = proj_operator(parsC, parsT, u_old - (parsC.beta^(kk)).*gradJ);
        kk_saver(iter) = kk;
    else % exit armijo step size loop with ascent ...
        u_next = u_old; % no update
        kk_saver(iter) = 1000; % set arbitrary high
        break; % break out of big loop
    end
    
    
    u_old = u_next; 
    
    %disp(iter)
end

% stack two periods dynamics (outbreak + intervention)
% dynamics
state_sol_intervention = practical_state_solver(parsM, parsC, parsT, initial_state, u_next);
state_sol_intervention_stack = [state_sol_outbreak; state_sol_intervention];
% controls
contact_rate_no_intervention = [control_sol_outbreak; control_sol];
contact_rate_intervention = [control_sol_outbreak; u_next];

contact_rate_no_intervention_duplicate = [contact_rate_no_intervention(:,1),...
    contact_rate_no_intervention(:,1),contact_rate_no_intervention(:,2),contact_rate_no_intervention(:,3)];
contact_rate_intervention_duplicate = [contact_rate_intervention(:,1),...
    contact_rate_intervention(:,1),contact_rate_intervention(:,2),contact_rate_intervention(:,3)];

% compute person-working fraction
workday_frac = person_workdays_fraction(u_next, state_sol_intervention, parsT, parsM);
% final time total deaths
total_death = state_sol_intervention(end, end);

total_deaths_list(i) = total_death;
work_frac_list(i) = workday_frac;

disp(i);
end

fid = fopen('s5.txt','wt');
fprintf(fid,'%10d %10d %10d\n',[log_weight_ratio_list; total_deaths_list'; work_frac_list']);
fclose(fid);

