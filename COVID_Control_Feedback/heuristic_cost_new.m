function cost = heuristic_cost_new(line_slope, cmin, cmax)
%% Parameters
% model parameters
parsM.Tinc = 4; % length of incubation period
parsM.Tinf = 6; % duration patient is infectious
parsM.etaI = 0.1; % transmission effectivenesss, note transmission rate beta ~ etaI*contactRate
parsM.mu = 1e-2; % case fatality ratio
parsM.cB = 5; % baseline contact rate
parsM.Ntot = 1e6; % total number of population (neglect death)
parsM.numSVar = 5; % number of state variabes
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
parsC.cmax = cmax; % constrain;
parsC.cmin = cmin; % constrain; if cmin is very low (0.1), decreasing I is enough

% optimization algorithm parameters (backtracking step size)
parsC.alpha = 0.1; % alpha is *very* important 
parsC.beta = 0.5;
% numerical solver parameters;
parsT.dt = 5e-2;

parsC.wr = 1; % weight ratio
parsC.Wd = parsC.wr*parsC.Wd;
parsC.WI = parsC.wr*parsC.WI;
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
parsT.tf = 360; % after 6 months intervention
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

u_heuristic(:,1) = parsC.cmin*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1,1); % initialize control sol
u_heuristic(:,2) = parsC.cmin*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1,1); % initialize control sol
u_heuristic(:,3) = parsC.cmax*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1,1);

cost = 0;

trajectory_count = 1;
sigma_fraction = 0.00;       % 2% of val is sigma       


for traj_count = 1:trajectory_count
    noise_sigma = sigma_fraction*initial_state;
    noise = noise_sigma.*randn(1,length(initial_state));
    initial_state_noisy = initial_state + noise;
    state_sol_heuristic(1,:) = [parsT.t0, initial_state_noisy]; 

for i = 2:round((parsT.tf - parsT.t0)/parsT.dt) + 1
    line_params =[line_slope  -0.001];
    u_heuristic(i-1,1) = bang_bang_control(line_params, ...
        state_sol_heuristic(i-1,4:5),parsC, parsM);
    state_sol_heuristic(i,2:end) = state_sol_heuristic(i - 1,2:end) + ...
        parsT.dt.*(state_f(parsM, parsC, state_sol_heuristic(i - 1,2:end), u_heuristic(i - 1, :)))';
    state_sol_heuristic(i,1) = parsT.t0 + (i - 1)*parsT.dt;
end
%% store data for extreme cases
total_deaths = state_sol_heuristic(end);
[working_days_frac, working_days]= compute_man_days_fraction(u_heuristic,state_sol_heuristic,parsT);

% % signed cost
% if total_deaths > optimal_deaths
%     cost = cost + (total_deaths - optimal_deaths)/optimal_deaths;
% end
% 
% if working_days_frac < optimal_working_frac
%     cost = cost + (optimal_working_frac - working_days_frac)/optimal_working_frac;
% end   

cost = cost+500*total_deaths/1e6 + (1-working_days_frac);
    
end
