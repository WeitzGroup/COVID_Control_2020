%% Parameters
% model parameters
parsM.Tinc = 4; % length of incubation period
parsM.Tinf = 6; % duration patient is infectious
parsM.etaI = 0.1; % *true* transmission effectivenesss, note transmission rate beta ~ etaI*contactRate
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

parsC.W2 = (0.1/parsM.cB^2).*eye(4);  % weight matrix
parsC.Ws = parsC.W2(1,1);
parsC.We = parsC.W2(2,2);
parsC.Wi = parsC.W2(3,3);
parsC.Wr = parsC.W2(4,4);


parsC.K = (7/(parsM.cB*parsM.Ntot)); % scale of exponential cost

% contacting rate range (policy efficiency)
parsC.cmax = parsM.cB*2; % constrain;
parsC.cmin = parsM.cB*0.5; % constrain; if cmin is very low (0.1), decreasing I is enough
% weight balancing
parsC.wr = 1; % weight ratio

% optimization algorithm parameters (backtracking step size)
parsC.alpha = 0.1; % alpha is *very* important 
parsC.beta = 0.5;
% numerical solver parameters;
parsT.dt = 5e-2;

parsC.Wd = parsC.wr*parsC.Wd;
parsC.WI = parsC.wr*parsC.WI;
parsC.WIf = 0;


% piecewise constant strategy or continouse open loop strategy
parsC.discrete = 1; % 0 - continous, 1 - discrete


%% Simulation of outbreak
% outbreak simulation with true model
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
lock_sol = (5*0.25).*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
initial_state = state_sol_outbreak(end, 2:end);

state_sol_no_intervention = state_solver(parsM, parsC, parsT, initial_state, control_sol);
state_sol_no_intervention_stack = [state_sol_outbreak; state_sol_no_intervention];

state_sol_lock = state_solver(parsM, parsC, parsT, initial_state, lock_sol);
state_sol_lock_stack = [state_sol_outbreak; state_sol_lock];


%% Projected GD implementation
% terminal condition of costate
terminal_costate = [0, parsC.WEf, parsC.WIf, 0, parsC.Wd];
u_init_1 = (parsC.cmin + 1).*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_init_2 = (parsC.cmax - 1).*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_random = parsC.cmin + (parsC.cmax - parsC.cmin).*rand(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_cB = control_sol;

u_old = 0.9.*u_cB; % initialize control sol
% projected GD iterations
num_iter = 80;
J_saver = zeros(num_iter, 1);
kk_saver = zeros(num_iter, 1);

parsM.etaI = 0.1; % parameter for assumed model and control (R0 = 3)

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
    
    disp(iter)
end

%% stack two periods dynamics (outbreak + intervention)

% dynamics
state_sol_intervention = practical_state_solver(parsM, parsC, parsT, initial_state, u_next);
state_sol_intervention_stack = [state_sol_outbreak; state_sol_intervention];

% controls
contact_rate_no_intervention = [control_sol_outbreak; control_sol];
contact_rate_intervention = [control_sol_outbreak; u_next];
contact_rate_lock = [control_sol_outbreak; lock_sol];

contact_rate_no_intervention_duplicate = [contact_rate_no_intervention(:,1),...
    contact_rate_no_intervention(:,1),contact_rate_no_intervention(:,2),contact_rate_no_intervention(:,3)];
contact_rate_intervention_duplicate = [contact_rate_intervention(:,1),...
    contact_rate_intervention(:,1),contact_rate_intervention(:,2),contact_rate_intervention(:,3)];
contact_rate_lock_duplicate = [contact_rate_lock(:,1),...
    contact_rate_lock(:,1),contact_rate_lock(:,2),contact_rate_lock(:,3)];

% compute person-working fraction
workday_frac = person_workdays_fraction(u_next, state_sol_intervention, parsT, parsM);
% final time total deaths
total_death = state_sol_intervention(end, end);


% compute Reff
RR_eff_c = R_eff_metric(contact_rate_intervention_duplicate, ...
    state_sol_intervention_stack, parsM);

RR_eff = R_eff_metric(contact_rate_no_intervention_duplicate, ...
    state_sol_no_intervention_stack, parsM);

RR_eff_l = R_eff_metric(contact_rate_lock_duplicate, ...
    state_sol_lock_stack, parsM);



%% Figures

set(0,'DefaultAxesLinewidth',2);
QQ = sum(state_sol_intervention_stack(:,2:5).*contact_rate_intervention_duplicate, 2);
opt_connection_loss = parsC.W1.*exp(parsC.K*(parsM.cB*parsM.Ntot - QQ)) + ...
    0.5.*sum((contact_rate_intervention_duplicate - parsM.cB).^2*parsC.W2, 2);

% without control socio-eco loss
QQ_wc = sum(state_sol_no_intervention_stack(:,2:5).*contact_rate_no_intervention_duplicate, 2);
wc_connection_loss = parsC.W1.*exp(parsC.K*(parsM.cB*parsM.Ntot - QQ_wc)) + ...
    0.5.*sum((contact_rate_no_intervention_duplicate - parsM.cB).^2*parsC.W2, 2);

QQ_lock = sum(state_sol_lock_stack(:,2:5).*contact_rate_lock_duplicate, 2);
lock_connection_loss = parsC.W1.*exp(parsC.K*(parsM.cB*parsM.Ntot - QQ_lock)) + ...
    0.5.*sum((contact_rate_lock_duplicate - parsM.cB).^2*parsC.W2, 2);


figure(1);
subplot(1, 3, 1);
plot1 = plot(state_sol_intervention_stack(:,1), (state_sol_no_intervention_stack(:,end)./parsM.Ntot).*1e5,...
    '-','Color', [0.4 0.4 0.4], 'lineWidth', 4); hold on;
plot1.Color(4) = 0.4;
plot2 = plot(state_sol_intervention_stack(:,1), ...
    (state_sol_intervention_stack(:,end)./parsM.Ntot + 1e-8).*1e5, '-m', 'LineWidth', 3); hold on;
plot3 = plot(state_sol_intervention_stack(:,1), ...
    (state_sol_lock_stack(:,end)./parsM.Ntot).*1e5,  'k', 'lineWidth', 2);
axis square;
h = legend('baseline interactions', 'optimal control', ...
    'fully lock down');
set(h,'FontName','Times New Roman','FontSize',12,'Interpreter','latex', 'Location','northeast');
legend boxoff;
set(gca,'FontSize',15);
yticks(0:200:1000); ylim([0 1000]); xlim([0 parsT.tf]); xticks(0:60:parsT.tf);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('days', 'FontName', 'Times New Roman','FontSize',15,...
    'Interpreter','latex', 'fontweight','bold'); 
ylabel('cummulative deaths per 100,000 ', 'FontName', 'Times New Roman','FontSize',15,...
    'Interpreter','latex', 'fontweight','bold');
axis square; %axis tight;


subplot(1, 3, 2);
plot1 = plot(state_sol_intervention_stack(:,1), wc_connection_loss, '-','Color', [0.4 0.4 0.4], 'lineWidth', 4);
plot1.Color(4) = 0.4;
hold on;
plot2 = plot(state_sol_intervention_stack(:,1), opt_connection_loss, '-m', 'LineWidth', 3); hold on;
plot3 = plot(state_sol_intervention_stack(:,1), lock_connection_loss, 'k', 'lineWidth', 2);
set(gca,'FontSize',15);
set(gca,'FontSize',16,'fontweight','bold');
h = legend('baseline interactions', 'optimal control', ...
    'fully lock down');
set(h,'FontName','Times New Roman','FontSize',12,'Interpreter','latex', 'Location','northeast');
legend boxoff;
%yticks(0:0.2:1); 
xticks(0:60:parsT.tf);
%ylim([0, 1]); 
xlim([0 parsT.tf]);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('days', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
ylabel('socioeconomic costs', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 


subplot(1, 3, 3);
plot1 = plot(state_sol_intervention_stack(:,1), RR_eff, '-','Color', [0.4 0.4 0.4], 'lineWidth', 4);
plot1.Color(4) = 0.4;
hold on;
plot2 = plot(state_sol_intervention_stack(:,1), RR_eff_c,  '-m', 'LineWidth', 3); hold on;
plot3 = plot(state_sol_intervention_stack(:,1), RR_eff_l, 'k', 'lineWidth', 2);
ylim([0, 4]); yticks(0:1:4);
xticks(0:60:parsT.tf); xlim([0 parsT.tf]);
xlabel('days', 'FontName', 'Times New Roman','FontSize',15,...
    'Interpreter','latex', 'fontweight','bold'); 
ylabel('$$\mathcal{R}_{eff}$$', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex', 'fontweight','bold');
axis square; %axis tight;
h = legend('baseline interactions', 'optimal control', ...
    'fully lock down');
set(h,'FontName','Times New Roman','FontSize',12,'Interpreter','latex', 'Location','northeast');
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
legend boxoff;


figure(2);
plot1 = plot(state_sol_intervention_stack(:,1), contact_rate_intervention(:,1), '-g', 'lineWidth', 5);
plot1.Color(4) = 0.4; hold on;
plot2 = plot(state_sol_intervention_stack(:,1), contact_rate_intervention(:,2), '-r', 'lineWidth', 3);
hold on;
plot3 = plot(state_sol_intervention_stack(:,1), contact_rate_intervention(:,3), '-b', 'lineWidth', 3); 
set(gca,'FontSize',15);
yticks(0:2:parsC.cmax); xticks(0:60:parsT.tf);
ylim([0, parsC.cmax]); xlim([0 parsT.tf]);
set(gca,'FontSize',15); set(gca,'TickLabelInterpreter', 'latex');
xlabel('days', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex', 'fontweight','bold'); 
ylabel('contact rate', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex', 'fontweight','bold');
h = legend('$$c_{S}$$',  ...
    '$$c_{I}$$', '$$c_{R}$$');
set(h,'FontName','Times New Roman','FontSize',12,'Interpreter','latex', 'Location','northeast');
legend boxoff;
axis square;


% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 30; Height = 11;
Name = 'fig_dynamics_compare';
ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');


% output figures in eps 
ff = figure(2);
ff.Units = 'inches';
Width = 10; Height = 10;
Name = 'fig_dynamics_control';
ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');

