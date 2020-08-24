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

% optimization algorithm parameters (backtracking step size)
parsC.alpha = 0.1; % alpha is *very* important 
parsC.beta = 0.5;
% numerical solver parameters;
parsT.dt = 5e-2;

% weight balancing
parsC.wr = 1; % weight ratio

parsC.Wd = parsC.wr*parsC.Wd;
parsC.WI = parsC.wr*parsC.WI;
parsC.WIf = 0;

% vary isolation eff from 0.25, 0.5, 0.75
parsC.cmin = parsM.cB*0.75; % constrain; if cmin is very low (0.1), decreasing I is enough
% plot name
Name = 'fig_dynamics_Tmiss_25';
% Name = 'fig_dynamics_oc_50';
% Name = 'fig_dynamics_oc_75';


% piecewise constant strategy or continouse open loop strategy
parsC.discrete = 1; % 0 - continous, 1 - discrete


%% Simulation of outbreak
% outbreak simulation (30 days)
parsT.t0 = 0;
parsT.tf = 30;
initial_state = [parsM.Ntot - 1, 0, 1, 0, 0];
control_sol_outbreak = parsM.cB.*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3); % baseline
state_sol_outbreak = state_solver(parsM, parsC, parsT, initial_state, control_sol_outbreak);

% outbreak simulation (60 days)
parsT1.t0 = 0;
parsT1.tf = 60;
parsT1.dt = 5e-2;
control_sol_outbreak1 = parsM.cB.*ones(round((parsT1.tf - parsT1.t0)/parsT.dt) + 1, 3); % baseline
state_sol_outbreak1 = state_solver(parsM, parsC, parsT1, initial_state, control_sol_outbreak1);


%% No intervention SEIRD dynamics
% reset time
parsT.t0 = 30;
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

parsT1.t0 = 60; parsT1.tf = 360; 
control_sol1 = parsM.cB.*ones(round((parsT1.tf - parsT1.t0)/parsT1.dt) + 1, 3);
initial_state1 = state_sol_outbreak1(end, 2:end);


state_sol_no_intervention = state_solver(parsM, parsC, parsT, initial_state, control_sol);
state_sol_no_intervention_stack = [state_sol_outbreak; state_sol_no_intervention];


%% Projected GD implementation
% terminal condition of costate
terminal_costate = [0, parsC.WEf, parsC.WIf, 0, parsC.Wd];
u_init_1 = (parsC.cmin + 1).*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_init_2 = (parsC.cmax - 1).*ones(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_random = parsC.cmin + (parsC.cmax - parsC.cmin).*rand(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 3);
u_cB = control_sol;

u_old = 0.8.*u_cB; % initialize control sol
% projected GD iterations
num_iter = 50;
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

u_next1 = u_next(1:length(control_sol1), :);
state_sol_intervention1 = ...
    practical_state_solver(parsM, parsC, parsT1, initial_state1,u_next1);

state_sol_intervention_stack = [state_sol_outbreak; state_sol_intervention];
state_sol_intervention_stack1 = [state_sol_outbreak1; state_sol_intervention1];

% controls
contact_rate_intervention = [control_sol_outbreak; u_next];
contact_rate_intervention1 = [control_sol_outbreak1; u_next1];

contact_rate_intervention_duplicate = [contact_rate_intervention(:,1),...
    contact_rate_intervention(:,1),contact_rate_intervention(:,2),contact_rate_intervention(:,3)];
contact_rate_intervention_duplicate1 = [contact_rate_intervention1(:,1),...
    contact_rate_intervention1(:,1),contact_rate_intervention1(:,2),contact_rate_intervention1(:,3)];

% compute person-working fraction
workday_frac = person_workdays_fraction(u_next, state_sol_intervention, parsT, parsM);
workday_frac1 = person_workdays_fraction(u_next1, state_sol_intervention1, parsT1, parsM);

% final time total deaths
total_death = state_sol_intervention(end, end);
total_death1 = state_sol_intervention1(end, end);


% compute Reff
RR_eff_c = R_eff_metric(contact_rate_intervention_duplicate, ...
    state_sol_intervention_stack, parsM);

RR_eff_c1 = R_eff_metric(contact_rate_intervention_duplicate1, ...
    state_sol_intervention_stack1, parsM);



%% Figures
%% Population Dynamics with/without control (assumed model)
set(0,'DefaultAxesLinewidth',2);
figure(1);
subplot(1,3,1);
semilogy(state_sol_intervention_stack(:,1), state_sol_intervention_stack1(:,2)./parsM.Ntot, 'g', 'lineWidth', 2); hold on;
plot(state_sol_intervention_stack(:,1), state_sol_intervention_stack1(:,3)./parsM.Ntot, 'r', 'lineWidth', 2); hold on;
plot(state_sol_intervention_stack(:,1), state_sol_intervention_stack1(:,4)./parsM.Ntot, 'b', 'lineWidth', 2); hold on;
plot(state_sol_intervention_stack(:,1), state_sol_intervention_stack1(:,5)./parsM.Ntot, 'm', 'lineWidth', 2); hold on;
plot(state_sol_intervention_stack(:,1), state_sol_intervention_stack1(:,6)./parsM.Ntot, 'k', 'lineWidth', 2); hold on;
axis square;
h = legend('susceptible, $$S$$', 'exposed, $$E$$', ...
    'infected, $$I$$',...
    'recovered, $$R$$', 'dead, $$D$$');
set(h,'FontName','Times New Roman','FontSize',12,'Interpreter','latex', 'Location','southeast');
legend boxoff;
% tl = title('dynamics with control'); set(tl,'FontName','Times New Roman','FontSize',15, 'Interpreter','latex');
set(gca,'FontSize',15);
yticks([1e-6 1e-4 1e-2 1]); ylim([1e-6 1]); xlim([0 parsT.tf]); xticks(0:60:parsT.tf);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('days', 'FontName', 'Times New Roman','FontSize',15,...
    'Interpreter','latex', 'fontweight','bold'); 
ylabel('population fraction', 'FontName', 'Times New Roman','FontSize',15,...
    'Interpreter','latex', 'fontweight','bold');
axis square; %axis tight;

subplot(1,3,2);
yyaxis left;
hold('on');
%plot(state_sol_intervention_stack(:,1), contact_rate_intervention(:,1), 'g', 'lineWidth', 3); hold on;
plot1 = plot(state_sol_intervention_stack(:,1), contact_rate_intervention1(:,1), '-g', 'lineWidth', 5);
plot1.Color(4) = 0.4;
hold on;
plot2 = plot(state_sol_intervention_stack(:,1), contact_rate_intervention1(:,2), '-r', 'lineWidth', 2);
hold on;
plot3 = plot(state_sol_intervention_stack(:,1), contact_rate_intervention1(:,3), '-b', 'lineWidth', 2); 
set(gca,'FontSize',15);
yticks(0:2:parsC.cmax); xticks(0:60:parsT.tf);
ylim([0, parsC.cmax]); xlim([0 parsT.tf]);
set(gca,'FontSize',15); set(gca,'TickLabelInterpreter', 'latex');
xlabel('days', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex', 'fontweight','bold'); 
ylabel('contact rate', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex', 'fontweight','bold');

yyaxis right;
plotX = plot(state_sol_intervention_stack(:,1), RR_eff_c1,...
    '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 2); % Reff
ylim([0, 5]); yticks(0:1:5);
plotX.Color(4) = 0.7;
ylabel('$$\mathcal{R}_{eff}$$', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex', 'fontweight','bold');
axis square; %axis tight;
box on;

h = legend('$$c_{S}$$',  ...
    '$$c_{I}$$', '$$c_{R}$$');
set(h,'FontName','Times New Roman','FontSize',12,'Interpreter','latex', 'Location','northeast');
legend boxoff;
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0.3 0.3 0.3];


subplot(1,3,3);
QQ = sum(state_sol_intervention_stack1(:,2:5).*contact_rate_intervention_duplicate1, 2);
opt_connection_loss = parsC.W1.*exp(parsC.K*(parsM.cB*parsM.Ntot - QQ)) + ...
    0.5.*sum((contact_rate_intervention_duplicate1 - parsM.cB).^2*parsC.W2, 2);
plot(state_sol_intervention_stack1(:,1), opt_connection_loss, 'k', 'LineWidth', 2);
% tl = title('socio-economic response'); set(tl,'FontName','Times New Roman','FontSize',15, 'Interpreter','latex');
set(gca,'FontSize',15);
set(gca,'FontSize',16,'fontweight','bold');
%yticks(0:0.2:1); 
xticks(0:60:parsT.tf);
%ylim([0, 1]); 
xlim([0 parsT.tf]);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('days', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
ylabel('socioeconomic costs', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 




% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 30; Height = 10;

ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];

print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');








