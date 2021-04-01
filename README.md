A1. Figures - Main Text 

Figure 2: 
Comparison of health and economic outcomes of COVID-19 given various interventions: baseline interactions (i.e., no intervention); optimal contact rate intervention (balance both health and economic outcomes) and fully lock down intervention (applied to all the subpopulations) with 75% isolation efficiency. 
Run “main_compare.m”

Figure 3: 
SEIR dynamics with contact rate interventions for various isolation efficiencies, (A) 25% isolation efficiency; (B) 50% isolation efficiency and (C) 75% isolation efficiency.
(A) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.75; parsC.wr = 1; Name = ‘fig_dynamics_w_control_lowisoeff_intermediate’;
(B) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.5; parsC.wr = 1; Name = ‘fig_dynamics_w_control_intermediate’;
(C) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.25; parsC.wr = 1; Name = ‘fig_dynamics_w_control_highisoeff_intermediate’;

Figure 4: 
Heuristic state feedback intervention policies varying with isolation efficiency: (A) 25% isolation efficiency; (B) 50% isolation efficiency and (C) 75% isolation efficiency.
Run plot_ga_reduced_grid.m for all 3 figures. A -> Fig. 1, B -> Fig. 2, C -> Fig. 3
Requires access to mat files to run.

A2. Figures - Supplementary Information

Figure S2: 
Population dynamics of SEIR model without control.
Run ‘main.m’ with setting: wo_control = 1;

Figure S4: 
Population dynamics in a SEIR model with controlled contact rate (25% isolation e_fficiency).
(Top) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.75; parsC.wr = 1e-4; Name = ‘fig_dynamics_w_control_lowisoeff_socioeco’;
(Middle) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.75; parsC.wr = 1; Name = ‘fig_dynamics_w_control_lowisoeff_intermediate’;
(Bottom) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.75; parsC.wr = 1e4; Name = ‘fig_dynamics_w_control_lowisoeff_infect’;

Figure S5: 
Population dynamics in a SEIR model with controlled contact rate (50% isolation e_fficiency).
(Top) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.5; parsC.wr = 1e-4; Name = ‘fig_dynamics_w_control_socioeco’;
(Middle) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.5; parsC.wr = 1; Name = ‘fig_dynamics_w_control_intermediate’;
(Bottom) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.5; parsC.wr = 1e4; Name = ‘fig_dynamics_w_control_infect’;

Figure S6: 
Population dynamics in a SEIR model with controlled contact rate (75% isolation e_fficiency).
(Top) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.25; parsC.wr = 1e-4; Name = ‘fig_dynamics_w_control_highisoeff_socioeco’;
(Middle) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.25; parsC.wr = 1; Name = ‘fig_dynamics_w_control_highisoeff_intermediate’;
(Bottom) Run “main.m” with setting: 
parsC.cmin = parsM.cB*0.25; parsC.wr = 1e4; Name = ‘fig_dynamics_w_control_highisoeff_infect’;

Figure S7: 
Population dynamics in a SEIR model with mis-timed control policy for various isolation e_fficiencies: (Top) 25%
isolation e_fficiency; (Middle) 50% isolation e_fficiency and (Bottom) 75% isolation e_fficiency.
(Top) Run “main_mistimed.m” with setting: 
parsC.cmin = parsM.cB*0.75; Name = ‘fig_dynamics_Tmiss25’;
(Middle) Run “main_mistimed.m” with setting: 
parsC.cmin = parsM.cB*0.5; parsC.wr = 1; Name = ‘fig_dynamics_Tmiss50’;
(Bottom) Run “main_mistimed.m” with setting: 
parsC.cmin = parsM.cB*0.25; parsC.wr = 1e4; Name = ‘fig_dynamics_Tmiss75’;


Figure S9: 
Heuristic state feedback intervention policies varying with isolation efficiency and shielding level on I-R phase plane. 11 figures with 5% increasing isolation efficiency in each row (from 25% to 75%). 4 columns with shielding levels of 2,3,4 and 5 resp.
Run plot_ga_full_grid.m to get 11 figures which are concatenated together to complete Fig. S9. Requires access to mat files to run.


B1. Description of key scripts - Optimal Control 

1. 	dLdx.m: implementation of dL/dx 
2. 	dLdu.m: implementation of dL/du 
3. 	dJdu.m: implementation of dJ/du
4. 	dfdx.m: implementation of df/dx
5.	dfdu.m: implementation of df/du
6. 	normdJdu: implementation of |dJ/du|
6. 	costate_f.m: costate system
7. 	costate_solver.m: numerical integration of costate system based on Euler method (backward direction)
8. 	state_f.m: state dynamical system
9. 	state_solver.m: numerical integration of state system based on Euler method (forward direction)
10. 	cost.m: cost functional
11. 	proj_operator: piecewise constant function approximation of continuous control input
12: 	L2innerProduct: inner product (in L2 space) of two signals 


B2. Description of key scripts - Feedback Control 

1.	bang_bang_control.m: Takes input of system state and provides the contact rate for susceptibles based on line params.
2.	compute_man_days_fraction.m: Based on contact rates and states, computes person-hours worked during lockdown
3.	ga_heuristic.m: Generates mat file with optimal line params for I-R plane to maximize cost function using inbuilt genetic algorithm function
4.	heuristic_applied_grid.m: Computes state trajectories based on optimal line parameters 
5.	heuristic_cost_new: Cost function to be maximized using GA
6.	plot_ga_full_grid: outputs 11 figures which make Fig. S9
7.	plot_ga_reduced_grid: outputs 3 figures which make Fig. 4
