function state_solution = practical_state_solver(parsM, parsC, parsT, initial_state, control_sol)
state_solution = zeros(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 6); % [t, S, E, I, R, D]
state_solution(1,:) = [parsT.t0, initial_state]; 

for i = 2:round((parsT.tf - parsT.t0)/parsT.dt) + 1
    ss = state_solution(i - 1,2:end);
    ss(ss < 0.8) = 0; % demographic reset
    state_solution(i,2:end) = state_solution(i - 1,2:end) + ...
        parsT.dt.*(state_f(parsM, parsC, ss, control_sol(i - 1, :)))';
    state_solution(i,1) = parsT.t0 + (i - 1)*parsT.dt;
end

end