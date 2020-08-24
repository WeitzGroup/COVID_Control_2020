function costate_solution = costate_solver(parsM, parsC, parsT, terminal_costate, state_sol, control_sol)
N = round((parsT.tf - parsT.t0)/parsT.dt) + 1;
costate_solution = zeros(N, parsM.numSVar + 1); % [t, lambda1, lambda2, lambda3, lambda4, lambda5]
costate_solution(end,:) = [parsT.tf, terminal_costate]; 

for i = 1:N-1
    f = costate_f(parsM, parsC, costate_solution(N + 1 - i, 2:end), ...
        state_sol(N + 1 - i, 2:end), control_sol(N + 1 - i, :));
    
    costate_solution(N - i, 2:end) = costate_solution(N + 1 - i, 2:end) - parsT.dt.*f;
    
    costate_solution(N - i, 1) = parsT.tf - i*parsT.dt;
end

end