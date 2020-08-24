function c = cost(parsM, parsC, parsT, state_sol, control_sol)
[numT, numC] = size(state_sol);
c = 0;

for i = 1:numT
    control_sol_duplicate = [control_sol(i, 1) control_sol(i, 1) control_sol(i, 2) control_sol(i, 3)];
    Q = control_sol_duplicate*(state_sol(i, 2:end - 1))'; % inner product
    E1 = parsC.W1*exp(parsC.K*(parsM.cB*parsM.Ntot - Q));
    E2 = 0.5*((parsM.cB - control_sol_duplicate)*parsC.W2*(parsM.cB - control_sol_duplicate)');
    c = c + (E1 + E2 + parsC.WI*state_sol(i, 4))*parsT.dt;
end

c = c + parsC.WEf*state_sol(end, 3) + parsC.WIf*state_sol(end, 4) + ...
    parsC.Wd*state_sol(end, end);

end