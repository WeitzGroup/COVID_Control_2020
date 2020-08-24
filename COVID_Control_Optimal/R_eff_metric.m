function R_eff = R_eff_metric(control_sol, state_sol, parsM)
ll = length(control_sol);
R_eff = zeros(1, ll);
for l = 1:ll
    Q = control_sol(l,:)*state_sol(l, 2:5)';
    R_eff(l) = parsM.Tinf*parsM.etaI*...
        control_sol(l,1)*control_sol(l,3)*state_sol(l,2)/Q;
end

end