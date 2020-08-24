function dJdu_v = dJdu(parsM, parsC, state_sol, costate_sol, control_sol)
% (dJ/du)(t) = lambda^T*Duf + DuL;
[numT, numS] = size(state_sol);
dJdu_v = zeros(numT, parsM.numCVar);
for i = 1:numT
    Duf = dfdu(parsM, parsC, state_sol(i, 2:end), control_sol(i, :));
    DuL = (dLdu(parsM, parsC, state_sol(i, 2:end), control_sol(i, :)))';
    dJdu_v(i, :) = costate_sol(i, 2:end)*Duf + DuL;
end

end