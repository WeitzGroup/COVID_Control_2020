function grad_dLdx = dLdx(parsM, parsC, state_vec, control_vec)
% dL/du - 5 by 1 vector
% syms cS cE cI cR S E I R cB K Ntot nu
% jacobian(1/(1 + exp(-K*(cB*Ntot*(1 - nu) - (cS*S + cE*E + cI*I + cR*R)))), [S, E, I, R])

S = state_vec(1);
E = state_vec(2);
I = state_vec(3);
R = state_vec(4);

cS = control_vec(1);
cI = control_vec(2);
cR = control_vec(3);

Q = E*cS + I*cI + R*cR + S*cS;

grad_dLdx = ...
[ -parsC.K*cS*exp(- parsC.K*(Q - parsM.cB*parsM.Ntot));
  -parsC.K*cS*exp(- parsC.K*(Q - parsM.cB*parsM.Ntot));
  -parsC.K*cI*exp(- parsC.K*(Q - parsM.cB*parsM.Ntot));
  -parsC.K*cR*exp(- parsC.K*(Q - parsM.cB*parsM.Ntot));
    0];

grad_dLdx = grad_dLdx.*parsC.W1 + [0; 0; 1; 0; 0].*parsC.WI; % second part for cost of I


end