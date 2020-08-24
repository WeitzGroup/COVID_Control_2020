function grad_dLdu = dLdu(parsM, parsC, state_vec, control_vec)
% dL/du - 4 by 1 vector
% syms cS cE cI cR S E I R cB K Ntot nu
% jacobian(0.5*((cS - cB)^2 + (cE - cB)^2 + (cI - cB)^2 + (cR - cB)^2) + ...
%    1/(1 + exp(-K*(cB*Ntot*(1 - nu) - (cS*S + cE*E + cI*I + cR*R)))), [cS, cE, cI, cR])

S = state_vec(1);
E = state_vec(2);
I = state_vec(3);
R = state_vec(4);

cS = control_vec(1);
cI = control_vec(2);
cR = control_vec(3);

Q = E*cS + I*cI + R*cR + S*cS;

grad_dLdu = ...
[parsC.Ws*(cS - parsM.cB) - ...
   parsC.W1*parsC.K*S*exp( - parsC.K*(Q - parsM.cB*parsM.Ntot)) + ...
   parsC.We*(cS - parsM.cB) - ...
   parsC.W1*parsC.K*E*exp( - parsC.K*(Q - parsM.cB*parsM.Ntot)); 
 parsC.Wi*(cI - parsM.cB) - ...
   parsC.W1*parsC.K*I*exp( - parsC.K*(Q - parsM.cB*parsM.Ntot));
 parsC.Wr*(cR - parsM.cB) - ...
   parsC.W1*parsC.K*R*exp( - parsC.K*(Q - parsM.cB*parsM.Ntot))];

end

