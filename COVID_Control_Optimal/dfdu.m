function grad_dfdu = dfdu(parsM, parsC, state_vec, control_vec)
%syms etaI cS cE cI cR S E I R cB 
%jacobian(etaI*cS*(cI*I*S/(cS*S + cE*E + cI*I + cR*R)), [cS, cE, cI, cR])

S = state_vec(1);
E = state_vec(2);
I = state_vec(3);
R = state_vec(4);

cS = control_vec(1);
cI = control_vec(2);
cR = control_vec(3);

Q = E*cS + I*cI + R*cR + S*cS;
dSdotdcS = - ((I*S*cI*parsM.etaI)/Q + (I*S*cI*cS*parsM.etaI*(E+S))/(Q)^2);
dSdotdcI = - ((I*S*cS*parsM.etaI)/Q + (I^2*S*cI*cS*parsM.etaI)/(Q)^2);
dSdotdcR = (I*R*S*cI*cS*parsM.etaI)/(Q)^2;

grad_dfdu = zeros(5, 3);
grad_dfdu(1, 1) = dSdotdcS;
grad_dfdu(1, 2) = dSdotdcI;
grad_dfdu(1, 3) = dSdotdcR;

grad_dfdu(2, :) =  - grad_dfdu(1, :);

% [ (I*S*cI*etaI)/(E*cS + I*cI + R*cR + S*cS) - (I*S*cI*cS*etaI*(E + S))/(E*cS + I*cI + R*cR + S*cS)^2, (I*S*cS*etaI)/(E*cS + I*cI + R*cR + S*cS) - (I^2*S*cI*cS*etaI)/(E*cS + I*cI + R*cR + S*cS)^2, -(I*R*S*cI*cS*etaI)/(E*cS + I*cI + R*cR + S*cS)^2]
end