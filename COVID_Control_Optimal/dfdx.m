function grad_dfdx = dfdx(parsM, parsC, state_vec, control_vec)
%syms etaI cS cE cI cR S E I R cB
%jacobian(etaI*cS*(cI*I*S/(cS*S + cE*E + cI*I + cR*R)), [S, E, I, R])

S = state_vec(1);
E = state_vec(2);
I = state_vec(3);
R = state_vec(4);

cS = control_vec(1);
cI = control_vec(2);
cR = control_vec(3);

Q = E*cS + I*cI + R*cR + S*cS;

dSdotdS = - ((I*cI*cS*parsM.etaI)/Q - (I*S*cI*cS^2*parsM.etaI)/(Q)^2);
dSdotdE = (I*S*cS*cI*cS*parsM.etaI)/(Q)^2;
dSdotdI = - ((S*cI*cS*parsM.etaI)/Q - (I*S*cI^2*cS*parsM.etaI)/(Q)^2);
dSdotdR = (I*S*cI*cR*cS*parsM.etaI)/(Q)^2;

dSdotdx = [dSdotdS, dSdotdE, dSdotdI, dSdotdR, 0];
dEdotdx = [- dSdotdS, - dSdotdE - 1/parsM.Tinc, -dSdotdI, -dSdotdR, 0];
dIdotdx = [0, 1/parsM.Tinc, - 1/parsM.Tinf, 0, 0];
dRdotdx = [0, 0, (1 - parsM.mu)/parsM.Tinf, 0, 0];
dDdotdx = [0, 0, parsM.mu/parsM.Tinf, 0, 0];

grad_dfdx = [dSdotdx; dEdotdx; dIdotdx; dRdotdx; dDdotdx];
end