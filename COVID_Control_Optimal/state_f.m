function f = state_f(parsM, parsC, state_vec, control_vec)
% vector field f(x) - 5 by 1 vector

% state vector
S = state_vec(1);
E = state_vec(2);
I = state_vec(3);
R = state_vec(4);
D = state_vec(5);

cS = control_vec(1);
cI = control_vec(2);
cR = control_vec(3);

Q = E*cS + I*cI + R*cR + S*cS;
F = parsM.etaI*cS*cI*I/Q;
f = [-F*S; F*S - E/parsM.Tinc; E/parsM.Tinc - I/parsM.Tinf;...
    (1 - parsM.mu)*I/parsM.Tinf; parsM.mu*I/parsM.Tinf];

end