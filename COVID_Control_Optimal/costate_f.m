function f = costate_f(parsM, parsC, costate_vec, state_vec, control_vec)

Dxf = dfdx(parsM, parsC, state_vec, control_vec);
DxL = dLdx(parsM, parsC, state_vec, control_vec);

f = - costate_vec*Dxf - DxL';

end