function normdJdu_v = normdJdu(parsT, dJdu_sol)
[numT, numC] = size(dJdu_sol);

normdJdu_v = 0;

for i = 1:numT
    normdJdu_v = normdJdu_v + norm(dJdu_sol(i, :))^2;
end

normdJdu_v = normdJdu_v*parsT.dt; % scalar

end