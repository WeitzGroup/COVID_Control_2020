function output = L2innerProduct(parsT, dJdu_sol, h)
[numT, numC] = size(dJdu_sol);

output = 0;

for i = 1:numT
    output = output + dJdu_sol(i,:)*h(i,:)'; 
end

output = output*parsT.dt; % scalar

end