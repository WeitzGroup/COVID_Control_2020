function U_proj = proj_operator(parsC, parsT, U_sol)
% project to hypercube
U_sol(U_sol > parsC.cmax) = parsC.cmax;
U_sol(U_sol < parsC.cmin) = parsC.cmin;
U_proj = U_sol;

if parsC.discrete % piecewise constant strategy
    U_piecewise = zeros(length(U_sol), 3);
    for k = 1:length(parsT.tpartition) - 1
        U_values = ...
            U_sol((parsT.tseq >= parsT.tpartition(k)) & (parsT.tseq < parsT.tpartition(k+1)),:);
        U_avg = mean(U_values, 1); % take mean on this piecewise part
        for col = 1:length(U_avg)
            U_piecewise((parsT.tseq >= parsT.tpartition(k)) &...
                (parsT.tseq < parsT.tpartition(k+1)),col) = U_avg(col);
        end
    end
    U_piecewise(end, :) = U_piecewise(end-1, :); 
    U_proj = U_piecewise;
end

end