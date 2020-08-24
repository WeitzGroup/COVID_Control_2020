function cS_heuristic = bang_bang_control(line_params,curr_state,parsC,parsM)
    % side of line with origin is go to work, other side is stay at home
    m = line_params(1);
    c = line_params(2);
    % line eqn is y - mx - c = 0
    I = curr_state(1);
    R = curr_state(2);
    if ((-c)*(R - m*I - c) > 0)     % point on same side of origin
        cS_heuristic = parsM.cB;    % control is go to work
    else
        cS_heuristic = parsC.cmin;  % lockdown
    end    
end

