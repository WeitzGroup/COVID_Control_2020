function [total_man_days_fraction, total_man_days]= ...
    compute_man_days_fraction(contact_rates, state_sol, parsT)
%COMPUTE_MAN_DAYS 
contact_duplicate = [contact_rates(:,1),contact_rates(:,1),...
    contact_rates(:,2),contact_rates(:,3)];
total_man_days_possible = 0;
total_man_days = 0;

for j = 1:4          
    for i = 1:length(contact_duplicate)
        total_man_days_possible = total_man_days_possible + ...
            state_sol(i,j+1)*parsT.dt;
        if contact_duplicate(i,j)>=5
            total_man_days = total_man_days ...
                + state_sol(i,j+1)*parsT.dt;
        else
            total_man_days = total_man_days ...
                + state_sol(i,j+1)*contact_duplicate(i,j)/5*parsT.dt;
        end    
    end
end
total_man_days_fraction = total_man_days/total_man_days_possible;
end
