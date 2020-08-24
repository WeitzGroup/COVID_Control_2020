close all 
clear all
idx = 1;

isolation_efficiency = 0.25:0.05:0.75;
%isolation_efficiency = 0.25:0.25:0.75;
%shield_level = 1;
shield_level = 2:1:5;
%shield_level = 2:1:2;
load('ga_line_params_constrained_no_intercept_latest.mat')
%load('ga_line_params_constrained_no_intercept_no_shield.mat')
parsM.cB = 5;
% x_f_test = x_f_grid;
% x_f_test(2,1,:) = (x_f_grid(1,1,:)+x_f_grid(3,1,:))/2;
% x_f_test(5,1,:) = (x_f_grid(4,1,:)+x_f_grid(7,1,:))/2;
% x_f_test(6,1,:) = (x_f_grid(4,1,:)+x_f_grid(7,1,:))/2;


set(0,'DefaultAxesTitleFontWeight','normal');

file_idx = [1,4,21,24,41,44];
%file_idx = [1,21,41,4,24,44];
%file_idx = [101,102,103];

for i = 1:5:length(isolation_efficiency)
    c_min = (1-isolation_efficiency(i))*parsM.cB;
    figure;
    subplot(1,3,1)
    %% draw dividing lines and color grey
    m_upper = x_f(i,1,1);
    m_lower = x_f(i,4,1);
    c = -0.001;    
    x_val = 0:0.001:1;
    y_val_line_upper = x_val*m_upper+c/1e6;
    y_val_line_lower = x_val*m_lower+c/1e6;
    y_val_top = ones(1,length(x_val)) - y_val_line_upper;
    y_val_mid = y_val_line_upper - y_val_line_lower;
    parsM.cB = 5;
    y_merged = [y_val_line_lower' y_val_mid' y_val_top'];
    h = area(x_val,y_merged);
    h(3).FaceColor = [1 1 1];
    h(2).FaceColor = [0.75 0.75 0.75];
    h(1).FaceColor = [0.5 0.5 0.5]; hold on;
%     h(1).HandleVisibility = 'off';
%     h(2).HandleVisibility = 'off';
%     h(3).HandleVisibility = 'off';

%% different shielding levels 
    for j = 1:3:length(shield_level)
        c_max = shield_level(j)*parsM.cB;
        
        matFileName = sprintf('heuristic_applied_grid_%d.mat', file_idx(idx));
        applied_result(j) = load(matFileName);
        if j == 1 
            R_max = max(applied_result(j).state_sol_heuristic(:,5));
            I_max = max(applied_result(j).state_sol_heuristic(:,4));
        end    
        
        xlim([0 min(1e6,1.2*I_max/1e6)]);
        ylim([0 min(1e6,1.2*R_max/1e6)]);
        plot(applied_result(j).state_sol_heuristic(:,4)/1e6,...
            applied_result(j).state_sol_heuristic(:,5)/1e6, 'LineWidth', 3);
        hold on;
        plot(applied_result(j).state_sol_heuristic(end,4)/1e6,...
            applied_result(j).state_sol_heuristic(end,5)/1e6,  'dg','MarkerSize',10,...
            'MarkerEdgeColor','k','MarkerFaceColor','g', 'HandleVisibility','off');
        idx = idx+1;
    end

            
    h = legend('$$\textbf{Lockdown~(both)}$$','$$\textbf{Lockdown~(200\%)}$$', ...
            '$$\textbf{Open~(both)}$$', '$$\textbf{200\%~Shielding}$$', '$$\textbf{500\%~Shielding}$$');
        set(h,'FontName','Times New Roman','FontSize',15,'Interpreter','latex', 'Location','northeast');
        legend boxoff;
        xlabel('Infected fraction, $$I$$', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex', 'fontweight','bold'); 
        ylabel('Recovered fraction, $$R$$', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex', 'fontweight','bold');
        %title_char = title(['Isolation - ',num2str(isolation_efficiency(i)*100),'%, Shielding - ',num2str(shield_level(j))]);
        %set(title_char,'FontName','Times New Roman','FontSize',16,'Interpreter','latex'
        %xticks(round(linspace(0, min(1e6,1.1*I_max/1e6), 5),4));
        %yticks(round(linspace(0, min(1e6,1.1*R_max/1e6), 5),2));

        set(gca,'TickLabelInterpreter', 'latex');
        set(gca,'FontSize',15, 'fontweight','bold'); 
        axis square; 
        
        for j = 1:2
            %% find switch time for lockdown to open
            for index = 1:length(applied_result((j-1)*3+1).u_heuristic(:,1))
                if applied_result((j-1)*3+1).u_heuristic(index,1) == 5
                    t_crit = applied_result((j-1)*3+1).state_sol_heuristic(index,1);
                    break;
                end
            end    
            %% plot this
            subplot (1,3,j+1)
            semilogy(applied_result((j-1)*3+1).state_sol_heuristic(:,1), applied_result((j-1)*3+1).state_sol_heuristic(:,4)/1e6, 'b', 'lineWidth', 2); hold on;
            semilogy(applied_result((j-1)*3+1).state_sol_heuristic(:,1), applied_result((j-1)*3+1).state_sol_heuristic(:,5)/1e6, 'm', 'lineWidth', 2); 
            xline(t_crit,'--k')
            axis square;
            %legend('infected, $$I$$','recovered, $$R$$');
            h = legend('infected, $$I$$','recovered, $$R$$');
            set(h,'FontName','Times New Roman','FontSize',12,'Interpreter','latex');
            legend boxoff;
            set(gca,'FontSize',15);
            %yticks([1e-6 1e-4 1e-2 1]); ylim([1e-6 1]); xlim([0 parsT.tf]);xticks(0:60:parsT.tf);
            set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
            xlabel('days', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
            ylabel('population fraction', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
            axis square; %axis tight;
        end
        
end

