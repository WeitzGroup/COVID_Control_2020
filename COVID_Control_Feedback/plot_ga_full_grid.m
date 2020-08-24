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

%file_idx = [1,4,21,24,41,44];
file_idx = 1:44;

for i = 1:length(isolation_efficiency)
    figure;
    c_min = (1-isolation_efficiency(i))*parsM.cB;
    
    for j = 1:length(shield_level)
        c_max = shield_level(j)*parsM.cB;
        subplot(1,4,j)
        m = x_f(i,j,1);
%        c = x_f(i,j,2);
        c = -0.001;    
%         m_test = x_f_test(i,j,1);
%         c_test = x_f_test(i,j,2);
        %heuristic_applied_grid([m c],c_min, c_max, idx);
%         cost_test(i) = heuristic_applied_grid([m_test c_test],c_min, c_max, idx);
        matFileName = sprintf('heuristic_applied_grid_%d.mat', file_idx(idx));
        applied_result = load(matFileName);
        if j == 1
            R_max = max(applied_result.state_sol_heuristic(:,5));
            I_max = max(applied_result.state_sol_heuristic(:,4));
        end
        x_val = 0:0.001:1;
        y_val_line = x_val*m+c/1e6;
        y_val_top = ones(1,length(x_val)) - y_val_line;
        parsM.cB = 5;
        y_merged = [y_val_line' y_val_top'];
            h = area(x_val,y_merged);
            h(2).FaceColor = [1 1 1]; 
            h(1).FaceColor = [0.7 0.7 0.7]; hold on;
            plot(applied_result.state_sol_heuristic(:,4)/1e6,...
                applied_result.state_sol_heuristic(:,5)/1e6, 'LineWidth', 3, 'Color', 'k');
            hold on;
            plot(applied_result.state_sol_heuristic(end,4)/1e6,...
                applied_result.state_sol_heuristic(end,5)/1e6,  'dg','MarkerSize',10,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g');
            h = legend('$$\textbf{Lockdown}$$', '$$\textbf{Open}$$');
            set(h,'FontName','Times New Roman','FontSize',15,'Interpreter','latex', 'Location','northeast');
            legend boxoff;
            xlabel('Infected fraction, $$I$$', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex', 'fontweight','bold'); 
            ylabel('Recovered fraction, $$R$$', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex', 'fontweight','bold');
            %title_char = title(['Isolation - ',num2str(isolation_efficiency(i)*100),'%, Shielding - ',num2str(shield_level(j))]);
            %set(title_char,'FontName','Times New Roman','FontSize',16,'Interpreter','latex'
            %xticks(round(linspace(0, min(1e6,1.1*I_max/1e6), 5),4));
            %yticks(round(linspace(0, min(1e6,1.1*R_max/1e6), 5),2));
            xlim([0 min(1e6,1.2*I_max/1e6)]);
            ylim([0 min(1e6,1.2*R_max/1e6)]);
            set(gca,'TickLabelInterpreter', 'latex');
            set(gca,'FontSize',15, 'fontweight','bold'); 
            axis square;
        
        idx = idx+1;
    end
end

% for i = 1:11
% Name = sprintf('plot_heuristic_ga_no_intercept_full_grid_%d', i);
% ff = figure(i);
% ff.Units = 'inches';
% Width = 150; Height = 15;
% 
% ff.PaperSize = [Width, Height];
% ff.PaperPosition = [0 0 Width, Height];
% ff.Position = [0 0 Width, Height];
% 
% print(ff, Name, '-depsc2','-r600');
% print(ff, Name, '-dpdf','-r600');
% end