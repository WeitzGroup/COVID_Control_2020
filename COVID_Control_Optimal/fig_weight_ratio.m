% load data
D2 = load('s2.txt');
D3 = load('s3.txt');
D4 = load('s4.txt');
D5 = load('s5.txt');

% multiple = 2
figure(1);
subplot(1,2,1);
plot(D2(:,1), D2(:,3),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue');
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylabel('$$\Psi$$, working-days fraction', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 

subplot(1,2,2);
semilogy(D2(:,1), D2(:,2),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','red'); hold on;
semilogy(D2(:,1), 43.*ones(length(D2(:,1)), 1), '--k', 'linewidth',3);
set(gca,'FontSize',15);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylim([1e1 1e4]);yticks([1e1, 1e2, 1e3, 1e4]); 
ylabel('$$D(t_{f})$$, total~deaths', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 

% multiple = 3
figure(2);
subplot(1,2,1);
plot(D3(:,1), D3(:,3),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue');
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylabel('$$\Psi$$, working-days fraction', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 

subplot(1,2,2);
semilogy(D3(:,1), D3(:,2),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','red');hold on;
semilogy(D2(:,1), 43.*ones(length(D2(:,1)), 1), '--k', 'linewidth',3);
set(gca,'FontSize',15);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylim([1e1 1e4]); yticks([1e1, 1e2, 1e3, 1e4]); 
ylabel('$$D(t_{f})$$, total~deaths', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 





% multiple = 4
figure(3);
subplot(1,2,1);
plot(D4(:,1), D4(:,3),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue');
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylabel('$$\Psi$$, working-days fraction', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 

subplot(1,2,2);
semilogy(D4(:,1), D4(:,2),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','red');hold on;
semilogy(D2(:,1), 43.*ones(length(D2(:,1)), 1), '--k', 'linewidth',3);
set(gca,'FontSize',15);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylim([1e1 1e4]); yticks([1e1, 1e2, 1e3, 1e4]); 
ylabel('$$D(t_{f})$$, total~deaths', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 


% multiple = 5
figure(4);
subplot(1,2,1);
plot(D5(:,1), D5(:,3),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue');
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylabel('$$\Psi$$, working-days fraction', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 

subplot(1,2,2);
semilogy(D5(:,1), D5(:,2),'-ok','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','red');hold on;
semilogy(D2(:,1), 43.*ones(length(D2(:,1)), 1), '--k', 'linewidth',3);
set(gca,'FontSize',15);
set(gca,'FontSize',15);set(gca,'TickLabelInterpreter', 'latex');
xlabel('log(weight~ratio)', 'FontName', 'Times New Roman','FontSize',15,'Interpreter','latex'); 
xticks(-4:2:4); xlim([-4 4]);
ylim([1e1 1e4]); yticks([1e1, 1e2, 1e3, 1e4]); 
ylabel('$$D(t_{f})$$, total~deaths', 'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
axis square; 





