green=[0.4660 0.6740 0.1880];
red=[0.6350 0.0780 0.1840];
blue=[0 0.4470 0.7410];
grey=[0.7 0.7 0.7];

% simtest10_3Dec.mat - error-vs-DO
% simtest11_4Dec.mat - error-vs-T

%% Gamma_e vs. Gamma_w for Non-stochastic


figure(1);hold on;
set(gca,'Box','on');
x=1:1:100;
y=((sims(2,2)-sims(1,2))*x);
y=y+(sims(1,2)-y(1,1));
% p2=plot(x,y,'LineWidth',1.5,'Color','k','LineStyle',':');

p1=plot(sims(:,1),sims(:,2).^2,'LineWidth',2,'Color',blue);
axis tight;
% ylim([0 70]);
% xlim([0 100]);
% legend([p1,p2],{'Non-stochastic exploration','Tangent at $\gamma_\mathrm{w}=1$'},'Interpreter','latex','FontSize',14,'Location','northwest');
xlabel('Disturbance energy-bound $\gamma_\mathrm{w}$','Interpreter','latex','FontSize',14);
ylabel('Input energy $\gamma_\mathrm{e}^2$','Interpreter','latex','FontSize',14);

%%

figure(1); hold on;
set(gca,'Box','on');
set(gca, 'YScale', 'log')
p1=plot(sims(1:5,4),sims(1:5,5),'LineWidth',2,'Color',blue);
p2=plot(sims(1:5,4),sims(1:5,6),'LineWidth',2,'Color',blue,'Linestyle',':');
p3=plot(sims(1:5,4),sims(1:5,7),'LineWidth',2,'Color',red);
p4=plot(sims(1:5,4),sims(1:5,8),'LineWidth',2,'Color',red,'Linestyle',':');
axis tight;
legend([p1,p2,p3,p4],{'$\|\theta_\mathrm{tr}-\hat{\theta}_{T,ns}\|$',...
    '$\|\theta_\mathrm{tr}-\hat{\theta}_{T,ns}\|_{D_\mathrm{des}}$',...
    '$\|\theta_\mathrm{tr}-\hat{\theta}_{T,s}\|$',...
    '$\|\theta_\mathrm{tr}-\hat{\theta}_{T,s}\|_{D_\mathrm{des}}$'},...
    'Interpreter','latex','FontSize',14,'Location','northwest');
xlabel('T (=L)','Interpreter','latex','FontSize',14);
ylabel('Error','Interpreter','latex','FontSize',14);


%% DO vs error

% simpars=[simpars,[6;5;4;3;2;1]];
temp1=["D0(1,1)","D_des_tilde(1,1)","gamma_w","T","err_ns","werr_ns","G_ns","normGP","nDdesinv"];

figure(1); hold on;
set(gca,'Box','on');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
p1=plot(simpars(3:6,1),simpars(3:6,9),'LineWidth',2,'Color',red);
p2=plot(simpars(3:6,1),simpars(3:6,8),'LineWidth',2,'Color',blue);
p3=plot(simpars(3:6,1),simpars(3:6,5).^2,'LineWidth',2,'Color',green);
axis tight;
legend([p1,p2,p3],{'Desired: $\|D_\mathrm{des}^{-1}\|$','Guaranteed: $\|GP\|$','Achieved: $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$'},...
    'Interpreter','latex','FontSize',12);
% ylabel('Squared error $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$','Interpreter','latex','FontSize',12);
ylabel('Squared error','Interpreter','latex','FontSize',12);
% xticks([1 2 3 4 5 6])
% xticklabels({'2.5\times 10^{2}','5\times 10^{2}','7.5\times 10^{2}','1\times 10^{3}','2.5\times 10^{3}','5\times 10^{3}'})
xlabel('$\|D_0\|$','Interpreter','latex','FontSize',12);

%% DO vs error - version 2

% simpars=[simpars,[6;5;4;3;2;1]];
temp1=["D0(1,1)","D_des_tilde(1,1)","gamma_w","T","err_ns","werr_ns","G_ns","normGP","nDdesinv"];

figure(1); hold on;
set(gca,'Box','on');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
p1=plot(simpars(3:6,1),simpars(3:6,9),'LineWidth',2,'Color',red);
p2=plot(simpars(3:6,1),simpars(3:6,8),'LineWidth',1.5,'Color',blue);
plot(simpars(3:6,1),simpars(3:6,8),'o','LineWidth',1,'Color',blue);
p3=plot(simpars(3:6,1),simpars(3:6,5).^2,'LineWidth',1.5,'Color',green);
plot(simpars(3:6,1),simpars(3:6,5).^2,'o','LineWidth',1,'Color',green);
% axis tight;
legend([p1,p2,p3],{'Desired: $\|D_\mathrm{des}^{-1}\|$','Guaranteed: $\|GP\|$','Achieved: $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$'},...
    'Interpreter','latex','FontSize',12);
% ylabel('Squared error $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$','Interpreter','latex','FontSize',12);
ylabel('Squared error','Interpreter','latex','FontSize',12);
% xticks([1 2 3 4 5 6])
% xticklabels({'2.5\times 10^{2}','5\times 10^{2}','7.5\times 10^{2}','1\times 10^{3}','2.5\times 10^{3}','5\times 10^{3}'})
xlabel('$\|D_0\|$','Interpreter','latex','FontSize',12);

%% T vs error
temp1=["D0(1,1)","D_des_tilde(1,1)","gamma_w","T","err_ns","werr_ns","G_ns","normGP","nDdesinv"];

figure(1); hold on;
set(gca,'Box','on');
set(gca, 'YScale', 'log');
p1=plot(simpars(:,4),simpars(:,9),'LineWidth',2,'Color',red);
p2=plot(simpars(:,4),simpars(:,8),'LineWidth',1.5,'Color',blue);
plot(simpars(:,4),simpars(:,8),'o','LineWidth',1,'Color',blue);
p3=plot(simpars(:,4),simpars(:,5).^2,'LineWidth',1.5,'Color',green);
plot(simpars(:,4),simpars(:,5).^2,'o','LineWidth',1,'Color',green);
axis tight;
legend([p1,p2,p3],{'Desired: $\|D_\mathrm{des}^{-1}\|$','Guaranteed: $\|GP\|$','Achieved: $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$'},...
    'Interpreter','latex','FontSize',12);
ylabel('Squared error $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$','Interpreter','latex','FontSize',12);
xlabel('$T$','Interpreter','latex','FontSize',12);


%% gamma_w vs error
temp1=["D0(1,1)","D_des_tilde(1,1)","gamma_w","T","err_ns","werr_ns","G_ns","normGP","nDdesinv"];

figure(1); hold on;
set(gca,'Box','on');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
p1=plot(simpars(:,3),simpars(:,9),'LineWidth',2,'Color',red);
p2=plot(simpars(:,3),simpars(:,8),'LineWidth',1.5,'Color',blue);
plot(simpars(:,3),simpars(:,8),'o','LineWidth',1,'Color',blue);
p3=plot(simpars(:,3),simpars(:,5).^2,'LineWidth',1.5,'Color',green);
plot(simpars(:,3),simpars(:,5).^2,'o','LineWidth',1,'Color',green);
axis tight;
legend([p1,p2,p3],{'Desired: $\|D_\mathrm{des}^{-1}\|$','Guaranteed: $\|GP\|$','Achieved: $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$'},...
    'Interpreter','latex','FontSize',12);
ylabel('Squared error $\|\theta_\mathrm{tr}-\hat{\theta}_T\|^2$','Interpreter','latex','FontSize',12);
xlabel('$T$','Interpreter','latex','FontSize',12);

%% 3D plot Gamma_w vs D0 vs error

ref1y=1e5*ones(6,1);
ref1x=simpars(1:6,3);
ref1z=1e-4*ones(6,1);

ref2y=[1e2;1e3;1e4;1e5];
ref2x=10*ones(4,1);
ref2z=1e-4*ones(4,1);
% a=simpars(1:6,3);b=simpars(1:6,1);c=simpars(1:6,8);
p1=plot3(simpars(1:6,3),simpars(1:6,1),simpars(1:6,8),'LineWidth',1.5,'Color',blue); hold on;
plot3(simpars(1:6,3),simpars(1:6,1),simpars(1:6,8),'o','LineWidth',1,'Color',blue);
% plot3(simpars(1:6,3),simpars(1:6,1),simpars(1:6,5).^2,'LineWidth',1.5,'Color',green);
plot3(simpars(7:12,3),simpars(7:12,1),simpars(7:12,8),'LineWidth',1.5,'Color',blue);
plot3(simpars(7:12,3),simpars(7:12,1),simpars(7:12,8),'o','LineWidth',1,'Color',blue);
% plot3(simpars(7:12,3),simpars(7:12,1),simpars(7:12,5).^2,'LineWidth',1.5,'Color',green);
plot3(simpars(13:18,3),simpars(13:18,1),simpars(13:18,8),'LineWidth',1.5,'Color',blue);
plot3(simpars(13:18,3),simpars(13:18,1),simpars(13:18,8),'o','LineWidth',1,'Color',blue);
% plot3(simpars(13:18,3),simpars(13:18,1),simpars(13:18,5).^2,'LineWidth',1.5,'Color',green);
plot3(simpars(19:24,3),simpars(19:24,1),simpars(19:24,8),'LineWidth',1.5,'Color',blue);
plot3(simpars(19:24,3),simpars(19:24,1),simpars(19:24,8),'o','LineWidth',1,'Color',blue);
% plot3(simpars(19:24,3),simpars(19:24,1),simpars(19:24,5).^2,'LineWidth',1.5,'Color',green);
p2=plot3(ref1x,ref1y,ref1z,'LineWidth',2,'Color',red);
plot3(ref2x,ref2y,ref2z,'LineWidth',2,'Color',red);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
set(gca, 'ZScale', 'log');
set(gca,'Box','on');
axis tight;
grid on;
legend([p2,p1],{'Desired: $\|D_\mathrm{des}^{-1}\|$','Guaranteed: $\|GP\|$'},...
    'Interpreter','latex','FontSize',20);
xlabel('$\gamma_\mathrm{w}$','Interpreter','latex','FontSize',20);
ylabel('$\|D_0\|$','Interpreter','latex','FontSize',20);
zlabel('Error bound','Interpreter','latex','FontSize',20);

%%
ref1y=min(simspars(:,3))*ones(6,1);
ref1x=simspars(:,1);
ref1z=1e-4*ones(6,1);

ref2y=simspars(:,3);
ref2x=max(simspars(:,1))*ones(6,1);
ref2z=1e-4*ones(6,1);

p1=plot3(simspars(:,3),simspars(:,1),simspars(:,4)); hold on;
plot3(simspars(:,3),simspars(:,1),simspars(:,4),'o');

% p2=plot3(ref1x,ref1y,ref1z,'LineWidth',2,'Color',red);
% plot3(ref2x,ref2y,ref2z,'LineWidth',2,'Color',red);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'ZScale', 'log');
set(gca,'Box','on');
grid on;
axis tight;
xlabel('$\gamma_\mathrm{e}$','Interpreter','latex','FontSize',20);
ylabel('$\|D_0\|$','Interpreter','latex','FontSize',20);
zlabel('Error bound','Interpreter','latex','FontSize',20);