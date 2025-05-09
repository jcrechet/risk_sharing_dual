%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Risk-sharing in a dual market
% Créchet (2020)
% Matlab script file
% file name: "figures.m"
% last updated: Oct 2023
%
% Description: script exporting figures from model's
% quantitative results (Figure 1, Panel a), b), )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load param. and equilibrium variables required for plots

% path to result folder
fig_path = '_results\figures\';

% workspace
load('workspaces\France.mat') % focus on French calibration

% parameters
bbeta = p.pval(p.ind.beta);
ssigma_x = p.pval(p.ind.sigma_x);
mu_x = - ssigma_x^2/2;
x_ub = mu_x + 2.5*ssigma_x;

% grids
x = eql.xgrid;
z = eql.zgrid;
Iz = length(z);

% distribution
G = @(x) logncdf(x, mu_x, ssigma_x);
xq = G(x);

% policy/value functions
D = eql.D;
w = eql.wgrid;
wr = eql.wr;
xp = eql.xp;
xt = eql.xt;
xhat = eql.xhat;
xtilde = eql.xtilde;
w_new = eql.w_new;

% option for graph prosition
pos = [680   678   500   350];


%% Panel a) surplus perm vs. temp (two dimensions, v0 and x)

% color
colormap winter

% variables
yy = (1-bbeta)*D*100;
[ww, xxq] = ndgrid(w, xq);

% plot
s = mesh(log(ww/wr), xxq, yy); 

% labels
zlabel('$\Delta (1-\beta) \mathcal{S}(\nu_0,x) \times 100$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\log(w_0(\nu_0) / \underline{w}_R)$', 'Interpreter', 'latex', 'Rotation', 40, 'Position', [-0.5 -0.55], 'FontSize', 14)
ylabel('$G_x(x)$', 'Interpreter', 'latex', 'FontSize', 14)

% axis
xlim([0 1])
ylim([0 1])
zlim([-3 3])
zticks([-1.5 0 1.5])
yticks([0 0.25 0.5 0.75 1])
xticks([0 0.5 1])
grid("minor")
box on

% additional options
view(-65, 30);
s.FaceColor = "none";

% position
set(gcf, 'position', pos)

% export
exportgraphics(gcf, [fig_path, 'panel_a.eps'])


%% Panel c) cross-section 1: surplus vs. x (i.e., fixing v0)

% fix v0
% (i): equilibrium value
[~, iv1] = min( abs(w_new - w) );

% dep variable
y1 = (1-bbeta)*D(iv1,:)*100;
% y2 = (1-bbeta)*D(iv2,:)*100;

% plot
plot(xq, y1, 'LineStyle', '-', 'LineWidth', 2, 'Color', '#4169E1');

% labels
ylabel('$\Delta (1-\beta) \mathcal{S}(\nu_0,x) \times 100$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Match-quality quantile, $G_x(x)$', 'Interpreter', 'latex', 'FontSize', 14)

% zero line
yline(0, 'LineWidth', 1, 'LineStyle', '--')

% cutoffs
xline(G(xt), 'LineWidth', 1, 'LineStyle', '--', 'Label', '$G_x(x_T)$', 'Interpreter', 'latex', 'LabelOrientation', 'aligned')
xline(G(xp), 'LineWidth', 1, 'LineStyle', '--', 'Label', '$G_x(x_P)$', 'Interpreter', 'latex', 'LabelOrientation', 'aligned')
xline(G(xhat), 'LineWidth', 1, 'LineStyle', '--', 'Label', '$G_x(\hat{x})$', 'Interpreter', 'latex', 'LabelOrientation', 'aligned')

% axis
xlim([0 1])
xticks([0 0.25 0.5 0.75 1])
ylim([-0.5 0.5])
yticks([-0.4 -0.2 0 0.2 0.4])
grid('on')

% position
set(gcf, 'position', pos)

% export
exportgraphics(gcf, [fig_path, 'panel_c.eps'])


%% Panel d) cross-section 2: surplus vs. v0 (i.e., fixing x0)

% fix x

% (i): midpoint between 
% [~,ix] = min( xhat - x);
ix = 36;

% dep variable
y = (1-bbeta)*D(:,ix)*100;

% plot
plot(log(w/wr), y, 'LineStyle', '-', 'LineWidth', 2, 'Color', '#4169E1');

% labels
% ylabel('              ', 'Interpreter', 'latex')
xlabel('Wage (rel.\ to reserv.\ wage), $\log(w_0(\nu_0) / \underline{w}_R)$', 'Interpreter', 'latex', 'FontSize', 12)

% zero line
yline(0, 'LineWidth', 1, 'LineStyle', '--')

% wage line
xline(log(w_new/wr), 'LineWidth', 1, 'LineStyle', '--', 'Label', '$\log(w_0(V_0))$', 'Interpreter', 'latex', 'LabelOrientation', 'aligned')

% axis
xlim([0 1])
xticks([0 0.25 0.5 0.75 1])
ylim([-0.5 0.5])
yticks([-0.4 -0.2 0 0.2 0.4])
grid('on')

% position
set(gcf, 'position', pos - [50 50 30 35])

% export
exportgraphics(gcf, [fig_path, 'panel_d.eps'])


%% Panel b) equilibrium distribution

clear gcf

% badwidth for kernel estimations
bw = 0.10;

% unconditional distribution
xsim = sim.xsim;
x = linspace(0, max(xsim), 500);
xq = G(x);
h = ksdensity(xsim, x, 'Kernel', 'normal', 'Support', [xt inf], 'BoundaryCorrection','reflection', 'Bandwidth', bw);

% permanent jobs
xsim = sim.xsim(sim.perm==1);
hp = ksdensity(xsim, x, 'Kernel', 'normal', 'Support', [xp inf], 'BoundaryCorrection','reflection', 'Bandwidth', bw);

% temporary
xsim = sim.xsim(sim.perm==0);
ht = ksdensity(xsim, x, 'Kernel', 'normal', 'Support', [xt inf], 'BoundaryCorrection','reflection', 'Bandwidth', bw);

% sampling
g = @(x) lognpdf(x, mu_x, ssigma_x);

% plot distributions
hold off
hold on
plot(log(x), hp, 'Color', '#87CEEB', 'LineWidth', 1.5)
plot(log(x), ht, 'Color', '#4169E1', 'LineWidth', 1.5)
plot(log(x), h, 'Color', '#808080', 'LineWidth', 1, 'LineStyle', '-.')
plot(log(x), g(x), 'Color', '#808080', 'LineWidth', 1, 'LineStyle', ':')
hold off

% vertical lines for hiring thresholds
xline(log(eql.xt), 'LineWidth', 1, 'LineStyle', '--', 'Label', '$\log(\underline{x}_T)$', 'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left')
xline(log(eql.xp), 'LineWidth', 1, 'LineStyle', '--', 'Label', '$\log(\underline{x}_P)$', 'Interpreter', 'latex')
xline(log(eql.xhat), 'LineWidth', 1, 'LineStyle', '--', 'Label', '$\log(\hat{x})$', 'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left')

% axis
yticks([0 0.5 1 1.5])
xticks([-1 -0.5 0 0.5 1])
xlabel('Log match quality, $\log(x)$', 'Interpreter', 'latex', 'FontSize', 14)
ylim([0 1.80])
xlim([-0.75 1.35])
grid("on")
box on

% legend
legend('Perm.\ jobs', 'Temp.\ jobs', 'All jobs', 'Sampling', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);

% position
set(gcf, 'position',  pos)

% export
exportgraphics(gcf, [fig_path, 'panel_b.eps'])

disp('Figures exported')

