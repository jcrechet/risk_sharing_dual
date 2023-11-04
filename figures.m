%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Risk-sharing in a dual market
% Créchet (2020)
% Matlab script file
% file name: "figures.m"
% last updated: Oct 2023
%
% Description: script exporting figures from model's
% quantitative results (Figure 1, Panel a) and b))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% path for figures

host = getenv('COMPUTERNAME');
if strcmp(host,'DESKTOP-LBOSF73')
    fig_path = 'D:\Dropbox\Applications\Overleaf\Risk Sharing in a Dual Labor Market\figures\';
end
if ( strcmp(host, 'DESKTOP-NSR4FR4') || strcmp(host, 'PC-SALON') || strcmp(host, 'DESKTOP-OV023AQ') )
    fig_path = ['C:\Users\',getenv('username'),'\Dropbox\Applications\Overleaf\Risk Sharing in a Dual Labor Market\figures\'];
end
clearvars host

%% Required variables

% workspace
load('workspaces\France.mat')

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
xp = eql.xp;
xt = eql.xt;
xhat = eql.xhat;
xtilde = eql.xtilde;
w_new = eql.w_new;


%% Panel a) surplus perm vs. temp

% color
colormap winter

% variables
yy = (1-bbeta)*D*100;
[ww, xxq] = ndgrid(w, xq);

% plot
s = mesh(log(ww/wr), xxq, yy); 

% labels
zlabel('$\Delta (1-\beta) \mathcal{S}(\nu_0,x) \times 100$', 'Interpreter', 'latex')
xlabel('$\log(w_0(\nu_0) / \underline{w}_R)$', 'Interpreter', 'latex', 'Rotation', 40, 'Position', [-0.5 -0.55])
ylabel('$G_x(x)$', 'Interpreter', 'latex')

% axis
xline(log(w_new/wr), 'LineWidth', 2, 'LineStyle', '--', 'Label', '$\log(w_0(V_0))$', 'Interpreter', 'latex', 'LabelOrientation', 'aligned')
xlim([0 1])
ylim([0 1])
zlim([-3 3])
zticks([-1.5 0 1.5])
yticks([0 0.25 0.5 0.75 1])
xticks([0 0.5 1])

% additional options
view(-65, 30);
grid("minor")
s.FaceColor = "none";

% export
exportgraphics(gcf, [fig_path, 'panel_a.eps'])


%% Panel b) equilibrium distribution

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
xlabel('Log match quality $\log(x)$', 'Interpreter', 'latex')
ylim([0 1.80])
xlim([-0.75 1.35])
grid("on")

% legend
legend('Perm.\ jobs', 'Temp.\ jobs', 'All jobs', 'Sampling', 'Location', 'northeast', 'Interpreter', 'latex');

% export
exportgraphics(gcf, [fig_path, 'panel_b.eps'])
disp('Figures exported')









