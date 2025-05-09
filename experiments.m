%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project: Risk-sharing in a dual market
% Créchet (2020)
% matlab script file
% file name: "experiment.m"
% created: 12-08-2021
% updated: 10-2023
% Description: quantitative experiments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Relaxing hiring regulation of temp contacts (early 1980-5.5% TC emp. share)

% France
load('workspaces\France.mat')
p.pval(p.ind.phi0) = 0.55;
p.equilibrium = 'general';
[eql, sim, agg_stat] = compute_equilibrium(p);
disp(agg_stat)
save('workspaces\counterfactuals\France_1980.mat', 'eql', 'sim', 'agg_stat', 'p') 

% Spain
load('workspaces\Spain.mat')
p.pval(p.ind.phi0) = 0.72;
p.equilibrium = 'general';
[eql, sim, agg_stat] = compute_equilibrium(p);
disp(agg_stat)
save('workspaces\counterfactuals\Spain_1980.mat', 'eql', 'sim', 'agg_stat', 'p')
disp('experiment: 1980 regulation on temp. contracts done.')


%% 2. Decomposition of the effect of temporary contracts on flows

run('flow_decomposition_TC.m')
disp('experiment: decomposition of the effect of temp. contracts, done.')


%% 3. Firing costs experiments, in baseline

% load baseline
load('workspaces\Baseline.mat', 'p', 'agg_stat')

% firing cost vector
F = [0,0.5,1,3,6,12]*agg_stat.Wmn;
clearvars agg_stat
nb_exp = length(F);

% preallocate cells
eql = cell(nb_exp,1);
agg_stat = cell(nb_exp,1);

% US, with Temp jobs (no restriction)
p.pval(p.ind.phi0) = 0;
p.pval(p.ind.phi) = 0;
p.equilibrium = "general";
p.compute_welfare = true;

for ii = 1:nb_exp

    % update firing costs
    p.pval(p.ind.F) = F(ii);

    % compute equilibrium
    [eql{ii}, ~, agg_stat{ii}] = compute_equilibrium(p);

    % update reservation wage guess for next iteration
    p.pval(p.ind.wr) = eql{ii}.wr;

end

% save results
save('workspaces\counterfactuals\US_F', 'eql', 'agg_stat', 'p', 'F')
disp('experiment: firing costs (TC allowed) in baseline economy, done.')
