function param_structure = setparameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project: Risk-sharing in a dual market
% Créchet (2020)
% matlab function
% file name: "setparameters.m"
% created: 28-10-2020
% last updated: sept 2023
%
% Description: set model parameters
% This function generates a structure with parameter
% values for numerical computation of the model (and simulations)

% out: param_structure. Contains:
% (i) model economic parameters; (ii) parameters
% specifying options for numerical algorithm; (iii)  options for
% simulations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% model parameters

% specifying whether numerical model is solved in general eql
% select 'partial' for calibrating the model; select 'general' for
% experiments
p.equilibrium = 'partial';
p.compute_welfare = true;
p.estimate_density = 'histogram';

% initialize parameter cell, vector, and index
pnames = cell(0);    % cell for name of param
pval = zeros(1,1);   % vector for initial values (for evaluating equilibrium)
ind = struct;        % structure with indexes
jj = 1;              % initialize integer scalar for assigning index to each param.

% specify parameter names, values, and indexes
% discount factor
ind.beta = jj;         % assign index
pval(jj) = 0.9957;     % value
pnames{jj} = '\\beta'; % name (latex)
jj = jj+1;             % update integer

% risk aversion coeff
ind.sigma = jj;
pval(jj) = 2.5;
pnames{jj} = '\\sigma';
jj = jj+1;

% reservation wage 
% (note: treated as a parameter when calibrated the benchmark model in partial equil
% /endogenous when model is solved in gen. equil)
ind.wr = jj;
pval(jj) = 0.80;
pnames{jj} = 'w_R';
jj = jj+1;

% labor market tightness 
% (note: normalized to one when bench. model is calibrated
% /endogenous otherwise)
ind.theta = jj;
pval(jj) = 1;
pnames{jj} = '\\theta';
jj = jj+1;

% unemployment benefit
% (note: 'backout' using reservation wage and tightness normalization when
% model is calibrated/treated as a parameter when conputing gen.
% equilibrium)
ind.b = jj;
pval(jj) = 0.37565;
pnames{jj} = 'b';
jj = jj+1;

% exogenous separations prob.
ind.delta = jj;
pval(jj) = 0.0024;
pnames{jj} = '\\delta';
jj = jj+1;

% elasticity of matching wrt unemployment/bargaining power
ind.eta = jj;
pval(jj) = 0.5;
pnames{jj} = '\\eta';
jj = jj+1;

% efficiency of matching
ind.A = jj;
pval(jj) =  0.50;
pnames{jj} = 'A';
jj = jj+1;

% vacancy posting cost
ind.kappa = jj;
pval(jj) = 2.0598;
pnames{jj} = '\\kappa';
jj = jj+1;

% std of log match quality
ind.sigma_x = jj;
pval(jj) = 0.40;
pnames{jj} = '\\sigma_x';
jj = jj+1;

% prob of idiosyncratic pdty shock
ind.lambda = jj;
pval(jj) = 0.10;
pnames{jj} = '\\lambda';
jj = jj+1;

% support of match-specific pdty 
% (stoch. component---uniform distr)

% lower bound 
ind.zlb = jj;
pval(jj) = 0;
pnames{jj} = '\\underline{z}';
jj = jj+1;

% upper bound
ind.zub = jj;
pval(jj) = 1;
pnames{jj} = '\\overline{z}';
jj = jj+1;

% firing costs
ind.F = jj;
pval(jj) = 0;
pnames{jj} = 'F';
jj = jj+1;

% TC duration/renewal restriction
ind.phi = jj;
pval(jj) = 0;
pnames{jj} = '\\phi';
jj = jj+1;

% TC hiring restriction
ind.phi0 = jj;
pval(jj) = 0;
pnames{jj} = '\\phi_0';

% transpose to get col. vectors
pval = pval'; pnames = pnames';

% stack in a structure
p.ind = ind;
p.pval = pval;

% build a table with parameter
p.ptable = table(pnames, pval);

%% computation algorithm parameters

% numerical computation
p.Ix = 50;              % match quality grid
p.Iz = 50;              % grid for productivity shock
p.vfi_tol = 10^(-12);   % vfi tolerance for computing PC and TC values
p.vfi_tol_U = (1-pval(ind.beta))^(-1)*10^(-6);   % vfi tolerance for unemployment value
p.max_iter = 1000;      % max number of vfi iterations


%% simulation parameters

% number of employment spells to simulate
sim.Ks = 1e4;  

% max length of a simulated employemnt history
sim.Lz = 50*12; % set to 50 years

% set random number gen
sim.rng_algo = 'combRecursive';
p.sim = sim;


%% calibration/minimization algo

% select parameters to be internally calibrated
cparam.names = {'w_R', 'A', '\\lambda', 'sigma_x'};      % names
cparam.ind =   [ind.wr; ind.A; ind.lambda; ind.sigma_x]; % indexes
cparam.values = pval(cparam.ind);

% select bounds for initial random grid search
lb = [0.7, 0.4, 0.05, 0.25]';
ub = [0.9, 0.6, 0.1, 0.5]';
cparam.bd = [lb,ub];
p.cparam = cparam;

% options for Matlab local minimization function 'fminsearch'
p.fminsearch_opt = ...
    optimset('Display', 'iter', 'MaxFunEvals', ...
    500*length(cparam.values), 'MaxIter', 500*length(cparam.values), 'TolFun', 1e-6);

% size of the parameter grid for the initial random search, preliminary to
% running minimizer fminsearch
p.random_grid_len = 500;

% number of points selected from the random grid and used as initial values
% for random search
p.initial_grid_len = 10;


%% final output: parameter structure
param_structure = p;


end

