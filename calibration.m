%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project: Risk-sharing in a dual market
% Créchet (2020)
% matlab script file
% file name: "calibration.m"
% created: 12-06-2021
% last updated: 10-2023
%
% Description: model calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% random grid

% seed for global stream
rng('default')

% bounds
lb = p.cparam.bd(:,1);
ub = p.cparam.bd(:,2);

% draw matrix of random param vector
K = length( p.cparam.names );
nb_draws = p.random_grid_len;
x_grid = zeros(K,nb_draws);
for k = 1:K    
    x_grid(k,:) = lb(k) + (ub(k)-lb(k))*rand(1,nb_draws)';
end

% evaluate model over random grid
obj_grid = zeros(nb_draws,1);
parfor jj = 1:nb_draws
    
    obj_grid(jj) = objfunction(x_grid(:,jj),p);
    disp(obj_grid(jj))
    
end

% select the best points for initial param values for the next step
nb_select = p.initial_grid_len;
[ obj_sorted_grid, iindexes ] = sort( obj_grid );
x_sorted_grid = x_grid(:,iindexes); 
x_initial = x_sorted_grid(:,1:nb_select);
obj_initial = obj_sorted_grid(1:nb_select);

% save intermediate results (random grid evaluation)
save('workspaces\calibration\step1.mat', 'x_grid', 'obj_grid', 'x_initial', 'obj_initial')
disp('random grid saved')


%% run minimization at best points
load('workspaces\calibration\step1.mat')

% preallocate arrays
obj_minimized = 0*obj_initial;
x_minimized = 0*x_initial;

% locally minimize with Nelder Mead algo
for jj = 1:length(obj_minimized)
    
   x0 = x_initial(:,jj); % initial param values
   [x_minimized(:,jj), obj_minimized(jj)] = fminsearch( @(x) objfunction(x,p), x0 , p.fminsearch_opt);
   
end

% save results
save('workspaces\calibration\step2.mat', 'x_minimized', 'obj_minimized')


%% Evaluate model at calibrated param values (baseline)

% load result workspace
load('workspaces\calibration\step2.mat')

% approximated global minimum
[ obj_star, ii ] = min(obj_minimized(obj_minimized>0));

% solution
x_star = x_minimized(:,ii); 

% update parameter structure
p_baseline = p;

% calibrated parameters
p_baseline.pval( p.cparam.ind ) = x_star; 
p_baseline.ptable.pval = p_baseline.pval;
    
% compute equilibrium
eql = compute_equilibrium(p_baseline);

% implied parameters
p_baseline.pval( p.ind.kappa ) = eql.kkappa;
p_baseline.pval( p.ind.b ) = eql.b;

% table
p_baseline.ptable = table;
% p_baseline.ptable.pname = p_baseline.pnames;
p_baseline.ptable.pval = p_baseline.pval;

% check equilibrium 
p_baseline.equilibrium = "general";
[eql, sim, agg_stat, ten_stat] = compute_equilibrium(p_baseline);

% results
disp( p_baseline.ptable )
disp( agg_stat )
disp( ten_stat )

% save outcomes
p = p_baseline;
p.estimate_density = 'histogram';
save('workspaces/Baseline.mat', 'p', 'eql', 'sim', 'agg_stat', 'ten_stat');


%% Calibration to French institutions
load('workspaces\Baseline.mat')

% unstack some variables
pval = p.pval;
ind = p.ind;

% ``French-style'' calibrated parameters
pval(ind.b) = 0.50;                          % (b/Ew = 0.55) Dolado et al. (2011)
pval(ind.F) = 1.37;                          % (F/Ew = 1.5) Dolado et al. (2011), Cahuc et al. (2016)
pval(ind.A) = p.pval(ind.A)/2;               % Jung Kuhn (2014)
pval(ind.phi0) = 0;                          % no hiring regulation on TC
pval(ind.phi) = 0.045;                       % 3.43% TP probability (Givord (2010))
pval(ind.delta) = 0.000825;                    % unemployment rate 9.20%
pval(ind.wr) = 0.796943023844906;            % initial guess for reservation wage (set to be equal to wr consistent with calibration)

p_France = p;
p_France.pval = pval;
p_France.ptable.pval = pval;
disp(p_France.ptable);      

% compute French general equilibrium
p_France.equilibrium = 'general';
[eql, sim, agg_stat, ten_stat] = compute_equilibrium(p_France);

% display outcomes
disp('French-style calibrated economy: targeted outcomes')

disp('urate')
disp(agg_stat.U);

disp('b/Ew')
disp(agg_stat.b_Wmn)

disp('F/E(w)')
disp(agg_stat.F_Wmn)

disp('TP rate')
disp(agg_stat.TP)

disp('tshare')
disp(agg_stat.T)

disp('Simulated stats:')
disp( agg_stat )

% save
p = p_France;
p.pval(ind.wr) = eql.wr;
save( 'workspaces\France.mat', 'p', 'eql', 'sim', 'agg_stat', 'ten_stat' )


%% Spanish institutions
load('workspaces\Baseline.mat')

% unstack some variables
pval = p.pval;
ind = p.ind;

% Spanish-style'' calibrated parameters
pval(ind.b) = 0.53;                       % (b/Ew = 0.58) Dolado et al. (2011)
pval(ind.F) = 1.65;                       % (F/Ew = 1.8) Dolado et al. (2011), Cahuc et al. (2016)
pval(ind.A) = p.pval(ind.A)/2;            % Jung Kuhn (2014)
pval(ind.phi0) = 0;                       % no hiring regulation on TC
pval(ind.phi) = 0.030;                    % 2.16% TP rate (Toledo et al.)
pval(ind.delta) = 0.0019;                 % unemployment rate 13.28%
pval(ind.wr) = 0.790194732141052;         % initial guess for reservation wage

p_Spain = p;
p_Spain.pval = pval;                
p_Spain.ptable.pval = pval;
disp(p_Spain.ptable);

% compute Spanish general equilibrium
disp('General equilibrium')
p_Spain.equilibrium = 'general';
[eql, sim, agg_stat, ten_stat] = compute_equilibrium(p_Spain);

% outcomes
disp('Spanish-style calibrated economy: targeted outcomes')

disp('urate')
disp(agg_stat.U);

disp('b/Ew')
disp(agg_stat.b_Wmn)

disp('F/E(w)')
disp(agg_stat.F_Wmn)

disp('TP rate')
disp(agg_stat.TP)

disp('tshare')
disp(agg_stat.T)

disp('Simulated stats:')
disp( agg_stat )

% update guess for reservation wage
p_Spain.pval(ind.wr) = eql.wr;

% save
p = p_Spain;
p.pval(ind.wr) = eql.wr;
save('workspaces\Spain.mat', 'p', 'eql', 'sim', 'agg_stat', 'ten_stat' )
