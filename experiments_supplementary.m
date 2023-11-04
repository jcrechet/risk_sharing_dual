%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project: Risk-sharing in a dual market
% Créchet (2020)
% matlab script file
% file name: "experiment.m"
% created: 4-10-2023
% Description: quantitative experiments - supplementary appendix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% A1. Sources of differences between France and the US

% load baseline outcomes
load('workspaces\France.mat', 'p', 'agg_stat')
agg_stat_France = agg_stat;
p_France = p;
clearvars agg_stat p

load('workspaces\Baseline.mat', 'p', 'agg_stat')
agg_stat_US = agg_stat;
p_US = p;
clearvars agg_stat p

% arrays for outcome pp. difference
U  = zeros(5,1) + agg_stat_France.U;
UE = zeros(5,1) + agg_stat_France.UE;
EU = zeros(5,1) + agg_stat_France.EU;

% pp. difference between US and French benchmark outcomes
ii = 1;
U(ii)  = U(ii)-agg_stat_US.U;
UE(ii) = UE(ii)-agg_stat_US.UE;
EU(ii) = EU(ii)-agg_stat_US.EU;

% structures for counterfactual parameters
% indexes for parameters
ind = p_US.ind;
p_US_1 = cell(4,1);
for ii = 1:4
    p_US_1{ii} = p_US;
end
% A
p_US_1{1}.pval(ind.A)     = p_France.pval(ind.A);
% delta
p_US_1{2}.pval(ind.delta) = p_France.pval(ind.delta);
% b
p_US_1{3}.pval(ind.b)     = p_France.pval(ind.b);
% F
p_US_1{4}.pval(ind.F)     = p_France.pval(ind.F);

%  benchmark pp. difference
for ii = 2:5
    
    % compute equilibrium
    p_US_1{ii-1}.equilibrium = 'general';
    [~, ~, agg, ~] = compute_equilibrium(p_US_1{ii-1});
    
    % fill in vectors
    U(ii)  = U(ii)-agg.U;
    UE(ii) = UE(ii)-agg.UE;
    EU(ii) = EU(ii)-agg.EU;

end

% save
save('workspaces\counterfactuals\France_vs_US.mat', 'U', 'UE', 'EU')
disp('experiment: France vs. the U.S. done')



%% A2. Sources of differences between France and Spain

% load baseline outcomes
load('workspaces\France.mat')
agg_stat_France = agg_stat;
p_France = p;
clearvars agg_stat p

load('workspaces\Spain.mat')
agg_stat_Spain = agg_stat;
p_Spain = p;
clearvars agg_stat p


% arrays for outcome pp. difference
U  = zeros(5,1) + agg_stat_Spain.U;
T  = zeros(5,1) + agg_stat_Spain.T;
UE = zeros(5,1) + agg_stat_Spain.UE;
EU = zeros(5,1) + agg_stat_Spain.EU;
UP = zeros(5,1) + agg_stat_Spain.UP;
TP = zeros(5,1) + agg_stat_Spain.TP;
PU = zeros(5,1) + agg_stat_Spain.PU;
TU = zeros(5,1) + agg_stat_Spain.TU;

% pp. difference between Spanish and French benchmark outcomes
ii = 1;
U(ii)  = U(ii)-agg_stat_France.U;
T(ii)  = T(ii)-agg_stat_France.T;
UE(ii) = UE(ii)-agg_stat_France.UE;
EU(ii) = EU(ii)-agg_stat_France.EU;
UP(ii) = UP(ii)-agg_stat_France.UP;    
TP(ii) = TP(ii)-agg_stat_France.TP;   
PU(ii) = PU(ii)-agg_stat_France.PU;
TU(ii) = TU(ii)-agg_stat_France.TU;

% structures for counterfactual parameters
% indexes for parameters
ind = p.ind;
p_France_1 = cell(4,1);
for ii = 1:4
    p_France_1{ii} = p_France;
end
% delta
p_France_1{1}.pval(ind.delta) = p_Spain.pval(ind.delta);
% b
p_France_1{2}.pval(ind.b)     = p_Spain.pval(ind.b);
% F
p_France_1{3}.pval(ind.F)     = p_Spain.pval(ind.F);
% phi
p_France_1{4}.pval(ind.phi)   = p_Spain.pval(ind.phi);

%  pp. difference between France and Spain
for ii = 2:5
    
    % compute equilibrium
    [~, ~, agg, ~] = compute_equilibrium(p_France_1{ii-1});
    
    % fill in vectors
    U(ii)  = U(ii)-agg.U;
    T(ii)  = T(ii)-agg.T;
    UE(ii) = UE(ii)-agg.UE;
    EU(ii) = EU(ii)-agg.EU;
    UP(ii) = UP(ii)-agg.UP;    
    TP(ii) = TP(ii)-agg.TP;   
    PU(ii) = PU(ii)-agg.PU;
    TU(ii) = TU(ii)-agg.TU;

end

% save
save('workspaces\counterfactuals\France_vs_Spain.mat', 'U', 'T', 'UE', 'EU', 'UP', 'TP', 'PU', 'TU')
disp('experiment: France vs. Spain done')
