%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Risk-sharing in a dual market
% Créchet (2020)
% Matlab script file
% file name: "flow_decomposition_TC.m"
% last updated: Oct 2023
%
% Description: Decomposition of UE, EU, and U rate change in response to liberalization
% of temporary contracts.

% Based on equations (44) and (45), and Supplementary Appendix C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for c = 1:2
    
    if c == 1
        ctry = 'France';
    else
        ctry = 'Spain';
    end

    % regime 1: pre-reform (TC restricted)
    load(['workspaces\counterfactuals\', ctry, '_1980.mat'])

    % statistics
    n = 1-agg_stat.T;
    ue = agg_stat.UE;
    eu = agg_stat.EU;

    % parameters
    phi = p.pval(p.ind.phi0);
    A = p.pval(p.ind.A);
    eeta = p.pval(p.ind.eta);
    sigma_x = p.pval(p.ind.sigma_x);
    mu_x = -sigma_x^2/2;
    G = @(x) logncdf(x, mu_x, sigma_x);

    % policy and distributions
    sp = eql.sp;
    st = eql.st;
    hp = eql.hp;
    ht = eql.ht;
    xt = eql.xt;
    xp = eql.xp;
    xhat = eql.xhat;

    % tightness
    theta = eql.theta;

    % regime 2: liberalization
    load(['workspaces\',ctry,'.mat'])

    % statistics
    phi(2) = p.pval(p.ind.phi0);
    n(2) = 1 - agg_stat.T;
    ue(2) = agg_stat.UE;
    eu(2) = agg_stat.EU;

    % policy and distributions
    sp(:,2) = eql.sp;
    st(:,2) = eql.st;
    hp(:,2) = eql.hp;
    ht(:,2) = eql.ht;
    xt(2) = eql.xt;
    xp(2) = eql.xp;
    xhat(2) = eql.xhat;

    % tightness
    theta(2) = eql.theta;

    % run EU decomposition
    decompositions = EU_decomposition(n, sp, st, hp, ht, eu);

    % UE decomposition
    ptheta = A*theta.^(1-eeta);
    decompositions.UE = UE_decomposition(phi, ptheta, xt, xp, xhat, G, ue);

    % U decomposition
    decompositions.U = U_decomposition(ue, eu);

    % save
    save(['workspaces\',ctry,'.mat'], 'decompositions', '-append')

end







