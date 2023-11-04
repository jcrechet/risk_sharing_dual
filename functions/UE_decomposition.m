function table_decomposition = UE_decomposition(phi, ptheta, xt, xp, xhat, G, ue)

%% UE decomposition
% decompose change in UE rate across policy regimes

%% in 
% phi: size-2 vector, hiring restriction param across policy regimes
% ptheta: size-2 vector, contact probability across regimes
% xt: size-2 vector, acceptance cutoff in match quality, temp contracts
% xp: size-2 vector, acceptance, perm contracts
% xhat: size-2 vector, acceptance, perm contracts when temp contracts are
% available
% G: function, sampling cdf of match quality upon matching

%% out
% table_decomposition: table with components of UE total change across
% regimes. Components: thightness, contribution of change in contact
% probability; acceptance, change in hiring probability conditional on a
% contacts; job creation: tightness + 'gross' change in job acceptance
% prob. 

%% perform

% initialize table
tb = table;

% acceptance probability
A = (1-phi).*(1 - G(xt)) + phi.*(1 - G(xp));

% total
%ue = ptheta.*A;
tb.total = (ue(2) - ue(1))*100;

% tightness
tb.tightness = (ptheta(2) - ptheta(1))*(A(1) + A(2))/2*100;

% selection
tb.selection = (A(2) - A(1))*(ptheta(1) + ptheta(2))/2*100;

% job creation
xhat_m = mean(xhat);
xp_m = mean(xp);
xt_m = mean(xt);
tb.temp_job_creation = - (phi(2) - phi(1))*(G(xhat_m) - G(xt_m))*(ptheta(1) + ptheta(2))/2*100;

% job substitution (approximated)
tb.perm_job_substitution = (phi(2) - phi(1))*(G(xhat_m) - G(xp_m))*(ptheta(1) + ptheta(2))/2*100;

% total computed (with approximation)
tb.total_computed = tb.tightness + tb.temp_job_creation + tb.perm_job_substitution;
tb.residual = tb.total - tb.total_computed;

% contribution shares
tmp = table2array(tb(1,:))/tb.total;
tmp = array2table(tmp);
tb(2,:) = tmp;

%% out table
table_decomposition = tb;


end