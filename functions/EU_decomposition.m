function tables_decomposition = EU_decomposition(n, sp, st, hp, ht, eu)

% Perform decomposition of equilibrium EU rate change across two policy
% regimes

%% in: n - vector of size 2, employment share of permanent jobs in the two
% regimes

% sp, st: matrices of size (I,2), where I is the size the discretized grid
% for match quality, separation policy functions in perm and temp jobs.

% hp, ht: matrices of size (I,2): discretized equilibrium distributions of
% match quality (simulated)

% table_U_decomposition: table with unemployment change components (UE and
% EU flows)

%% out: tables_decomposition - structure with 3 tables:
% (i) and (ii): PU and TU, components of changes in PU and TU rates across
% regimes
% (iii): EU: components of EU change.


%% initialize structure
tb = struct;

%% Component 1: difference between average PU and TU
tb.bw = table;

% total, implied by discretization
pu_mean = (hp(:,2)'*sp(:,2) + hp(:,1)'*sp(:,1))/2;
tu_mean = (ht(:,2)'*st(:,2) + ht(:,1)'*st(:,1))/2;
tb.bw.total = (pu_mean - tu_mean)*100;

% initialize
I = length(hp);
stilde = zeros(I,2);
htilde = zeros(I,2);

tb.bw.distribution = 0;
tb.bw.policy = 0;

% loop over policy regimes
for ii = 1:2
    
    % weights for distribution comp
    stilde(:,ii) = (hp(:,ii).*ht(:,ii)>0).*(sp(:,ii)+st(:,ii))/2 ...
        + (hp(:,ii).*ht(:,ii)==0).*( (hp(:,ii)>0).*sp(:,ii) + (ht(:,ii)>0).*st(:,ii) );
    
    % weights for separation policy comp. 
    htilde(:,ii) = (hp(:,ii).*ht(:,ii)>0).*(hp(:,ii)+ht(:,ii))/2;

    % compute component in each regime and sum
    tb.bw.distribution = tb.bw.distribution + (1/2)*(hp(:,ii)-ht(:,ii))'*stilde(:,ii)*100;
    tb.bw.policy = tb.bw.policy + (1/2)*(sp(:,ii)-st(:,ii))'*htilde(:,ii)*100;

end

tb.bw.total_computed = tb.bw.distribution + tb.bw.policy;
tb.bw.residual = tb.bw.total - tb.bw.total_computed;


%% Component 2 & 3: change in PU and TU
% initialize table
tb.PU = table;
tb.TU = table;

for ii = 1:2 
        
    tb_tmp = table;

    if ii == 1
        s = sp; h = hp;
    else
        s = st; h = ht;
    end

% compute total implied by discretization of x
tb_tmp.total = (h(:,2)'*s(:,2) - h(:,1)'*s(:,1))*100;

% (i) separation policy
tmp = (h(:,1).*h(:,2)>0).*(s(:,2)-s(:,1)).*(h(:,1)+h(:,2))/2;
tb_tmp.policy = sum(tmp)*100;

% (ii) distribution
s_tilde = (h(:,1).*h(:,2)>0).*(s(:,1)+s(:,2))/2 ...
           + (h(:,1)>0 & h(:,2)==0).*s(:,1) ...
           + (h(:,1)==0 & h(:,2)>0).*s(:,2);
tmp = s_tilde.*(h(:,2) - h(:,1));
tb_tmp.distribution = sum(tmp)*100;


% total, computed
tb_tmp.total_computed = tb_tmp.policy + tb_tmp.distribution;

% residual
tb_tmp.residual = tb_tmp.total - tb_tmp.total_computed;

% stack in table
if ii == 1
   tb.PU = tb_tmp;
else
   tb.TU = tb_tmp;
end

end

%% between + within table

% total change
tb.EU = table;
tb.EU.total = (eu(2) - eu(1))*100;

% total change implied by discretization
tb.EU.total_disc = (n(2)*(hp(:,2)'*sp(:,2)) + (1-n(2))*(ht(:,2)'*st(:,2)) ...
                          - n(1)*(hp(:,1)'*sp(:,1)) - (1-n(1))*(ht(:,1)'*st(:,1)))*100;


% PU change contribution
n_mean =  (n(1)+n(2))/2;
tb.EU.PU = n_mean*tb.PU.total_computed;
tb.EU.PU_distribution = n_mean*tb.PU.distribution;

% TU change contribution
tb.EU.TU = (1-n_mean)*tb.TU.total_computed;
tb.EU.TU_distribution = (1-n_mean)*tb.TU.distribution;

% reallocation contributions
tb.EU.bw_distribution = (n(2)-n(1))*tb.bw.distribution;
tb.EU.bw_policy = (n(2)-n(1))*tb.bw.policy;

% distribution (between + within)
tb.EU.tot_distribution = tb.EU.bw_distribution + tb.EU.PU_distribution + tb.EU.TU_distribution;

% total, computed and residual
tb.EU.total_computed = tb.EU.PU + tb.EU.TU + tb.EU.bw_distribution + tb.EU.bw_policy;
tb.EU.residual = tb.EU.total - tb.EU.total_computed;

% contribution shares
% between
tb.bw = compute_shares(tb.bw);

% PU
tb.PU = compute_shares(tb.PU);

% TU
tb.TU = compute_shares(tb.TU);

% EU
tb.EU = compute_shares(tb.EU);


%% output
tables_decomposition = tb;


end


function table = compute_shares(table)

% compute contribution share of each element of employment outlflow rate in
% table

tmp = table2array(table(1,:)) / table.total;
tmp = array2table(tmp);
table(2,:) = tmp;

end