function table_U = U_decomposition(ue, eu)


%% contribution of aggregate transitions

% initialize table
tb = table;

% compute unemployment (as implied by discretization)
u = eu./(ue + eu);

% total unemployment change
tb.total = (u(2) - u(1))*100;

% contribution of ue rate change
tb.ue = (1/2)*( eu(2) / (eu(2) + ue(2)) -  eu(2) / (eu(2) + ue(1)) ...
               + eu(1) / (eu(1) + ue(2)) - eu(1) / (eu(1) + ue(1)) )*100;

% contribution of EU rate change
tb.eu = (1/2)*( eu(2) / (eu(2) + ue(2)) -  eu(1) / (eu(1) + ue(2)) ...
               + eu(2) / (eu(2) + ue(1)) - eu(1) / (eu(1) + ue(1)) )*100;

% computed and residual
tb.total_computed = tb.eu + tb.ue;
tb.residual = tb.total - tb.total_computed;


%% contribution of subcomponents

% % through UE 
% % temp job creation
% tb.temp_job_creation = table_UE.temp_job_creation(2)*tb.ue;
% 
% % job substitution
% tb.perm_job_substitution = table_UE.perm_job_substitution(2)*tb.ue;
% 
% 
% % through EU
% % within reallocation
% tb.PU_distribution = table_EU.PU_distribution(2)*tb.eu;
% tb.TU_distribution = table_EU.TU_distribution(2)*tb.eu;
% 
% % between reallocation
% tb.bw_distribution = table_EU.bw_distribution(2)*tb.eu;

%% output table
table_U = tb;

end