%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Risk-sharing in a dual market
% Créchet (2020)
% Matlab script file
% file name: "tables.m"
% last updated: Oct 2023
%
% Description: script exporting tables from model's
% quantitative results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Path to table folder
table_path = '_results\tables\';


%% Table 1: Benchmark parameter values

% Parameter names and description
parameters = {'\beta', '\sigma', '\eta', '\gamma', '\delta', 'F', ...
    'b', 'A', '\lambda', '\sigma_x^2', '\mu_x', '\kappa'};

descriptions = {'Discount factor', 'Relative risk aversion', ...
    'Elasticity of matching', 'Worker''s bargaining power', ...
    'Exogenous separation probability', 'Firing costs', ...
    'Non-work income', 'Matching efficiency', ...
    'Probability of idiosyncratic shock', ...
    'Log match-quality variance', 'Log match-quality mean', ...
    'Vacancy posting cost'};

% Parameter values
load('workspaces\Baseline.mat', 'p')
x = p.ptable.pval;
ind = p.ind;
mu_x = - x(ind.sigma_x)^2/2;

% Function to convert num. values in string
s = @(x, fmt) num2str( x, fmt);

values = {s(x(ind.beta), '%.3f'), s(x(ind.sigma), '%.1f'), s(x(ind.eta), '%.1f'), s(x(ind.eta), '%.1f'), s(x(ind.delta), '%.4f'), s(x(ind.F), '%.0f'), ...
    s(x(ind.b), '%.3f'), s(x(ind.A), '%.3f'), s(x(ind.lambda), '%.3f'), s(x(ind.sigma_x)^2, '%.3f'), s(mu_x, '%.3f'), s(x(ind.kappa), '%.3f')};

clearvars p

% Open file for writing
file = [table_path, 'benchmark_parameters.tex'];
fid = fopen(file, 'w');

% Start table environment
fprintf(fid, '\\begin{table}[!h]\n');
fprintf(fid, '\t\\centering\n');
fprintf(fid, '\t\\captionof{table}{Benchmark parameter values}\n');
fprintf(fid, '\t\\label{tab:param:benchmark}\n');
fprintf(fid, '\t\\begin{tabular}{l l c}\n');
fprintf(fid, '\t\t\\hline \\hline\n');

% Iterate through parameters and write to file
for i = 1:length(parameters)
    if i == 1
        fprintf(fid, '\t\t\\addlinespace\n');
        fprintf(fid, '\t\t\\textit{Preset} & & \\\\ \n');
    end

    if i == 7
        fprintf(fid, '\t\t\\addlinespace\n');
        fprintf(fid, '\t\t\\textit{Internal} & & \\\\ \n');
    end

    fprintf(fid, '\t\t$%s$ & %s & %s \\\\ \n', parameters{i}, descriptions{i}, values{i});
end

% Table notes
notes = ['\t\\caption*{\\footnotesize Notes: Parameter values in the benchmark U.S.\\ calibrated economy. ' ...
    'The exogenous probability of separation $\\delta$ is set following \\cite{jung_kuhn:JEEA:2019}. ' ...
    'The other preset parameters are standard. The targets for the internally calibrated parameters ' ...
    'are an unemployment rate equal to 5.6\\%%, an aggregate EU probability of 2\\%%, non-work income $b=0.4$ ' ...
    '(\\cite{shimer:2005:AER}), and relative separation rates by tenure estimated from the Job Tenure Supplement ' ...
    'of the CPS (less than one year of tenure to 1-15 years). Moreover, $\\mu_x = -\\sigma_x^2/2$. The value of $\\kappa$ ' ...
    'is set to be consistent with the normalization that labor-market tightness $\\theta=1$ (\\cite{shimer:2005:AER}). ' ...
    'In addition, the model is calibrated in partial equilibrium with the value of unemployment $U$ exogenous and the ' ...
    'parameter $b$ is computed accordingly (see the text for details). The model fit to the calibration targets is ' ...
    'shown in Table \\ref{tab:targets}, Panel a).}'];

% Close table environment
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\hline \\hline\n');
fprintf(fid, '\t\\end{tabular}\n');
fprintf(fid,  notes);
fprintf(fid, '\\end{table}\n');

% Close file
fclose(fid);

% Cleanup
clearvars -EXCEPT table_path

disp('Table 1 exported')


%% Table 2: Model outcomes and data targets

% Open file for writing
file = [table_path, 'targeted_moments.tex'];
fid = fopen(file, 'w');

% Write table beginning structure
fprintf(fid, '\\begin{table}[!h]\n');
fprintf(fid, '\t\\centering\n');
fprintf(fid, '\t\\captionof{table}{Data and model targeted statistics}\n');
fprintf(fid, '\t\\label{tab:targets}\n');
fprintf(fid, '\t\\begin{tabular}{l c c}\n');
fprintf(fid, '\t\t\\hline \\hline\n');

% Define header
headers = {'\hspace{250pt}', '\textit{Target}', '\textit{Model}'};
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t%s & %s & %s \\\\ \n', headers{:});
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\addlinespace\n');

% Define sections, descriptions, and target values
sections = {'\textit{a) U.S.\ benchmark}', '\textit{b) ``French'''' counterfactual}', '\textit{c) ``Spanish'''' counterfactual}'};
descriptions = {
    '\hspace{10pt} Unem.\ rate (\%)', '\hspace{10pt} EU rate  (\%)', '\hspace{10pt} Sep.\ rate tenure ratio', '\hspace{10pt} Non-work income, $b$', ...
    '\hspace{10pt} Unem.\ rate', '\hspace{10pt} TP rate', '\hspace{10pt} $b/E(w)$', '\hspace{10pt} $F/E(w)$', ...
    '\hspace{10pt} Unem.\ rate', '\hspace{10pt} TP rate', '\hspace{10pt} $b/E(w)$', '\hspace{10pt} $F/E(w)$'
    };
targets = {'5.61', '2.00', '3.13', '0.40', '9.20', '3.43', '0.55', '1.50', '13.28', '2.16', '0.58', '1.80'};

% Model outcomes
% Baseline
load('workspaces\Baseline.mat', 'agg_stat')
model = [agg_stat.U*100, agg_stat.EU*100, agg_stat.sr_tn, agg_stat.b];
clearvars agg_stat ten_stat

% French
load('workspaces\France.mat', 'agg_stat')
model = [model, agg_stat.U*100, agg_stat.TP*100, agg_stat.b_Wmn, agg_stat.F_Wmn];
clearvars agg_stat

% Spanish
load('workspaces\Spain.mat', 'agg_stat')
model = [model, agg_stat.U*100, agg_stat.TP*100, agg_stat.b_Wmn, agg_stat.F_Wmn];
clearvars agg_stat

% Function to convert num. values in string
s = @(x) num2str( x, '%.2f');

% Write rows of table
for i = 1:length(descriptions)
    if i == 1 || i == 5 || i == 9
        fprintf(fid, '\t\t\\addlinespace\n');
        fprintf(fid, '\t\t%s & & \\\\ \n', sections{(i-1)/4 + 1});
    end
    fprintf(fid, '\t\t%s & %s & %s \\\\ \n', descriptions{i}, targets{i}, s(model(i)) );
end

fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\addlinespace\n');

% Table notes
notes = ['\t\\caption*{\\footnotesize Notes: targeted and model statistics for the U.S.\\ benchmark economy (Panel a), the counterfactual model economies '...
'with ``French-type'''' institutions (b) and with ``Spanish-type'''' institutions (c). The U.S.\\ job tenure separation ratio ' ...
'(sep.\\ rate of jobs with less than one year to jobs with one to fifteen years) is estimated using CPS Job supplement ' ...
'data (1996-2010). The U.S.\\ value of non-work income of $b=0.4$ is from \\cite{shimer:2005:AER}; replacement ratios ' ...
'($b/E(w)$) are from \\cite{bentolila_cahuc_dolado_lebarbanchon:2012:EJ}, and firing cost (relative to the average ' ...
'equilibrium wage $F/E(w)$) are from \\cite{cahuc_charlot_malherbet:2016:IER}. Transition probabilities are estimates '...
'from \\cite{jung_kuhn:2014:EJ} (U.S.), \\cite{gw2015} (France), and \\cite{Silva_Vazquez_Grenno:2013:LabourEcon} (Spain).'...
 '}\n' ];

% Write table end structure
fprintf(fid, '\t\t\\hline \\hline\n');
fprintf(fid, '\t\\end{tabular}\n');
fprintf(fid, notes);
fprintf(fid, '\\end{table}\n');

% Close file
fclose(fid);

% cleanup
clearvars -EXCEPT table_path

disp('Table 2 exported')

%% Table 3: Counterfactual parameter values

% Open file for writing
file = [table_path, 'counterfactual_parameters.tex'];
fid = fopen(file, 'w');

% Write table beginning structure
fprintf(fid, '\\begin{table}[!h]\n');
fprintf(fid, '\t\\centering\n');
fprintf(fid, '\t\\captionof{table}{Parameter values for the French and Spanish counterfactual model economies}\n');
fprintf(fid, '\t\\label{tab:param:counterfactual}\n');
fprintf(fid, '\t\\begin{tabular}{l l c c c}\n');
fprintf(fid, '\t\t\\hline \\hline\n');

% Define headers
fprintf(fid, '\t\t  & & \\multicolumn{2}{c}{\\textit{Counterfactual}} & \\textit{Benchmark} \\\\ \n');
fprintf(fid, '\t\t& & \\textit{France} & \\textit{Spain} \\\\ \n');
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\addlinespace\n');

% Define parameter names and description
parameters = {'$A$', '$b$', '$F$', '$\phi$', '$\delta$'};
descriptions = {
    'Matching efficiency    \hspace{80pt}', ...
    'Non-work income', ...
    'Firing costs', ...
    'TC duration regulation', ...
    'Exogenous separation'
    };

% Parameter values
% format
s = @(x, fmt) num2str(x, fmt);

% function to extract values for each calibrated economy
values = @(x, ind) { s(x(ind.A), '%.3f'), s(x(ind.b), '%.3f'), s(x(ind.F), '%.2f'), s(x(ind.phi), '%.3f'), s(x(ind.delta), '%.4f') } ;

% France
load('workspaces\France.mat', 'p')
franceValues = values(p.pval, p.ind);
clearvars p;

% Spain
load('workspaces\Spain.mat', 'p')
spainValues = values(p.pval, p.ind);
clearvars p;

% Baseline
load('workspaces\Baseline.mat', 'p')
benchmarkValues = values(p.pval, p.ind);
benchmarkValues{3} = '0';
benchmarkValues{4} = '-';
clearvars p;

% Write rows of table
for i = 1:length(parameters)
    fprintf(fid, '\t\t%s & %s & %s & %s & %s \\\\ \n', parameters{i}, descriptions{i}, franceValues{i}, spainValues{i}, benchmarkValues{i});
end

% Table notes
notes = ['\\caption*{\\footnotesize Notes: parameter values for the ``French`'' and ``Spanish'' counterfactual ' ...
     'economies and the benchmark. The matching efficiency $A$ and TC hiring regulation strictness $\\phi_0$ are preset; ' ...
     '$A$ is equal to half of its benchmark and $\\phi_0=0$ (no hiring regulatory constraints on TC). ' ...
     'The parameters $b$, $F$, $\\phi$, and $\\delta$ are set to match the statistics in Table \\ref{tab:targets}.}'];

% Write table end structure
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\hline \\hline\n');
fprintf(fid, '\t\\end{tabular}\n');
fprintf(fid, notes);
fprintf(fid, '\\end{table}\n');

% Close file
fclose(fid);

% cleanup
clearvars -EXCEPT table_path

disp('Table 3 exported')

%% Table 4: Untargeted moments

% Open file for writing
file = [table_path, 'untargeted_moments.tex'];
fid = fopen(file, 'w');

% Write table beginning structure
fprintf(fid, '\\begin{table}[!h]\n');
fprintf(fid, '\t\\centering\n');
fprintf(fid, '\t\\captionof{table}{Counterfactual stock and flows (non-targeted)}\n');
fprintf(fid, '\t\\label{tab:outcomes}\n');
fprintf(fid, '\t\\begin{tabular}{l c c c c}\n');
fprintf(fid, '\t\t\\hline \\hline\n');

% Define headers
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\hspace{100pt} & \\multicolumn{2}{c}{\\textit{France}} & \\multicolumn{2}{c}{\\textit{Spain}} \\\\ \n');
fprintf(fid, '\t\t                & \\textit{Data} & \\textit{Model} & \\textit{Data} & \\textit{Model} \\\\ \n');
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\textit{Transition rates (\\%%)} \\\\ \n');

% Define descriptions
descriptions = {
    '\hspace{5pt} Unem.\ to Emp.\ (UE)', ...
    '\hspace{5pt} Emp.\ to Unem.\ (EU)', ...
    '\hspace{5pt} Unem.\ to Perm.\ (UP)', ...
    '\hspace{5pt} Perm.\ to Unem.\ (PU)', ...
    'Employment stock, temp.\ contracts (TC) share (\%)  \hspace{20pt}', ...
    'Unemployment to employment flows, TC share', ...
    'Employment exit flows, TC share'
    };

% Values
franceData = {'8.43', '1.09', '1.83', '0.33', '11.88', '73.18', '53.93'};
load('workspaces\France.mat', 'agg_stat')
franceModel = [agg_stat.UE, agg_stat.EU, agg_stat.UP, agg_stat.PU, agg_stat.T, agg_stat.T_inflow, agg_stat.T_outflow]*100;
clearvars agg_stat

spainData = {'5.32', '0.76', '0.41', '0.15', '24.26',  '92.29', '67.36'};
load('workspaces\Spain.mat', 'agg_stat')
spainModel = [agg_stat.UE, agg_stat.EU, agg_stat.UP, agg_stat.PU, agg_stat.T, agg_stat.T_inflow, agg_stat.T_outflow]*100;
clearvars agg_stat

% Format
s = @(x) num2str(x, '%.2f');

% Write rows of table
for i = 1:length(descriptions)
    fprintf(fid, '\t\t%s & %s & %s & %s & %s \\\\ \n', descriptions{i}, franceData{i}, s(franceModel(i)), spainData{i}, s(spainModel(i)));
    if i == 4
        fprintf(fid, '\t\t\\addlinespace\n');
    end
end

fprintf(fid, '\t\t\\addlinespace\n');
fprintf(fid, '\t\t\\addlinespace\n');

% Notes
notes = [' \\caption*{\\footnotesize Notes: select statistics for the counterfactual ``French'''' and ``Spanish'''' model economies. ' ...
 'The transition  probabilities are taken from \\cite{gw2015} ' ...
 'for France (private sector excluding self-employment) and \\cite{Silva_Vazquez_Grenno:2013:LabourEcon} for Spain. ' ...
 'Monthly transition probabilities for France are computed by applying the time-aggregation adjustment in ' ...
 '\\cite{shimer:2012:RED} to quarterly transition probabilities in \\cite{gw2015}. The empirical employment ' ...
 'share of temporary jobs is computed from steady-state stocks implied by transition probabilities in \\cite{gw2015} ' ...
 'and \\cite{Silva_Vazquez_Grenno:2013:LabourEcon}. Empirical employment exit flows are defined as flows into unemployment ' ...
 'and non-participation.}'];

% Write table end structure
fprintf(fid, '\t\t\\hline \\hline\n');
fprintf(fid, '\t\\end{tabular}\n');
fprintf(fid, notes);
fprintf(fid, '\\end{table}\n');

% Close file
fclose(fid);

% cleanup
clearvars -EXCEPT table_path

disp('Table 4 exported')

%% Table 5: The employment effect of temporary contract reforms, France and Spain

% preallocate vectors for statistics
U = zeros(4,1);
T = U;
UE = U;
EU = U;
UP = U;
TP = U;
PU = U;
TU = U;

Yagg = U;
Ymn = U;
Wmn = U;

cm = U;
cu = U;
cpop = U;

theta = U;
ptheta = U;
Gp = U;
Gt = U;
Acc = U;


for i = 1:length(U)

    % France, reference
    if i==1
        load('workspaces\France.mat', 'agg_stat', 'p')
        pphi = p.pval(p.ind.phi0);
        eeta = p.pval(p.ind.eta);
        A = p.pval(p.ind.A);
        sigma_x = p.pval(p.ind.sigma_x);
        mu_x = - sigma_x^2/2;
        clearvars p
    end

    % France, pre-reforms
    if i==2
        load('workspaces\counterfactuals\France_1980.mat', 'agg_stat', 'p')
        pphi = p.pval(p.ind.phi0);
        clearvars p
    end

    % Spain, reference
    if i==3
        load('workspaces\Spain.mat', 'agg_stat', 'p')
        pphi = p.pval(p.ind.phi0);
        clearvars p
    end

    % Spain, pre-reforms
    if i==4
        load('workspaces\counterfactuals\Spain_1980.mat', 'agg_stat', 'p')
        pphi = p.pval(p.ind.phi0);
        clearvars p
    end


    % fill in vector
    U(i) = agg_stat.U;
    T(i) = agg_stat.T;
    UE(i) = agg_stat.UE;
    EU(i) = agg_stat.EU;
    UP(i) = agg_stat.UP;
    TP(i) = agg_stat.TP;
    PU(i) = agg_stat.PU;
    TU(i) = agg_stat.TU;

    Yagg(i) = agg_stat.Ytot;
    Ymn(i)  = agg_stat.Ymn;
    Wmn(i)  = agg_stat.Wmn;

    cu(i) = agg_stat.cU;
    cm(i)  = agg_stat.cV + agg_stat.cJ;
    cpop(i) = agg_stat.cTot;

    theta(i)  = agg_stat.theta;
    ptheta(i) = A*theta(i)^(1-eeta);
    Gp   = 1-normcdf(log( agg_stat.xp ), mu_x, sigma_x);
    Gt   = 1-normcdf(log( agg_stat.xt ), mu_x, sigma_x);
    Acc(i) = (1-pphi)*Gt + pphi*Gp;

    clearvars agg_stat

end


% function for printing line in table
s = @(v) [ num2str(v(1)*100, '%.2f'),' & ',num2str(v(2)*100, '%.2f'), ...
    ' & ',num2str(v(3)*100, '%.2f'),' & ', num2str(v(4)*100, '%.2f')];

% table
TC_liberalization = [ '\\begin{table}[!h]'   '\n' ...
    '\\centering'                    '\n' ...
    '\\captionof{table}{The employment effect of temporary contracts}' '\n' ...
    '\\begin{tabular}{l c c c c}'                           '\n' ...
    '\\hline \\hline'                                       '\n' ...
    '\\addlinespace'                                        '\n' ...
    ' & \\multicolumn{2}{c}{ \\textit{French institutions} } & \\multicolumn{2}{c}{ \\textit{Spanish institutions} } ' '\\\\' '\n' ...
    ' & \\textit{No restrictions}  & \\textit{Pre-reform} &  \\textit{No restrictions}  & \\textit{Pre-reform} ' '\\\\' '\n' ...
    ' & $\\phi_0=0$ & $\\phi_0=0.55$ & $\\phi_0=0$ & $\\phi_0=0.72$ ' '\\\\' '\n' ...
    ' \\addlinespace '  '\n' ...
    ' Unemp.\\ rate (\\%%)                                &' s(U)  '\\\\' '\n' ...
    ' TC emp.\\ share                              &' s(T)  '\\\\' '\n' ...
    ' \\addlinespace'                        '\n' ...
    'UE rate (\\%%)           &' s(UE)  '\\\\' '\n' ...
    'UP rate                  &' s(UP)  '\\\\' '\n' ...
    ' \\addlinespace' '\n' ...
    'EU rate (\\%%)           &' s(EU)  '\\\\' '\n' ...
    'PU rate                  &' s(PU)  '\\\\' '\n' ...
    '\\addlinespace' '\n' ...
    '\\hline \\hline'                                       '\n' ...
    '\\end{tabular}'  '\n' ...
    '\\label{tab:effect_EPL}'        '\n' ...
    '\\caption*{ \\footnotesize Notes: evaluation of the effect of temporary contracts on employment stocks and flows in the model economies with ' ...
    'French and Spanish-type institutions. The columns ``No restrictions'''' present equilibrium outcomes (in percentage) in reference economies without hiring restrictions ' ...
    'on temporary contracts, i.e., $\\phi_0=0$; the columns ``Pre-reforms'''' present outcomes in counterfactual economies with $\\phi_0$ that is set to generate ' ...
    'an employment share of temporary jobs equal to 5.5\\%%, equal to the average in Continental Europe in 1983 (\\cite{faccini:2014:EJ}).}' ...
    '\\end{table}' ];
% Export table
% open (or create) tex file
fileID = fopen([table_path, 'TC_liberalization.tex'], 'w');
% write table
fprintf(fileID, TC_liberalization);
% close file
fclose(fileID);

% cleanup
% clearvars -except table_path
disp('Table 5 exported')


    %' \\hspace{5pt} \\textit{Contact rate, $p(\\theta)$}         &' s(ptheta)  '\\\\' '\n' ...
    %' \\hspace{5pt} \\textit{Acceptance rate, $\\mathcal{A}$}    &' s(Acc)  '\\\\' '\n' ...
    %' \\hspace{5pt} \\textit{Temp.\\ to unemp.\\ (TU)}           &' s(TU)  '\\\\' '\n' ...


%% Table 6: Decomposition of the effect of liberalization of TC on employment flows

d = cell(2,1);
load('workspaces\France.mat', 'decompositions')
d{1} = decompositions;
clearvars decompositions

load('workspaces\Spain.mat', 'decompositions')
d{2} = decompositions;
clearvars decompositions

% EU component
eu = [d{1}.EU.total(1), d{1}.U.eu, d{2}.EU.total(1), d{2}.U.eu]; % total

retention = [d{1}.EU.bw_policy(1), d{1}.EU.bw_policy(2)*d{1}.U.eu, ...
    d{2}.EU.bw_policy(1), d{2}.EU.bw_policy(2)*d{2}.U.eu ]; % policy

distr = [d{1}.EU.tot_distribution(1), d{1}.EU.tot_distribution(2)*d{1}.U.eu, ...
    d{2}.EU.tot_distribution(1), d{2}.EU.tot_distribution(2)*d{2}.U.eu ]; % total distribution

bw_distr = [d{1}.EU.bw_distribution(1), d{1}.EU.bw_distribution(2)*d{1}.U.eu, ...
    d{2}.EU.bw_distribution(1), d{2}.EU.bw_distribution(2)*d{2}.U.eu ]; % total distribution

tu_distr = [d{1}.EU.TU_distribution(1), d{1}.EU.TU_distribution(2)*d{1}.U.eu, ...
    d{2}.EU.TU_distribution(1), d{2}.EU.TU_distribution(2)*d{2}.U.eu ]; % TU (distribution)

pu_distr = [d{1}.EU.PU_distribution(1), d{1}.EU.PU_distribution(2)*d{1}.U.eu, ...
    d{2}.EU.PU_distribution(1), d{2}.EU.PU_distribution(2)*d{2}.U.eu ]; % PU (distribution)


% UE component
ue = [d{1}.UE.total(1), d{1}.U.ue, d{2}.UE.total(1), d{2}.U.ue]; % total

tightness = [d{1}.UE.tightness(1), d{1}.UE.tightness(2)*d{1}.U.ue, ...
    d{2}.UE.tightness(1), d{2}.UE.tightness(2)*d{2}.U.ue]; % tightness

selection = [d{1}.UE.selection(1), d{1}.UE.selection(2)*d{1}.U.ue, ...
    d{2}.UE.selection(1), d{2}.UE.selection(2)*d{2}.U.ue]; % selection

creation = [d{1}.UE.temp_job_creation(1), d{1}.UE.temp_job_creation(2)*d{1}.U.ue, ...
    d{2}.UE.temp_job_creation(1), d{2}.UE.temp_job_creation(2)*d{2}.U.ue]; % creation of temp contracts

substitution = [d{1}.UE.perm_job_substitution(1), d{1}.UE.perm_job_substitution(2)*d{1}.U.ue, ...
    d{2}.UE.perm_job_substitution(1), d{2}.UE.perm_job_substitution(2)*d{2}.U.ue]; % substitution of perm contracts


% format
s = @(v) ...
    [ num2str(v(1), '%.2f'),' & ',num2str(v(2), '%.2f'),' & ',num2str(v(3), '%.2f'),' & ',...
    num2str(v(4), '%.2f') ];

% table
tab = ['\\begin{table}[!h]' '\n' ...
    '\\centering' '\n' ...
    '\\captionof{table}{Decomposition of the effect of temporary contracts on employment steady-state flows}'   '\n' ...
    '\\label{tab:flow_decomposition} '          '\n' ...
    '\\begin{tabular}{l c c c c}' 		 '\n' ...
    '\\hline \\hline' '\n' ...
    '\\addlinespace'  '\n' ...
    ' \\hspace{100pt} &  \\multicolumn{2}{c}{\\textit{France}}  & \\multicolumn{2}{c}{\\textit{Spain}}  \\\\ ' '\n' ...
    '    &  \\textit{Absolute} &  \\textit{Unem.\\ change} & \\textit{Absolute} & \\textit{Unem.\\ change} \\\\ ' '\n' ...
    '\\addlinespace'  '\n' ...
    '\\textit{UE rate change (p.p.)}   &' s(ue)   '\\\\' '\n' ...
    '\\hspace{2.5pt} (i) tightness  &' s(tightness)  '\\\\' '\n' ...
    '\\hspace{2.5pt} (ii) selection &'  s(selection)  '\\\\' '\n' ...
    '\\hspace{10pt} (ii-a) temp.\\ job creation &' s(creation)  '\\\\' '\n' ...
    '\\hspace{10pt} (ii-b) perm.\\ job substitution &' s(substitution)  '\\\\' '\n' ...
    '\\addlinespace'  '\n' ...
    '\\textit{EU rate change (p.p.)}   &' s(eu)   '\\\\' '\n' ...
    '\\hspace{2.5pt} (i) retention &' s(retention) '\\\\' '\n' ...
    '\\hspace{2.5pt} (ii) reallocation, total &' s(distr) '\\\\' '\n' ...
    '\\hspace{10pt} (ii-a) between perm.\\ and temp.\\ jobs &' s(bw_distr) '\\\\' '\n' ...
    '\\hspace{10pt} (ii-b) within perm.\\ jobs &' s(pu_distr) '\\\\' '\n' ...
    '\\hspace{10pt} (ii-c) within temp.\\ jobs &' s(tu_distr) '\\\\' '\n' ...
    '\\addlinespace' '\n' ...
    '\\hline \\hline' '\n' ...
    '\\end{tabular}'  '\n' ...
    '\\caption*{ \\footnotesize Notes: decomposition of the effect of lowering the temporary contract hiring restriction parameter' '\n' ...
    '$\\phi_0$ on steady-state equilibrium unemployment flows. ``France'''' refers to the model economy with French-type institutions, and '  '\n' ...
    '``Spain'''' refers to Spanish-type institutions. The table reports changes in equilibrium outcomes implied by a decline in the parameter'  '\n' ...
    '$\\phi_0$ from 0.54 (France) and 0.71 (Spain) to 0. The second and fourth columns show changes in the percentage-point value of components in '  '\n' ...
    'equations \\ref{UE_change} and \\ref{EU_change}; the third and fifth columns show the percentage-point contribution of these components to '  '\n' ...
    'the change in the unemployment rate. }'  '\n' ...
    '\\end{table}' ];

% Export table
% open (or create) tex file
fileID = fopen([table_path, 'flow_decomposition.tex'], 'w');
% write table
fprintf(fileID, tab);
% close file
fclose(fileID);

% cleanup
clearvars -except table_path

disp('Table 6 exported')

%% Table 7: Welfare

% load
load('workspaces\counterfactuals\US_F.mat', 'p', 'agg_stat', 'eql', 'F')

bbeta = p.pval(p.ind.beta);

% initialize vectors for outcomes
L = length(F);

% LT equivalent consumption
cu = zeros(L,1);
cw = zeros(L,1);
c_new = zeros(L,1);

% average profits
pi = zeros(L,1);

% matching profits
pi_new = zeros(L,1);

% aggregate firing and hiring costs
FC = zeros(L,1);
HC = zeros(L,1);

% employment, productivity, output, wage
E = zeros(L,1);
Ymn = zeros(L,1);
Ytot = zeros(L,1);
Wmn = zeros(L,1);
Wnew = zeros(L,1);

% welfare (consumption equivalent net of firing and hiring costs)
% and net output
Wagg = zeros(L,1);
Yagg = zeros(L,1);

% stack out required parameters
ddelta = p.pval(p.ind.delta);   % separation
kkappa = p.pval(p.ind.kappa);   % vacancy posting costs
bb = p.pval(p.ind.b);           % unemployment benefits

% fill in outcome vectors
for j = 1:L

    % equiv. consumption
    cu(j) = agg_stat{j}.cU;
    cw(j) = agg_stat{j}.U*agg_stat{j}.cU + (1-agg_stat{j}.U)*agg_stat{j}.cV;
    c_new(j) = agg_stat{j}.U*agg_stat{j}.cU + (1-agg_stat{j}.U)*agg_stat{j}.cV;

    % average profits
    pi(j) = agg_stat{j}.cJ;

    % aggregate hiring/firing costs
    HC(j) = (agg_stat{j}.theta*agg_stat{j}.U)*kkappa;
    FC(j) = (1-agg_stat{j}.U)*(1-agg_stat{j}.T)*(agg_stat{j}.PU-ddelta)*F(j);

    % employment, output, prdty
    E(j) = 1-agg_stat{j}.U;
    Ytot(j) = agg_stat{j}.Ytot;
    Ymn(j) = agg_stat{j}.Ymn;

    % wage
    Wmn(j) = agg_stat{j}.logwmn;

    % welfare: LT consumption + profit annuity net of firms' aggregate costs
    Wagg(j) = cw(j) + (1-agg_stat{j}.U)*pi(j) - HC(j) - FC(j);

    % alternative welfare: total income net of aggregate costs
    Yagg(j) = agg_stat{j}.U*bb + agg_stat{j}.Ytot - HC(j) - FC(j);

end

% firing costs relative to total output
FC = FC./Ytot;
HC = HC./Ytot;

% outcomes relative to baseline
cu = (cu - cu(1))/cu(1);
cw = (cw - cw(1))/cw(1);
pi = (pi - pi(1));
HC = (HC - HC(1));
FC = (FC - FC(1));
Yagg = (Yagg - Yagg(1))/Yagg(1); 
Wagg = Wagg - Wagg(1);
Ytot = (Ytot - Ytot(1))/Ytot(1);
Ymn = (Ymn - Ymn(1))/Ymn(1);
Wmn = (Wmn - Wmn(1));
E = (E - E(1))/E(1);
  
% string format
s = @(v) ...F
    [ num2str(v(3)*100, '%.2f'),' & ',num2str(v(4)*100, '%.2f'),' & ', ...
    num2str(v(5)*100, '%.2f'),' & ', num2str(v(6)*100, '%.2f') ];

% table
tab = ['\\begin{table}[!h]' '\n' ...
    '\\centering'           '\n' ...
    '\\captionof{table}{Baseline model outcomes under alternative firing-cost regimes (relative to $F=0$, in \\%%)}'   '\n' ...
    '\\label{tab:welfare} '          '\n' ...
    '\\begin{tabular}{l c c c c}' 		        '\n' ...
    '\\hline \\hline'                           '\n' ...
    '\\addlinespace'                            '\n' ...
    '\\textit{Firing costs/benchmark wage}  & \\textit{one month} & \\textit{three months}  & \\textit{six months}  & \\textit{twelve months}' '\\\\' '\n' ...
    '\\addlinespace'                            '\n' ...
    '\\addlinespace' '\n' ...
    'Total output (\\%% change) & '             s(Ytot)        '\\\\'    '\n' ...
    'Employment (\\%% change) &'                 s(E)           '\\\\'    '\n' ...
    'Labor productivity (\\%% change) & '       s(Ymn)         '\\\\'    '\n' ...
    '\\addlinespace' '\n' ...
    '\\textit{Profits and labor costs (p.p.\\ change)}' '\\\\'  '\n' ...
    '\\hspace{5pt} Average profits $\\overline{\\pi}$ & '     s(pi)         '\\\\' '\n' ...
    '\\hspace{5pt} Hiring costs/output & '                    s(HC)         '\\\\' '\n' ...
    '\\hspace{5pt} Average log wage & '                       s(Wmn)        '\\\\' '\n' ...
    '\\addlinespace' '\n' ...
    '\\textit{Welfare (\\%% change)}' '\\\\'  '\n' ...
    '\\hspace{5pt} Unemployed workers, ${c}_u$ \\hspace{100pt}& ' s(cu)   '\\\\' '\n' ...
    '\\hspace{5pt} All workers, $\\overline{c}_e$ & '             s(cw)   '\\\\'              '\n' ...
    '\\hspace{5pt} Social welfare, $\\mathcal{W}(F)$& '           s(Wagg) '\\\\'    '\n' ...
    '\\addlinespace'    '\n' ...
    '\\addlinespace'    '\n' ...
    '\\hline \\hline'   '\n' ...
    '\\end{tabular}'    '\n' ...
    '\\caption*{ \\footnotesize Notes: equilibrium outcomes associated with an increase in firing costs $F$' ...
    ' in the baseline economy, assuming no regulation on temporary jobs: $\\phi_0 = \\phi = 0$. Profits are measured in lifetime ' ...
    'annuity income terms and worker''s welfare are measured in lifetime certainty equivalent consumption terms. The aggregate ' ...
    'social welfare statistics is given by \\eqref{welfare}; aggregate income is defined as total output added to aggregate non-work ' ...
    'income, net of aggregate hiring costs of firms. }' '\n' ...
    '\\end{table}' ];

% Export table
% open (or create) tex file
fileID = fopen([table_path, 'welfare.tex'], 'w');
% write table
fprintf(fileID, tab);
% close file
fclose(fileID);

% cleanup
%clearvars -except table_path
disp('Table 7 exported')





%% APPENDIX

%% Table A1: differences between France and the U.S.

load('workspaces/counterfactuals/France_vs_US.mat');

% parameter difference
load('workspaces\Baseline.mat', 'p')
ii = [p.ind.A; p.ind.delta; p.ind.b; p.ind.F];
pval0 = p.pval(ii);
clearvars p

load('workspaces\France.mat', 'p')
pval1 = p.pval(ii);
pdiff = pval1 - pval0;
clearvars p

% functions to convert number in string format
s = @(nb, fmt) num2str(nb, fmt);

% table: France vs. Spain
France_vs_US = [ '\\begin{table}[!h]' '\n' ...
    '\\centering' '\n' ...
    '\\captionof{table}{Decomposition of differences between the French and U.S.\\ economies}'   '\n' ...
    '\\label{tab:France_vs_Spain} '          '\n' ...
    '\\begin{tabular}{l c c c c c}' 		 '\n' ...
    '\\hline \\hline' '\n' ...
    '\\addlinespace'  '\n' ...
    ' \\hspace{100pt} &    Total  & $A$ & $\\delta$ & $b$ & $F$   \\\\'      '\n' ...
    ' \\addlinespace'  '\n' ...
    ' Reference param.\\ value (U.S.) & &' s(pval0(1), '%.2f') '&' s(pval0(2), '%.4f')  '&' s(pval0(3), '%.2f') '&' s(pval0(4), '%.0f')   '\\\\'      '\n' ...
    ' Comparison param.\\ (France) & &'  s(pval1(1), '%.2f') '&' s(pval1(2), '%.4f')  '&' s(pval1(3), '%.2f') '&' s(pval1(4), '%.2f')   '\\\\'      '\n' ...
    ' Param.\\ differential &  &' s(pdiff(1), '%.2f') '&' s(pdiff(2), '%.4f')  '&' s(pdiff(3), '%.2f') '&' s(pdiff(4), '%.2f')   '\\\\'      '\n' ...
    ' \\addlinespace'  '\n' ...
    ' \\multicolumn{6}{l}{\\textit{Labor-market stocks differential (percentage-point)}}  \\\\' '\n' ...
    ' \\hspace{5pt} Unem.\\ rate &'   s(U(1,1)*100, '%.2f')   '&' s(U(2,1)*100, '%.2f')    '&' s(U(3,1)*100, '%.2f')  '&' s(U(4,1)*100, '%.2f')   '&' s(U(5,1)*100, '%.2f')   '\\\\'  '\n' ...
    ' \\addlinespace'  '\n' ...
    ' \\textit{Transition rates (p.p.)}  &   &   &   &   &  \\\\' '\n' ...
    ' \\hspace{5pt} UE           &'   s(UE(1,1)*100, '%.2f') '&' s(UE(2,1)*100, '%.2f') '&' s(UE(3,1)*100, '%.2f')  '&' s(UE(4,1)*100, '%.2f') '&' s(UE(5,1)*100, '%.2f')   '\\\\'  '\n' ...
    ' \\hspace{5pt} EU           &'   s(EU(1,1)*100, '%.2f') '&' s(EU(2,1)*100, '%.2f') '&' s(EU(3,1)*100, '%.2f')  '&' s(EU(4,1)*100, '%.2f') '&' s(EU(5,1)*100, '%.2f')   '\\\\'  '\n' ...
    '\\addlinespace'  '\n' ...
    '\\hline \\hline' '\n' ...
    '\\end{tabular}'  '\n' ...
    '\\caption*{ \\footnotesize Notes: Decomposition of differences in aggregate equilibrium outcomes ' ...
    'between the ``French'''' and U.S.\\ baseline model economy. The first column reports the total difference ' ...
    'in equilibrium outcomes between these two economies. The subsequent columns report the difference in ' ...
    'outcomes after counterfactually imposing a given parameter in the U.S.\\ baseline to be equal to ' ...
    'its counterpart in the French economy, i.e., a measure of the total difference in equilibrium outcomes ' ...
    'attributable to this parameter. The parameters under focus are the matching efficiency $A$ , the ' ...
    'exogenous probability of separation $\\delta$, the non-work income level $b$, and firing costs $F$. ' ...
    'The first row indicates parameter values in the baseline (U.S.) reference economy. ' ...
    'The second row indicates the difference in parameter values between the baseline and French model economies. }' '\n' ...
    '\\end{table} ' ];

% Export table
% open (or create) tex file
fileID = fopen([table_path, 'France_vs_US.tex'], 'w');
% write table
fprintf(fileID, France_vs_US);
% close file
fclose(fileID);

% cleanup
clearvars -EXCEPT table_path
disp('Table A1 exported')


%% Table A2: differences between France and Spain

load('workspaces/counterfactuals/France_vs_Spain.mat');

% parameter difference
load('workspaces\France.mat', 'p')
ii = [p.ind.delta; p.ind.b; p.ind.F; p.ind.phi];
pval0 = p.pval(ii);
clearvars p

load('workspaces\Spain.mat', 'p')
pval1 = p.pval(ii);
pdiff = pval1 - pval0;
clearvars p

% functions to convert number in string format
s = @(nb, fmt) num2str(nb, fmt);

% table: France vs. Spain
France_vs_Spain = [ '\\begin{table}[!h]' '\n' ...
    '\\centering' '\n' ...
    '\\captionof{table}{Decomposition of differences between the French and Spanish economies}'   '\n' ...
    '\\label{tab:France_vs_Spain} '          '\n' ...
    '\\begin{tabular}{l c c c c c}' 		 '\n' ...
    '\\hline \\hline' '\n' ...
    '\\addlinespace'  '\n' ...
    ' \\hspace{100pt} &    Total  & $\\delta$ & $b$ & $F$ & $\\phi$   \\\\'      '\n' ...
    ' \\addlinespace'  '\n' ...
    ' Reference param.\\ value (French) & &' s(pval0(1), '%.4f') '&' s(pval0(2), '%.2f')  '&' s(pval0(3), '%.2f') '&' s(pval0(4), '%.3f')   '\\\\'      '\n' ...
    ' Comparison param.\\ (Spanish) & &'  s(pval1(1), '%.4f') '&' s(pval1(2), '%.2f')  '&' s(pval1(3), '%.2f') '&' s(pval1(4), '%.3f')   '\\\\'      '\n' ...
    ' Param.\\ differential &  &' s(pdiff(1), '%.4f') '&' s(pdiff(2), '%.2f')  '&' s(pdiff(3), '%.2f') '&' s(pdiff(4), '%.3f')   '\\\\'      '\n' ...
    ' \\addlinespace'  '\n' ...
    ' \\multicolumn{6}{l}{\\textit{Labor-market stocks differential (percentage-point)}}  \\\\' '\n' ...
    ' \\hspace{5pt} Unem.\\ rate &'   s(U(1,1)*100, '%.2f')   '&' s(U(2,1)*100, '%.2f')    '&' s(U(3,1)*100, '%.2f')  '&' s(U(4,1)*100, '%.2f')   '&' s(U(5,1)*100, '%.2f')   '\\\\'  '\n' ...
    ' \\hspace{5pt} TC emp.\\ share           &'   s(T(1,1)*100, '%.2f')  '&' s(T(2,1)*100, '%.2f')  '&' s(T(3,1)*100, '%.2f')   '&' s(T(4,1)*100, '%.2f')  '&' s(T(5,1)*100, '%.2f')    '\\\\'  '\n' ...
    ' \\addlinespace'  '\n' ...
    ' \\textit{Transition rates (p.p.)}  &   &   &   &   &  \\\\' '\n' ...
    ' \\hspace{5pt} UE           &'   s(UE(1,1)*100, '%.2f') '&' s(UE(2,1)*100, '%.2f') '&' s(UE(3,1)*100, '%.2f')  '&' s(UE(4,1)*100, '%.2f') '&' s(UE(5,1)*100, '%.2f')   '\\\\'  '\n' ...
    ' \\hspace{5pt} EU           &'   s(EU(1,1)*100, '%.2f') '&' s(EU(2,1)*100, '%.2f') '&' s(EU(3,1)*100, '%.2f')  '&' s(EU(4,1)*100, '%.2f') '&' s(EU(5,1)*100, '%.2f')   '\\\\'  '\n' ...
    '\\addlinespace'  '\n' ...
    '\\hline \\hline' '\n' ...
    '\\end{tabular}'  '\n' ...
    '\\caption*{ \\footnotesize Notes: Decomposition of differences in aggregate equilibrium outcomes ' ...
    'between the ``French'''' and ``Spanish'''' model economy. The first column reports the total difference ' ...
    'in equilibrium outcomes between these two economies. The subsequent columns report the difference in ' ...
    'outcomes after counterfactually imposing a given parameter in the French model economy to be equal to ' ...
    'its counterpart in the Spanish economy, i.e., a measure of the total difference in equilibrium outcomes ' ...
    'attributable to this parameter. The parameters under focus are the exogenous probability of separation ' ...
    '$\\delta$, the non-work income level $b$, firing costs $F$, and the regulatory restriction on the duration ' ...
    'of temporary contracts, $\\phi$. The first row indicates parameter values in the French reference economy. ' ...
    'The second row indicates the difference in parameter values between the Spanish and French model economies. }' '\n' ...
    '\\end{table} ' ];

% Export table
% open (or create) tex file
fileID = fopen([table_path, 'France_vs_Spain.tex'], 'w');
% write table
fprintf(fileID, France_vs_Spain);
% close file
fclose(fileID);

% cleanup
clearvars -EXCEPT table_path
disp('Table A2 exported')
