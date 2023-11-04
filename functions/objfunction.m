function obj = objfunction(x, param_structure)

%% Objective function for calibrating the model

% update the vector of parameters with the calibrated values
indexes = param_structure.cparam.ind;
param_structure.pval( indexes ) = x;

% do not compute equilibrium or welfare
param_structure.equilibrium = "partial";
param_structure.compute_welfare = false;

% Compute equilibrium evaluated at calibrated param values
[~,~,agg] = compute_equilibrium(param_structure);

% targets       
b             = .4;               % Shimer  (2005)
U             = .0561;            % U rate in the US, 1996-2010
EU            = .02;              % EU rate in the US
EU_tn         = 3.17;             % separation rates ratio

% vecotr
vector = abs( [ b - agg.b; U - agg.U; EU - agg.EU; EU_tn - agg.EU_tn ] );
vector = vector./[ b; U; EU; EU_tn ];

obj = mean( vector );

% nan
if isnan(obj)
   obj = inf; 
end

end