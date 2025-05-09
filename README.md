Replication package for:

Jonathan Créchet,
A model of risk sharing in a dual labor market,
Journal of Monetary Economics,
Volume 147,
2024,
103591,
ISSN 0304-3932,
https://doi.org/10.1016/j.jmoneco.2024.103591.
(https://www.sciencedirect.com/science/article/pii/S0304393224000448)


# Instructions

- Download repository content and drop in desired working directory;
- in MATLAB script 'main.m', specify the working directory path;
- run 'main.m';
- run latex file '_result/main.tex'.


# Specs

- Processor	Intel(R) Xeon(R) W-2265 CPU @ 3.50GHz   3.50 GHz
- Installed RAM	64.0 GB (63.7 GB usable)
- System type	64-bit operating system, x64-based processor
- Windows 11 Pro for Workstations
- MATLAB R2022a
- Toolboxes: Optimization, Statistics and Machine Learning Toolbox.


# Description of MATLAB files

Scripts:

- 'main': main script.
- 'calibration': run calibration algorithm, produces results for Tables 1 and 2.
- 'experiments' and 'flow_decomposition_TC': run quantitative experiments; results for Tables 4 to 7.
- 'experiments_supplementary': results for Tables A1 and A2 in supplementary paper's appendix.
- 'tables' and 'figures': produce latex output exported in folder '_results'.

Functions:

- 'setparameters': specify numerical settings and model's economic parameters.
- 'compute_equilibrum': solve model's equilibrium and compute moments.
- 'n_discretized': vector and pmf for discretized match quality distribution (to compute equilibrium).
- 'objfunction': compute distance between model and empirical moments (used in calibration).
- 'EU_decomposition', 'UE_decomposition', 'U_decomposition': conduct decomposition in Table 6.
