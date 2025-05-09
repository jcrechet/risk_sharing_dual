%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Risk-sharing in a dual market
% Créchet (2020)
% Matlab script file
% file name: "main.m"
% created: 28-10-2020
% last updated: 2024
%
% Description: main script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Preamble
clear, clc
clear global

% specify directory and result path
host = getenv('COMPUTERNAME');

% SPECIFY COMPUTER NAME AND WORKING DIRECTORY PATH HERE

%if strcmp(host,'[MYCOMPNAME]')
    
%    path = '[MYPATH]';
    
%end

clearvars host

% set working dir path
cd(path);
                                                                                                      
% folder with functions
addpath('functions');

% set model and algorithms parameters
p = setparameters;

clearvars host path

%% 1.1. Evaluate model at preset parameter values

% start parallel pool using local profile
parpool('local')

% solve model's equil at user's preset parameter values
[eql, sim, agg_stat] = compute_equilibrium(p);

% save outcomes
save('workspaces/Preset.mat','eql','sim','agg_stat', 'p');

% show timer
disp('Computation time (in sec.):')
disp( eql.time ) 

% show aggregate statistics
disp('Simulated stats:')
disp( agg_stat )


%% 1.2. Run calibration algo
run('calibration')


%% 1.3.Solve baseline calibrated model equilibrium

% load calibrated model's outcome (from calibration algo)
load('workspaces\Baseline.mat')

% show paramater values
disp('Parameter values:')
disp( p.ptable )

% solve equilibrium at calibrated parameters
eql = compute_equilibrium(p);

% show timer and equil outcomes
disp('Computation time (in sec.):')
disp( eql.time )

disp('Simulated stats:')
disp( agg_stat )

% save outcomes
save('workspaces\Baseline.mat', 'p', 'eql', 'sim', 'agg_stat')


%% 2. Experiments
run('experiments')
run('experiments_supplementary')


%% Tables
run('tables')


%% Figures
run('figures')

