%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Risk-sharing in a dual market
% Créchet (2020)
% Matlab script file
% file name: "main.m"
% created: 28-10-2020
% last updated: Oct 2023
%
% Description: main script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
clear, clc
clear global

% specify directory path
host = getenv('COMPUTERNAME');
if strcmp(host,'DESKTOP-LBOSF73') 
    path = 'D:\Dropbox\Research\Projects\Contracts\Matlab';
elseif strcmp(host,'DESKTOP-OV023AQ')
    path = 'C:\Users\jcrechet\Dropbox\Research\Projects\Contracts\Matlab';
end

% set path
cd(path);
                                                                                                      
% folder with functions
addpath('functions');

% parameters
p = setparameters;
clearvars host path


%% Evaluate model at preset parameter values

p.pval(p.ind.F) = 0;
[eql, sim, agg_stat] = compute_equilibrium(p);

% save outcomes
save('workspaces/Preset.mat','eql','sim','agg_stat', 'p');      

% timer
disp('Computation time (in sec.):')
disp( eql.time ) 

% show aggregate statistics
disp('Simulated stats:')
disp( agg_stat )


%% Calibration
run('calibration')


%% Baseline calibrated model

% load calibrated model
load('workspaces\Baseline.mat')

disp('Parameter values:')
disp( p.ptable )

% run model
eql = compute_equilibrium(p);

% outcomes
disp('Computation time (in sec.):')
disp( eql.time )

disp('Simulated stats:')
disp( agg_stat )

% save
save('workspaces\Baseline.mat', 'p', 'eql', 'sim', 'agg_stat')


%% Experiments
run('experiments')
run('experiments_supplementary')


%% Tables
run('tables')


%% Figures
run('figures')
