% main
% Ethan LoCicero
% 01/06/2023

% Adds subfolder paths, then runs three major scripts to
% 1. Generate a randomized stable LTI plant
% 2. Initialize a feasible controller
% 3. Apply iterative convex overbounding to reduce the H-infinity norm
%    while satisfying dissipativity constraints

clear
close all
clc

%% Add Paths
addpath('./Main Functions')
addpath('./Weight Update Functions')
addpath('./Auxiliary Functions')

%% Run Scripts
disp('Generating Plant')
plant
disp('Initializing Controller')
initialize
disp('Applying Iterative Convex Overbounding')
ICO

%% Display Results
disp(['Initial H-inf norm:  ',num2str(gamma(1))])
disp(['Final H-inf norm:    ',num2str(gamma(k))])
disp(['Percent Improvement: ',num2str(percent_improvement)])
