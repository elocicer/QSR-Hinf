% plant
% Ethan LoCicero
% 01/06/2023

% This script randomly generates a plant and saves the data in
% plantdata.mat

clear

%% make random stable plant
state     = randi(10); % plant state dimension
ctrlstate = randi(10); % desired controller state dimension
in        = randi(10); % input dimension
out       = randi(10); % measured output dimension
dist      = randi(10); % disturbance signal dimension
signal    = randi(10); % controlled output dimension
A = rand(state,state); A = A-(1+rand)*max(real(eig(A))).*eye(size(A));
B2 = rand(state,in);
C2 = rand(out,state);
D22 = zeros(out,in); % must be zero for the Hinf LMI
D21 = rand(out,dist);
B1 = rand(state,dist);
C1 = rand(signal,state);
D11 = zeros(signal,dist);
D12 = rand(signal,in);


%% save data
save('plantdata.mat')
