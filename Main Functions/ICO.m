% main
% Ethan LoCicero
% July 13, 2022

% This script runs the iterations of iterative convex overbounding
% (Algorithm 1).

close all
clear
clc

%% load initial feasible point
load('initialization.mat')

%% set hyperparameters
eps_main = 10^(-4);             % definiteness tolerance in main iteration LMIs
eps_update = .01;               % definiteness tolerance in weight update LMIs
check = true;                   % double check validity of results?
Nmax = 1001;                    % maximum number of iterations
epsilon = 10^-2;                % iteration convergence tolerance
gamma = [gamma, nan(1,Nmax-1)]; % initialize gamma vector
soltime = nan(Nmax);
soltimeL1 = nan(Nmax); soltimeL2 = nan(Nmax); soltimeL4 = nan(Nmax); 
soltimeLAM1 = nan(Nmax); soltimeLAM2 = nan(Nmax); soltimeLAM4 = nan(Nmax);
dJ = inf;                       % initialize chenge in objective
k = 1;                          % initialize iteration counter

%% iterative convex overbounding
while (dJ > epsilon) && (k < Nmax)
    [Ac, Bc, Cc, Dc, P, X, gamma(k+1), soltime(k+1), change] = HinfQSRsimple(A,B2,C2,D21,B1,C1,D12,Qcp,Qcm,Sc,Rc,P,X,Ac,Bc,Cc,Dc,L1,L2,L4,LAM1,LAM2,LAM4,eps_main);
    [L1,soltimeL1(k+1)] = UpdateL1(change.dAc,change.dBc,change.dP,eps_update);
    [L2,soltimeL2(k+1)] = UpdateL2(change.dCc,change.dDc,Qcp,eps_update);
    [L4,soltimeL4(k+1)] = UpdateL4(B2,C2,D21,change.dAc,change.dBc,change.dCc,change.dDc,change.dX,eps_update);
    [LAM1,soltimeLAM1(k+1)] = UpdateLambda1(change.dAc,change.dBc,change.dP,L1,eps_update);
    [LAM2,soltimeLAM2(k+1)] = UpdateLambda2(change.dCc,change.dDc,Qcp,L2,eps_update);
    [LAM4,soltimeLAM4(k+1)] = UpdateLambda4(B2,C2,D21,change.dAc,change.dBc,change.dCc,change.dDc,change.dX,L4,eps_update);
    dJ = gamma(k) - gamma(k+1);
    disp([num2str(k),' | ',num2str(dJ)])
    k = k+1;
end
Time_calc = nanmean(soltime) + nanmean(soltimeL1) + nanmean(soltimeL2) + nanmean(soltimeL4) + nanmean(soltimeLAM1) + nanmean(soltimeLAM2) + nanmean(soltimeLAM4);

%% double check
if check == true
    plantIsQSR = checkQSR(A,B2,C2,D22,Qp,Sp,Rp)
    controllerIsQSR = checkQSR(Ac,Bc,Cc,Dc,Qcp+Qcm,Sc,Rc)
    isStable = checkQSRtheorem(Qp,Sp,Rp,Qcp+Qcm,Sc,Rc)
    Acl = [A-B2*Dc*C2 , -B2*Cc ; Bc*C2 , Ac];
    Bcl = [B1-B2*Dc*D21; Bc*D21];
    Ccl = [C1-D12*Dc*C2 , -D12*Cc];
    Dcl = -D12*Dc*D21;
    isGammaBounded_new = checkHinf(Acl,Bcl,Ccl,Dcl,gamma(k))
    format bank
    percent_improvement = round((gamma(1)-gamma(k))/gamma(1)*10000)/100
    format long
end

%% save data
save('data.mat')
