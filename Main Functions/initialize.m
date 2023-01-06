% initialize
% Ethan LoCicero
% July 13, 2023

% This script loads a stable LTI plant from plantdata.mat and generates an
% initial feasible point for Algorithm 1. The data is saved in
% initialization.mat. The method for initialization is to identify the L2
% gain of the plant, then find a controller that satisfies the small gain
% theorem. Thi contorller must result in a finite H-infinity norm and must
% have QSR properties satisfying the Dissipativity Theorem.


%% load plant
load('plantdata.mat')

%% set dimensions
sL1       = ctrlstate;
sL2       = in;
sL4       = ctrlstate + in;
sLAM1     = ctrlstate;
sLAM2     = in;
sLAM4     = ctrlstate + in;

%% find feasible initial controller and QSR bounds for stability
% find gain of stable plant
[plant_gamma,~] = findHinf(A,B2,C2,D22);
% satisfy small gain theorem with 5% error
ctrl_gamma = .95/plant_gamma;
% make random stable controller
Ac = rand(ctrlstate,ctrlstate); Ac = Ac-(1+rand)*max(real(eig(Ac))).*eye(size(Ac));
Bc = rand(ctrlstate,out);
Cc = rand(in,ctrlstate);
Dc = rand(in,out);
% make controller have appropriate gain
[Cc,Dc] = HinfCD(Ac,Bc,Cc,Dc,ctrl_gamma);
% construct closed-loop
Acl = [A-B2*Dc*C2 , -B2*Cc ; Bc*C2 , Ac];
Bcl = [B1-B2*Dc*D21; Bc*D21];
Ccl = [C1-D12*Dc*C2 , -D12*Cc];
Dcl = -D12*Dc*D21;
% find closed-loop Hinf gain
[gamma,X] = findHinf(Acl,Bcl,Ccl,Dcl);
% find plant and controller QSR bounds
[Qp,Sp,Rp,Qc,Sc,Rc,Pp,P] = findQSRboth(A,B2,C2,D22,Ac,Bc,Cc,Dc);
% split Q0 into positive and negative semidefinite portions
[V,D] = eig(Qc);
Dp = D; Dp(Dp<=0)=0;
Dm = D; Dm(Dm>=0)=0;
Qcp = V*Dp*V'; Qcp = .5.*(Qcp+Qcp');
Qcm  = V*Dm*V'; Qcm = .5.*(Qcm+Qcm');
% make Q0p, Q0m definite so that Q0m is invertable
eps = 10^(-1); % changed from 10^-4 to reduce numerical errors
Qcp = Qcp + eps.*eye(in);
Qcm = Qcm - eps.*eye(in);
% double check QSR bounds
Qcnew = Qcp+Qcm;

%% Initialize L's and Lambda's
eps = .01;
[L1,~] = UpdateL1(Ac,Bc,P,eps);
[L2,~] = UpdateL2(Cc,Dc,Qcp,eps);
[L4,~] = UpdateL4(B2,C2,D21,Ac,Bc,Cc,Dc,X,eps);
[LAM1,~] = UpdateLambda1(Ac,Bc,P,L1,eps);
[LAM2,~] = UpdateLambda2(Cc,Dc,Qcp,L2,eps);
[LAM4,~] = UpdateLambda4(B2,C2,D21,Ac,Bc,Cc,Dc,X,L4,eps);

%% save data
save('initialization.mat')
