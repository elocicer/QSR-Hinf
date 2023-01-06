function [Lambda,soltime] = UpdateLambda4(B2,C2,D2,Ac,Bc,Cc,Dc,X,L,eps)
% Ethan LoCicero
% 2022
% Solve SDP to initialize or update weight matrix Lambda4 between iterations

%% dimensions
sL2 = size(L,2);
sC2_2 = size(C2,2);
sB2_2 = size(B2,2);
sAc_2 = size(Ac,2);
sD2_2 = size(D2,2);
sA   = size(B2,1);

%% vairables
Lambda = sdpvar( sL2,sL2 );

%% objective
X1 = X(1:sA,1:sA);
X2 = X(1:sA,sA+1:end);
X3 = X(sA+1:end,sA+1:end);
[sX2_1,sX2_2] = size(X2);
sX3_1 = size(X3,1);
invL = inv(L); invL = .5.*(invL+invL');
I = eye(sL2);
blk1 = [X2 , -X1*B2 ; X3 , -X2'*B2];
blk1_padded = [blk1 ; zeros(sC2_2+sAc_2+sD2_2-sX2_1-sX3_1,sX2_2+sB2_2)];
blk2 = [C2'*Bc' , C2'*Dc' ; Ac' , Cc' ; D2'*Bc' , D2'*Dc'];
gamma = trace(blk1*L*Lambda*L'*blk1');
chi = trace(blk2*invL'*Lambda*invL*blk2');
eta = trace(blk1*L*(I-Lambda)*L'*blk1') + trace(blk2*invL'*(I-Lambda)*invL*blk2') - 2.*trace(blk2*invL'*(I-Lambda)*L'*blk1_padded');
objective = gamma + chi + eta;

%% constraints
c1 = Lambda >= eps.*eye(sL2);
c2 = Lambda <= (1-eps).*eye(sL2);
constraints = [c1, c2];

%% solve
options   = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
tic
sol       = optimize( constraints , objective , options );
soltime = toc;
Lambda = value(Lambda);

end
