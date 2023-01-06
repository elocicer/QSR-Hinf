function [L,soltime] = UpdateL1(A,B,P,eps)
% Ethan LoCicero
% 2022
% Solve SDP to initialize or update weight matrix L1 between iterations

%% dimensions
sA = size(A,1);
sB2 = size(B,2);
sblk = 2*sA + sB2;

%% variables
Z1 = sdpvar( sA,sA );
Z2 = sdpvar( sA,sB2,'full' );
Z3 = sdpvar( sB2,sB2 );
W = sdpvar( sA,sA );
Z = [Z1 , Z2 ; Z2' , Z3];

%% objective
objective = trace(Z);

%% constraints
constraint = [Z1-P*W*P , Z2 , A' ; Z2' , Z3 , B' ; A , B , W ] >= eps*eye(sblk);

%% solve
options   = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
tic
sol       = optimize( constraint , objective , options );
soltime = toc;
W = value(W);
L = sqrtm(W);

end