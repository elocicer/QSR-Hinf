function [L,soltime] = UpdateL2(C,D,Qp,eps)
% Ethan LoCicero
% 2022
% Solve SDP to initialize or update weight matrix L2 between iterations

%% dimensions
[sC1,sC2] = size(C);
sD2 = size(D,2);
sZ = sC2+sD2;
sblk = sZ+sC1;

%% variables
Z = sdpvar( sZ,sZ );
W = sdpvar( sC1,sC1 );

%% objective
objective = trace(Z);

%% constraints
CtW = C'*W;
DtW = D'*W;
block = [CtW*C , CtW*D ; DtW*C , DtW*D]; block = .5.*(block+block');
edge = [Qp*C , Qp*D];
constraint = [Z - .25.*block , edge' ; edge , W ] >= eps*eye(sblk);

%% solve
options   = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
tic
sol       = optimize( constraint , objective , options );
soltime = toc;
W = value(W);
L = sqrtm(W);

end