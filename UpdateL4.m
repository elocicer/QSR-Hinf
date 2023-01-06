function [L,soltime] = UpdateL4(B2,C2,D2,Ac,Bc,Cc,Dc,X,eps)
% Ethan LoCicero
% 2022
% Solve SDP to initialize or update weight matrix L4 between iterations

%% dimensions
[sA,sB2_2] = size(B2);
sD2_2 = size(D2,2);
sBc_1 = size(Bc,1);
sDc_1 = size(Dc,1);
X1 = X(1:sA,1:sA);
X2 = X(1:sA,sA+1:end);
X3 = X(sA+1:end,sA+1:end);
[sX2_1,sX2_2] = size(X2);
sX3_1 = size(X3,1);
sblock = sX2_1+sX3_1+sD2_2+sBc_1+sDc_1;
sZ = sX2_1 + sX3_1 + sD2_2;
sW = sX2_2+sB2_2;

%% variables
Z = sdpvar( sZ,sZ );  
W = sdpvar( sW,sW );

%% constraints
blk1 = [X2 , -X1*B2 ; X3 , -X2'*B2];   
sblk1_1 = size(blk1,1);
block1 = [blk1*W*blk1' , zeros(sblk1_1, sD2_2) ; zeros(sD2_2, sblk1_1+sD2_2) ];
block2 = [Bc*C2, Ac, Bc*D2; Dc*C2, Cc, Dc*D2];     
constraint = [Z-block1 , block2' ; block2 , W] >= eps*eye(sblock);

%% objective
objective = trace(Z);

%% solve
options   = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
tic
sol       = optimize( constraint , objective , options );
soltime = toc;
W = value(W);
L = sqrtm(W);
if ~isreal(L)
    disp('Tolerance error in L4 update') % this happens sometimes with sqrtm
end
 
end