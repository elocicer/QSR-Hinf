function [Lambda,soltime] = UpdateLambda1(A,B,P,L,eps)
% Ethan LoCicero
% 2022
% Solve SDP to initialize or update weight matrix Lambda1 between iterations

%% dimensions
sL2 = size(L,2);
Lambda = sdpvar( sL2,sL2 );

%% objective
invL = inv(L); invL = .5.*(invL+invL');
I = eye(sL2);
PL = P*L;
invLA = invL*A;
invLB = invL*B;
gamma = trace(PL*Lambda*PL');
chi = trace(invLA'*Lambda*invLA) + trace(invLB'*Lambda*invLB);
eta = trace(PL*(I-Lambda)*PL') + trace(invLA'*(I-Lambda)*invLA) + trace(invLB'*(I-Lambda)*invLB) - trace(invLA'*(I-Lambda)*PL') - trace(PL*(I-Lambda)*invLA);
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
