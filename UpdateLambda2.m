function [Lambda,soltime] = UpdateLambda2(C,D,Qp,L,eps)
% Ethan LoCicero
% 2022
% Solve SDP to initialize or update weight matrix Lambda2 between iterations

%% dimensions
sL2 = size(L,2);
Lambda = sdpvar( sL2,sL2 );

%% objective
invL = inv(L); invL = .5.*(invL+invL');
I = eye(sL2);
CtL = C'*L;
DtL = D'*L;
invLQp = invL*Qp;
invLQpC = invLQp*C;
invLQpD = invLQp*D;
gamma = .25*( trace(CtL*Lambda*CtL') + trace(DtL*Lambda*DtL') );
chi = trace(invLQpC'*Lambda*invLQpC) + trace(invLQpD'*Lambda*invLQpD);
eta = .25*( trace(CtL*(I-Lambda)*CtL') + trace(DtL*(I-Lambda)*DtL') ) + .5*trace( invLQpC'*(I-Lambda)*CtL' + CtL*(I-Lambda)*invLQpC ) + .5*trace( invLQpD'*(I-Lambda)*DtL' + DtL*(I-Lambda)*invLQpD );
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