function [Ac,Bc,Cc,Dc,P,X,gamma,soltime,change] = HinfQSRsimple(A,B2,C2,D21,B1,C1,D12,Qp,Qm,S,R,P0,X0,Ac0,Bc0,Cc0,Dc0,L1,L2,L4,LAM1,LAM2,LAM4,eps)
% Ethan LoCicero
% 2022
% Semidefinite program to locally minimize the closed-loop H-infinity norm
% subject to an open loop dissipativity constraint on the controller.

%% Dimensions
sA            = size(A,1);
sAc           = size(Ac0,1);
[sBc1,sBc2]   = size(Bc0);
[sCc1,sCc2]   = size(Cc0);
[sDc1,sDc2]   = size(Dc0);
sP        = size(P0,1);
sX1   = sA;
sX2_1 = sA; 
sX2_2 = sAc;
sX3   = sAc;
X10 = X0(1:sX1,1:sX1);
X20 = X0(1:sX1,sX1+1:end);
X30 = X0(sX1+1:end,sX1+1:end);

%% Variables
dAc = sdpvar( sAc, sAc, 'full' );
dBc = sdpvar( sBc1,sBc2,'full' );
dCc = sdpvar( sCc1,sCc2,'full' );
dDc = sdpvar( sDc1,sDc2,'full' );
dP = sdpvar( sP,sP );
dX1 = sdpvar( sX1,sX1 );
dX2 = sdpvar( sX2_1,sX2_2,'full' );
dX3 = sdpvar( sX3,sX3 );
gamma = sdpvar;

%% Pre-multiply to eliminate redundancies for speed
invQm = inv(Qm); invQm = .5.*(invQm+invQm');
P0A0 = P0*Ac0;
P0B0 = P0*Bc0;
S0tC0 = S'*Cc0;
S0tD0 = S'*Dc0;
P0dA = P0*dAc;
P0dB = P0*dBc;
S0tdC = S'*dCc;
S0tdD = S'*dDc;
dPA0 = dP*Ac0;
dPB0 = dP*Bc0;
QC0 = Qp*Cc0;
QD0 = Qp*Dc0;
QdC = Qp*dCc;
QdD = Qp*dDc;
C0tQC0 = Cc0'*QC0; C0tQC0 = .5.*(C0tQC0 + C0tQC0');
C0tQD0 = Cc0'*QD0;
D0tQD0 = Dc0'*QD0; D0tQD0 = .5.*(D0tQD0 + D0tQD0');
C0tQdC = QC0'*dCc;
C0tQdD = QC0'*dDc;
D0tQdD = QD0'*dDc;
dCtQD0 = dCc'*QD0;
L1inv = inv(L1); L1inv = .5.*(L1inv + L1inv');
L2inv = inv(L2); L2inv = .5.*(L2inv + L2inv');
L1L1t = L1*L1';
L2L2t = L2*L2';
LAM1inv = inv(LAM1); LAM1inv = .5.*(LAM1inv + LAM1inv');
LAM2inv = inv(LAM2); LAM2inv = .5.*(LAM2inv + LAM2inv');
LAMI1inv = inv(eye(size(LAM1))-LAM1); LAMI1inv = .5.*(LAMI1inv + LAMI1inv');
LAMI2inv = inv(eye(size(LAM2))-LAM2); LAMI2inv = .5.*(LAMI2inv + LAMI2inv');

X1 = X10+dX1;
X2 = X20+dX2;
Bc0C2 = Bc0*C2;
Bc0D21 = Bc0*D21;
B2Dc0D21 = B2*Dc0*D21;
dBcD21 = dBc*D21;
B2dDcD21 = B2*dDc*D21;
B2Cc0 = B2*Cc0;
B2Dc0C2 = B2*Dc0*C2;
dBcC2 = dBc*C2;
B2tX20 = B2'*X20;
X10B2 = X10*B2;
dDcC2 = dDc*C2;
C2tDc0t = C2'*Dc0';
L4inv = inv(L4);
L4L4t = L4*L4';
LAM4inv = inv(LAM4); LAM4inv = .5.*(LAM4inv+LAM4inv');
LAMI4inv = inv(eye(size(LAM4))-LAM4);

%% QSR LMI
phi = [P0A0 + P0dA + dPA0    ,  P0B0 + P0dB + dPB0   ;...
      -S0tC0 - S0tdC , -S0tD0 - S0tdD];
xi_off = C0tQD0 + C0tQdD + dCtQD0;
xi = [C0tQC0 + C0tQdC + C0tQdC' , xi_off ;...
     xi_off' , R + D0tQD0 + D0tQdD + D0tQdD' ];

Pi11 = phi + phi' - xi;
Pi21 = -.5.*[dCc , dDc];
Pi31 = [QdC , QdD];
Pi41 = (-.5.*L2L2t)*[dCc , dDc] + [QdC , QdD];
Pi51 = [dP, zeros(sBc1,sBc2)];
Pi61 = [dAc , dBc];
Pi71 = [L1L1t*dP + dAc , dBc ];
Pi81 = [Cc0+dCc , Dc0+dDc];

Pi22 = -L2inv'*LAM2inv*L2inv; 
Pi33 = -L2*LAM2inv*L2'; 
Pi44 = -L2*LAMI2inv*L2';
Pi55 = -L1inv'*LAM1inv*L1inv;
Pi66 = -L1*LAM1inv*L1';
Pi77 = -L1*LAMI1inv*L1';
Pi88 = invQm;

Pi22 = .5.*(Pi22+Pi22');
Pi33 = .5.*(Pi33+Pi33');
Pi44 = .5.*(Pi44+Pi44');
Pi55 = .5.*(Pi55+Pi55');
Pi66 = .5.*(Pi66+Pi66');
Pi77 = .5.*(Pi77+Pi77');

sPi11_1 = size(Pi11,1);
sPi21_1 = size(Pi21,1);
sPi31_1 = size(Pi31,1);
sPi41_1 = size(Pi41,1);
sPi51_1 = size(Pi51,1);
sPi61_1 = size(Pi61,1);
sPi71_1 = size(Pi71,1);
sPi81_1 = size(Pi81,1);
emptyside = sPi21_1+sPi31_1+sPi41_1+sPi51_1+sPi61_1+sPi71_1+sPi81_1;

column1 = [zeros(sPi11_1,sPi11_1); Pi21; Pi31; Pi41; Pi51; Pi61; Pi71; Pi81;];
sc1 = size(column1,1);
PI1sub = [column1, zeros(sc1,emptyside)];
PI1diag = blkdiag(Pi11,Pi22,Pi33,Pi44,Pi55,Pi66,Pi77,Pi88);
PI1 = PI1diag + PI1sub + PI1sub';

%% Hinf LMI
X1A         = X1*A;
X2BcC2      = X20*Bc0C2 + X20*dBcC2 + dX2*Bc0C2;
X1B2DcC2    = X10*B2Dc0C2 + dX1*B2Dc0C2 + X10B2*dDcC2;
fi11_1      = X1A + X2BcC2 - X1B2DcC2 +(X1A + X2BcC2 - X1B2DcC2)';
AtX2        = A'*X2;
X2Ac        = X20*Ac0 + dX2*Ac0 + X20*dAc;
C2tBctX3    = Bc0C2'*X30 + dBcC2'*X30 + Bc0C2'*dX3;
X1B2Cc      = X10*B2Cc0 + dX1*B2Cc0 + X10B2*dCc;
C2tDctB2tX2 = C2tDc0t*B2tX20 + dDcC2'*B2tX20 + C2tDc0t*B2'*dX2;
fi11_2      = AtX2 + X2Ac + C2tBctX3 - X1B2Cc - C2tDctB2tX2;
X3Ac        = X30*Ac0 + dX3*Ac0 + X30*dAc;
X2tB2Cc     = X20'*B2Cc0 + dX2'*B2Cc0 + B2tX20'*dCc;
fi11_3      = X3Ac - X2tB2Cc + (X3Ac - X2tB2Cc)';
fi11        = [fi11_1, fi11_2; fi11_2', fi11_3];
X2BcD2    = X20*Bc0D21 + dX2*Bc0D21 + X20*dBcD21;
X1B2DcD2  = X10*B2Dc0D21 + dX1*B2Dc0D21 + X10*B2dDcD21; 
fi12_1    = X1*B1 + X2BcD2 - X1B2DcD2;
X3BcD2    = X30*Bc0D21 + dX3*Bc0D21 + X30*dBcD21;
X2tB2DcD2 = X20'*B2Dc0D21 + dX2'*B2Dc0D21 + X20'*B2dDcD21;
fi12_2    = X2'*B1 + X3BcD2 - X2tB2DcD2;
Fi12 = [fi12_1 ; fi12_2];
Fi13 = [C1' - C2'*(Dc0 + dDc)'*D12'; -(Cc0 + dCc)'*D12'];
Fi14 = [dX2 , -dX1*B2 ; dX3 , -dX2'*B2];
Fi15 = [C2'*dBc' , C2'*dDc' ; dAc' , dCc'];
Fi16 = Fi14*L4L4t + Fi15;

sFi11_2 = size(fi11,2);
sFi12_2 = size(Fi12,2);
sFi13_2 = size(Fi13,2);
sFi14_2 = size(Fi14,2);
sFi15_2 = size(Fi15,2);
sFi16_2 = size(Fi16,2);

Fi23 = -(D12*(Dc0 + dDc)*D21)';
Fi256 = [D21'*dBc' , D21'*dDc'];   

Fi22 = -gamma*eye(sFi12_2);
Fi33 = -gamma*eye(sFi13_2);
Fi44 = -L4inv'*LAM4inv*L4inv; Fi44 = .5.*(Fi44 + Fi44');
Fi55 = -L4*LAM4inv*L4'; Fi55 = .5.*(Fi55 + Fi55');
Fi66 = -L4*LAMI4inv*L4'; Fi66 = .5.*(Fi66 + Fi66');

toprow = [zeros(sFi11_2,sFi11_2), Fi12, Fi13, Fi14, Fi15, Fi16];
secondrow = [zeros(sFi12_2,sFi11_2+sFi12_2), Fi23, zeros(sFi12_2,sFi14_2), Fi256, Fi256];
emptylength = sFi13_2+sFi14_2+sFi15_2+sFi16_2;
fulllength = sFi11_2+sFi12_2+emptylength;

FI1sub = [toprow; secondrow; zeros(emptylength,fulllength)];
FI1diag = blkdiag(fi11,Fi22,Fi33,Fi44,Fi55,Fi66);
FI1 = FI1diag + FI1sub + FI1sub';

bigX = [X10 + dX1 , X20 + dX2 ; X20' + dX2' , X30 + dX3];

%% Constraints
c1a = PI1 <= 0;
c1b = P0+dP >= eps*eye(sP);
c2a = FI1 <= 0.*eye(size(FI1));
c2b = bigX >= eps.*eye(sX1+sX3);
constraints = [c1a, c1b, c2a, c2b];

%% Objective
objective = gamma;

%% Solve
options   = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
tic
sol       = optimize( constraints , objective , options );
soltime = toc;
change.dAc = value(dAc);
change.dBc = value(dBc);
change.dCc = value(dCc);
change.dDc = value(dDc);
change.dP = value(dP);
change.dX = value([dX1 , dX2 ; dX2' , dX3]);
Ac = Ac0 + change.dAc;
Bc = Bc0 + change.dBc;
Cc = Cc0 + change.dCc;
Dc = Dc0 + change.dDc;
P = P0 + change.dP;
X = value(bigX);
gamma = value(gamma);

end