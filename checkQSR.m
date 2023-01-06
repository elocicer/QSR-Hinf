function [isQSR,P] = checkQSR(A,B,C,D,Q,S,R,varargin)
% Ethan LoCicero
% 2022
% Checks that the system satisfies the QSR bounds

%% Determine output information level
if isempty(varargin)
    quiet = false;
elseif strcmp(varargin{1},'quiet')
    quiet = true;
else
    quiet = false;
end

%% Check symmetry
if ~issymmetric(Q) || ~issymmetric(R)
    disp('Q and R must be symmetric')
    isQSR = [];
    return
end

%% Parameters
epsilon = 10^-6; % definiteness tolerance
sA = size(A,1);

%% Variables
P = sdpvar( sA,sA );

%% Premultiply for speed
PA = P*A;
QD = Q*D;
DtS = D'*S;
CtQC = C'*Q*C; CtQC = .5.*(CtQC + CtQC');
DtQD = D'*Q*D; DtQD = .5.*(DtQD + DtQD');

%% Constraints
offdiag = P*B - C'*QD - C'*S;
BigMat =[PA + PA' - CtQC , offdiag ; offdiag' , -R - DtS - DtS' - DtQD ];
c1 = BigMat <= 0;
c2 = P >= epsilon*eye(sA);
constraints = [c1 c2];

%% Solve
options = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
sol = optimize( constraints , [] , options );

%% Check
if sol.problem == 0
    isQSR = true; P = value(P);
elseif sol.problem == 4
    if min(eig(value(P))) > 0 && max(eig(value(BigMat))) <= 0
        if ~quiet
            disp([sol.info,', but constraints satisfied'])
        end
        isQSR = true; P = value(P);
    else
        if ~quiet
            disp(sol.info)
        end
        isQSR = false; P = [];
    end
else
    if ~quiet
        disp(sol.info)
    end
    isQSR = false; P = [];
end

end