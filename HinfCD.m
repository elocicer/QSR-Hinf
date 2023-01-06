function [Cnew,Dnew] = HinfCD(A,B,C,D,gamma)
% Ethan LoCicero
% 2022
% Changes the output matrices C and D of the system to satisfy the
% H-infinity gain gamma.

    %% Parameters
    epsilon = 10^-6; % definiteness tolerance
    sA = size(A,1);
    [sC1,sC2] = size(C);
    sB2 = size(B,2);
    [sD1,sD2] = size(D);

    %% Variables
    P_ = sdpvar( sA,sA );
    C_ = sdpvar( sC1,sC2,'full' );
    D_ = sdpvar( sD1,sD2,'full' );

    %% Premultiply to eliminate redundancies for speed
    PA = P_*A;
    PB = P_*B;
    
    %% Constraints
    mainblock = [PA + PA' , PB , C_' ; PB' , -gamma*eye(sB2) , D_' ; C_ , D_ , -gamma*eye(sC1) ];
    c1 = mainblock <= 0;
    c2 = P_ >= epsilon*eye(sA);
    constraints = [c1 c2];

    %% Objective 
    dC = C-C_;
    dD = D-D_;
    objective = trace(dC*dC') + trace(dD*dD');

    %% Solve
    options = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
    sol = optimize( constraints , objective , options );    

    %% Check
    if sol.problem == 0
        Cnew = value(C_); Dnew = value(D_);
    elseif sol.problem == 4
        if min(eig(value(P_))) > 0 && max(eig(value(mainblock))) <= 0
            disp([sol.info,', but constraints satisfied'])
            Cnew = value(C_); Dnew = value(D_);
        else
            disp(sol.info)
            Cnew = []; Dnew = [];
        end
    else
        disp(sol.info)
        Cnew = []; Dnew = [];
    end
    
end