function [gamma,P] = findHinf(A,B,C,D)
% Ethan LoCicero
% 2022
% Finds the H-infininty norm of the input system.

    %% Parameters
    epsilon = 10^-4; % definiteness tolerance
    sA = size(A,1);
    sB2 = size(B,2);
    sC1 = size(C,1);
    %% Variables
    gamma_ = sdpvar;
    P_     = sdpvar( sA,sA );
    
    %% Premultiply to eliminiate redundances for speed
    PA = P_*A;
    PB = P_*B;

    %% Constraints
    mainblock = [PA + PA' , PB , C' ; PB' , -gamma_*eye(sB2) , D' ; C , D , -gamma_*eye(sC1) ];
    c1 = mainblock <= 0;
    c2 = P_ >= epsilon*eye(sA);
    c3 = gamma_ >= 10^(-6);
    constraints = [c1 c2 c3];

    %% Solve
    options = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
    sol = optimize( constraints , gamma_ , options );
    
    %% Check
    if sol.problem == 0
        gamma = value(gamma_); P = value(P_);
    elseif sol.problem == 4
        if min(eig(value(P_))) > 0 && max(eig(value(mainblock))) <= 0
            disp([sol.info,', but constraints satisfied'])
            gamma = value(gamma_); P = value(P_);
        else
            disp(sol.info)
            gamma = []; P = [];
        end
    else
        disp(sol.info)
        gamma = []; P = [];
    end
    
end
