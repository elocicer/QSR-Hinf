function isGammaBounded = checkHinf(A,B,C,D,gamma)
% Ethan LoCicero
% 2022
% Determines whether the LTI system has an H-infinity norm below gamma

    %% parameters
    epsilon = 10^-6; % definiteness tolerance
    sA = size(A,1);
    sB2 = size(B,2);

    %% variables
    P_ = sdpvar( sA,sA );

    %% constraints
    PA = P_*A;
    CtC = C'*C; CtC = .5.*(CtC + CtC');
    DtD = D'*D; DtD = .5.*(DtD + DtD');
    offdiag = P_*B + C'*D;
    mainblock = [PA + PA' + CtC , offdiag ; offdiag' , -gamma^2*eye(sB2)  + DtD ];
    c1 = mainblock <= 0;
    c2 = P_ >= epsilon*eye(sA);
    constraints = [c1 c2];

    %% Solve
    options = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
    sol = optimize( constraints , [] , options );
    
    %% Check
    if sol.problem == 0
        isGammaBounded = true;
    elseif sol.problem == 4
        if min(eig(value(P_))) > 0 && max(eig(value(mainblock))) <= 0
            disp([sol.info,', but constraints satisfied'])
            isGammaBounded = true;
        else
            disp(sol.info)
            isGammaBounded = false; disp('Here 3')
        end
    else
        disp(sol.info)
        isGammaBounded = false; disp('Here 4')
    end
    
end