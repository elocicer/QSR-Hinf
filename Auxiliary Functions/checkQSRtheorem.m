function isStable = checkQSRtheorem(Qg,Sg,Rg,Qc,Sc,Rc,varargin)
% Ethan LoCicero
% 2022
% Checks that two sets of QSR-dissipativity parameters satisfy the
% Dissipativity Theorem. Varargin can be empty, leaving alpha as a
% variable, or can provide a fixed value for alpha.

    %% Check symmetry
    if ~issymmetric(Qg) || ~issymmetric(Rg) || ~issymmetric(Qc) || ~issymmetric(Rc)
        disp('Q and R must be symmetric')
        isStable = [];
        return
    end
    
    %% parameters
    epsilon = 10^-6; % definiteness tolerance
    if isempty(varargin)
        %% variable alpha
        % variables
        alpha = sdpvar;
        % constraints
        c2 = alpha >= epsilon;
        Theorem = [alpha.*Qg + Rc , -alpha.*Sg + Sc' ; -alpha.*Sg' + Sc , alpha.*Rg + Qc];
        sThm = size(Theorem,1);
        c1 = Theorem <= -eps*eye(sThm);
        constraints = [c1 c2];
        % solve
        options = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
        sol = optimize( constraints , [] , options );
        % check
        if sol.problem == 0
            isStable = true;
        elseif sol.problem == 4
            if value(alpha) > 0 && max(eig(value(Theorem))) <= 10^-6
                disp([sol.info,', but constraints satisfied'])
                isStable = true;
            else
                disp(sol.info)
                isStable = false;
            end
        else
            disp(sol.info)
            isStable = false;
        end
    else
        %% fixed alpha provided
        alpha = varargin{1};
        % manually check constraint satisfaction
        Theorem = [alpha.*Qg + Rc , -alpha.*Sg + Sc' ; -alpha.*Sg' + Sc , alpha.*Rg + Qc];
        if ~issymmetric(Theorem) || alpha <= 0; disp('Invalid choices of Q,S,R, or alpha'); return; end
        if max(real(eig(Theorem))) < 0
            isStable = true;
        else 
            isStable = false;
        end
    end
    
end
