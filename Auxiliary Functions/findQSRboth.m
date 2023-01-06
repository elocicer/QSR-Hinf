function [Q,S,R,Qc,Sc,Rc,Pp,Pc] = findQSRboth(A,B,C,D,Ac,Bc,Cc,Dc)
%{
Ethan LoCicero
2022
Finds a plant QSR and controller (Qc,Sc,Rc) = (-R,S',-Q).
I.e. finds a pair that satisfies the stability theorem.
(Qc,Sc,Rc) = (-R,S',-Q) is the best choice for QSR, but QSR not necessarily
optimized for the given plant or controller
%}

    %% Parameters
    epsilon = 10^-6; % definiteness tolerance
    sA = size(A,1);
    sAc = size(Ac,1);
    sC1 = size(C,1);
    sB2 = size(B,2);

    %% Variables
    Q_ = sdpvar( sC1,sC1 );
    S_ = sdpvar( sC1,sB2,'full' );
    R_ = sdpvar( sB2,sB2 );
    Qc = -R_;
    Sc = S_';
    Rc = -Q_;
    P_ = sdpvar( sA,sA );
    Pc_ = sdpvar( sAc,sAc );
    
    %% Constraint 1
    PA = P_*A;
    QD = Q_*D;
    DtS = D'*S_;
    CtQC = C'*Q_*C; CtQC = .5.*(CtQC + CtQC');
    DtQD = D'*Q_*D; DtQD = .5.*(DtQD + DtQD');
    offdiag = P_*B - C'*QD - C'*S_;
    Block1 = [PA + PA' - CtQC , offdiag ; offdiag' , -R_ - DtS - DtS' - DtQD ];
    c1 = Block1 <= 0;
    c2 = P_ >= epsilon*eye(sA);
    
    %% Constraint 2
    PAc = Pc_*Ac;
    QDc = Qc*Dc;
    DtSc = Dc'*Sc;
    CtQCc = Cc'*Qc*Cc; CtQCc = .5.*(CtQCc + CtQCc');
    DtQDc = Dc'*Qc*Dc; DtQDc = .5.*(DtQDc + DtQDc');
    offdiagc = Pc_*Bc - Cc'*QDc - Cc'*Sc;
    Block2 = [PAc + PAc' - CtQCc , offdiagc ; offdiagc' , -Rc - DtSc - DtSc' - DtQDc ];
    c3 = Block2 <= 0;
    c4 = Pc_ >= epsilon*eye(sAc);
    
    %% Solve
    constraints = [c1 c2 c3 c4];
    options = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
    sol = optimize( constraints , [] , options );
    
    %% Check
    if sol.problem == 0
        Q = value(Q_); S = value(S_); R = value(R_); Pp = value(P_);
        Qc = value(Qc); Sc = value(Sc); Rc = value(Rc); Pc = value(Pc_);
    elseif sol.problem == 4
        if min(eig(value(P_))) > 0 && max(eig(value(Block1))) < 10^(-5) && min(eig(value(Pc_))) > 0 && max(eig(value(Block2))) < 10^(-5)  
            disp([sol.info,', but constraints satisfied'])
            Q = value(Q_); S = value(S_); R = value(R_); Pp = value(P_);
            Qc = value(Qc); Sc = value(Sc); Rc = value(Rc); Pc = value(Pc_);
        else
            disp(sol.info)
            Q = []; S = []; R = []; Pp = [];
            Qc = []; Sc = []; Rc = []; Pc = [];
        end
    else
        disp(sol.info)
        Q = []; S = []; R = []; Pp = [];
        Qc = []; Sc = []; Rc = []; Pc = [];
    end
    
end
