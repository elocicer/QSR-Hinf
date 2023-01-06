function [L,soltime] = UpdateL(B,R,U,V,S,C)

sB = size(B,1);
sU = size(U,2);
Z = sdpvar( sB,sB );
W = sdpvar( sU,sU );
objective = trace(Z);
constraint = [Z-B*R*U*W*U'*R'*B' , (V*S*C)' ; V*S*C , W] > 0;
options   = sdpsettings('verbose' , 0 , 'solver' , 'mosek' , 'cachesolvers' , 1);
tic
sol       = optimize( constraint , objective , options );
soltime = toc;
L = sqrt(value(W));

end
