function [constraints,BRL,X] = LMIs_BRL(A,B,C,D,X,gamma)

m = size(B,2);
p = size(C,1);

BRL         = [A'*X+X*A  X*B            C'; 
               B'*X      -gamma*eye(m)  D'; 
               C         D              -gamma*eye(p)]; 
constraints = [BRL <= -eye(size(BRL))*1e-6, X>=eye(size(X))*1e-6];
