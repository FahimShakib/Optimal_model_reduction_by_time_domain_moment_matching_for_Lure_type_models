function BRL = BRL_matrix(A,B,C,D,X,gamma)

m = size(B,2);
p = size(C,1);

if D == 0
    D = zeros(p,m);
end

BRL = [X*A+A'*X, X*B,           C';
       B'*X,     -gamma*eye(m), D';
       C,        D,             -gamma*eye(p)];