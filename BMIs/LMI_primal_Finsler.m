function LMI_out = LMI_primal_Finsler(A,Bu,Bw,Ccal,Hzeta,G,X1,X2,N,S,L,gamma,form)
% 2024 - 01 - 17

% X1 = randn(n+nu);
% X2 = randn(nu);

n  = size(A,1);
nu = size(G,1);
q  = size(Bw,2);
m  = size(Bu,2);
p  = size(Ccal,1)-q;


Ncal_gamma  = BRL_matrix(blkdiag(A,S),[Bu Bw; zeros(nu,m+q)],Ccal,0,X1,gamma);
Ncal        = BRL_matrix(S,zeros(nu,q),Hzeta,0,X2,1);
M1          = blkdiag(-X1,-X2,Ncal_gamma,Ncal);

M2 = blkdiag(zeros(n+2*nu),[X1;zeros(m+p+2*q,n+nu)],[X2;zeros(2*q,nu)]);

Gbar    = [zeros(n,n+nu+m+q+p+q);
           zeros(nu,n) -G*L G zeros(nu,p+q)];
Glambda = G*[zeros(m,q);eye(q)];
M3      = blkdiag(zeros(n+2*nu),Gbar',[-(G*L)';Glambda';zeros(q,nu)]);

if form == 1
    LMI_primal = M1 + M2*M3' + M3*M2';
    LMI_out = LMI_primal;
else
    M = [M1 M2;M2' zeros(2*n+4*nu)];

    LMI_Finsler = M + [M3;-eye(2*n+4*nu)]*N' + N*[M3;-eye(2*n+4*nu)]';
    LMI_out = LMI_Finsler;
end
