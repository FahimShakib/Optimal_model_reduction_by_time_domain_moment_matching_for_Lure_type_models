function LMI_out = LMI_primal_Finsler_LTI(A,B,Ccal,G,X,N,S,L,gamma,form)
% 2024 - 01 - 17

n  = size(A,1);
nu = size(G,1);
m  = 1;
p  = 1;

Ncal_gamma  = BRL_matrix(blkdiag(A,S),[B; zeros(nu,m)],Ccal,0,X,gamma);
M1          = blkdiag(-X,Ncal_gamma);

% X11 = X(1:n,1:n);
% X12 = X(1:n,n+1:end);
% X22 = X(n+1:end,n+1:end);

M2 = blkdiag(zeros(n+nu),[X;zeros(m+p,n+nu)]);

Gbar    = [zeros(n,n+nu+m+p);
           zeros(nu,n) -G*L G zeros(nu,p)];
M3      = blkdiag(zeros(n+nu),Gbar');

if form == 1
    LMI_primal = M1 + M2*M3' + M3*M2';
    LMI_out = LMI_primal;
else
    M = [M1 M2;M2' zeros(2*n+2*nu)];

    LMI_Finsler = M + [M3;-eye(2*n+2*nu)]*N' + N*[M3;-eye(2*n+2*nu)]';
    LMI_out = LMI_Finsler;
end
