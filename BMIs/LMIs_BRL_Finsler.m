function [constraints,C1,X] = LMIs_BRL_Finsler(S,L,sys0,H,v,n,gamma,X1,X2,X3,G,N)
A = sys0.A;
B = sys0.B;
C = sys0.C;
D = 0;

M0 = [A'*X1+X1*A   A'*X3+X3*S   X1*B   C';
      X3'*A+S'*X3' S'*X2+X2*S   X3'*B  -H';
      B'*X1        B'*X3        -gamma 0;
      C            -H           0      -gamma];
  
M1      = [X3' X2 zeros(v,2)];
Rcomp   = [zeros(v,n) -G*L G zeros(v,1) -eye(v)];

C1 = [M0 M1'; M1 zeros(v)]+N*Rcomp+Rcomp'*N';

X  = [X1  X3;
      X3' X2];

 constraints = [C1<=-eye(n+2*v+2)*1e-6, X>=eye(n+v)*1e-6];