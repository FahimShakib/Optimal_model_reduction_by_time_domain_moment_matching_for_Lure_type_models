function lure = sys2lure(sys,nonlin,q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Lur'e-type system:     %%%
%%%                             %%%
%%%   x_dot = A*x + B*u + L*w   %%%
%%%   y     = C*x       + D*w   %%%
%%%   z     = F*x + G*u + H*w   %%%
%%%   u     = -phi(y)           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size(sys.B,2)-q;
p = size(sys.C,1)-q;

lure.A = sys.A;
lure.B = sys.B(:,m+1:end);
lure.L = sys.B(:,1:m);
lure.C = sys.C(p+1:end,:);
lure.D = 0;
lure.F = sys.C(1:p,:);
lure.G = 0;
lure.H = 0;
lure.nonlin = nonlin;



