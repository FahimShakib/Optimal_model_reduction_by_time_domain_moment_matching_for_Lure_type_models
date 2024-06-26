function [IO] = MTF_Solution(lure, sim_data, mtf_pars)

%% System definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Lur'e-type system:     %%%
%%%                             %%%
%%%   x_dot = A*x + B*u + L*w   %%%
%%%   y     = C*x       + D*w   %%%
%%%   z     = F*x + G*u + H*w   %%%
%%%   u     = -phi(y)           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define linear transfers
G_yw = ss(lure.A, lure.L, lure.C, lure.D);
G_yu = ss(lure.A, lure.B, lure.C, 0);
G_zw = ss(lure.A, lure.L, lure.F, lure.H);
G_zu = ss(lure.A, lure.B, lure.F, lure.G);

q = 2;

%% MTF algorithm setup

% Extract parameters
n = mtf_pars.n;
T = mtf_pars.T;
f = (0 : 1 : n-1) * (1/T);                                                 % Corresponding frequency vector with frequency steps 1/T

% Evaluate linear transfers
G_yu_vec = squeeze(freqresp(G_yu,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));       % Notes: 1) These spectrum repeat themselves
G_yw_vec = squeeze(freqresp(G_yw,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi)).';       %        2) The negative frequencies are stored in fft order behind the positive ones)
G_zu_vec = squeeze(freqresp(G_zu,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi)).';       %        3) Use squeeze to make it a vector
G_zw_vec = squeeze(freqresp(G_zw,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%   Calculate initial response of linear dynamics     %%%
%%%     1) Subject to excitation w(t)                   %%%
%%%     2) Assume u(t) = 0,                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W  = fft(sim_data.w);                                                         % Fourier coefficients W of excitation w(t). 
Y0 = G_yw_vec.*W;
Z0 = G_zw_vec.*W;

% MTF algorithm execution
y0      = zeros(n,q);                                                      % initial guess for steady state response for e
y       = y0;                                                              % output signal (error e in this case)
Y       = zeros(n,q);                                                      % fft of error signal
Y_old   = zeros(n,q);                                                      % used for determining convergence in MTF algorithm
y_old   = zeros(n,q);                                                      % used for determining convergence in MTF algorithm
U       = zeros(n,q);                                                      % fft of output of the nonlinearity

yerror  = 1;                                                               % set initial error > tol
k       = 0;
while k < mtf_pars.max_iter && yerror > mtf_pars.convtol
    if any(isnan(Y))
        disp('Response diverges')
        break;
    end
    u     = lure.nonlin(y.').';
    U     = fft(u);                                                        % fft of output of nonlinearity
    Y     = Y0;
    for i = 1:q
        Y = Y + (squeeze(G_yu_vec(:,i,:)).*(U(:,i).')).';                                              % linear dynamics in frequency domain , superposition ?!
    end
    y     = real(ifft(Y));                                                 % use real value because there is a very small imaginary part present (1e-21)
    Z     = Z0 + sum(G_zu_vec.*U,2);
    z     = real(ifft(Z));
    
    % Error based on Fourier coefficients to check convergence
    yerror   = norm(Y-Y_old)/norm(Y_old);
    %yerror   = norm((y-y_old)./y_old);
    nrm(k+1) = yerror;
    if k < mtf_pars.max_iter - 1
    Y_old    = Y;
    y_old    = y;
    end
    
    k     = k + 1;                                                         % update the iteration index
end


% Formulate steady state solution
IO.w        = sim_data.w;
IO.t        = sim_data.t;
IO.y        = y;
IO.u        = lure.nonlin(IO.y.').';
IO.z        = z;

end