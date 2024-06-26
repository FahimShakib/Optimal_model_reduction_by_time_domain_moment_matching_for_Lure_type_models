% This code corresponds to the results in [Scarciotti & Astolfi, 2017]
clear all; clc
addpath('BMIs')
addpath('Functions')

%% BMI solver settings
kmax   = 10;        % Maximum number of iterations inside BMI iterations
gtol   = 1e-8;      % Tolerance for the objective
eigtol = 1e-6;      % Tolerance for checking negativity of eigenvalues

%% Load beam model
load Data/Beam_Model.mat
q = 2;              % Number of nonlinear functions
p_locations = 1;    % Outputs
m = 1;              % Inputs

%% Harmonic oscillator dot{s} = Ws, u = Lw
f0  = [10.2 64.1];                  % Base frequencies
nu   = 2*length(f0) + 1;
S   = [0];
for k = 1:length(f0)
    S = blkdiag(S,[0 1;-1 0]*2*pi*f0(k));
end

L_ = [1 1 0 1 0];
L  = [L_;L_;L_];
L  = L(:,1:nu);
w0  = L_';

SG = ss(S,[],L,[]);

%% Bring into balanced form
[sys,hsing] = balreal(sys);       % Create a balanced realization and compute Hankel singular values
display(['Lower bound for gamma = ' num2str(hsing(nu+1))])


%% Compute CPi
Pi   = sylvester(sys.A,-S,-sys.B*L);
CPi  = sys.C*Pi;
Hzeta = Cz*Pi;

% Definte the ROM as a function of G
sys_red = @(G) ss(S-G*L,G,CPi,0);

%% See whether we achieve tangential moment matching
% Take a random G
G_random = randn(nu,m+q);
sys_red_random = sys_red(G_random);

plot_bode_tangential = 0;
% As all tangential directions are the same (ell = [1 1 1]'), there should 
% be a match at 0, 10.2 and 64.1 Hz
if plot_bode_tangential
    figure
    bode(sys*L(:,1))
    hold all
    bode(sys_red_random*L(:,1))
end

%% Initialize problem (Theorem 5)
Q1 = sdpvar(nu,m);
Q2 = sdpvar(nu,q);
X_init  = sdpvar(nu,nu);

eigtol = 1e-9;

[~, LMI_BRL_init] = LMIs_BRL(X_init*S-[Q1 Q2]*L,Q2,Hzeta,zeros(2,2),eye(nu),1);

LMI_init = [X_init >= eigtol*eye(nu); LMI_BRL_init <= -eigtol*eye(size(LMI_BRL_init))];

optimize(LMI_init)

G_init       = inv(double(X_init))*double([Q1 Q2]);
sys_red_init = sys_red(G_init);
gamma_init   = hinfnorm(sys-sys_red_init)*1.1;

sys_red_lemma = sys_red_init;

display('Theorem 5')
display(['||Gamma_{zeta,lambda}||_\infty = ' num2str(hinfnorm(sys_red_init(2:3,2:3)))])
display('which should be smaller than one')


%% Solve BMI problem to find the optimal G
% Options
opts_BMI.kmax    = 1000;
opts_BMI.gtol    = 1e-6;
opts_BMI.eigtol  = 1e-9;
opts_BMI.imax    = 10;
plot_figure      = 1;

% Initialisation
GG    = {{G_init}};
gamma = {gamma_init};

gammatr = {};
X1 = {};
X2 = {};

%% Solve BMI problem
% To stop the iterations manually, pause Matlab and commond stopp = 1 in
% the command window
tic
for i = 1:100
    stopp = 0;
    display(' ')
    display(['Lower bound for gamma = ' num2str(hsing(nu+1))])
    if mod(i,2)
        % Iterate in primal form
        display(' ')
        display('Start iterations in primal form')
        opts_BMI.kmax    = 100;
        
        form = 1;
        [~,~,gamma{end+1},X1{end+1},X2{end+1},gammatr{end+1},GG] = ...
            BMI_solver(sys,n,nu,GG,S,L,gamma{end}(end),sys_red,Pi,plot_figure,form,opts_BMI);
    else
        % Iterate in Finsler's form
        display(' ')
        display('Start iterations Finslers form')
        opts_BMI.kmax    = 100;

        form = 0;
        [~,~,gamma{end+1},X1{end+1},X2{end+1},gammatr{end+1},GG] = ...
            BMI_solver(sys,n,nu,GG,S,L,gamma{end}(end),sys_red,Pi,plot_figure,form,opts_BMI);
    end

    if stopp
        display('Stopped')
        break;
    end
end
toc


%% Extract Lur'-type models
G_final = GG{1}{end};
sys_red_final = sys_red(G_final);