% This code corresponds to the results in
clear all; clc
addpath('BMIs')
addpath('Functions')

%% BMI solver settings
kmax   = 10;        % Maximum number of iterations inside BMI iterations
gtol   = 1e-8;      % Tolerance for the objective
eigtol = 1e-6;      % Tolerance for checking negativity of eigenvalues

%% LTI model
p = 1;
m = 1;
n = 6;

sys = ss(zpk([1 2 3 4 5],[-1 -2 -3 -4 -5 -6],1));

A   = sys.A;
B   = sys.B;
C   = sys.C;

%% Harmonic oscillator dot{s} = Ws, u = Lw
f0  = 1;                  % Interpolation frequency [Hz]
S   = [0 1;-1 0]*2*pi*f0;
nu  = size(S,1);
L   = [1 0];

SG = ss(S,[],L,[]);

%% Bring into balanced form
[sys,hsing] = balreal(sys);       % Create a balanced realization and compute Hankel singular values

display(['Lower bound for gamma = ' num2str(hsing(nu+1))])

%% Compute CPi
Pi = sylvester(sys.A,-S,-sys.B*L);
CPi = sys.C*Pi;

% Definte the ROM as a function of G
sys_red = @(G) ss(S-G*L,G,CPi,0);

%% figure
clear Hnorm min_G min_error
min_error = inf;
N = 100;
GG2 = linspace(-2,2,N);
GG1 = logspace(-6,2,N+1);
Hnorm = inf(length(GG1),length(GG2));
for i = 1:size(GG1,2)
    for k = 1:size(GG2,2)
        tmp = hinfnorm(sys-sys_red([GG1(i);GG2(k)]));
        if tmp < min_error
            min_G = [GG1(i);GG2(k)];
            min_error = tmp;
        end
        Hnorm(i,k) = tmp;
    end
end

min_error
hinfnorm(sys-sys_red(min_G))
min(min(Hnorm))

% Plot figure
figure
surf((GG1),GG2,db(Hnorm'))
xlabel('$G_{(1)}$')
ylabel('$G_{(2)}$')
zlabel('Error')
set(gca,'xscale','log')


%% See whether we achieve moment matching
% Take a random G
G_random        = randn(nu,1);
sys_red_random  = sys_red(G_random);
gamma_random    = hinfnorm(sys_red_random-sys);

plot_bode = 0;
if plot_bode
    figure
    bode(sys)
    hold all
    bode(sys_red_random)
end

display('Pole placement model')
display(['||Error||_\infty = ' num2str(hinfnorm(sys_red_random-sys))])

%% Solve BMI problem to find the optimal G
% Options
opts_BMI.kmax    = 1000;
opts_BMI.gtol    = 1e-6;
opts_BMI.eigtol  = 1e-6;
opts_BMI.imax    = 10;

% Initial points
G_init = [0.1 0.1 1 1 10; 1 -1 1 -1 25];
G1 = [0.1; 1];
G2 = [0.1; -1];
G3 = [1; 1];
G4 = [1; -1];
G5 = [10; 25];

for k = 1:5
    GG      = {{G_init(:,k)}};
    gamma   = {hinfnorm(sys-sys_red(GG{end}{end}))*1.05};
    
    plot_figure = 0;
    
    gammatr = {};
    X = {};
    
    %% Solve BMI problem
    tic
    for i = 1:100
        stopp = 0;
        display(' ')
        display(['Lower bound for gamma = ' num2str(hsing(nu+1))])
        if mod(i,2)
            % Iterate in primal form
            display(' ')
            display('Start iterations in primal form')
            opts_BMI.kmax    = 1000;
            
            form = 1;
            [~,~,gamma{end+1},X{end+1},gammatr{end+1},GG] = ...
                BMI_solver_LTI(sys,n,nu,GG,S,L,gamma{end}(end),sys_red,Pi,plot_figure,form,opts_BMI);
        else
            % Iterate in Finsler's form
            display(' ')
            display('Start iterations Finslers form')
            opts_BMI.kmax    = 1000;
    
            form = 0;
            [~,~,gamma{end+1},X{end+1},gammatr{end+1},GG] = ...
                BMI_solver_LTI(sys,n,nu,GG,S,L,gamma{end}(end),sys_red,Pi,plot_figure,form,opts_BMI);
        end
    
        if stopp
            display('Stopped')
            break;
        end
    end
    elapsed_time = toc
end