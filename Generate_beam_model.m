clear all; close all; clc

% Select the number of beam elements
Ne  = 12;

%% Beam parameters
Le  = 2;            % Beam length [m] 
h   = 3e-2;         % Beam heigth [m]
w   = 5e-2;         % Beam width [m]

E   = 200e9;        % Young modulus [Pa]
rho = 7746;         % Density [kg/m^3]

A   = w*h;          % Cross section [m^2]
I   = 1/12*w^3*h;   % Second moment of inertia [m^4]

%% Collect beam matrices 
[K,M] = composeMatrices(Ne, Le, E, I, rho, A);
% Clamp the cantilever beam at one side
K = K(3:end,3:end);
M = M(3:end,3:end);
% Solve eigenvalue problem to obtain modes  
[U,Ld] = eig(full(K),full(M));
U = real(U);
% Add 0.1 percent modal damping for stability
Bd = (U.')\(0.1*sqrt(Ld))/U; 

%% Construct state-space model
% Retrieve size of x
n   = size(M,1);

% Select input location
q_input          = n/2+1;          % External loading at the middle of the beam
q_spring_middle  = n/2+1;          % Spring force acting on the middle of the beam
q_spring_end     = n;              % Spring force acting on the end of the beam
Q                       = zeros(2,n);
Q(1,q_input)            = 1;
Q(2,q_spring_middle)    = 1;
Q(3,q_spring_end)       = 1;

% Select sensor location
q_sensor = n;

Ass = [zeros(n) eye(n); -inv(M)*K -inv(M)*Bd];
Bss = [zeros(n); inv(M)];
Css = zeros(3,2*n);
Css(1,q_sensor)         = 1;
Css(2,q_spring_middle)  = 1;
Css(3,q_spring_end)     = 1;
Dss = 0;

sys = ss(Ass,Bss*Q',Css,Dss);
n   = size(Ass,1);

loop_gain = hinfnorm(sys(2:3,2:3));

% One-sided spring has a spring constant of
K_spring = 3e3*eye(2);

phi = @(z) K_spring*max(0,z);

%% Simulate full-order LTI model (with fixed spring forces)
% Time vector definition
ts  = 0.01;
t   = 0:ts:10-ts;

% Input signal definition
freq   = 1;
u(1,:) = sin(2*pi*freq*t)*0;
u(2,:) = sin(2*pi*freq*t);
u(3,:) = sin(2*pi*freq*t)*0;

[y,t,x] = lsim(sys,u',t);

%% Plot deflection
plot_animation = 0;
if plot_animation
    l = linspace(0,Le,Ne);
    P = 100;
    ind = round(length(t)/P);
    figure
    for i = 1:ind:length(t)
        subplot(211)
        plot(l,x(i,1:2:end/2))
        title(['Time ' num2str(t(i)) ' sec'])
        ylabel('Beam deflection [m]')
        ylim([min(min(x(:,1:2:end/2))) max(max(x(:,1:2:end/2)))]*1.1)
        subplot(212)
        plot(l,x(i,2:2:end/2))
        ylim([min(min(x(:,2:2:end/2))) max(max(x(:,2:2:end/2)))]*1.1)
        xlabel('Spatial coordinate [m]')
        ylabel('Beam inclination [m]')
        pause(0.1)
    end
end

% return

%% Plot deflection at sensor location
plot_deflection = 0;
if plot_deflection
    figure
    plot(t,y)
    xlabel('Time [sec]')
    ylabel('Beam end deflection [m]')
end

%% Plot bode
plot_bode = 0;
if plot_bode
    figure
    bodemag(sys)
    xlim([1 1e4])
    grid on
end

%% Extract matrices for the Lur'e-type model
A  = sys.A;
Bu = sys.B(:,1);
Bw = sys.B(:,[2 3]);
Cz = sys.C([2 3],:);
Cy = sys.C(1,:);

%% Transform the Lur'e-type model such that the resulting nonlinearity satisfies (3)
A_tf  = A + 0.5*Bw*K_spring*Cz;
Bw_tf = 0.5*Bw*K_spring;

Gzw = ss(A_tf,Bw_tf,Cz,0);

phi_tf = @(z)(inv(0.5*K_spring)*(phi(z)-0.5*K_spring*z));

% View plot of transformed nonlinearity - it satisfies (3) in the paper
zz = [linspace(-10,10,100);linspace(-10,10,100)];
ww = phi_tf(zz);

figure
plot(zz,ww(1,:))
clear zz ww

%% Perform loop-transformation on the linear part - define A, B, etc. from here
sys_tf = ss(A_tf,[Bu Bw_tf],[Cy;Cz],0);

%% Rename variables before saving
sys = sys_tf;
A   = sys_tf.A;
Bu  = sys_tf.B(:,1);
Bw  = sys_tf.B(:,[2 3]);
Cz  = sys_tf.C([2 3],:);
Cy  = sys_tf.C(1,:);
phi_sym = phi_tf;

display(['||Gamma_{z,w}||_\infty = ' num2str(hinfnorm(sys_tf(2:3,2:3)))])
display('which should be smaller than one')

% save 'Data\Beam_Model.mat' sys A Bu Bw Cz Cy n phi_sym K_spring
