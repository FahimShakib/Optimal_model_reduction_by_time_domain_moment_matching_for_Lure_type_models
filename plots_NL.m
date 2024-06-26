clear all; clc
addpath('Functions')

load Data/'20240129 - NL Example.mat'

%% Manual Bode plot (Figure 3)

% Take the best G and best reduced model
G_final         = GG{1}{end};
sys_red_final   = sys_red(G_final);

% Compute FRF
[FRF_sys_mag,~,f_vector_sys]            = bode(sys);
[FRF_sys_final_mag,~,f_vector_final]    = bode(sys_red(GG{end}{end})-sys);
[FRF_sys_pp_mag,~,f_vector_pp]          = bode(sys_red_pp-sys);
[FRF_sys_init_mag,~,f_vector_init]      = bode(sys_red(GG{end}{1})-sys);

hb = figure;
hb.Position = [300 600 1200 375];
for k = 1:3
    for i = 1:3
        subplot(3,3,(k-1)*3+i)
        p1 = semilogx(f_vector_final/2/pi,db(squeeze(FRF_sys_final_mag(k,i,:))),'b');
        hold all
        semilogx(f_vector_final/2/pi,f_vector_final*0+max(db(squeeze(FRF_sys_final_mag(k,i,:)))),'b--')
        p2 = semilogx(f_vector_init/2/pi,db(squeeze(FRF_sys_init_mag(k,i,:))),'r');
        semilogx(f_vector_init/2/pi,f_vector_init*0+max(db(squeeze(FRF_sys_init_mag(k,i,:)))),'r--')
        xlim([1e0 1e3])
        if k == 3
            xlabel('Frequency [Hz]')
        end
        if k == 2
            ylim([-100 0])
        end
        set(gca,'fontsize', 12)
        set(findall(gcf,'type','line'),'linewidth',2)
        if i == 1 && k == 2
            ylabel(['$z_{(1)}-\zeta_{(1)}$ [dB]' ])
        end
        if i == 1 && k == 3
            ylabel(['$z_{(2)}-\zeta_{(2)}$ [dB]' ])
        end
        if i == 1 && k == 1
            ylabel(['$y-\psi$ [dB]' ])
        end
        if k == 3 && i == 3
            legend([p1,p2],'$\Sigma - \hat\Sigma$','$\Sigma - \hat\Sigma^\circ$','location','SE')
            ylim([-50 0])
        end
        if i == 1 && k == 1
            title(['$u$'])
            ylim([-160 -65])
        end
        if i == 2 && k == 1
            title(['$w_{(1)}$'])
            ylim([-155 -75])
        end
        if i == 3 && k == 1
            title(['$w_{(2)}$'])
            ylim([-120 -72])
        end
    end
end

%% Time-domain simulations with block-wave inputs (Figure 4)
% Generate Lur'e-type models
clc
lure_full  = sys2lure(sys,phi_sym,q);
lure_final = sys2lure(sys_red(GG{end}{end}),phi_sym,q);
lure_init  = sys2lure(sys_red(GG{end}{1}),phi_sym,q);

% Time-domain simulation using the MTF algorithm
% Square wave simulation
N    = 1000;
fsim = 0.5;

mtf_pars.n = N;
mtf_pars.T = 1/fsim;
mtf_pars.max_iter = 1000;
mtf_pars.convtol  = 1e-10;

simdata.t = linspace(0,1/fsim,N);
simdata.w = 1e4*square(simdata.t*2*pi*fsim)';
w_low     = simdata.w;

IO_full   = MTF_Solution(lure_full, simdata, mtf_pars);
IO_final  = MTF_Solution(lure_final, simdata, mtf_pars);
IO_init   = MTF_Solution(lure_init, simdata, mtf_pars);

display(['Norm full ' num2str(signal_norm(IO_full.z,1))])
display(['Error full - final ' num2str(signal_norm(IO_full.z-IO_final.z,1))])
display(['Error full - init ' num2str(signal_norm(IO_full.z-IO_init.z,1))])

h = figure;
h.Position = [100 100 600 300];
subplot(221)
plot(IO_full.t,IO_full.z,'k--')
hold all
plot(IO_final.t,IO_final.z,'b')
plot(IO_final.t,IO_init.z,'r')

ylabel('Output')
set(gca,'fontsize', 12)

subplot(223)
plot(IO_full.t,IO_full.z-IO_final.z,'b')
hold all
plot(IO_full.t,IO_full.z-IO_init.z,'r')
xlabel('Time [sec]')
ylabel('Output error')
set(gca,'fontsize', 12)
ylim([-1 1]*0.25)

% Square wave simulation (with a different frequency now)
N = 1000;
fsim = 10;

mtf_pars.n = N;
mtf_pars.T = 1/fsim;
mtf_pars.max_iter = 1000;
mtf_pars.convtol  = 1e-10;

simdata.t = linspace(0,1/fsim,N);
simdata.w = 1e4*square(simdata.t*2*pi*fsim)';
w_high    = simdata.w;

IO_full   = MTF_Solution(lure_full, simdata, mtf_pars);
IO_final  = MTF_Solution(lure_final, simdata, mtf_pars);
IO_init   = MTF_Solution(lure_init, simdata, mtf_pars);

subplot(222)
plot(IO_full.t,IO_full.z,'k--')
hold all
plot(IO_final.t,IO_final.z,'b')
plot(IO_final.t,IO_init.z,'r')
set(gca,'fontsize', 12)
ylim([-1 1]*3.3)

legend('$\bar y_u$','$\bar \psi_u$','$\bar \psi_u^\circ$')

subplot(224)
plot(IO_full.t,IO_full.z-IO_final.z,'b')
hold all
plot(IO_full.t,IO_full.z-IO_init.z,'r')
legend('$\bar y_u - \bar \psi_u$','$\bar y_u - \bar \psi_u^\circ$')
xlabel('Time [sec]')
ylim([-1 1]*2.5)

set(gca,'fontsize', 12)
set(findall(gcf,'type','line'),'linewidth',2)

% Display values useful for Table 1
display(['Norm full ' num2str(signal_norm(IO_full.z,1))])
display(['Error full - final ' num2str(signal_norm(IO_full.z-IO_final.z,1))])
display(['Error full - init ' num2str(signal_norm(IO_full.z-IO_init.z,1))])

%% Compute error bounds (for Table 1)
% Gains of the full-order model
g_yw = hinfnorm(sys(1,2:3));
g_zw = hinfnorm(sys(2:3,2:3));

% General formula (eqn (12) in the paper)
error = @(sys,sys_red,norm_u)hinfnorm(sys-sys_red)*(1+g_yw/(1-g_zw))*(1+hinfnorm(sys_red(2:3,1))/(1-hinfnorm(sys_red(2:3,2:3))))*norm_u;

display(['Error bound init wlow ' num2str(error(sys,sys_red_init,signal_norm(w_low,1)))])
display(['Error bound init whigh ' num2str(error(sys,sys_red_init,signal_norm(w_low,1)))])

display(['Error bound final wlow ' num2str(error(sys,sys_red_final,signal_norm(w_low,1)))])
display(['Error bound final whigh ' num2str(error(sys,sys_red_final,signal_norm(w_low,1)))])

%% Simulate with Simulink to find computation time for numerical forward integration (for Table 1)
fsim = 0.5;
dur = 1/fsim*10; % run all simulation for 10 periods

% Average over 5 simulations
tic
for k = 1:5
    sim Simulink_Beam_full.slx
end
t_full = toc/5;
tic
for k = 1:5
    sim Simulink_Beam_reduced
end
t_reduced = toc/5;
tic
for k = 1:5
    sim Simulink_Beam_init.slx
end
t_init = toc/5;

display('Computation times')
display(['Full ' num2str(t_full)])
display(['Reduced ' num2str(t_reduced)])
display(['Init ' num2str(t_init)])
