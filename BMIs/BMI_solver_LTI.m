function [sysr,G,gamma_,X,gammatr,G_array] = BMI_solver_LTI(sys,n,nu,G_array,S,L,gamma_prev,sys_red,Pi,plot_figure,form,opts_BMI)


% Optimize in the primal form for X and G
% Iteratively: fix G and solve for X, then X and solve G (etc) 


% imax is the number of iterations in bi-section
imax = opts_BMI.imax;

options = sdpsettings('solver','mosek','verbose',0);

gamma_  = [];
gammatr = [];
k       = 0;
flg     = 1;

m = 1;
p = 1;

A  = sys.A;
B = sys.B;
C = sys.C;

H     = C*Pi;

f_plot = figure(1);
try
    f_plot_p = get(0, "MonitorPositions");
    f_plot.Position = f_plot_p(2, :); % second display
    f_plot.WindowState = "maximized";
end

% N is required as a dummy variable
N = zeros(4*n+4*nu+p+m,2*n+2*nu);

while flg
    k = k +1;

    if mod(k,2) % Fix G and solve for X (or X and N)
        % Compute the Hinf-norm
        G       = double(G_array{end}{end});

        % Compute a balanced realisation for numerical conditioning
        [~,~,T] = balreal(sys_red(G));
        if abs(T(1)) > 1e3
            a = 9;
        end
        Sb      = T*S/T;
        Gb      = T*G;
        Lb      = L/T;
        Hb      = H/T;
        Ccalb   = [sys.C -Hb];

        X = sdpvar(n+nu);
        if form == 0
            % Finsler's form, N is an additional variable
            N  = sdpvar(4*n+4*nu+p+m,2*n+2*nu);
        end


    else % Fix X12 and X22 (or fix only N) and solve for G
        if form
            % Primal form
            X = double(X);
            X11 = sdpvar(n);
            X12 = X(1:n,n+1:end);
            X22 = X(n+1:end,n+1:end);
            X = [X11 X12;X12' X22];
        else
            % Finsler's form: only fix N whereas X remains free
            N  = double(N);
            X  = sdpvar(n+nu);
        end
        Gb      = sdpvar(nu,m);
    end

    LMI_             = @(g)LMI_primal_Finsler_LTI(A,B,Ccalb,Gb,X,N,Sb,Lb,g,form); 
    constrfnc        = @(g)[LMI_(g) <= -opts_BMI.eigtol*eye(size(LMI_(1)))];
    if mod(k,2) % Computation of the Hinfnorm
        gamma_tmp = hinfnorm(sys-sys_red(double(G)));
        constrfnc        = [LMI_(gamma_tmp*1.05) <= -opts_BMI.eigtol*eye(size(LMI_(1)))];
        % tic
        sol              = optimize(constrfnc,[],options);
        % toc
    else
        [gamma_array,ff] = bisection_BMIs(constrfnc, options, gamma_prev*1.05, 1e-4, opts_BMI.gtol, imax);
        gamma_prev       = gamma_array(end);
    end
    
    if gamma_prev == 0
        display('Bisection Infeasible')
        return
    end
    
    if ~mod(k,2)
        Gb       = double(Gb);
        G        = T\Gb;
        % Fb       = Sb-Gb*Lb;
        % sysr     = ss(Fb,Gb,Hb,0);
        sysr     = sys_red(double(G));

        gammatr(end+1) = hinfnorm(sys-sysr);

        if plot_figure
            figure(f_plot)
            subplot(121)
            bodemag(sys,sysr)
            grid on
            subplot(122)
            bodemag(sys-sysr)
            grid on
        end
        
        display(['Norm error system ' num2str(gammatr(end))])
        display(['gamma ' num2str(double(gamma_prev))])
        gamma_  = [gamma_ double(gamma_prev)];
        if k > 2
            display(['Stop crt: ' num2str((gamma_(end-1)-gamma_(end))/gamma_(end))])
            display(['Stop crt true: ' num2str((gammatr(end-1)-gammatr(end))/gammatr(end))])
        end
        if ((k > 2) && ((gamma_(end-1)-gamma_(end))/gamma_(end)<opts_BMI.gtol) && ((gammatr(end-1)-gammatr(end))/gammatr(end)<opts_BMI.gtol)) || (k >= opts_BMI.kmax)
            flg = 0;
        end

        G_array{end}{end+1} = double(G);
    end
end

X = double(X);
G  = double(G);