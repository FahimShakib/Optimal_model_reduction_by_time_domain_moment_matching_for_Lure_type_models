function [sysr,G,gamma_,X1,X2,gammatr,G_array] = BMI_solver(sys,n,nu,G_array,S,L,gamma_prev,sys_red,Pi,plot_figure,form,opts_BMI)

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
q = 2;

A  = sys.A;
Bu = sys.B(:,1:m);
Bw = sys.B(:,m+1:end);
Cy = sys.C(1:p,:);
Cz = sys.C(p+1:end,:);

H     = sys.C*Pi;
Hzeta = Cz*Pi;

f_plot = figure(1);
try
    f_plot_p = get(0, "MonitorPositions");
    f_plot.Position = f_plot_p(2, :); % second display
    f_plot.WindowState = "maximized";
end

% N is required as a dummy variable
N = zeros(4*n+8*nu+4*q+p+m,2*n+4*nu);

while flg
    k = k +1;

    if mod(k,2) % Fix G and solve for X1 and X2 (or X1, X2, and N)
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
        Hzetab  = Hzeta/T;
        Ccalb   = [sys.C -Hb];

        X1 = sdpvar(n+nu);
        X2 = sdpvar(nu);
        if form == 0
            % Finsler's form, N is an additional variable
            N  = sdpvar(4*n+8*nu+4*q+p+m,2*n+4*nu);
        end
        

    else % Fix X1 and X2 (or fix only N) and solve for G
        if form
            % Primal form
            X1 = double(X1);
            X2 = double(X2);
        else
            % Finsler's form: only fix N whereas X1 and X2 remain free
            N  = double(N);
            X1 = sdpvar(n+nu);
            X2 = sdpvar(nu);
        end
        Gb      = sdpvar(nu,m+q);
    end

    LMI_             = @(g)LMI_primal_Finsler(A,Bu,Bw,Ccalb,Hzetab,Gb,X1,X2,N,Sb,Lb,g,form); 
    constrfnc        = @(g)[LMI_(g) <= -opts_BMI.eigtol*eye(size(LMI_(1)))];
    if mod(k,2) % Computation of the Hinfnorm
        gamma_tmp        = hinfnorm(sys-sys_red(double(G)));
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

X1 = double(X1);
X2 = double(X2);
G  = double(G);