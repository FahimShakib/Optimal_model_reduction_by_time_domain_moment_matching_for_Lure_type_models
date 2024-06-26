function [gammanew,feas] = bisection_BMIs(constrfnc, options, gammau, gammal,tol,imax)
% Bisection Method in MATLAB

% gammau: feasible gamma
% gammal: infeasible gamma
 
% Check feasibility of upper bound
constraints = constrfnc(gammau);
sol = optimize([constraints], [],options);

if sol.problem ~= 0
    gammanew    = 0;
    feas        = 0;
    return
end


for i = 1:imax
    gammar              = (gammau+gammal)/2;                                % Select new gamma as the mean
    constraints         = constrfnc(gammar);                                % Get constraints (LMIs)
    sol                 = optimize([constraints], [],options);              % Check feasibility
    
    if double(constrfnc(gammar))
       constr_satis = 1;
    else
        constr_satis = 0;
    end
    
    if constr_satis
        if sol.problem ~= 0
            % Infeasible or other problems --> Update the lower bound
            gammal  = gammar;
            feas(i) = 0;
        else
            % Feasible -> Decrease the upper bound
            gammau  = gammar;
            feas(i) = 1;
        end
    else
        % Constraints not satisfied --> Update the lower bound
        gammal  = gammar;
        feas(i) = 0;
    end

    gammanew(i) = gammar;
    if i > 1
        if abs((gammanew(i)-gammanew(i-1))/gammanew(i))<tol,break,end
    end
end

% Return the last feasible value
gammanew(end+1) = gammau;
feas(end+1)     = 1;
constraints     = constrfnc(gammau);                                % Get constraints (LMIs)
sol             = optimize([constraints], [],options);              % Check feasibility

