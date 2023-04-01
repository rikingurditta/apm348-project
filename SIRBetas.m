function [t, x] = SIRBetas(x0, t_betas, betas, gamma, T)
%SIR
%   SIR model simulation with evolving beta values
    opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-1);
    [t, x] = ode15s(@(t, x) SIR_betas_rhs(t, x, t_betas, betas, gamma), T, x0, opts);
end

function dxdt = SIR_betas_rhs(t, x, t_betas, betas, gamma)
    S = x(1); I = x(2); R = x(3);

    % find beta for current t
    beta = interp1(t_betas, betas, t, 'previous');
    
    dxdt = [ -beta * S * I;             ... dS/dt
              beta * S * I - gamma * I; ... dI/dt
                             gamma * I  ... dR/dt
           ];
end