function [t, x] = SIRBetas(x0, t_betas, betas, gamma, T)
%SIR standard SIR model simulation
    opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-2);
    [t, x] = ode15s(@(t, x) SIR_betas_rhs(t, x, t_betas, betas, gamma), T, x0, opts);
end

function dxdt = SIR_betas_rhs(t, x, t_betas, betas, gamma)
    S = x(1); I = x(2); R = x(3);

    beta = interp1(t_betas, betas, t);
    
    dxdt = [ -beta * S * I;             ... dS/dt
              beta * S * I - gamma * I; ... dI/dt
                             gamma * I  ... dR/dt
           ];
end