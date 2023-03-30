function [t, x] = SIR(x0, beta, gamma, T)
%SIR standard SIR model simulation
%   params = [beta, gamma]
    % TODO: make this a global variable or take in argument
    opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-2);
    [t, x] = ode15s(@(t, x) SIR_rhs(t, x, beta, gamma), T, x0, opts);
end

function dxdt = SIR_rhs(t, x, beta, gamma)
%     disp(t)
    S = x(1); I = x(2); R = x(3);
    
    dxdt = [ -beta * S * I;             ... dS/dt
              beta * S * I - gamma * I; ... dI/dt
                             gamma * I  ... dR/dt
           ];
end