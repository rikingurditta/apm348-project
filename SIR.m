function [t, x] = SIR(x0, params, T)
%SIR standard SIR model simulation
%   params = [beta, gamma]
    param.beta = params(1);
    param.gamma = params(2);
    [t, x] = ode45(@(t, x) SIR_rhs(t, x, param), T, x0);
end

function dxdt = SIR_rhs(t, x, param)
    S = x(1); I = x(2); R = x(3);
    
    dxdt = [ -param.beta * S * I;                   ... dS/dt
              param.beta * S * I - param.gamma * I; ... dI/dt
                                   param.gamma * I  ... dR/dt
           ];
end