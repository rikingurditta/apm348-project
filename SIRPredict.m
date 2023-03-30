function x = SIRPredict(x0, beta, gamma, T, t)
%SIRPREDICT use SIR model to predict x = [S, I, R] at time t
%   use linear interpolation to make prediction (good if T has higher
%   resolution than t)
%   params = [beta, gamma]
    [t_sir, x_sir] = SIR(x0, beta, gamma, T);
    x = interp1(t_sir, x_sir, t);
%     disp(beta)
end

