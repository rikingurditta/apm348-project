function [betas, grad_x, d_gamma] = FindBetas(t_data, x_data, gamma, beta0)
%FindBetas
%   given SIR data, find beta values that take each timestep to the next
    betas = zeros(size(t_data));
    grad_x = zeros(length(t_data), 3);
    d_gamma = zeros(size(t_data));
    h = 0.0001;
    for i=2:length(t_data)
        if mod(i, floor(length(t_data) / 50)) == 0
            fprintf('*');
        end
        xprev = x_data(i-1,:)';
        xcurr = x_data(i,:);
        tprev = t_data(i-1);
        tcurr = t_data(i);
        T_sim = linspace(tprev, tcurr, 2);
        f = @(beta) SIRSquareError(xprev, beta, gamma, T_sim, tcurr, xcurr);
        sol = fminsearch(f, beta0);
        betas(i-1) = sol;
        beta_of_x1 = @(x1) fminsearch(@(beta) SIRSquareError(xprev, beta, gamma, T_sim, tcurr, x1), beta0);
        grad_x(i-1, :) = FiniteDiffGradient(beta_of_x1, xcurr, h);
        beta_of_gamma = @(y) fminsearch(@(beta) SIRSquareError(xprev, beta, y, T_sim, tcurr, xcurr), beta0);
        d_gamma(i-1) = FiniteDiffGradient(beta_of_gamma, gamma, h);
%         disp([beta_of_gamma(gamma + h) beta_of_gamma(gamma - h)])
    end
    fprintf('\n')
end


