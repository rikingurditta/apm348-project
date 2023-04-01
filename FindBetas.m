function [betas, grads] = FindBetas(t_data, x_data, gamma, beta0)
%FindBetas
%   given SIR data, find beta values that take each timestep to the next
    betas = zeros(size(t_data));
    grads = zeros(length(t_data), 3);
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
        g = @(x1) SIRSquareError(xprev, sol, gamma, T_sim, tcurr, x1);
        grads(i-1, :) = FiniteDiffGradient(g, xcurr);
    end
    fprintf('\n')
end

