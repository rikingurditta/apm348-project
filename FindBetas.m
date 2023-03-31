function betas = FindBetas(t_data, x_data, gamma, beta0)
%FindBetas
%   given SIR data, find beta values that take each timestep to the next
    betas = zeros(size(t_data));
    for i=2:length(t_data)
        i
        xprev = x_data(i-1,:)';
        tprev = t_data(i-1);
        tcurr = t_data(i);
        T_sim = linspace(tprev, tcurr, 2);
        f = @(beta) SIRSquareError(xprev, beta, gamma, T_sim, tcurr, x_data(i,:));
        sol = fminsearch(f, beta0);
        betas(i-1) = sol;
    end
end

