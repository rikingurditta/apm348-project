function L = SIRSquareError(x0, beta, gamma, T, t_data, x_data)
%SIRSQUAREERROR
%   squared error between simulated SIR model and given data
    % initialize error (L for Loss)
    L = 0;
    % predict values for given t values
    x_predict = SIRPredict(x0, beta, gamma, T, t_data);
    % for each t, compute squared error between given x
    % aka 0.5 * ||x_data - x_predicted||^2
    parfor i=1:length(t_data)
%         i
        % d is a row vector
        d = x_data(i, 1) - x_predict(i, 1);
        L = L + d * d';
    end
%     disp(num2str(beta))
%     disp(num2str(L))
end

