function L = SIRSquareError(x0, params, T, t_data, x_data)
%SIRSQUAREERROR squared error between simulated SIR model and given data
    % initialize error (L for Loss)
    L = 0;
    % predict values for given t values
    x_predict = SIRPredict(x0, params, T, t_data);
    % for each t, compute squared error between given x
    % aka 0.5 * ||x_data - x_predicted||^2
    for i=1:length(t_data)
        % d is a row vector
        d = x_data(i, :) - x_predict(i, :);
        L = L + d * d' / 2;
    end
end

