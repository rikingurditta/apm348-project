clear;
close all;

x0 = [850, 150, 0];
tf = 50;
T = linspace(0, tf, 100);

% generate data points for certain params
params_data = [0.001; 0.1];
t_data = T(:)';
x_data = SIRPredict(x0, params_data, T, t_data);

% initial guess for parameters
p0 = [0.01; 0];
% the function that we want to minimize
% f = @(params) SIRSquareError(x0, params, T, t_data, x_data);
% do gradient descent
% p = GradientDescent(f, p0, 0.000000001, 100)

betas = zeros(size(t_data));
for i=2:length(t_data)
    xprev = x_data(i-1,:)';
    tprev = t_data(i-1);
    tcurr = t_data(i);
    T_sim = linspace(tprev, tcurr, 50);
    f = @(params) SIRSquareError(xprev, params, T_sim, [tprev; tcurr], x_data(i-1:i,:));
    sol = fminsearch(f, p0);
    betas(i-1) = sol(1);
end
betas

subplot(2, 1, 1)
plot(t_data, x_data(:, 1), '-r', t_data, x_data(:, 2), '-g', t_data, x_data(:, 3), '-b')
subplot(2, 1, 2)
plot(t_data, betas, '-k')
