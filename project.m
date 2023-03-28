clear;
close all;

x0 = [850, 150, 0];
tf = 50;
T = linspace(0, tf, 100);

% generate data points for certain params
params_data = [0.001; 0.1];
t_data = linspace(0, tf, 100)';
x_data = SIRPredict(x0, params_data, T, t_data);

% initial guess for parameters
p0 = [0.001; 0.1];
% the function that we want to minimize
f = @(params) SIRSquareError(x0, params, T, t_data, x_data);
% do gradient descent
p = GradientDescent(f, p0, 0.00001, 10)
