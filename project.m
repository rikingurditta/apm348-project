clear;
close all;

T_FACTOR = 1e-5;

GAMMA = 0.2 / T_FACTOR;

opts = detectImportOptions("covidtesting.csv")
opts.VariableTypes{8} = 'double';
opts.SelectedVariableNames = opts.SelectedVariableNames([1, 5, 6, 7, 8]);
all_data = readtable("covidtesting.csv", opts);
all_data.ReportedDate = datenum(all_data.ReportedDate);
all_data.Resolved(isnan(all_data.Resolved)) = 0;
all_data.Deaths(isnan(all_data.Deaths)) = 0;
all_data.Deaths_New_Methodology(isnan(all_data.Deaths_New_Methodology)) = 0;
all_data.Deaths = all_data.Deaths + all_data.Deaths_New_Methodology;
all_data = removevars(all_data, "Deaths_New_Methodology");
all_data = all_data(3:end, :)

t_data = all_data.ReportedDate * T_FACTOR;
I_data = all_data.ConfirmedPositive;
R_data = all_data.Resolved + all_data.Deaths;
S_data = ones(length(t_data), 1) * 14.57e6 - I_data - R_data;
x_data = [S_data, I_data, R_data];


% initial guess for parameters
beta0 = 40;

betas = zeros(size(t_data));
for i=2:length(t_data)
    i
    xprev = x_data(i-1,:)';
    tprev = t_data(i-1);
    tcurr = t_data(i);
    T_sim = linspace(tprev, tcurr, 2);
    f = @(beta) SIRSquareError(xprev, beta, GAMMA, T_sim, [tprev; tcurr], x_data(i-1:i,:));
    sol = fminsearch(f, beta0);
%     sol = GradientDescent(f, beta0, 0.000001, 10)
    betas(i-1) = sol(1);
end
% betas

subplot(3, 1, 1)
plot(t_data, x_data(:, 1), '-r', t_data, x_data(:, 2), '-g', t_data, x_data(:, 3), '-b')
subplot(3, 1, 2)
plot(t_data, betas, '-k')

% simulated with betas
[~, x_betas] = SIRBetas(x_data(1, :)', t_data, betas, GAMMA, t_data);
size(x_betas)
subplot(3, 1, 3)
% figure()
plot(t_data, x_betas(:, 1), '-r', t_data, x_betas(:, 2), '-g', t_data, x_betas(:, 3), '-b')
