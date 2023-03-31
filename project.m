clear;
close all;

%% parameters
% factor to modify time in simulation
T_FACTOR = 1;
% T_FACTOR = 1e-5;
% recovery/death rate
GAMMA = 1 / 30 / T_FACTOR;
% how many times to repeat beta finding procedure
NUM_PERTURBATIONS = 1;

%% import data
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
all_data = all_data(3:end, :);
all_data = all_data(600:900, :);

t_data = (all_data.ReportedDate - all_data.ReportedDate(1)) * T_FACTOR;
I_data = all_data.ConfirmedPositive;
R_data = all_data.Resolved + all_data.Deaths;
S_data = 14.57e6 - I_data - R_data;
% S_data = 2e6 - I_data - R_data;
x_data = [S_data, I_data, R_data];

%% find betas
% initial guess for parameters
beta0 = 1e-7;
betas = zeros(size(t_data));

for i=1:NUM_PERTURBATIONS
    i
%     x_perturbed = x_data + rand(size(x_data)) * 1000;
    x_perturbed = x_data;
    betas = betas + FindBetas(t_data, x_perturbed, GAMMA, beta0) / NUM_PERTURBATIONS;
end
betas = betas / T_FACTOR;

%% simulate with discovered betas
[~, x_betas] = SIRBetas(x_data(1, :)', t_data, betas, GAMMA, t_data);

%% plot
subplot(3, 1, 1)
plot(t_data, x_data(:, 1), '-r', t_data, x_data(:, 2), '-g', t_data, x_data(:, 3), '-b')

subplot(3, 1, 2)
plot(t_data, betas, '-k')

subplot(3, 1, 3)
% figure()
plot(t_data, x_betas(:, 1), '-r', t_data, x_betas(:, 2), '-g', t_data, x_betas(:, 3), '-b')
