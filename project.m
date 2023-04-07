clear;
close all;

%% parameters

% factor to modify time in simulation
T_FACTOR = 1;

% recovery/death rate
GAMMA = 1 / 9.32 / T_FACTOR;

%% import data
opts = detectImportOptions("covidtesting.csv");
opts.VariableTypes{8} = 'double';
opts.SelectedVariableNames = opts.SelectedVariableNames([1, 5, 6, 7, 8]);
all_data = readtable("covidtesting.csv", opts);
all_data.Resolved(isnan(all_data.Resolved)) = 0;
all_data.Deaths(isnan(all_data.Deaths)) = 0;
all_data.Deaths_New_Methodology(isnan(all_data.Deaths_New_Methodology)) = 0;
all_data.Deaths = all_data.Deaths + all_data.Deaths_New_Methodology;
all_data = removevars(all_data, "Deaths_New_Methodology");
all_data = all_data(150:1030, :);
% all_data = all_data(400:900, :);

dates = all_data.ReportedDate;

t_data = datenum(dates);
t_data = t_data - t_data(1);
t_data = t_data * T_FACTOR;
I_data = all_data.ConfirmedPositive;
R_data = all_data.Resolved + all_data.Deaths;
S_data = 14223942 - I_data - R_data;
x_data = [S_data, I_data, R_data];

%% find betas
% initial guess for parameters
beta0 = 1e-7;

[betas, grad_x, d_gamma] = FindBetas(t_data, x_data, GAMMA, beta0);
betas = smoothdata(betas);
betas = smoothdata(betas);

%% simulate with discovered betas
[~, x_betas] = SIRBetas(x_data(1, :)', t_data, betas, GAMMA, t_data);

%% plot1

fig = figure();
compartments = {'S', 'I', 'R'};
for i=1:3
    subplot(3, 1, i)
    p = plot(dates, x_data(:, i), '-k', dates, x_betas(:, i), '-k');
    v = [0, 0, 0];
    v(i) = 1;
    set(p(1), 'Color', v * 0.5);
    set(p(2), 'Color', v);
    title(strcat(compartments{i}, ' compartment vs time'))
    xlabel('time (days)'), ylabel('people')
    legend({'data', 'simulated'})
end

exportgraphics(fig, 'SIR.png', 'Resolution', 300)

%% plot2

fig = figure();
subplot(3, 1, 1)
plot(dates, betas, '-k')
title('β vs time')
xlabel('time (days)'), ylabel('β (1/people days)')

subplot(3, 1, 2)

L_beta = zeros(size(betas));
for i=1:length(betas)-1
    T_sim = linspace(t_data(i), t_data(i+1), 2);
    L_beta(i) = SIRSquareError(x_data(i, :)', betas(i), GAMMA, T_sim, T_sim, x_data(i:i+1, :));
end
plot(dates, L_beta, '-k');
title('squared error of best fit β vs time')
xlabel('time (days)'), ylabel('L (people^2)')

% plot(dates, d_gamma * GAMMA / betas, '-k')
% title('sensitivity of β to γ vs time')
% xlabel('time (days)'), ylabel('S(β, γ)')

subplot(3, 1, 3)
p = plot(dates, x_data(:, 2), '-k', dates, x_betas(:, 2), '-k');
set(p(1), 'Color', [0, 0.5, 0]);
set(p(2), 'Color', [0, 1, 0]);
title('I compartment vs time')
xlabel('time (days)'), ylabel('I (people)')
legend({'data', 'simulated'})

exportgraphics(fig, 'beta.png', 'Resolution', 300)
