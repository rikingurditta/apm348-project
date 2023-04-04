clear;
close all;

%% parameters

% factor to modify time in simulation
T_FACTOR = 1;
% T_FACTOR = 1e-5;
% recovery/death rate
GAMMA = 1 / 9.32 / T_FACTOR;
% GAMMA = 1 / 30 / T_FACTOR;

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
% S_data = 2e6 - I_data - R_data;
x_data = [S_data, I_data, R_data];

%% find betas
% initial guess for parameters
beta0 = 1e-7;

[betas, grad_x, d_gamma] = FindBetas(t_data, x_data, GAMMA, beta0);
betas = smoothdata(betas);
betas = smoothdata(betas);

%% calculate betas over weeks
% skip = 7;
% dates_week = dates(1:skip:end);
% t_week = t_data(1:skip:end);
% x_data_week = x_data;
% for i=1:3
% %     x_data_week(:, i) = smoothdata(x_data_week(:, i));
% end
% x_data_week = x_data_week(1:skip:end, :);
% betas = FindBetas(t_week, x_data_week, GAMMA, beta0);
% betas = betas / 1.1;

%% simulate with discovered betas
[~, x_betas] = SIRBetas(x_data(1, :)', t_data, betas, GAMMA, t_data);
% [~, x_betas_week] = SIRBetas(x_data_week(1, :)', t_week, betas, GAMMA, t_week);

%% plot1
% subplot(3, 1, 1)
% plot(dates, x_data(:, 1), '-r', dates, x_data(:, 2), '-g', dates, x_data(:, 3), '-b')
% % plot(t_data, x_data(:, 2), '-g', t_data, x_data(:, 3), '-b')
% % plot(t_data, x_data(:, 2), '-g')
% 
% subplot(3, 1, 2)
% plot(dates, betas, '-k')
% % plot(dates_week, betas, '-k')
% 
% subplot(3, 1, 3)
% % figure()
% plot(dates, x_betas(:, 1), '-r', dates, x_betas(:, 2), '-g', dates, x_betas(:, 3), '-b')
% % plot(dates_week, x_betas_week(:, 1), '-r', dates_week, x_betas_week(:, 2), '-g', dates_week, x_betas_week(:, 3), '-b')
% % plot(t_data, x_betas(:, 2), '-g', t_data, x_betas(:, 3), '-b')
% % plot(t_data, x_betas(:, 2), '-g')

%% plot2

fig = figure();
compartments = {'S', 'I', 'R'};
for i=1:3
    subplot(3, 1, i)
    p = plot(dates, x_data(:, i), '-k', dates, x_betas(:, i), '-k');
%     p = plot(dates, x_data(:, i), '-k', dates_week, x_betas_week(:, i), '-k');
    v = [0, 0, 0];
    v(i) = 1;
    set(p(1), 'Color', v * 0.5);
    set(p(2), 'Color', v);
    title(strcat(compartments{i}, ' compartment vs time'))
    xlabel('time (days)'), ylabel('people')
    legend({'data', 'simulated'})
end

exportgraphics(fig, 'SIR.png', 'Resolution', 300)

%% plot3

fig = figure();
subplot(3, 1, 1)
plot(dates, betas, '-k')
title('β vs time')
xlabel('time (days)'), ylabel('β (1/people days)')
% plot(dates_week, betas, '-k')

subplot(3, 1, 2)
plot(dates, d_gamma * GAMMA / betas, '-k')
title('sensitivity of β to γ vs time')
xlabel('time (days)'), ylabel('S(β, γ)')

subplot(3, 1, 3)
p = plot(dates, x_data(:, 2), '-k', dates, x_betas(:, 2), '-k');
% p = plot(dates, x_data(:, 2), '-k', dates_week, x_betas_week(:, 2), '-k');
set(p(1), 'Color', [0, 0.5, 0]);
set(p(2), 'Color', [0, 1, 0]);
title('I compartment vs time')
xlabel('time (days)'), ylabel('I (people)')
legend({'data', 'simulated'})

exportgraphics(fig, 'beta.png', 'Resolution', 300)


%% fit both params
% p0 = [betas(1); GAMMA];
% [params, dp_dx] = FindParams(t_data, x_data, p0);
% [~, x_params] = SIRParams(x_data(1, :)', t_data, params, t_data);
% 
% figure()
% subplot(3, 1, 1)
% plot(dates, x_data(:, 1), '-r', dates, x_data(:, 2), '-g', dates, x_data(:, 3), '-b')
% 
% subplot(3, 1, 2)
% yyaxis left
% plot(dates, params(:, 1), '-k')
% yyaxis right
% plot(dates, params(:, 2), '-r')
% 
% subplot(3, 1, 3)
% % figure()
% plot(dates, x_params(:, 1), '-r', dates, x_params(:, 2), '-g', dates, x_params(:, 3), '-b')
% % plot(dates_week, x_betas_week(:, 1), '-r', dates_week, x_betas_week(:, 2), '-g', dates_week, x_betas_week(:, 3), '-b')
% % plot(t_data, x_betas(:, 2), '-g', t_data, x_betas(:, 3), '-b')
% % plot(t_data, x_betas(:, 2), '-g')
