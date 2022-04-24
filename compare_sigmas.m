% Authors: Nathaniel Ruhl and Shinya Kondo
% This script compares the singular values

sigma_moving = readmatrix("20sigmas_moving.txt");
sigma_stationary = readmatrix("20sigmas_stationary.txt");

figure(1)
semilogy(sigma_stationary, 'o')
hold on
semilogy(sigma_moving, 'o')

legend(["Stationary Cylinder", "Oscillating Cylinder"])
xlabel("Mode number (k)")
ylabel('$\lambda_k/\Sigma_n\lambda_n$', 'Interpreter', 'latex')
title("POD Eigenvalues")

xlim([0,20])

%% Read in data

dt = 0.005; % Time steps (s)
nt = 100;  % number of time steps
total_time = 0.5;   % sec, totaltime duration of experiment
xmin = -10; xmax = 25; ymin = -15; ymax = 15; 
nx = 500; ny = 500;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
t = linspace(0, total_time, nt);

% Relevant information for reading the Ansys Data for both stationary and
% oscillating cylinders (ratios 1, 0.5, and 1.5)
stationary_arg_list = {x, y, t, "velocityfield5/FFF-2-01500-1", 600};
moving_1_args = {x, y, t, "movingfield1/FFF-Setup-Output-0", 200};
moving_1p5_args = {x, y, t, "movingfield2/FFF-Setup-Output-0", 141};
moving_p5_args = {x, y, t, "movingfield3/FFF-Setup-Output-0", 165};

[XX, YY, stationary_v_matrix] = readData(stationary_arg_list);
[~, ~, ratio_1_v_matrix] = readData(moving_1_args);
[~, ~, ratio_1p5_v_matrix] = readData(moving_1p5_args);
[~, ~, ratio_p5_v_matrix] = readData(moving_p5_args);

%% Calculate DMD eigenvalues for each case

[~, S_stationary, ~] = svd(stationary_v_matrix, 'econ');
[~, S_1, ~] = svd(ratio_1_v_matrix, 'econ');
[~, S_p5, ~] = svd(ratio_1p5_v_matrix, 'econ');
[~, S_1p5, ~] = svd(ratio_p5_v_matrix, 'econ');

writematrix(diag(S_stationary),"sigmas_staionary.txt")
writematrix(diag(S_1),"sigmas_ratio1-0.txt")
writematrix(diag(S_1p5),"sigmas_ratio1-5.txt")
writematrix(diag(S_p5),"sigmas_ratio0-5.txt")

%% Plot the Singular Values

figure(1)
semilogy(diag(S_stationary)/sum(diag(S_stationary)), 'o')
hold on
semilogy(diag(S_p5)/sum(diag(S_p5)), 'o')
hold on
semilogy(diag(S_1)/sum(diag(S_p5)), 'o')
hold on
semilogy(diag(S_1p5)/sum(diag(S_p5)), 'o')

legend(["Stationary Cylinder", "$f_e/f_0=0.5$", "$f_e/f_0=1.0$", "$f_e/f_0=1.5$"], 'Interpreter', 'latex')
xlabel("Mode number (k)")
ylabel('$\sigma_k/\Sigma_n\sigma_n$', 'Interpreter', 'latex')
title("SVD Singular Values")

xlim([0,20])
