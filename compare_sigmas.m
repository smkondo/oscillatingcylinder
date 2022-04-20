% Authors: Nathaniel Ruhl and Shinya Kondo
% This script compares the first 20 SVD signular values
% (the values for sigma have already beean normalized)

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

xlim([0,14])