%force.m
%extracting drag and lift coefficient from ANSYS data to compute the
%shedding frequency of each case

% Authors: Nathaniel Ruhl and Shinya Kondo
% This file does dmd analysis on the flow past a stationary cylinder

clc 
close all 

% Define relevant info for reading in data files
args1 = 'FinalData/stationaryforce.txt';  % STATIONARY CYLINDER
args2 = 'FinalData/moving1force.txt';  % OSCILLATIONG CYLINDER (frequency ratio R=0.5)
args3 = 'FinalData/moving2force.txt';% OSCILLATING CYLINDER (R=1.0)
args4 = 'FinalData/moving3force.txt'; % OSCILLATING CYLINDER (R=1.5)

data = readmatrix(args2);

t = data(:,2); cd = data(:,3); cl = data(:,4);
dt = t(2)-t(1); Lt = (t(end) - t(1));
nt = length(t);
om = 1/Lt.*[-nt/2:1:nt/2-1];

fcl = fftshift(fft(cl));

figure(1)
plot(t,cl./max(cl))
title('life coefficient vs time (f_e/f_0 = 1.0)')
xlabel('time')
ylabel('C_L')
xlim([2 2.5])
ylim([-1 1])

figure(2)
plot(om./100,abs(fcl),'o')
title('Frequency vs Fourier coefficient magnitude (f_e/f_0 = 1.0)')
xlabel('St')
ylabel('magnitude')
xlim([-0.4 0.4])
