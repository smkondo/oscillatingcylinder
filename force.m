%drag and lift coefficient
%Computing the Strohaul number from the velocity output file
clc 
close all 

% Define relevant info for reading in data files
args1 = 'FinalData/stationaryforce.txt';  % STATIONARY CYLINDER
args2 = 'FinalData/moving1force.txt';  % OSCILLATIONG CYLINDER (frequency ratio R=0.5)
args3 = 'FinalData/moving1force.txt';% OSCILLATING CYLINDER (R=1.0)
args4 = 'FinalData/moving1force.txt'; % OSCILLATING CYLINDER (R=1.5)

data = readmatrix(args2);
%%
t = data(:,2); cd = data(:,3); cl = data(:,4);
dt = t(2)-t(1); Lt = (t(end) - t(1));
nt = length(t);
om = 1/Lt.*[-nt/2:1:nt/2-1];

fcl = fftshift(fft(cl));

figure(1)
plot(t,cl./max(cl))
title('life coefficient vs time')
xlabel('time')
ylabel('C_L')

%% 
figure(2)
plot(om,abs(fcl),'o')
title('Frequency vs Fourier coefficient magnitude')
xlabel('Frequency')
ylabel('magnitude')