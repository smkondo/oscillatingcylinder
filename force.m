%drag and lift coefficient
%Computing the Strohaul number from the velocity output file
clc 
close all 

% Define relevant info for reading in data files
args1 = 'FinalData/stationaryforce.txt';  % STATIONARY CYLINDER
args2 = 'FinalData/moving1force.txt';  % OSCILLATIONG CYLINDER (frequency ratio R=0.5)
args3 = 'FinalData/moving2force.txt';% OSCILLATING CYLINDER (R=1.0)
args4 = 'FinalData/moving3force.txt'; % OSCILLATING CYLINDER (R=1.5)

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

figure(3)
plot(t,cd)

%% 
figure(2)
plot(om./100,abs(fcl),'o')
title('Frequency vs Fourier coefficient magnitude')
xlabel('St')
ylabel('magnitude')
xlim([-0.4 0.4])

%%
frequency1 = [0,1,1,2,2,3,3].*0.14286;
frequency2 = [0,0.5,0.5,1,1,1.5,1.5].*0.14286;

figure(1)
plot(frequency1,'o')
title('First 7 Dominant Frequencies')
xlabel('mode')
ylabel('Strouhal Number')

figure(2)
plot(frequency1,'o')
hold on
plot(frequency2,'o')
hold off
title('First 7 Dominant Frequencies')
xlabel('mode')
ylabel('Strouhal Number')
legend('Stationary','f_e/f_0 = 1','Location','northwest')
