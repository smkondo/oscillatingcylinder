%Computing the Strohaul number from the velocity output file
clc 
clear all

data = readmatrix('vel_pt1.txt');
%%
t = data(:,2); vel = data(:,3);
v = vel/max(vel);
dt = t(2)-t(1); Lt = t(end) - t(1);
nt = length(t);
om = 2*pi/Lt*[-nt/2:1:nt/2-1];

ft = fftshift(fft(v));



figure(1)
plot(t,vel) 
%% 
figure(2)
plot(om,abs(ft),'o')