%drag and lift coefficient
%Computing the Strohaul number from the velocity output file
clc 
close all 

data = readmatrix('force3.txt');
%%
t = data(:,2); cd = data(:,3); cl = data(:,4);
dt = t(2)-t(1); Lt = (t(end) - t(1));
nt = length(t);
om = 1/Lt*[-nt/2:1:nt/2-1];

% fcd = fftshift(fft(cd));
fcl = fftshift(fft(cl));

figure(2)
plot(t,cl./max(cl))
% hold on
% plot(t,sin(2*pi*14.4.*t))
% hold off

%%
plot(cl./max(cl),'o')
%% 
figure(3)
plot(om,abs(fcl),'o')