%PODmodes.m
%computes the spatial POD modes of each flow past a cylinder (only done for
%the stationary cylinder case)

% Authors: Nathaniel Ruhl and Shinya Kondo


%% Read in CFD
dt = 0.005; % Time steps (s)
nt = 100;  % number of time steps
total_time = 0.5;   % sec, totaltime duration of experiment
xmin = -2.5; xmax = 17.5; ymin = -4; ymax = 4; 
nx = 500; ny = 500;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
t = linspace(0, total_time, nt);


% Define relevant info for reading in data files
args1 = {"FinalData/StationaryCylinderFinal/FFF-0",200};  % STATIONARY CYLINDER

[XX, YY, v_matrix] = readData(x,y,t,args1{1},args1{2});

% Subtract temporal mean from data matrix
vx_mean = mean(v_matrix(1:nx*ny,:), 2);
vy_mean = mean(v_matrix(nx*ny+1:2*nx*ny,:), 2);
v_mean = [vx_mean; vy_mean];

v_matrix_fluc = v_matrix - v_mean;


[U,S,V] = svds(v_matrix_fluc,20);

figure(1)
semilogy(diag(S)./sum(S,'all'),'o')
title('20 Largest Singular Values')
xlabel('j','FontSize',14)
ylabel('\sigma_j','FontSize',14)

velmag = sqrt(reshape(U(1:nx*ny,5),[nx,ny]).^2 + reshape(U(nx*ny+1:2*nx*ny,5),[nx,ny]).^2);
figure(2)
contourf(XX,YY,reshape(velmag,[nx,ny]),'LineStyle','none')
title('Velocity Magnitude of Fluctuations of Fifth POD Mode')
xlabel('x')
ylabel('y')
axis([-2.5 17.5 -4 4])
daspect([1 1 1])
colorbar()

