%DMDmodes.m
%computes the spatial DMD modes as well as the frequency and growth rates
%of each flow past a cylinder. Also, computes the reduced order model of
%the flow using enough modes to capture ~95% of the energy

% Authors: Nathaniel Ruhl and Shinya Kondo
% This file does dmd analysis on the flow past a stationary cylinder

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
args2 = {"FinalData/moving1final/FFF--0",300};   % OSCILLATIONG CYLINDER (frequency ratio R=1)
args3 = {"FinalData/moving2final/FFF--0", 300}; % OSCILLATING CYLINDER (R=1.5)
args4 = {"FinalData/moving3final/FFF--0", 300}; % OSCILLATING CYLINDER (R=0.5)

[XX, YY, v_matrix] = readData(x,y,t,args4{1},args4{2});


%% DMD analysis
X = v_matrix(:,1:nt-1);
X2 = v_matrix(:,2:nt);

r = 17; %use the first 7 modes for stationary cylinder and 17 for oscillating case
[Ured, Sred, Vred] = svds(X,r);

% Build the best-fit linear model that shows how POD modes evolve in time
Atilde = Ured'*X2*Vred/Sred;  % project A onto U bases vectors
[W, eigs] = eig(Atilde);  % compute e-vecs and evals of Atilde
Phi = X2*Vred/Sred*W;  % recover full-dimensional eigenflow

lambda = log(eigs)./2./pi./dt; %real component contains the growth rate. imaginary component contains frequency

figure(1)
plot(abs(imag(diag(lambda))),'o')
title('First 7 Dominant Frequencies')
xlabel('mode','FontSize',14)
ylabel('St','FontSize',14)

%%
plot(abs(imag(diag(lambda))),'o')
%% 
figure(2)
velmag = sqrt(reshape(Ured(1:nx*ny,4),[nx,ny]).^2 + reshape(Ured(nx*ny+1:2*nx*ny,4),[nx,ny]).^2);
% contourf(XX, YY, reshape(real(Phi(1:nx*ny,12)) , [nx, ny]),'edgecolor','none')
contourf(XX, YY, reshape(velmag , [nx, ny]),'edgecolor','none')
title('Velocity Magnitude of Fourth DMD mode')
axis([-2.5 17.5 -4 4])
daspect([1 1 1])
colorbar()

%% Optimization problem to compute the weights of each mode
%solving the optimization problem min alpha J(alpha) = ||S*V' -
%eigs*D(alpha)*Vandermonde_matrix||_F

%compute the matrices for optimization
Vand = fliplr(vander(diag(eigs)));

Dalpha = W\Sred*Vred(1:r,:)'/Vand;

Dalpha = diag(diag(Dalpha));
Dalpha2 = Dalpha./sum(Dalpha,'all');

%computing the evolution by adding the first r modes
dt = 0.005; tend = 0.5;
t = dt:dt:tend;

mode_evolve = zeros(nx*ny,length(t));

for i = 1:length(t)
    evol = zeros(nx*ny,1);
    for j = 1:r
        evol = evol + real(Dalpha(j,j).*Phi(1:nx*ny,j)).*real(exp(2*pi*lambda(j,j)*t(i)));
    end 
    mode_evolve(:,i) = evol;
end 

%% create the animation for the reduced-order model of the flow
cmax = max(mode_evolve(1:nx*ny,1));
cmin = min(mode_evolve(1:nx*ny,1));

%writing the animation video
video = VideoWriter('ROM.avi');
video.FrameRate = 10; %set frames per second
open(video);


for i = 1:length(t)
    contourf(XX,YY,real(reshape(mode_evolve(:,i),[nx,ny])),'edgecolor','none')
    colorbar()
    axis([-2.5 17.5 -4 4])
    daspect([1 1 1])
    title(strcat('ROM Streamwise velocity contour at t=',num2str(dt*i,'%.3f'),' (f_e/f_0 = 1)'))
    xlabel('x')
    ylabel('y')
%     caxis([cmin cmax])
    pause(0.5)
    frame = getframe(gcf);
    writeVideo(video,frame)
end

close(video)
