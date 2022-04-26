%DMDmodes.m
%computes the frequency and growth rates of the reduced SVD corresponding
%to approximately 95% of the total energy content


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
args2 = {"FinalData/moving1final/FFF--0",300};   % OSCILLATIONG CYLINDER (frequency ratio R=0.5)
args3 = {"FinalData/moving2final/FFF--0", 300}; % OSCILLATING CYLINDER (R=1.0)
args4 = {"FinalData/moving3final/FFF--0", 300}; % OSCILLATING CYLINDER (R=1.5)

[XX, YY, v_matrix] = readData(x,y,t,args2{1},args2{2});


%% DMD analysis
X = v_matrix(:,1:nt-1);
X2 = v_matrix(:,2:nt);

%% compute the first 20 singular values
[~,ss,~] = svds(X,20);

figure(1)
semilogy(diag(ss)./sum(ss,'all'),'o')
title('Largest 20 modes')
xlabel('j')
ylabel('\sigma_j')

%remove the first mode corresponding to the mean flow
singval = diag(ss(2:end,2:end));

%compute the ratio between the sum of the first n modes to the total from
%all 19 modes
ratio = zeros(19,1);

for i = 1:19
    ratio(i) = sum(singval(1:i))./sum(singval);
end

figure(2)
semilogy(diag(ss)./sum(diag(ss)),'o')
title('first 19 singular values')

figure(3)
plot(ratio,'o')
title('the importance of the first n modes')

%%
r = 17; %using the first 15 modes (note this may change)
[Ured, Sred, Vred] = svds(X,r);

% Build the best-fit linear model that shows how POD modes evolve in time
Atilde = Ured'*X2*Vred/Sred;  % project A onto U bases vectors
[W, eigs] = eig(Atilde);  % compute e-vecs and evals of Atilde
Phi = X2*Vred/Sred*W;  % recover full-dimensional eigenflow

lambda = log(eigs)./2./pi./dt;

%% 
figure(4)
plot(real(diag(lambda)),'o')
title('lambda real')

%% 
lambda2 = flip(diag(lambda(6:13,6:13)));
figure(5)
plot(abs(imag(lambda2)),'o')
% plot(abs(imag(diag(lambda))),'o')
title('First 7 Dominant Frequencies')
xlabel('mode')
ylabel('Frequency (Hz)')
xlim([1 7])

%% 
figure(5)
velmag = sqrt(reshape(Ured(1:nx*ny,4),[nx,ny]).^2 + reshape(Ured(nx*ny+1:2*nx*ny,4),[nx,ny]).^2);
% contourf(XX, YY, reshape(real(Phi(1:nx*ny,12)) , [nx, ny]),'edgecolor','none')
contourf(XX, YY, reshape(velmag , [nx, ny]),'edgecolor','none')
title('Velocity Magnitude of Second DMD mode')
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

%%
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

%% create the animation for the flow
cmax = max(mode_evolve(1:nx*ny,1));
cmin = min(mode_evolve(1:nx*ny,1));
%% 
%writing the animation video
video = VideoWriter('uvel12.avi');
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

%% 
amp = vecnorm(Phi(1:nx*ny,:),1);

%%
plotlambda = abs(diag(imag(lambda)));
% plotlambda2 = [plotlambda(1:39);plotlambda(42:end)]./100;
plotamp = amp;
% plotamp = [amp(1,1:39), amp(1,42:end)];
figure(7)
% plot(plotlambda,plotamp./sum(plotamp),'o')
plot(plotlambda,plotamp,'o')
title('normalized magnitudes of modes at each frequency')