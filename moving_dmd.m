% Authors: Nathaniel Ruhl and Shinya Kondo
% This script does DMD analysis on an oscillating cylinder

%% Read in the CFD data

dt = 0.005; % Time steps (s)
nt = 100;  % number of time steps
total_time = 0.5;   % sec, totaltime duration of experiment
xmin = -10; xmax = 25; ymin = -15; ymax = 15; 
nx = 500; ny = 500;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
t = linspace(0, total_time, nt);

% Relevant information for reading the oscillating cylinder data
moving_1_args = {x, y, t, "movingfield1/FFF-Setup-Output-0", 200};
moving_1p5_args = {x, y, t, "movingfield3/FFF-Setup-Output-0", 141};
moving_p5_args = {x, y, t, "movingfield3/FFF-Setup-Output-0", 165};

% [XX, YY, v_matrix] = readData(x, y, nt, 'movingvelocity/FFF-0', 100);
[XX, YY, v_matrix] = readData(moving_1_args);

% Subtract temporal mean from data matrix (for POD analysis)
vx_mean = mean(v_matrix(1:nx*ny,:), 2);
vy_mean = mean(v_matrix(nx*ny+1:2*nx*ny,:), 2);
v_mean = [vx_mean; vy_mean];

v_matrix_fluc = v_matrix - v_mean;

%% Test that reading in CFD data went successfully
contourf(XX,YY,reshape(v_matrix(1:nx*ny,50),[nx,ny]))
axis([-2.5 17.5 -3 3])
daspect([1 1 1])
colorbar()

%% Write animation video
video = VideoWriter('uvel3.avi');
video.FrameRate = 10; %set frames per second
open(video);

for i = 1:100
    contourf(XX,YY,reshape(v_matrix(1:nx*ny,i),[nx,ny]),'edgecolor','none')
    colorbar()
    ylim([-5 5])
%     daspect([1 1 1])
    title(int2str(i))
    frame = getframe(gcf);
    writeVideo(video, frame)
end

close(video)

%% DMD analysis
X = v_matrix(:,1:nt-1);
X2 = v_matrix(:,2:nt);

% Find POD modes of X
[U, S, V] = svd(X, 'econ');

% Compute DMD using "reduced order matrices" (only use r modes)
% Phi are evecs of A
r = 7;
% Ured = U(:,1:r);
% Sred = S(1:r,1:r);
% Vred = V(1:nt-1,1:r);
Ured = U;
Sred = S;
Vred = V;

% Build the best-fit linear model that shows how POD modes evolve in time
Atilde = Ured'*X2*Vred/Sred;  % project A onto U bases vectors
[W, eigs] = eig(Atilde);  % compute e-vecs and evals of Atilde
Phi = X2*Vred/Sred*W;  % recover full-dimensional eigenflow
%% 
dt = 0.01;

lambda = 1/dt.*eigs;

% Plot the eigenflows

%% 
figure(4)
plot(real(diag(lambda)),'o')
title('lambda real')
%% 
figure(5)
plot(abs(imag(diag(lambda))),'o')
title('lambda imag')

%% 
amp = diag(abs(lambda))'.*vecnorm(Phi,1);
% amp = vecnorm(Phi,1);
% amp = diag(abs(lambda));
%%
figure(6)
plot(abs(diag(imag(lambda))),amp/sum(amp),'o')
xlim([1 100])

%%
plot(diag(real(lambda)),diag(imag(lambda)),'o')
daspect([1 1 1])
%% 
figure(3)
title("Contour Plot of First Eigenflow of Reduced Operator");
velmagPhi1 = real(sqrt(Phi(1:nx*ny, 6).^2+Phi(nx*ny+1:2*nx*ny,6).^2));
contourf(XX, YY, reshape(velmagPhi1, [nx, ny]))
title('third DMD mode')
axis([-2.5 17.5 -3 3])
daspect([1 1 1])
colorbar()