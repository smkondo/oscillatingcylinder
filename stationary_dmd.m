
% Authors: Nathaniel Ruhl and Shinya Kondo
% This file does dmd analysis on the flow past a stationary cylinder

% The script readData.md must be in the active working directory

%% Read in the CFD data
fn_string_format = "velocityfield3/FFF-1-00450-0";
start_int = 200;
nt = 100; % Number of time steps
total_time = 1;   % sec, totaltime duration of experiment
xmin = -10; xmax = 25; ymin = -15; ymax = 15; 
nx = 500; ny = 500;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);

[XX, YY, v_matrix] = readData(x, y, nt, fn_string_format, start_int);

% Subtract temporal mean from data matrix (for POD analysis)
vx_mean = mean(v_matrix(1:nx*ny,:), 2);
vy_mean = mean(v_matrix(nx*ny+1:2*nx*ny,:), 2);
v_mean = [vx_mean; vy_mean];

v_matrix_fluc = v_matrix - v_mean;

%% Test that reading in CFD data went successfully
contourf(XX,YY,reshape(v_matrix(1:nx*ny,50),[nx,ny]))
axis([-2.5 17.5 -3 3])
daspect([1 1 1])

%% Determine the dominant frequncy of the vortex street via DMD

% Full DMD analysis
X = v_matrix(:,1:nt-1);
X2 = v_matrix(:,2:nt);

% Find POD modes of X
[U, S, V] = svd(X, 'econ');

% Compute DMD (use all svd modes)
% Phi are evecs of A
Ured = U;
Sred = S;
Vred = V;

% Build the best-fit linear model that shows how POD modes evolve in time
Atilde = Ured'*X2*Vred/Sred;  % project A onto U bases vectors
[W, eigs] = eig(Atilde);  % compute e-vecs and evals of Atilde
Phi = X2*Vred/Sred*W;  % recover full-dimensional eigenflow
 
dt = 0.01;  % total_time/nt

lambda = 1/(2*pi*dt).*eigs;

frequency = imag(diag(lambda));
growth_rate = log(abs(diag(lambda)));

% Determine the dominant frequency of the vortex stree
[~, f0_index] = min(abs(growth_rate));

f0 = frequency(f0_index);

disp("The dominant frequency of the vortex street is: ")
disp(f0)

% Plot the eigenflows

figure(2)
plot(real(diag(lambda)), imag(diag(lambda)), "o")
xlabel("Re(\lambda_j)", "Interpreter", "latex")
ylabel("Im(\lambda_j)", "Interpreter", "latex")

%% Write animation video
video = VideoWriter('uvel.avi');
video.FrameRate = 10; %set frames per second
open(video);

for i = 1:100
    contourf(XX,YY,reshape(v_matrix(1:nx*ny,i),[nx,ny]),'edgecolor','none')
    frame = getframe(gcf);
    writeVideo(video, frame)
end

close(video)
%%
% contourf(XX,YY,reshape(v_matrix(1:nx*ny,50),[nx,ny]),'LineStyle','none')
%% SVD
[U,S,V] = svds(v_matrix_fluc,20);
%% 
figure(1)
semilogy(diag(S)./sum(S,'all'),'o')
title('Largest 20 modes')
xlabel('j')
ylabel('\sigma_j')
%%
velmag = sqrt(U(1:nx*ny,5).^2 + U(nx*ny+1:2*nx*ny,5).^2);
% velmag = sqrt(reshape(U(1:nx*ny,1),[nx,ny]).^2 + reshape(U(nx*ny+1:2*nx*ny,1),[nx,ny]).^2);

figure(2)
contourf(XX,YY,reshape(velmag,[nx,ny]),'LineStyle','none')
% contourf(XX,YY,reshape(U(nx*ny+1:2*nx*ny,1),[nx,ny]),'LineStyle','none')
% contourf(XX,YY,reshape(U(1:nx*ny,1),[nx,ny]),'LineStyle','none')
title('Velocity Magnitude of the fifth mode')
xlabel('x')
ylabel('y')
axis([-2.5 17.5 -3 3])
daspect([1 1 1])
colorbar()


%% DMD analysis
X = v_matrix(:,1:nt-1);
X2 = v_matrix(:,2:nt);

% Find POD modes of X
[U, S, V] = svd(X, 'econ');

% Compute DMD using "reduced order matrices" (only use r modes)
% Phi are evecs of A
r = 9;
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

lambda = 1/(2*pi*dt).*eigs;

% Plot the eigenflows

%% 
figure(4)
plot(real(diag(lambda)),'o')
title('lambda real')

figure(5)
plot(abs(imag(diag(lambda))),'o')
title('lambda imag')

figure(6)
plot(real(diag(lambda)), imag(diag(lambda)), "o")

%% 
figure(3)
title("Contour Plot of First Eigenflow of Reduced Operator");
velmagPhi1 = real(sqrt(Phi(1:nx*ny, 2).^2+Phi(nx*ny+1:2*nx*ny,2).^2));
contourf(XX, YY, reshape(velmagPhi1, [nx, ny]),'edgecolor','none')
title('third DMD mode')
axis([-2.5 17.5 -3 3])
daspect([1 1 1])
colorbar()

    
    


    