% Authors: Nathaniel Ruhl and Shinya Kondo
% This script creates an animation of the leading eigenflows

%% Read in data for the stationary cylinder

fn_string_format = "velocityfield3/FFF-1-00450-0";
start_int = 200;
nt = 100; % Number of time steps
total_time = 1;   % sec, totaltime duration of experiment
t = linspace(0, total_time, nt);
xmin = -10; xmax = 25; ymin = -15; ymax = 15; 
nx = 500; ny = 500;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);

stationary_arg_list = {x, y, t, fn_string_format, start_int};

[XX, YY, v_matrix] = readData(stationary_arg_list);

%% Calculate non-reduced DMD
[Phi, Lambda] = calcDMD(v_matrix, t);

%% Identify n'th leading Eigenflow and eigenvalue for plotting

[phi1, lambda1] = getNthEigenflow(Phi, Lambda, 1);

%% Plot the eigenflow

% Velocity components (may be complex)
phi_matrix_x = reshape(phi1(1:nx*ny), [nx, ny]);
phi_matrix_y = reshape(phi1(nx*ny+1:2*nx*ny), [nx, ny]);

% Velocity magnitude field (may be complex)
phi_matrix_mag = sqrt(phi_matrix_x.^2+phi_matrix_y.^2);

figure()
contourf(XX, YY, real(phi_matrix))

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

%% Functions

% Function to identify the N'th eigenflow
% Returns 2 matrices, evecs and evals of linear operator
function [Phi, Lambda] = calcDMD(v_matrix, t)
dt = t(2) - t(1);
nt = length(t);
%DMD Analysis (Not reduced)
X = v_matrix(:,1:nt-1);
X2 = v_matrix(:,2:nt);

% Find POD modes of X
[U, S, V] = svd(X, 'econ');

% Compute DMD using "reduced order matrices" (only use r modes)
% Phi are evecs of A
% r = 9;
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

Lambda = 1/(2*pi*dt).*eigs;

end

% Function to identify the N'th eigenflow and complex eigenvector
function [phi_n, lambda_n] = getNthEigenflow(Phi, Lambda, n)
lambda = abs(diag(Lambda));
lambda_imag = diag(Lambda);

% identify the largest eigenvalues
[~, sortIdx] = sort(lambda,'descend');

phi_n = Phi(:,sortIdx(n));
lambda_n = lambda_imag(sortIdx(n));

end

