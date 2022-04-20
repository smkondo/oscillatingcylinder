

%% 
% fn = append('movingvelocity/FFF-',int2str(101));
% Create new meshgrids
nx = 500; ny = 500; % number of discrete points

nt = 100;

v_matrix = zeros(2*nx*ny, nt);

xmin = -10; xmax = 25; ymin = -15; ymax = 15; 
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);

%create a mesh
[XX,YY] = meshgrid(x,y);

for i = 1:nt
    data = readmatrix(append('movingvelocity/FFF-0',int2str(100+i)));
    x_data = data(:,2); 
    y_data = data(:,3); 
    vx_data = data(:,4);
    vy_data = data(:,5);
    
    % Create interpolating function
    vx_func = scatteredInterpolant(x_data, y_data, vx_data);
    vy_func = scatteredInterpolant(x_data, y_data, vy_data);
    
    % Apply interpolating function to the meshgrid
    vx = vx_func(XX, YY);
    vy = vy_func(XX, YY);
    
    % Create data column
    vx_col = reshape(vx, [nx*ny, 1]);
    vy_col = reshape(vy, [nx*ny, 1]);
    v_col = [vx_col; vy_col];
    
    v_matrix(:, i) = v_col;
end


%% 
contourf(XX,YY,reshape(v_matrix(1:nx*ny,50),[nx,ny]))
axis([-2.5 17.5 -3 3])
daspect([1 1 1])

%% 
% Subtract temporal mean from data matrix
vx_mean = mean(v_matrix(1:nx*ny,:), 2);
vy_mean = mean(v_matrix(nx*ny+1:2*nx*ny,:), 2);
v_mean = [vx_mean; vy_mean];

v_matrix_fluc = v_matrix - v_mean;

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

%% SVD - write singular values to a text file
[U,S,V] = svds(v_matrix_fluc,20);

%% DMD analysis
X = v_matrix(:,1:nt-1);
X2 = v_matrix(:,2:nt);

% Find POD modes of X
[U, S, V] = svd(X, 'econ');

% Compute DMD using "reduced order matrices" (only use r modes)
% Phi are evecs of A
r = 7;
Ured = U(:,1:r);
Sred = S(1:r,1:r);
Vred = V(1:nt-1,1:r);
% Ured = U;
% Sred = S;
% Vred = V;

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