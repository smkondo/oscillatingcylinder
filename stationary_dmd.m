
% Authors: Nathaniel Ruhl and Shinya Kondo
% This file does dmd analysis on the flow past a stationary cylinder

%% Read in the CFD data
% This list is an output of the python script "create_fn_list.py"
% fn_list = {'FFF-1-00250-0351', 'FFF-1-00250-0352', 'FFF-1-00250-0353', 'FFF-1-00250-0354', 'FFF-1-00250-0355', 'FFF-1-00250-0356', 'FFF-1-00250-0357', 'FFF-1-00250-0358', 'FFF-1-00250-0359', 'FFF-1-00250-0360', 'FFF-1-00250-0361', 'FFF-1-00250-0362', 'FFF-1-00250-0363', 'FFF-1-00250-0364', 'FFF-1-00250-0365', 'FFF-1-00250-0366', 'FFF-1-00250-0367', 'FFF-1-00250-0368', 'FFF-1-00250-0369', 'FFF-1-00250-0370', 'FFF-1-00250-0371', 'FFF-1-00250-0372', 'FFF-1-00250-0373', 'FFF-1-00250-0374', 'FFF-1-00250-0375', 'FFF-1-00250-0376', 'FFF-1-00250-0377', 'FFF-1-00250-0378', 'FFF-1-00250-0379', 'FFF-1-00250-0380', 'FFF-1-00250-0381', 'FFF-1-00250-0382', 'FFF-1-00250-0383', 'FFF-1-00250-0384', 'FFF-1-00250-0385', 'FFF-1-00250-0386', 'FFF-1-00250-0387', 'FFF-1-00250-0388', 'FFF-1-00250-0389', 'FFF-1-00250-0390', 'FFF-1-00250-0391', 'FFF-1-00250-0392', 'FFF-1-00250-0393', 'FFF-1-00250-0394', 'FFF-1-00250-0395', 'FFF-1-00250-0396', 'FFF-1-00250-0397', 'FFF-1-00250-0398', 'FFF-1-00250-0399', 'FFF-1-00250-0400', 'FFF-1-00250-0401', 'FFF-1-00250-0402', 'FFF-1-00250-0403', 'FFF-1-00250-0404', 'FFF-1-00250-0405', 'FFF-1-00250-0406', 'FFF-1-00250-0407', 'FFF-1-00250-0408', 'FFF-1-00250-0409', 'FFF-1-00250-0410', 'FFF-1-00250-0411', 'FFF-1-00250-0412', 'FFF-1-00250-0413', 'FFF-1-00250-0414', 'FFF-1-00250-0415', 'FFF-1-00250-0416', 'FFF-1-00250-0417', 'FFF-1-00250-0418', 'FFF-1-00250-0419', 'FFF-1-00250-0420', 'FFF-1-00250-0421', 'FFF-1-00250-0422', 'FFF-1-00250-0423', 'FFF-1-00250-0424', 'FFF-1-00250-0425', 'FFF-1-00250-0426', 'FFF-1-00250-0427', 'FFF-1-00250-0428', 'FFF-1-00250-0429', 'FFF-1-00250-0430', 'FFF-1-00250-0431', 'FFF-1-00250-0432', 'FFF-1-00250-0433', 'FFF-1-00250-0434', 'FFF-1-00250-0435', 'FFF-1-00250-0436', 'FFF-1-00250-0437', 'FFF-1-00250-0438', 'FFF-1-00250-0439', 'FFF-1-00250-0440', 'FFF-1-00250-0441', 'FFF-1-00250-0442', 'FFF-1-00250-0443', 'FFF-1-00250-0444', 'FFF-1-00250-0445', 'FFF-1-00250-0446', 'FFF-1-00250-0447', 'FFF-1-00250-0448', 'FFF-1-00250-0449', 'FFF-1-00250-0450'};

% Create new meshgrids
nx = 500; ny = 500; % number of discrete points

nt = 100;

v_matrix = zeros(2*nx*ny, nt);

xmin = -10; xmax = 25; ymin = -15; ymax = 15; 
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);

%create a mesh
[XX,YY] = meshgrid(x,y);

% Assemble data matrix
for i = 1:nt
    data = readmatrix(append("velocityfield3/FFF-1-00450-0",int2str(200+i)));
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

figure(5)
plot(abs(imag(diag(lambda))),'o')
title('lambda imag')
%% 
figure(3)
title("Contour Plot of First Eigenflow of Reduced Operator");
velmagPhi1 = real(sqrt(Phi(1:nx*ny, 2).^2+Phi(nx*ny+1:2*nx*ny,2).^2));
contourf(XX, YY, reshape(velmagPhi1, [nx, ny]),'edgecolor','none')
title('third DMD mode')
axis([-2.5 17.5 -3 3])
daspect([1 1 1])
colorbar()

    
    


    