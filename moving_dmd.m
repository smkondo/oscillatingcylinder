

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
