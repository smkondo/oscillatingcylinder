%vorticity
% Authors: Nathaniel Ruhl and Shinya Kondo
% This script creates the animation of the voriticty contour

%% Read in CFD data
dt = 0.005; % Time steps (s)
nt = 100;  % number of time steps
total_time = 0.5;   % sec, totaltime duration of experiment
xmin = -2.5; xmax = 17.5; ymin = -10; ymax = 10; 
nx = 500; ny = 500;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
t = linspace(0, total_time, nt);

% Define relevant info for reading in data files
args1 = {"FinalData/StationaryCylinderFinal/FFF-0",200};  % STATIONARY CYLINDER
args2 = {"FinalData/moving1final/FFF--0",300};   % OSCILLATIONG CYLINDER (frequency ratio R=0.5)
args3 = {"FinalData/moving2final/FFF--0", 300}; % OSCILLATING CYLINDER (R=1.0)
args4 = {"FinalData/moving3final/FFF--0", 300}; % OSCILLATING CYLINDER (R=1.5)

% Define the argument list for the flow that you want to use in this
% notebook
% arg_list = {x, y, t, args2{1}, args2{2}};

[XX, YY, om_matrix] = readvorticityData(x,y,t,args2{1},args2{2});

%% create the animation for the flow
cmax = max(om_matrix(:,1));
cmin = min(om_matrix(:,1));
%% Write animation video
video = VideoWriter('vorticitymoving1slow.avi');
video.FrameRate = 2; %set frames per second
open(video);

for i = 1:100
%     contourf(XX,YY,reshape(om_matrix(:,i),[nx,ny]),'edgecolor','none')
    contourf(XX,YY,reshape(om_matrix(:,i),[nx,ny]))
    colorbar()
    caxis([cmin cmax])
%     axis([-2.5 17.5 -4 4])
    axis([-2.5 5 -2 2])
    daspect([1 1 1])
%     colormap bone
%     title(strcat('Vorticity contour at t=',num2str(dt*i,'%.3f'),' (f_e/f_0 = 1)'))
    title(strcat('Vorticity contour at t=',num2str(dt*i,'%.3f'),' (stationary)'))
    xlabel('x')
    ylabel('y')
    pause(0.5)
    frame = getframe(gcf);
    writeVideo(video, frame)
end

close(video)

%%
contourf(XX,YY,reshape(om_matrix(:,i),[nx,ny]),'edgecolor','none')
colorbar()
caxis([-50 50])
axis([-2.5 17.5 -4 4])
daspect([1 1 1])