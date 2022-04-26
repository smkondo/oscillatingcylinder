%animation.m
% Authors: Nathaniel Ruhl and Shinya Kondo
% This script creates the animation of the streamwise velocity contours

%% Read in CFD data
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

[XX, YY, v_matrix] = readData(x,y,t,args1{1},args1{2});

%% create the animation for the flow
cmax = max(v_matrix(1:nx*ny,1));
cmin = min(v_matrix(1:nx*ny,1));
%% Write animation video
video = VideoWriter('moving3final.avi');
video.FrameRate = 10; %set frames per second
open(video);

for i = 1:100
    contourf(XX,YY,reshape(v_matrix(1:nx*ny,i),[nx,ny]),'edgecolor','none')
    colorbar()
    caxis([cmin cmax])
    axis([-2.5 17.5 -4 4])
    daspect([1 1 1])
%     colormap bone
    title(strcat('Streamwise velocity contour at t=',num2str(dt*i,'%.3f'),' (f_e/f_0 = 0.5)'))
    xlabel('x')
    ylabel('y')
    pause(0.5)
    frame = getframe(gcf);
    writeVideo(video, frame)
end

close(video)


