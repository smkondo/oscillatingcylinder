%% Optimization problem to compute the weights of each mode
%solving the optimization problem min alpha J(alpha) = ||S*V' -
%eigs*D(alpha)*Vandermonde_matrix||_F

%compute the matrices for optimization
Vand = fliplr(vander(diag(eigs)));

Dalpha = W\Sred*Vred(1:7,:)'/Vand;

Dalpha = diag(diag(Dalpha));
Dalpha2 = Dalpha./sum(Dalpha,'all');

%%
%computing the evolution by adding the first 7 modes
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

%% 
video = VideoWriter('uvel10.avi');
video.FrameRate = 10; %set frames per second
open(video);


for i = 1:length(t)
    contourf(XX,YY,real(reshape(mode_evolve(:,i),[nx,ny])),'edgecolor','none')
    colorbar()
    axis([-2.5 17.5 -3 3])
    daspect([1 1 1])
    title(int2str(i))
    pause(0.5)
    frame = getframe(gcf);
    writeVideo(video,frame)
end

close(video)