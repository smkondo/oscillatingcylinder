% This is a function to read the vorticity fields

% Function input: experiment_arg_list = {} contains the following 5 fields
% x_list: xmin:xmax:nx
% y_list: ymin:ymax:ny
% t_list: timelist
% fn_string_format: Format string of data file names before the last 3
% digits (Ex: "velocityfield3/FFF-1-00450-0")
% start_int: example: 400, if the first data file name ends with "451"


function [XX, YY, om_matrix] = readvorticityData(x_list, y_list, t_list, fn_string_format, start_int)
% [x_list, y_list, t_list, fn_string_format, start_int] = experiment_arg_list{:};
nx = length(x_list);
ny = length(y_list);
nt = length(t_list);
om_matrix = zeros(nx*ny, nt);

%create a mesh
[XX,YY] = meshgrid(x_list,y_list);

for i = 1:nt
    
    data = readmatrix(append(fn_string_format,int2str(start_int+i)));
    x_data = data(:,2); 
    y_data = data(:,3); 
    dvdx_data = data(:,6);
    dudy_data = data(:,7);
    
    % Create interpolating function
    dvdx_func = scatteredInterpolant(x_data, y_data, dvdx_data);
    dudy_func = scatteredInterpolant(x_data, y_data, dudy_data);
    
    % Apply interpolating function to the meshgrid
    dvdx= dvdx_func(XX, YY);
    dudy = dudy_func(XX, YY);
    
    % Create data column
    om_matrix(:,i) = reshape(dvdx, [nx*ny, 1]) - reshape(dudy, [nx*ny, 1]);
end

end