% This is a function to read the velocity fields

% Function input: experiment_arg_list = {} contains the following 5 fields
% x_list: xmin:xmax:nx
% y_list: ymin:ymax:ny
% t_list: timelist
% fn_string_format: Format string of data file names before the last 3
% digits (Ex: "velocityfield3/FFF-1-00450-0")
% start_int: example: 400, if the first data file name ends with "451"


function [XX, YY, v_matrix] = readData(x_list, y_list, t_list, fn_string_format, start_int)
% [x_list, y_list, t_list, fn_string_format, start_int] = experiment_arg_list{:};
nx = length(x_list);
ny = length(y_list);
nt = length(t_list);
v_matrix = zeros(2*nx*ny, nt);

%create a mesh
[XX,YY] = meshgrid(x_list,y_list);

for i = 1:nt
    
    data = readmatrix(append(fn_string_format,int2str(start_int+i)));
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

v_matrix = v_matrix./max(v_matrix,[],'all');

end

