clear
clc

% Fourier cut-offs
lambda_x_plus = 2000*pi; % Stream-wise cut-off
lambda_z_plus = 250*pi; % Span-wise cut-off

% Corresponding cut-off wavenumbers
kx = round(Lx/lambda_x_plus); 
kz = round(Lz/lambda_z_plus);

% Load the list of names of the input data
Files_input = dir([dir_input,'vel_wall_*.mat']);

% Apply the large-scale filter to the input data
for n = 1:size(Files_input,1)

    load([dir_input,Files_input(n).name]);

    nx = size(u,1); % Number of grid points in x direction
    nz = size(u,2); % Number of grid points in z direction

    u_hat = fft2(u); % 2D Fast Fourier Transform, stream-wise velocity
    w_hat = fft2(w); % 2D Fast Fourier Transform, span-wise velocity

    u_hat((kx+2):(nx-kx),:) = 0; % Zero out the specified modes
    u_hat(:,(kz+2):(nz-kz)) = 0; 

    w_hat((kx+2):(nx-kx),:) = 0;
    w_hat(:,(kz+2):(nz-kz)) = 0; 

    uL = ifft2(u_hat); % 2D Inverse FFT
    wL = ifft2(w_hat); 

    duLdy = uL/y; % Large-scale wall shear stresses in x and z directions
    dwLdy = wL/y;

    tauL = sqrt(duLdy.^2+dwLdy.^2); % Resultant wall shear stress
    u_tauL_plus = sqrt(tauL); % Large-scale friction velocity

    % Direction of the large-scale friction velocity
    theta = -atan(dwLdy./duLdy); 

    % Name of files of the filtered field
    newname = strrep(Files_input(n).name,'vel_wall','u_tauL'); 

    % Store the output data
    save([dir_output,newname],'u_tauL_plus','theta','-v7.3'); 

end

% Load the list of names of the output data
Files_output = dir([dir_output,'u_tauL_*.mat']);

% Mean value of u_tauL_plus and theta

% Pre-allocate the matrix 
u_tauL_mean_loc = zeros([size(Files_output,1) 1]);
theta_mean_loc = zeros([size(Files_output,1) 1]);

for n = 1:size(Files_output,1)

    load([dir_output,Files_output(n).name]);

    u_tauL_mean_loc(n,:) = mean(u_tauL_plus,'all');

    theta_mean_loc(n,:) = mean(theta,'all');

end

u_tauL_plus_mean = mean(u_tauL_mean_loc,'all');
theta_mean = mean(theta_mean_loc,'all');

% Second moments of u_tauL_plus and theta

% Pre-allocate the matrix 
u_tauL_var_loc = zeros([size(Files_output,1) 1]);
theta_var_loc = zeros([size(Files_output,1) 1]);
theta_var_include_theta_mean_loc = zeros([size(Files_output,1) 1]);

for n = 1:size(Files_output,1)

    load([dir_output,Files_output(n).name]);

    u_tauL_var_loc(n,:) = ...
        mean(((u_tauL_plus/u_tauL_plus_mean - 1).^2),'all');

    theta_var_loc(n,:) = mean(((theta).^2),'all');

    theta_var_include_theta_mean_loc(n,:) = ...
        mean(((theta - theta_mean).^2),'all');

end

% Variances of u_tauL and theta
u_tauL_var = mean(u_tauL_var_loc,'all');
theta_var = mean(theta_var_loc,'all');
theta_var_include_theta_mean = ...
    mean(theta_var_include_theta_mean_loc,'all');

% Store the statistics of the output data
save([dir_output,'u_tauL_stats.mat'], ...
    'u_tauL_plus_mean','theta_mean','u_tauL_var','theta_var',...
    'theta_var_include_theta_mean');