%%
clear
clc

% User input parameters

lines = readlines('user_input_1D.txt');
InputParameters = extractAfter(lines(:),',');

% Size of the input friction field in wall unit
% (The test input has the size of Lx = 8000pi and Lz = 3000pi.)
Lx = str2double(InputParameters(1)); % Length of the stream-wise domain
Lz = str2double(InputParameters(2)); % Length of the span-wise domain

% Choose the maximum spatial ranges in wall unit
% (xmax=0 & zmax=0 for 1D output field)
xmax = str2double(InputParameters(3)); %xmax = {0, 204}
zmax = str2double(InputParameters(4)); %zmax = {0, 120}

ymin = str2double(InputParameters(5)); %y = {0, 100}
ymax = str2double(InputParameters(6)); 

% Is this run a test? (true/false)
test = str2num(InputParameters(7));

% Starting loop (starting from 1 for a new run)
starting_loop = str2double(InputParameters(8));
ending_loop = str2double(InputParameters(9));

% The following information of directories should remain unchanged for the
% test, unless the different file separator "/" is used instead of "\". 

% The directory of the tool (current folder)
workdir = pwd;

% The directory of the input friction field
dir_input = [workdir,convertStringsToChars(InputParameters(10))];

% The directory of the reference wall-normal grid
dir_y_ref = [workdir,convertStringsToChars(InputParameters(11))];

% The directory of the 1D universal velocity field
dir_vel_tilde = [workdir,convertStringsToChars(InputParameters(12))];

% The folder that will store the 1D output of the synthetic model
dir_output = [workdir,convertStringsToChars(InputParameters(13))];

%% Data converter

% Fourier cut-offs
lambda_x_plus = 2000*pi; % Stream-wise cut-off
lambda_z_plus = 250*pi; % Span-wise cut-off

% Corresponding cut-off wavenumbers
kx = round(Lx/lambda_x_plus); 
kz = round(Lz/lambda_z_plus);

% Load the list of names of the input data
Files_input = dir([dir_input,'tau_*.mat']);

% Apply the large-scale filter to the input data
for n = 1:size(Files_input,1)

    load([dir_input,Files_input(n).name]);

    nx = size(taux,1); % Number of grid points in x direction
    nz = size(taux,2); % Number of grid points in z direction

    taux_hat = fft2(taux); % 2D Fast Fourier Transform, stream-wise velocity
    tauz_hat = fft2(tauz); % 2D Fast Fourier Transform, span-wise velocity

    taux_hat((kx+2):(nx-kx),:) = 0; % Zero out the specified modes
    taux_hat(:,(kz+2):(nz-kz)) = 0; 

    tauz_hat((kx+2):(nx-kx),:) = 0;
    tauz_hat(:,(kz+2):(nz-kz)) = 0; 

    tauxL = ifft2(taux_hat); % 2D Inverse FFT
    tauzL = ifft2(tauz_hat); 

    tauL = sqrt(tauxL.^2+tauzL.^2); % Resultant wall shear stress
    u_tauL_plus = sqrt(tauL); % Large-scale friction velocity

    % Direction of the large-scale friction velocity
    theta = -atan(tauzL./tauxL); 

    % Name of files of the filtered field
    newname = strrep(Files_input(n).name,'tau','u_tauL'); 

    % Store the output data
    save([dir_input,newname],'u_tauL_plus','theta','-v7.3'); 

end

% Load the list of names of the output data
Files_output = dir([dir_input,'u_tauL_*.mat']);

% Mean value of u_tauL_plus and theta

% Pre-allocate the matrix 
u_tauL_mean_loc = zeros([size(Files_output,1) 1]);
theta_mean_loc = zeros([size(Files_output,1) 1]);

for n = 1:size(Files_output,1)

    load([dir_input,Files_output(n).name]);

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

    load([dir_input,Files_output(n).name]);

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
save([dir_input,'utauL_stats.mat'], ...
    'u_tauL_plus_mean','theta_mean','u_tauL_var','theta_var',...
    'theta_var_include_theta_mean');

%% Is this a test?
if test == true
    outputending = '_test';
elseif test == false
    outputending = '';
end

% Init progress bar
progressbar
%% Load the input data
load([dir_vel_tilde,'vel_tilde_1D.mat']);
load([dir_y_ref,'y_ref.mat'])
load([dir_input,'utauL_stats.mat']);

%% Set up the rid
[valmin,idmin] = min(abs(ymin - y));
[valmax,idmax] = min(abs(ymax - y));
y_plus_ref = y(idmin:idmax).';

y_plus = y_plus_ref;

save([dir_output,'y_plus.mat'],'y_plus','-v7.3');

clear y_plus;

%% QSQH model

Files = dir([dir_input,'u_tauL_*.mat']);

load([dir_input,Files(1).name]);

X = size(u_tauL_plus,1);
Z = size(u_tauL_plus,2);

pw_x = ceil(log2(X));

pw_z = ceil(log2(Z));

pw_t = ceil(log2(size(Files,1)));

Nx0 = 2; Nz0 = 2; T0 = 2; pw = starting_loop-1;

while pw < pw_x || pw < pw_z || pw < pw_t

    pw = pw+1;

    if pw > ending_loop
        break
    end

    if pw > pw_x || pw == pw_x

        Nx = X;

    elseif pw < pw_x

        Nx = Nx0^pw;

    end

    if pw > pw_z || pw == pw_z

        Nz = Z;

    elseif pw < pw_z

        Nz = Nz0^pw;

    end

    if pw > pw_t || pw == pw_t

        T = size(Files,1);

    elseif pw < pw_t

        T = T0^pw;

    end
    
    nx = round(linspace(1,X,Nx));
    nz = round(linspace(1,Z,Nz));

    tstep = round(linspace(1,size(Files,1),T));

    u_plus_qsqh = zeros(T,(Nx*Nz),numel(y_plus_ref));

    v_plus_qsqh = zeros(T,(Nx*Nz),numel(y_plus_ref));

    w_plus_qsqh = zeros(T,(Nx*Nz),numel(y_plus_ref));

    for n = 1:size(tstep,2)

        LS_wall = load([dir_input,Files(tstep(n)).name]);

        u_tauL = reshape(LS_wall.u_tauL_plus(nx,nz),[],1)/u_tauL_plus_mean;

        theta = reshape(LS_wall.theta(nx,nz),[],1);

        theta_fluc = theta - theta_mean;

        [u_tilde_sample,idx] = datasample(u_tilde,numel(u_tauL));

        v_tilde_sample = v_tilde(idx,:);

        w_tilde_sample = w_tilde(idx,:);

        y_tilde_sample = y_tilde(idx,:);

        for num = 1:numel(u_tauL)
        
            u_plus = (u_tilde_sample(num,:)*cos(theta_fluc(num,:))+...
                w_tilde_sample(num,:)*sin(theta_fluc(num,:)))*...
                u_tauL(num)*u_tauL_plus_mean;

            v_plus = v_tilde_sample(num,:)*u_tauL(num,:)*u_tauL_plus_mean;

            w_plus = (-u_tilde_sample(num,:)*sin(theta_fluc(num,:))+...
                w_tilde_sample(num,:)*cos(theta_fluc(num,:)))...
                *u_tauL(num,:)*u_tauL_plus_mean;

            y_plus = (y_tilde_sample(num,:)/u_tauL(num,:))/u_tauL_plus_mean;

            u_plus_qsqh(n,num,:) = csapi(y_plus,u_plus,y_plus_ref);

            v_plus_qsqh(n,num,:) = csapi(y_plus,v_plus,y_plus_ref);

            w_plus_qsqh(n,num,:) = csapi(y_plus,w_plus,y_plus_ref);

            progressbar((n*num)/(size(tstep,2)*numel(u_tauL)));

        end

    end

    save([dir_output,'u_plus_qsqh_loop_',num2str(pw,'%03d'),outputending, ...
        '.mat'],'u_plus_qsqh','-v7.3');

    save([dir_output,'v_plus_qsqh_loop_',num2str(pw,'%03d'),outputending, ...
        '.mat'],'v_plus_qsqh','-v7.3');

    save([dir_output,'w_plus_qsqh_loop_',num2str(pw,'%03d'),outputending, ...
        '.mat'],'w_plus_qsqh','-v7.3');

    disp(['Loop ',num2str(pw,'%03d'),' finished.']);

end
