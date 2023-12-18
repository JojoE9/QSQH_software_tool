%% Interface for the data converter

clear
clc

% User input parameters

% Size of the input friction field in wall unit
% (The test input has the size of Lx = 8000pi and Lz = 3000pi.)
Lx = 8*1000*pi; % Length of the stream-wise domain
Lz = 3*1000*pi; % Length of the span-wise domain

% Choose the dimension of the output synthetic field (1D/3D)
dim = '1D';

% Choose the maximum spatial ranges in wall unit
% (xmax=0 & zmax=0 for 1D output field)
xmax = 0; %xmax = {0, 204}
zmax = 0; %zmax = {0, 120}

ymin = 0; %y = {0, 100}
ymax = 50; 

% Is this run a test? (true/false)
test = true;

% Starting loop (starting from 1 for a new run)
starting_loop = 8;

% The following information of directories should remain unchanged for the
% test, unless the different file separator "/" is used instead of "\". 

% The directory of the tool (current folder)
workdir = pwd;

% The directory of the input friction field
dir_input = [workdir,'\input_friction_field_test\'];

% The directory of the universal velocity field
dir_vel_tilde = [workdir,'\universal_velocity_field_',dim,'_test\'];

% The directory of the reference wall-normal grid
dir_y_ref = [workdir,'\'];

% The folder that will store the output of the synthetic model
dir_output = [workdir,'\output_',dim,'_test\'];



% Excuting code file
Data_converter

% Excuting code file
if (xmax == 0) && (zmax == 0)
    QSQH_model_1D
else 
    QSQH_model_3D
end
