%% Interface for the data converter

clear
clc

% Locate the directory of the data of the velocity field at the wall
dir_input = '';

% Locate the directory of the output of the converted data
dir_output = '';

% Size of the target flow field (in wall unit)
Lx = []; % Length of the stream-wise domain
Lz = []; % Length of the span-wise domain

% Wall-normal position of the input data (in wall unit)
y = [];

% Excuting code file
Data_converter

%% Interface for the QSQH synthetic model

clear
clc

% Locate the directory of the universal velocity field
dir_vel_tilde = '';

% Locate the directory of reference wall-normal grid
dir_y_ref = '';

% Locate the directory of u_tauL and theta
dir_u_tauL = '';

% Locate the directory of statistics of u_tauL and theta
dir_u_tauL_stats  = '';

% Locate the directory of the output of the synthetic model
dir_save = '';

% Maximum wall-normal position (choose between 2 and 104)
ymax = 104; 

% Excuting code file
QSQH_model
