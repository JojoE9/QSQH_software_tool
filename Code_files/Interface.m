%% Interface for the data converter

clear
clc

% Locate the directory of the data of the velocity field at the wall
dir_input = '';

% Locate the directory of the output of the converted data
dir_save = '';

% Size of the target flow field (in wall unit)
Lx = 8*1000*pi; % Length of the stream-wise domain, the pre-allocated value is the testing value
Lz = 3*1000*pi; % Length of the span-wise domain

% Excuting code file
Data_converter

%% Interface for the QSQH synthetic model

clear
clc

% Locate the directory of the universal velocity field
dir_vel_tilde = 'Z:\MK_Channel\MK5200\vel_tilde_X_2000pi_Z_250pi_box\';

% Locate the directory of reference wall-normal grid
dir_y_ref = 'Z:\MK_Channel\MK5200\';

% Locate the directory of u_tauL and theta
dir_u_tauL = 'Z:\MK_Channel\MK5200\u_tauL_X_2000pi_Z_250pi\';

% Locate the directory of statistics of u_tauL and theta
dir_tau_stats  = 'Z:\MK_Channel\MK5200\u_tauL_X_2000pi_Z_250pi\';

% Locate the directory of the output of the synthetic model
dir_output = 'Z:\MK_Channel\MK5200\vel_qsqh_X_2000pi_Z_250pi_box\';

% Maximum spatial positions in wall unit
% xmax = [];
ymax = 5; 
% zmax = [];

% Excuting code file
QSQH_model
