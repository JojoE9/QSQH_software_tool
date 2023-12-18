%% Interface for the data converter

clear
clc

% User input parameters

lines = readlines('user_input_3D.txt');
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

%% Excuting Data converter
Data_converter

%% Excuting QSQH model
if xmax == 0 && zmax == 0
    QSQH_model_1D
else
    QSQH_model_3D
end
