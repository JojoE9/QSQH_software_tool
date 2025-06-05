%% Interface of QSQH tool
clear
clc

% User input parameters
lines = readlines('InputParameters.txt');

% Excute the code

InputParameters = extractAfter(lines(:),',');

xmax = str2double(InputParameters(3));
zmax = str2double(InputParameters(4));

workdir = pwd;
file_separator = convertStringsToChars(InputParameters(10));

if xmax == 0 && zmax == 0
    run([workdir,file_separator,'QSQH_model',file_separator,...
        'QSQH_model_1D'])
else
    run([workdir,file_separator,'QSQH_model',file_separator,...
        'QSQH_model_3D'])
end