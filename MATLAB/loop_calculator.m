clear
clc

workdir = pwd;
lines = readlines('InputParameters.txt');
InputParameters = extractAfter(lines(:),',');

test = str2num(InputParameters(7));
file_separator = convertStringsToChars(InputParameters(10));

% Is this a test?
if test == true
    folderending = '_test';
elseif test == false
    folderending = '';
end


% The directory of the input friction field
dir_input = [workdir,file_separator,'input_friction_field',folderending,...
    file_separator];

Files_input = dir([dir_input,'tau_*.mat']);

load([dir_input,Files_input(1).name]);


%%
X = size(taux,1);
Z = size(taux,2);

pw_x = ceil(log2(X));

pw_z = ceil(log2(Z));

pw_t = ceil(log2(size(Files_input,1)));

pw_max = max([pw_x,pw_z,pw_t]);
n_total = size(taux,1)*size(taux,2)*size(Files_input,1);

disp(['There are up to ',num2str(pw_max),' loops, ','the maximum' ...
    ' number of output instances is ', num2str(n_total),'.'])
%%
pw = 0; Nx0 = 2; Nz0 = 2; T0 =2;

while pw < pw_x || pw < pw_z || pw < pw_t

    pw = pw+1;

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

        T = size(Files_input,1);

    elseif pw < pw_t

        T = T0^pw;

    end
    
    nx = round(linspace(1,X,Nx));
    nz = round(linspace(1,Z,Nz));

    tstep = round(linspace(1,size(Files_input,1),T));

    n = numel(nx)*numel(nz)*numel(tstep);

    disp(['Loop ',num2str(pw),' generates ',num2str(n),' output instances.']);
end