% Is this a test?
if test == true
    outputending = '_test';
elseif test == false
    outputending = '';
end

% Init progress bar
progressbar
%% Load the all the input data
% load(dir_vel_tilde);
load([dir_y_ref,'y_ref.mat']);
load([dir_input,'utauL_stats.mat']);

%% Set up the uniform grid

% 1D grid in x,z,y directions
x_plus_ref = 0:12:xmax;
z_plus_ref = 0:6:zmax;

[valmin,idmin] = min(abs(ymin - y));
[valmax,idmax] = min(abs(ymax - y));
y_plus_ref = y(idmin:idmax).';

% 3D grid
[Xref,Zref,Yref] = ndgrid(x_plus_ref,z_plus_ref,y_plus_ref);

x_plus = x_plus_ref;
y_plus = y_plus_ref;
z_plus = z_plus_ref;

save([dir_output,'boxsize.mat'],'x_plus','z_plus','y_plus');

clear x_plus y_plus z_plus;

%% QSQH model

Files = dir([dir_input,'u_tauL_*.mat']);

load([dir_input,Files(1).name]);

load([dir_vel_tilde,'xz_tilde.mat']);

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

    u_plus_qsqh = zeros(T,(Nx*Nz),...
        numel(x_plus_ref),numel(z_plus_ref),numel(y_plus_ref));

    v_plus_qsqh = zeros(T,(Nx*Nz),...
        numel(x_plus_ref),numel(z_plus_ref),numel(y_plus_ref));

    w_plus_qsqh = zeros(T,(Nx*Nz),...
        numel(x_plus_ref),numel(z_plus_ref),numel(y_plus_ref));

    for n = 1:T

        LS_wall = load([dir_input,Files(tstep(n)).name]);

        u_tauL = reshape(LS_wall.u_tauL_plus(nx,nz),[],1)/u_tauL_plus_mean;

        theta = reshape(LS_wall.theta(nx,nz),[],1);

        theta_fluc = theta - theta_mean;

        TildeFileList = dir(fullfile(dir_vel_tilde,'vel_tilde_*.mat'));

        index = randi(numel(TildeFileList),1,numel(u_tauL));

        for k = 1:numel(u_tauL)

            vel_tilde = load([dir_vel_tilde,TildeFileList(index(k)).name]);

            u_tilde = vel_tilde.u_tilde;

            v_tilde = vel_tilde.v_tilde;

            w_tilde = vel_tilde.w_tilde;

            y_tilde = vel_tilde.y_tilde;

            u_plus = (u_tilde*cos(theta_fluc(k))+...
                w_tilde*sin(theta_fluc(k)))*...
                u_tauL(k)*u_tauL_plus_mean;

            v_plus = v_tilde*u_tauL(k)*u_tauL_plus_mean;

            w_plus = (-u_tilde*sin(theta_fluc(k))+...
                w_tilde*cos(theta_fluc(k)))...
                *u_tauL(k)*u_tauL_plus_mean;

            u_plus = double(u_plus); 
            v_plus = double(v_plus); 
            w_plus = double(w_plus);

            y_plus = (y_tilde/u_tauL(k))/u_tauL_plus_mean;

            x_plus_r = (x_tilde/u_tauL(k))/u_tauL_plus_mean;

            z_plus_r = (z_tilde/u_tauL(k))/u_tauL_plus_mean;

            [Xq,Zq,Yq] = ndgrid(x_plus_r,z_plus_r,y_plus);

            Xrot = (Xq-max(u_plus)/2)*cosd(theta(k)) +...
                (Zq-max(w_plus)/2)*sind(theta(k)) + max(u_plus)/2;
            Zrot = -(Xq-max(u_plus)/2)*sind(theta(k)) +...
                (Zq-max(w_plus)/2)*cosd(theta(k)) + max(w_plus)/2;
            Xrot = double(Xrot); Zrot = double(Zrot); Yq = double(Yq);

            u_plus_qsqh(n,k,:,:,:) = griddata(Xrot,Zrot,Yq,u_plus, ...
                Xref,Zref,Yref);
            v_plus_qsqh(n,k,:,:,:) = griddata(Xrot,Zrot,Yq,v_plus, ...
                Xref,Zref,Yref);
            w_plus_qsqh(n,k,:,:,:) = griddata(Xrot,Zrot,Yq,w_plus, ...
                Xref,Zref,Yref);
            
            progressbar((n*k)/(size(tstep,2)*numel(u_tauL)));
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
