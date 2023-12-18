%%
% Is this a test?
if test == true
    outputending = '_test';
elseif test == false
    outputending = '';
end

% Init progress bar
progressbar
%% Extrapolation


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
