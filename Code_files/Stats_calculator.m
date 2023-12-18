clear
clc

% Make sure the current folder is the one containing this code file
workdir = pwd;

% Input parameters

% Directory of the 1D output
diroutput_1D = [workdir,'\output_1D_test\'];

% The names of the files of the output and the grid of 1D output
File1D_vel_name = 'u_plus_qsqh_loop_008_test.mat';
File1D_y_name = 'y_plus.mat';

% Directory of the 3D output
diroutput_3D = [workdir,'\output_3D_test\'];

% The names of the files of the output and the grid of 3D output
File3D_vel_name = 'u_plus_qsqh_loop_004_test.mat';
File3D_xyz_name = 'boxsize.mat';

%% Mean profile and urms square of 1D output

load([diroutput_1D,File1D_vel_name]);
load([diroutput_1D,File1D_y_name]);

U_plus_qsqh = mean(u_plus_qsqh,[1 2]);

u_fluc_plus_qsqh = u_plus_qsqh - U_plus_qsqh;

urms2_plus_qsqh = mean((u_fluc_plus_qsqh.^2),[1 2]);
%% Plot the mean profile of the 1D output 
figure
hold on 

plot(y_plus,squeeze(U_plus_qsqh),'LineStyle','-','Marker','none', ...
    'MarkerIndices',1:10:length(y_plus),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'User Generated Output')

load([diroutput_1D,'u_plus_qsqh_loop_008_test_ref.mat'],'U_plus_qsqh_ref');
load([diroutput_1D,'y_plus_ref.mat'])

plot(y_plus_ref,squeeze(U_plus_qsqh_ref),'LineStyle','none','Marker','+', ...
    'MarkerIndices',1:10:length(y_plus_ref),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'Reference Ouput')

legend('FontSize',16,'Interpreter','latex','Location','Southeast', ...
    'NumColumns',1);
legend('boxoff')

xlim([0 50])
% ylim([-0.2 0.02])
ax = gca; 
ax.FontSize = 16; 
xlabel('$y^+$','FontSize',16,'Interpreter','latex');
ylabel('$\langle U \rangle^+$','FontSize',16,'Interpreter','latex');

box on
hold off

%% Plot the urms square of the 1D output 
figure
hold on 

plot(y_plus,squeeze(urms2_plus_qsqh),'LineStyle','-','Marker','none', ...
    'MarkerIndices',1:10:length(y_plus),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'User Generated Output')

load([diroutput_1D,'u_plus_qsqh_loop_008_test_ref.mat'],'urms2_plus_qsqh_ref');

plot(y_plus_ref,squeeze(urms2_plus_qsqh_ref),'LineStyle','none','Marker','+', ...
    'MarkerIndices',1:10:length(y_plus_ref),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'Reference Ouput')
legend('off')

xlim([0 50])
% ylim([-0.2 0.02])
ax = gca; 
ax.FontSize = 16; 
xlabel('$y^+$','FontSize',16,'Interpreter','latex');
ylabel('$\langle u''^2 \rangle^+$','FontSize',16,'Interpreter','latex');

box on
%% Two-point correlation coefficients in wall-parallel directions of 3D output 

load([diroutput_3D,File3D_vel_name]);
load([diroutput_3D,File3D_xyz_name]);

ny = 2;

u = double(u_plus_qsqh(:,:,:,:,ny));

Dx_max = size(u,3);
Dz_max = size(u,4);

Dx_qsqh = linspace(0,(Dx_max-1),Dx_max)*12;
Dz_qsqh = linspace(0,(Dz_max-1),Dz_max)*6;

Rux_qsqh = zeros([1 Dx_max]);
Ruz_qsqh = zeros([1 Dz_max]);

for i = 1:size(Dx_qsqh,2)

    u_reshape = u;

    u_xshift = circshift(u_reshape,[0,0,-(i-1),0]);

    if i > 1
        u_xshift(:,:,((Dx_max-(i-1)):Dx_max),:) = [];

        u_reshape(:,:,((Dx_max-(i-1)):Dx_max),:) = [];
    end

    U_reshape = mean(u_reshape,'all','omitnan');

    U_xshift = mean(u_xshift,'all','omitnan');

    up_reshape = u_reshape - U_reshape;

    up_xshift = u_xshift - U_xshift;

    correlation = mean((up_reshape.*up_xshift),'all','omitnan');

    std_up_reshape = sqrt(mean((up_reshape.^2),'all','omitnan'));

    std_up_xshift = sqrt(mean((up_xshift.^2),'all','omitnan'));

    Rux_qsqh(i) = correlation/(std_up_reshape*std_up_xshift);
end

for i = 1:size(Dz_qsqh,2)
    u_reshape = u;

    u_zshift = circshift(u_reshape,[0,0,0,-(i-1)]);

    if i > 1
        u_zshift(:,:,:,((Dz_max-(i-1)):Dz_max)) = [];

        u_reshape(:,:,:,((Dz_max-(i-1)):Dz_max)) = [];
    end

    U_reshape = mean(u_reshape,'all','omitnan');

    U_zshift = mean(u_zshift,'all','omitnan');

    up_reshape = u_reshape - U_reshape;

    up_zshift = u_zshift - U_zshift;

    correlation = mean((up_reshape.*up_zshift),'all','omitnan');

    std_up_reshape = sqrt(mean((up_reshape.^2),'all','omitnan'));

    std_up_zshift = sqrt(mean((up_zshift.^2),'all','omitnan'));

    Ruz_qsqh(i) = correlation/(std_up_reshape*std_up_zshift);

end

%% Plot two-point correlation coefficient in x direction of the 3D output 
figure
hold on 

plot(Dx_qsqh,Rux_qsqh,'LineStyle','-','Marker','none', ...
    'MarkerIndices',1:10:length(Dx_qsqh),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'User Generated Ouput')

load([diroutput_3D,'u_plus_qsqh_loop_005_test_ref.mat'],'Dx_qsqh_ref', ...
    'Rux_qsqh_ref');

plot(Dx_qsqh_ref,Rux_qsqh_ref,'LineStyle','none','Marker','+', ...
    'MarkerIndices',1:2:length(Dx_qsqh_ref),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'Reference Ouput')

legend('FontSize',16,'Interpreter','latex','Location','Southeast', ...
    'NumColumns',1);
legend('boxoff')

xlim([0 102])
ylim([0 1])
ax = gca; 
ax.FontSize = 16; 
xlabel('$y^+$','FontSize',16,'Interpreter','latex');
ylabel('$$R_{u''}(\Delta x^+)$$','FontSize',16,'Interpreter','latex');

box on
hold off

%% Plot the urms square of the 1D output 
figure
hold on 

plot(Dz_qsqh,Ruz_qsqh,'LineStyle','-','Marker','none', ...
    'MarkerIndices',1:2:length(Dz_qsqh),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'User Generated Ouput')

load([diroutput_3D,'u_plus_qsqh_loop_005_test_ref.mat'],'Dz_qsqh_ref', ...
    'Ruz_qsqh_ref');

plot(Dz_qsqh_ref,Ruz_qsqh_ref,'LineStyle','none','Marker','+', ...
    'MarkerIndices',1:2:length(Dz_qsqh_ref),'LineWidth',1, ...
    'Color','k','DisplayName', ...
    'Reference Ouput')

legend('FontSize',16,'Interpreter','latex','Location','Southeast', ...
    'NumColumns',1);
legend('off')

xlim([0 60])
ylim([0 1])
ax = gca; 
ax.FontSize = 16; 
xlabel('$y^+$','FontSize',16,'Interpreter','latex');
ylabel('$$R_{u''}(\Delta z^+)$$','FontSize',16,'Interpreter','latex');

box on
hold off
