%% Drift and Rotation correction
% Load data
% Manual select 1st and 2nd lattice vectors, and input number of NP along the vectors.
% Rough drift and rotation correction to the 1st frame
% Fine drift and rotation correction
% Rotation so that lattice vector 2 is along x-axis.

clear; close all; clc;
addpath('functions/');

%%
Input_path = 'Data raw/';
Output_path = 'Data drift rotation corrected/';
Plot_path = 'Plots/';

activate_sets = [4];
% Input format
%   Input matrix name
%   3 indices to determine 2 lattice vectors:
%       p1-p2 as vector 1
%       p1-p3 as vector 2
%   2 numbers indicating number of NP along each vector
% If the input only include 
Input_params{1} = {'Cube_22mM_raw.mat',     1,13,79,        13,7};
Input_params{2} = {'Cube_27mM_raw.mat',     1,8,57,         8,8};
Input_params{3} = {'Cube_110mM_raw.mat',    7,1,56,         7,8};
Input_params{4} = {'Rod_raw.mat',           32,135,201,     10,10};
Input_params{5} = {'Prism_raw.mat',         1,46,41,        6,6};

if ~exist(Output_path, 'dir'); mkdir(Output_path); end
if ~exist(Plot_path, 'dir'); mkdir(Plot_path); end

%% Const
Rot_mat = @(x) [cos(x), sin(x); -sin(x), cos(x)];

%%
figure(1)
set(gcf,'Position',[122.3333 414.3333 1304 420.0000])
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for activate_set = activate_sets
    In_name = Input_params{activate_set}{1};
    load([Input_path, In_name]); % Position matrix

    % Inspection:
    nexttile(1)
    hold off
    PlotPositionMatrix(pos)
    hold on
    temp = mean(pos,3);
    for i = 1:size(pos,1)
        text(temp(i,1),temp(i,2),num2str(i));
    end
    axis equal

    % If input not ready, skip processing.
    if length(Input_params{activate_set})<6; continue; end

    %% Rough correction
    % Drift
    pos = pos - mean(pos,1);
    
    % % Rotation to align to 1st frame, by averaged rotation angle between
    % % two frames.
    % % First select particles that are away from the center, for accuracy
    % % and speed.
    % d = sqrt(sum(pos(:,:,1).^2,2));
    % select = find(d > max(d)/2);
    % for i = 2:size(pos,3)
    %     angle_temp = [];
    %     for j = select'
    %         angle_temp = [angle_temp, AngleCalc([0,0], pos(j,:,i), pos(j,:,1))];
    %     end
    %     pos(:,:,i) = pos(:,:,i) * Rot_mat(mean(angle_temp)/180*pi);
    % end

    %% Fine correction
    R0 = zeros(size(pos,3),3); % Initial guess
    Loss = @(x) FrameCompare_TransRot(x,pos); % Loss function
    
    % Fitting without display
    R = fminsearch(Loss, R0);
    % Fitting with display
    % options = optimset('Display','iter','PlotFcns',@optimplotfval);
    % R = fminsearch(Loss, R0);

    R(1,:) = 0;
    pos_corrected = zeros(size(pos));
    for i = 1:size(pos,3)
        pos_corrected(:,:,i) = (pos(:,:,i) - R(i,2:3)) * Rot_mat(R(i,1));
    end
    CoM = mean(mean(pos_corrected,3),1);
    pos_corrected = pos_corrected - CoM;

    figure(1)
    nexttile(2)
    PlotPositionMatrix(pos_corrected)
    axis equal
    title('After drift/rotation correction')

    %% Rotation
    [id1, id2, id3] = deal(Input_params{activate_set}{2:4});
    [n1, n2] = deal(Input_params{activate_set}{5:6});

    temp = mean(pos_corrected, 3);
    angle_temp = AngleCalc(temp(id1,:), temp(id3,:));
    pos_rotate = pos_corrected;
    for i = 1:size(pos,3)
        pos_rotate(:,:,i) = pos_rotate(:,:,i) * Rot_mat(angle_temp/180*pi);
    end

    nexttile(3)
    PlotPositionMatrix(pos_rotate)
    axis equal
    title('After rotation')

    %% Save output
    pos = pos_rotate;
    
    temp = mean(pos,3);
    v1 = (temp(id2,:) - temp(id1,:))/(n1-1);
    v2 = (temp(id3,:) - temp(id1,:))/(n2-1);

    % Prepare for lattice fitting
    x0 = [v2,v1, -(n2-1)/2, (n2-1)/2, -(n1-1)/2, (n1-1)/2];

    Out_name = strrep(In_name,'_raw','');
    save([Output_path, Out_name], 'pos','v1','v2','n1','n2','x0');
    saveFigures([Plot_path,strrep(Out_name,'.mat','_pos')])

end

%% Functions











