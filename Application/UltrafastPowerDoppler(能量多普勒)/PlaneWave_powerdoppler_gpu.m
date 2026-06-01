% author:seu
% date:2025-07-30

% This is the second step for Power Doppler (Energy Doppler); step one is data acquisition.

% Data processing flow:
% Beamforming with demodulation -> SVD filtering -> Calculate Energy

% Can save beamformed IQ data and power maps
cd(fileparts(mfilename('fullpath')));
clear all
clc
close all

currentPath = pwd;
parentDir = fileparts(fileparts(currentPath));
addpath(genpath(parentDir));

%% Necessary Parameters
filefolder = 'D:\software_matlab\exampledata\doppler\20251020190356';  % Directory where data is located

% 
[fs,prf,sampleNum,scanLine,imagedepth,focus_depth,cstartoffset,frame_nums,numperfile,steering_deg,scaninfo] = read_adc_para(strcat(filefolder,'\Param.txt'),'plane wave');
fc=7.5e6; % Transmit frequency
BF_SampleN = sampleNum;

                            % fs               Sampling rate
                            % sampleNum        Number of sampling points
                            % numperfile       Number of frames contained in each saved data bin
                            % steering_deg     Plane wave steering angles
                            % ImageDepth       Image depth
                            % prf              Pulse Repetition Frequency used for data acquisition 

load_file_num = 20;         % Number of files to read and process

Imagestart = 0.002;         % B-mode image start depth
ImageDepth = 0.0360;        % B-mode image depth
BeamN = 256;                % Beams (number of beams in x-direction for beamforming)       

fft_period = 50;            % Number of frames used for one blood flow calculation (Window Size)

x1_loc_real = -0.006;       % ROI top-left physical coordinate x1
z1_loc_real = 0.003;        % ROI top-left physical coordinate z1
x2_loc_real = 0.006;        % ROI bottom-right physical coordinate x2
z2_loc_real = 0.015;        % ROI bottom-right physical coordinate z2

svd_auto = 1;               % Whether to use adaptive threshold for SVD filtering. If 1, svd_ord1 and svd_ord2 are ignored
svd_ord1 = 60;              % SVD filter start order
svd_ord2 = 90;              % SVD filter end order

is_save_BF = 1;             % Whether to save beamformed IQ data. 1 for save, 0 for do not save
Bmode_save_index = 1:(numperfile/numel(steering_deg));     % Frame indices to save for each data bin
                            % Only effective when is_save_BF = 1. Example: 1 means save only the 1st frame, 1:numperfile/numel(steering_deg) means save from 1st to last frame
is_save_power = 1;          % Whether to save blood flow data. 1 for save, 0 for do not save

drange_B = 60;              % B-mode dynamic range
drange_power = 30;          % Power Doppler dynamic range

% Default parameters
probe_name = 'L5-10';       % Probe name
TxChannel = 128;            % Number of transmit channels
RxChannel = 128;            % Number of receive channels
sos = 1540;                 % Speed of sound



%% Beamforming Parameter Calculation

probe = Probe_para(probe_name);
ch_map = probe.rx_ele_map(1:RxChannel);
elex = probe.element_pos.x;
elez = probe.element_pos.z;

x_axis = linspace(elex(1),elex(end),BeamN);
step_z = sos/fs/2;
z_axis = (0:step_z:(BF_SampleN-1) * step_z);
[x_grid, z_grid] = meshgrid(x_axis,z_axis);
scan.xx  =  x_grid;
scan.zz  =  z_grid;

allbeamx = reshape(scan.xx,[BF_SampleN*BeamN,1]);
allbeamz = reshape(scan.zz,[BF_SampleN*BeamN,1]);

probe_type_ptr = libpointer('cstring', "linear");

single_Steering = single(steering_deg);
single_Steering_ptr = libpointer('singlePtr', single_Steering);
single_Steering_len = length(single_Steering);

single_allbeamx = single(allbeamx);
single_allbeamx_ptr = libpointer('singlePtr', single_allbeamx);
single_allbeamx_len = length(single_allbeamx);

single_allbeamz = single(allbeamz);
single_allbeamz_ptr = libpointer('singlePtr', single_allbeamz);
single_allbeamz_len = length(single_allbeamz);

single_cstartoffset = single(cstartoffset);
single_cstartoffset_ptr = libpointer('singlePtr', single_cstartoffset);
single_cstartoffset_len = length(single_cstartoffset);

single_ch_map = single(ch_map);
single_ch_map_ptr = libpointer('singlePtr', single_ch_map);
single_ch_map_len = length(single_ch_map);

elexz = cat(2, elex, elez)';
single_elexz = single(elexz);
single_elexz_ptr = libpointer('singlePtr', single_elexz);
single_elexz_len = length(single_elexz);

lambda = sos/fc;
probe.element_width = probe.element_pitch/2;
f_number = est_fNumber(probe.element_width,lambda,0.71);

f_mask = zeros(BF_SampleN*BeamN,probe.element_num,'single');
for i = 1:probe.element_num
    f_mask(:,i) = allbeamz./abs(allbeamx -  probe.element_pos.x(i))/2 > f_number;
end
rx_apod = f_mask;
xm = bsxfun(@minus, probe.element_pos.x,allbeamx);
zm = bsxfun(@minus,probe.element_pos.z,allbeamz);
rx_delay = sqrt(xm.^2+zm.^2)/sos*fs;
SteeringNum = numel(steering_deg);
tx_delay = zeros(BF_SampleN*BeamN,SteeringNum);
steer = deg2rad(steering_deg);
for i = 1:SteeringNum
    if steer(i) >= 0
        tx_delay(:,i) = ((allbeamx - min(probe.element_pos.x))*sin(steer(i)) + allbeamz*cos(steer(i)))/sos*fs;
    else
        tx_delay(:,i) = ((max(probe.element_pos.x) - allbeamx)*sin(-steer(i)) + allbeamz*cos(steer(i)))/sos*fs;
    end
end

single_rx_apod = reshape(rx_apod,[],1);
single_rx_apod_ptr = libpointer('singlePtr', single_rx_apod);
single_rx_apod_len = length(single_rx_apod);

single_rx_delay = reshape(rx_delay,[],1);
single_rx_delay_ptr = libpointer('singlePtr', single_rx_delay);
single_rx_delay_len = length(single_rx_delay);

single_tx_delay = reshape(tx_delay,[],1);
single_tx_delay_ptr = libpointer('singlePtr', single_tx_delay);
single_tx_delay_len = length(single_tx_delay);


%% Parameter Calculation

% Demodulation filter
bandwidth = 80; 
Wn = (fc*bandwidth/100)/(fs/2);
% Determine filter length (empirical formula)
M = ceil(6.64 * fs/2 / (fc*bandwidth/100));   % Formula specific to Hamming window
% Ensure odd length for linear phase
if mod(M, 2) == 0
    M = M + 1;
end
% Design low-pass FIR filter using Hamming window
b_fir = fir1(M - 1, Wn, 'low', hamming(M));
% Pass to GPU
single_filter = single(b_fir);
single_filter_ptr = libpointer('singlePtr', single_filter);
single_filter_len = length(single_filter);

% 
buffer_num = floor(fft_period/(numperfile/SteeringNum))+2; % Number of buffers
t_idx = fft_period; % This idx is the current index in the image, initially fft_period 
t_buffer_idx = fft_period; % This idx is the index in the buffer, initially fft_period 


%% Coordinates of Pixels Used for Blood Flow Calculation


x_dif = abs(x_axis - x1_loc_real);
z_dif = abs(z_axis - z1_loc_real);
[~, x_left] = min(x_dif);
[~, z_up] = min(z_dif);

x_dif = abs(x_axis - x2_loc_real);
z_dif = abs(z_axis - z2_loc_real);
[~, x_right] = min(x_dif);
[~, z_down] = min(z_dif);

rec_x_num = x_right - x_left + 1;
rec_z_num = z_down - z_up + 1;

x_loc_all = [];
z_loc_all = [];
for i = x_left:x_right
    for j = z_up:z_down
        z_loc_all = [z_loc_all,j];
        x_loc_all = [x_loc_all,i];
    end
end

loc_num = numel(x_loc_all);

single_x_loc_all = single(x_loc_all);
single_x_loc_all_ptr = libpointer('singlePtr', single_x_loc_all);
single_z_loc_all = single(z_loc_all);
single_z_loc_all_ptr = libpointer('singlePtr', single_z_loc_all);
single_points_len = loc_num;

%% Loading

% Load DLL
if ~libisloaded('US_APP')
    loadlibrary('US_APP.dll', 'ApplicationMatlabInterface.h');
end

bprocesspara.Demod_AFE_Dynamic = -9;
gpu_handle = calllib('US_APP', 'initializepowerGPU', ...
    prf, ...
    numperfile, ...
    buffer_num, ...
    fft_period, ...
    2, ... %lag
    2, ... %t_axis_span
    t_idx, ...
    t_buffer_idx, ...
    probe_type_ptr, ...
    bprocesspara.Demod_AFE_Dynamic, ...
    1, ... %AcqConfig.Tx.FsNum
    TxChannel, ... %AcqConfig.Tx.Channel
    RxChannel, ... %RxChannel
    probe.element_num, ... %AcqConfig.Probe.element_num
    probe.element_pitch, ... %AcqConfig.Probe.element_pitch
    32, ...
    BF_SampleN, ...
    BeamN, ...
    sos, ...
    fs, ...
    fc, ...
    svd_auto,... svd
    svd_ord1,... svd 
    svd_ord2,... svd
    single_Steering_ptr,single_Steering_len,...
    single_allbeamx_ptr,single_allbeamx_len,...
    single_allbeamz_ptr,single_allbeamz_len,...
    single_cstartoffset_ptr,single_cstartoffset_len,...
    single_ch_map_ptr,single_ch_map_len,...
    single_elexz_ptr,single_elexz_len,...
    single_filter_ptr,single_filter_len,...
    single_rx_apod_ptr,single_rx_apod_len,...
    single_rx_delay_ptr,single_rx_delay_len,...
    single_tx_delay_ptr,single_tx_delay_len, ...
    single_x_loc_all_ptr,single_z_loc_all_ptr,single_points_len, ...
    rec_x_num, rec_z_num);


%% Get All Valid Data Files
files =dir (fullfile (filefolder ,'**' ,'*bin' ));
num_files =length (files);
file_numbers =zeros (num_files, 1);

for i =1 :num_files
    [~,file_name ]=fileparts (files(i).name );
    try
        file_numbers(i)=str2double(file_name);
    catch
        file_numbers(i)=inf ;
    end
end

valid_indices =~isinf(file_numbers);
valid_files =files (valid_indices);
valid_numbers =file_numbers (valid_indices);
[~,sort_indices ]=sort (valid_numbers);
sorted_files =valid_files (sort_indices);


%% Define Image

% B-mode image to display (beamformed image)
bfdata = zeros(1, 2* BF_SampleN * BeamN * numperfile / SteeringNum);
bfdata = single(bfdata);
bfdata_ptr = libpointer('singlePtr', bfdata);

valid_indices = find(z_axis >= Imagestart & z_axis <= ImageDepth);
zz_cut = z_axis(valid_indices);
Nz_cut = numel(valid_indices);

% Get screen dimensions
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);
            
real_width = x_grid(end)- x_grid(1);
real_height = z_grid(end) - z_grid(1);

% Define window proportion relative to screen
fig_width_ratio = 0.3;
fig_height_ratio = real_height/real_width*fig_width_ratio*1.5;
% Set window size
fig_width = screen_width * fig_width_ratio;
fig_height = screen_height * fig_height_ratio;
% Set window position
fig_left = screen_width * (0.5-fig_width_ratio/2);
fig_bottom = screen_height * (1-fig_height_ratio)/2;

hFig1 = figure('Name',"Bmode",'Position', [fig_left, fig_bottom, fig_width, fig_height]);
hIm1 = imagesc(x_axis, zz_cut, zeros(Nz_cut,BeamN)-100,[-60 0]);
colormap(gray);title("Power Doppler");
axis equal;axis tight

power_min = -drange_power;  % Min power dB
power_max = 0;   % Max power dB
colormap(hot); 
cb = colorbar; % Add colorbar
cb.Label.String = 'power (dB)'; % Set colorbar label
% Set colorbar limits to correspond with blood flow data range
caxis([power_min, power_max]);

hold on
rectangle('Position', [x_axis(x_left), z_axis(z_up), x_axis(x_right) - x_axis(x_left), z_axis(z_down) - z_axis(z_up)], 'EdgeColor', 'r', 'LineWidth', 2);

% Blood flow to display
power_matrix = single(zeros(rec_z_num, rec_x_num));
power_matrix_ptr = libpointer('singlePtr', power_matrix);

% Create directory if saving is enabled
if is_save_BF
    if ~exist(fullfile(filefolder, 'bfdata'), 'dir')
        mkdir(fullfile(filefolder, 'bfdata'));
    end
end
if is_save_power
    if ~exist(fullfile(filefolder, 'power'), 'dir')
        mkdir(fullfile(filefolder, 'power'));
    end
    if ~exist(fullfile(filefolder, 'B+power'), 'dir')
        mkdir(fullfile(filefolder, 'B+power'));
    end
end


%% Read Data and Process
for file_i = 1:load_file_num
    disp(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name))
    fileID = fopen(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name), 'rb'); 

    fseek(fileID,0,1);
    nFileLen = ftell(fileID);
    fseek(fileID,0,-1);
    % Note type control here
    alldata = fread(fileID, nFileLen,'int8=>int8');
    
    sid = 0;
    alldata_len = SteeringNum*(2*BF_SampleN*TxChannel+32);
    cur_data = alldata(1:(numperfile/SteeringNum)*alldata_len);
    
    alldata_ptr = libpointer('int8Ptr', cur_data);
    bag_idx = file_i-1;  % This bag_idx must start from 0
    ret = calllib('US_APP', 'processPowerDataBeamformingandPostGPU', gpu_handle, alldata_ptr, ...
        alldata_len, bfdata_ptr, power_matrix_ptr, bag_idx); 
    
    rev = bfdata_ptr.Value;
    bfdata_real = reshape(rev(1:2:end),BF_SampleN,BeamN,[]);
    bfdata_imag = reshape(rev(2:2:end),BF_SampleN,BeamN,[]);
    bfdata_iq = bfdata_real + 1i* bfdata_imag;
    B_image = log_compressed(abs(bfdata_iq(valid_indices,:,end)));
    B_image(B_image<-drange_B) = -drange_B;
    B_image(B_image>0) = 0;
    B_image = B_image + drange_B;
    B_image = B_image/drange_B; % Now B_image range is [0, 1]
    B_image_rgb = repmat(B_image, [1, 1, 3]); % Replicate 3 times to create [rows, cols, 3] RGB matrix

    if is_save_BF
        bfdata_save = bfdata_iq(:,:,Bmode_save_index);
        save(fullfile(filefolder,'bfdata',num2str(bag_idx)+".mat"), 'bfdata_save',"x_axis","z_axis");
    end

    % Get power data and map back to original image position 
    power_matrix_update = reshape(power_matrix_ptr.Value, rec_z_num, rec_x_num);
    power_matrix_update = log_compressed(power_matrix_update);
    power_full = zeros(BF_SampleN,BeamN);
    power_full(z_up:z_down,x_left:x_right) = power_matrix_update;
    power_full = power_full(valid_indices,:);

    % Get area mask and map back to original image position
    flow_mask_full = zeros(BF_SampleN,BeamN);
    flow_mask_full(z_up:z_down,x_left:x_right) = 1;
    flow_mask_full = logical(flow_mask_full); 
    flow_mask_full = flow_mask_full(valid_indices,:);
    % Expand mask to 3 color channels
    mask_3D = repmat(flow_mask_full, [1, 1, 3]);

    % Normalize data to [1, 256] integer indices
    power_indices = round( (power_full - power_min) / (power_max - power_min) * 255 + 1 );
    % Handle out-of-range values
    power_indices(power_indices < 1) = 1;
    power_indices(power_indices > 256) = 256;
    % Use ind2rgb to convert index matrix and colormap to RGB image
    power_image_rgb = ind2rgb(power_indices, hot(256));

    % In masked areas, replace original grayscale data with blood flow data
    B_image_rgb(mask_3D) = power_image_rgb(mask_3D);

    hIm1.CData = B_image_rgb;
    pause(0.001)

    if is_save_power
        save(fullfile(filefolder,'power',num2str(bag_idx)+".mat"), 'power_matrix_update');
        save(fullfile(filefolder,'B+power',num2str(bag_idx)+".mat"), "B_image_rgb","x_axis","zz_cut");
    end

    fclose(fileID);
end


%% Unload

% Release GPU memory and RAM
calllib('US_APP', 'deleteBeamformingGPUHandle', gpu_handle);

% Unload DLL
unloadlibrary('US_APP');







