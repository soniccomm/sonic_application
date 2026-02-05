% author:seu
% date:2025-07-30

% This is the second step of Spectral Doppler; step one is data acquisition (single angle).

% Data processing flow:
% Beamforming with demodulation -> Wall filtering on fft_period frames -> Calculate spectrum for each point and sum them up -> (Convert frequency shift to flow velocity and) Display

% Can save beamformed IQ data and spectrograms

clear all
clc
close all

currentPath = pwd;
parentDir = fileparts(fileparts(currentPath));
addpath(genpath(parentDir));

%% Necessary Parameters
filefolder = 'D:\software_matlab\exampledata\doppler\20251020190426';  % Directory where data is located

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


load_file_num = 40;         % Number of files to read and process

Imagestart = 0.002;         % B-mode image start depth
ImageDepth = 0.040;         % B-mode image depth
BeamN = 256;                % Beams (number of beams in x-direction for beamforming)       

t_axis_span = 2;            % Seconds. Defines the time span of the cyclically refreshed spectrogram
fft_period = 50;            % Number of frames used to calculate one spectrum (Spectral Window Size)
lag = 5;                    % Number of frames between adjacent spectra (Spectral Update Interval/Hop Size)
cutoff_fre = 150;           % Hz Wall filter cutoff frequency

x_loc_real = 0.00635;       % Sample line center point physical coordinate x
z_loc_real = 0.011183;      % Sample line center point physical coordinate z


is_save_BF = 1;             % Whether to save beamformed IQ data. 1 for save, 0 for do not save
Bmode_save_index = 1:5;     % Frame indices to save for each data bin
                            % Only effective when is_save_BF = 1. Example: 1 means save only the 1st frame, 1:numperfile means save from 1st to last frame
is_save_spectral = 1;       % Whether to save spectral data. 1 for save, 0 for do not save

is_fftshift = 1;            % Whether to shift the zero-frequency component to the center of the spectrum
is_velocity = 1;            % Spectrogram y-axis display. 1 to display flow velocity, 0 to display frequency shift
reference_angle = 60;       % Plane wave and blood flow reference angle

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

% Wall filtering
fre_max = prf/2; %nyquist
filter_order = 4;
[b3, a3] = butter(filter_order,cutoff_fre/fre_max,"high");
% initial state
b3 = b3(:);a3 = a3(:);
A = eye(filter_order) - [-a3(2:end), [eye(filter_order-1); zeros(1, filter_order-1)]];
zi = A \ (b3(2:end) - b3(1)*a3(2:end));
% Pass to GPU
single_iir_b = single(b3);
single_iir_a = single(a3);
single_iir_zi = single(zi);
single_iir_b_ptr = libpointer('singlePtr', single_iir_b);
single_iir_a_ptr = libpointer('singlePtr', single_iir_a);
single_iir_zi_ptr = libpointer('singlePtr', single_iir_zi);
single_iir_len = filter_order;

% 
buffer_num = floor(fft_period/(numperfile/SteeringNum))+2; % Number of buffers
t_idx = fft_period; % This idx is the current index in the spectrogram, initially fft_period 
t_buffer_idx = fft_period; % This idx is the index in the buffer, initially fft_period 


%% Coordinates of Pixels Used for Spectrum Calculation

% Plane wave transmit angle
tx_angle_deg = steering_deg;
tx_angle_rad = deg2rad(tx_angle_deg);
ratio_xz = tan(tx_angle_rad);
x_step = mean(diff(scan.xx(1,:)));
z_step = mean(diff(scan.zz(:,1)));

x_dif = abs(x_axis - x_loc_real);
z_dif = abs(z_axis - z_loc_real);
[~, x_middle] = min(x_dif);
[~, z_middle] = min(z_dif);

z_pixel_shift = -30:30;
x_pixel_shift = []; 
for z_range = z_pixel_shift
    x_range = round(z_range*z_step*ratio_xz/x_step);
    x_pixel_shift = [x_pixel_shift x_range];
end

x_loc_all = [];
z_loc_all = [];
for x_width = -2:2
    for shift_idx = 1:numel(x_pixel_shift)
        x_loc = x_middle + x_width + x_pixel_shift(shift_idx);
        x_loc_all = [x_loc_all x_loc];
        z_loc = z_middle + z_pixel_shift(shift_idx);
        z_loc_all = [z_loc_all z_loc];
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

bprocesspara.Demod_AFE_Dynamic = -5;
gpu_handle = calllib('US_APP', 'initializespectralGPU', ...
    prf, ...
    numperfile, ...
    buffer_num, ...
    fft_period, ...
    lag, ...
    t_axis_span, ...
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
    single_iir_b_ptr, single_iir_a_ptr,single_iir_zi_ptr, single_iir_len);


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
fig_left = screen_width * 0.1;
fig_bottom = screen_height * (1-fig_height_ratio)/2;

hFig1 = figure('Name',"Bmode",'Position', [fig_left, fig_bottom, fig_width, fig_height]);
hIm1 = imagesc(x_axis, zz_cut, zeros(Nz_cut,BeamN)-100,[-60 0]);
colormap(gray);title("Plane Wave Imaging");
axis equal;axis tight

hold on
for loc_idx = 1:loc_num
    scatter(x_axis(x_loc_all(loc_idx)), z_axis(z_loc_all(loc_idx)), 2, 'g', 'filled');
end

% Spectrogram to display
if is_fftshift
    y_axis = -prf/2: prf/fft_period : prf/2-prf/fft_period;
    if is_velocity
        y_axis = y_axis.*(sos/(2*fc*cos(deg2rad(reference_angle))))*100;
    end
else
    y_axis = 0: prf/fft_period : prf-prf/fft_period;
    if is_velocity
        y_axis = y_axis.*(sos/(2*fc*cos(deg2rad(reference_angle))))*100;
    end
end
y_axis_num = fft_period;
t_axis_num = t_axis_span*prf/lag;
t_axis = (1:t_axis_num)*(lag/prf);
spectral_matrix = single(zeros(y_axis_num, t_axis_num));

% Set window size
fig_width_ratio = 0.4;
fig_height_ratio = 0.3;
fig_width = screen_width * fig_width_ratio;
fig_height = screen_height * fig_height_ratio;
% Set window position
fig_left = screen_width * 0.5;
fig_bottom = screen_height * (1-fig_height_ratio)/2;

hFig2 = figure('Name',"spectral doppler",'Position', [fig_left, fig_bottom, fig_width, fig_height]);
hIm2 = imagesc(t_axis, y_axis, spectral_matrix);
set(gca, 'YDir', 'normal')
colormap(gray);title("Spectral Doppler Imaging");
xlabel('Time (s)');
if is_velocity
    ylabel('Velocity (cm/s)');
else
    ylabel('Frequency Shift (Hz)');
end
spectral_matrix_ptr = libpointer('singlePtr', spectral_matrix);


% Create directory if saving is enabled
if is_save_BF
    if ~exist(fullfile(filefolder, 'bfdata'), 'dir')
        mkdir(fullfile(filefolder, 'bfdata'));
    end
end
if is_save_spectral
    if ~exist(fullfile(filefolder, 'spectral'), 'dir')
        mkdir(fullfile(filefolder, 'spectral'));
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
    ret = calllib('US_APP', 'processSpectralDataBeamformingandPostGPU', gpu_handle, alldata_ptr, ...
        alldata_len, bfdata_ptr, spectral_matrix_ptr, bag_idx); 
    
    rev = bfdata_ptr.Value;
    bfdata_real = reshape(rev(1:2:end),BF_SampleN,BeamN,[]);
    bfdata_imag = reshape(rev(2:2:end),BF_SampleN,BeamN,[]);
    bfdata_iq = bfdata_real + 1i* bfdata_imag;
    hIm1.CData = log_compressed(abs(bfdata_iq(valid_indices,:,end)));
    pause(0.001)

    if is_save_BF
        bfdata_save = bfdata_iq(:,:,Bmode_save_index);
        save(fullfile(filefolder,'bfdata',num2str(bag_idx)+".mat"), 'bfdata_save','x_axis','z_axis');
    end

    spectral_matrix_update = reshape(spectral_matrix_ptr.Value, y_axis_num, t_axis_num);
    if is_fftshift
        mid_row = ceil(fft_period / 2);
        spectral_matrix_update = [spectral_matrix_update(mid_row+1:end, :); 
                                  spectral_matrix_update(1:mid_row, :)];
    end    
    hIm2.CData = log_compressed(spectral_matrix_update);
    pause(0.001)

    if is_save_spectral
        spectral_matrix = flipud(spectral_matrix_update);
        save(fullfile(filefolder,'spectral',num2str(bag_idx)+".mat"), 'spectral_matrix','t_axis','y_axis');
    end

    fclose(fileID);
end


%% Unload

% Release GPU memory and RAM
calllib('US_APP', 'deleteBeamformingGPUHandle', gpu_handle);

% Unload DLL
unloadlibrary('US_APP');







