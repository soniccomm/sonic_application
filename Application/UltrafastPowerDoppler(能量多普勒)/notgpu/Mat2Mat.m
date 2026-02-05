% Functionality:
% Read beamformed IQ data and assemble it into a matrix with a specified number of frames

clear all
clc
close all
% Load current environment variables
currentPath = pwd;
parentDir = fileparts(fileparts(fileparts(currentPath)));
addpath(genpath(parentDir));


%% Beamformed data path (this path already contains beamformed data)

data_filepath = 'D:\software_matlab\exampledata\doppler\20251020190356\bfiq';

%% Get data file list
[load_file_start_idx,min_num,max_num,sorted_files] = getfiles_mat(data_filepath);

% Get dimensions
load(fullfile(sorted_files(1).folder, sorted_files(1).name))
[H,W,frameperfile] = size(bfdata_iq);

%% Read and Concatenate
% How many frames are needed to calculate Doppler
framenum = 200;

% Calculate how many .mat files need to be read
need_filenum = 1;
while (need_filenum*frameperfile<framenum)
    need_filenum = need_filenum + 1;
end

% Pre-allocation
bfiq_com = single(zeros(H, W, need_filenum*frameperfile));

% Check
if (load_file_start_idx+need_filenum-1)>max_num
    error("Insufficient data")
end

disp("Reading started...")
for file_i = load_file_start_idx+1-min_num:load_file_start_idx+1-min_num+need_filenum-1
    disp(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name))
    load(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name))
    idx = file_i - load_file_start_idx + min_num;
    bfiq_com(:,:,(idx-1)*frameperfile+1:(idx)*frameperfile) = bfdata_iq;
end
disp("Reading completed...")

% Crop
bfiq_com = bfiq_com(:,:,1:framenum);

% Save
save(fullfile(data_filepath, "bfiq_com.mat"),"bfiq_com","x_axis","z_axis")
disp("Saving completed...")

