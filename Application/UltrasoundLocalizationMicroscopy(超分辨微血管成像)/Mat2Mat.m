% Function Functionality:
% This script aims to take the original segmented saved beamformed IQ data (usually with few or variable frames per file), and according to the specified 
% frame length (framenum) required by subsequent algorithms (like Doppler calculation, ULM imaging), 
% re-stitch, crop, and repackage them into a new data file sequence.

clear all
clc
close all
% Load current environment variables
currentPath = pwd;
parentDir = fileparts(fileparts(fileparts(currentPath)));
addpath(genpath(parentDir));


% Beamformed data path (path for reading data)
data_filepath = 'D:\software_matlab\exampledata\ULM\20250703165301\bfiq';

% Beamformed data path (path for saving data)
data_save_filepath = 'D:\software_matlab\exampledata\ULM\20250703165301\bfiq_com';
if ~exist(data_save_filepath,'dir')
        mkdir(data_save_filepath);
end


%% Get data file list
[load_file_start_idx,min_num,max_num,sorted_files] = getfiles_mat(data_filepath);

% Get dimensions
load(fullfile(sorted_files(1).folder, sorted_files(1).name))
[H,W,frameperfile] = size(bfdata_iq);

%% Read and Stitch
% Number of frames needed for Doppler calculation
framenum = 100;

% Calculate total number of files from start index to end
total_available_files = max_num - load_file_start_idx + 1;

% Calculate total frames
total_frames_all = total_available_files * frameperfile;

% Calculate number of complete packages (round down, discard incomplete packages)
num_packages = floor(total_frames_all / framenum);

fprintf('Total available frames: %d\n', total_frames_all);
fprintf('Frames per package: %d\n', framenum);
fprintf('Estimated packages to generate: %d\n', num_packages);

if num_packages == 0
    error("Insufficient total data to assemble a complete package (%d frames)", framenum);
end

%% Loop for packaging process
disp("---------------------------------")
disp("Starting package processing...")

% Outer loop: Corresponds to each package to be generated
for pkg_idx = 1 : num_packages
    
    % 1. Pre-allocate memory for the current package
    IQ = (zeros(H, W, framenum));
    
    % 2. Calculate Start Frame and End Frame on global timeline for current package (absolute index)
    global_req_start = (pkg_idx - 1) * framenum + 1;
    global_req_end   = pkg_idx * framenum;
    
    % 3. Calculate which files this package spans (relative index to sorted_files)
    % Index starts from 0 for easier modulo and division operations
    file_rel_idx_start = floor((global_req_start - 1) / frameperfile);
    file_rel_idx_end   = floor((global_req_end   - 1) / frameperfile);
    
    % Inner loop: Iterate through all source files covering the current package
    for f_rel = file_rel_idx_start : file_rel_idx_end
        
        % Locate actual file
        % load_file_start_idx is the position of the starting file in the list returned by getfiles_mat
        % Note: Assuming sorted_files are arranged in order
        current_file_list_idx = load_file_start_idx + f_rel; 
        
        % File path
        file_path = fullfile(sorted_files(current_file_list_idx-min_num+1).folder, sorted_files(current_file_list_idx-min_num+1).name);
        
        % Load source file
        % disp(['  -> Reading source file segment: ', sorted_files(current_file_list_idx).name]);
        load(file_path, 'bfdata_iq');
        
        % --- Core Logic: Calculate cropping and stitching indices ---
        
        % Global frame range contained in the current file
        file_global_start = f_rel * frameperfile + 1;
        file_global_end   = (f_rel + 1) * frameperfile;
        
        % Calculate intersection range between current file and current package requirements
        overlap_start = max(global_req_start, file_global_start);
        overlap_end   = min(global_req_end, file_global_end);
        
        % Calculate intersection length
        len = overlap_end - overlap_start + 1;
        
        if len > 0
            % Source data indices (position in the currently read bfdata_iq)
            src_idx_start = overlap_start - file_global_start + 1;
            %disp("src_idx_start "+src_idx_start)
            src_idx_end   = src_idx_start + len - 1;
            %disp("src_idx_end "+src_idx_end)
            
            % Destination data indices (position in the current package IQ)
            dst_idx_start = overlap_start - global_req_start + 1;
            %disp("dst_idx_start "+dst_idx_start)
            dst_idx_end   = dst_idx_start + len - 1;
            %disp("dst_idx_end "+dst_idx_end)
            
            % Fill data
            IQ(:, :, dst_idx_start:dst_idx_end) = bfdata_iq(:, :, src_idx_start:src_idx_end);
        end
    end
    
    % 4. Save current package
    % Naming format: bfiq_com_1.mat, bfiq_com_2.mat ...
    save_filename = sprintf('bfiq_com_%03d.mat', pkg_idx);
    save_fullpath = fullfile(data_save_filepath, save_filename);
    
    % Save variables, note that the variable saved here is IQ
    % If x_axis and z_axis exist, save them together, assuming they are invariant
    if exist('x_axis', 'var') && exist('z_axis', 'var')
        save(save_fullpath, "IQ", "x_axis", "z_axis");
    else
        save(save_fullpath, "IQ");
    end
    
    fprintf('Saved package %d: %s (contains %d frames)\n', pkg_idx, save_filename, framenum);
    
end

disp("---------------------------------")
disp("All complete data packages processed.");