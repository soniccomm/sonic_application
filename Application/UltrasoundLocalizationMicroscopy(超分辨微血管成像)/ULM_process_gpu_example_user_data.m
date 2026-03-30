%% IMPORTANT NOTE
% Before running this script:
% 1 Run PlaneWave_Acq_BF_IQ.m​ -> Acquires raw plane-wave channel data.
% 2 Run PlaneWave_ADC_BF_IQ_Analysis.m​ -> Performs beamforming on the raw 
% data to generate beamformed IQ (In-phase/Quadrature) data.
% 3 Run Mat2Mat.m​ -> Consolidates multiple smaller .mat files (each 
% containing a limited number of frames) into larger .mat files with
% a sufficient number of frames. 

%% ========================================================================
% ENVIRONMENT SETUP
% =========================================================================
clear; clc; close all;

% Add paths for the PALA_SRUS toolbox and parent directories
currentPath = pwd;
parentDir = fileparts(fileparts(currentPath));
addpath(genpath(parentDir));

%% ========================================================================
% RENDERING FUNCTION CONTROL FLAGS
% =========================================================================
% Set to 1 to calculate and render specific velocity maps
flag_vel_norm = 1; % Absolute velocity magnitude map
flag_vel_z    = 1; % Axial directional velocity map (Z-axis)

%% ========================================================================
% LOAD DATA & DIRECTORY SETUP
% =========================================================================
dataFolder = 'D:\software_matlab\exampledata\S256-ULM\bfiq';
resultFolder = fullfile(dataFolder, 'results');

% Create the result directory if it does not exist
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);
end

% Create a list of all .mat files in the data folder
IQdataList = {};
IQnum = 1;
fileList = dir([dataFolder, filesep, '*.mat']);
for i = 1:numel(fileList)
    IQdataList{IQnum} = [fileList(i).folder, filesep, fileList(i).name];
    IQnum = IQnum + 1;
end
IQnum = IQnum - 1; % Total number of data blocks


%% ========================================================================
% PARAMETERS CONFIGURATION
% =========================================================================
% Load the first file to extract dimensions
load(IQdataList{1});
bfiq = IQ;

% if micro bubbles are too much, choose ROI to get better tracking
ROI_z = 100:400;
ROI_x = 1:256;

NbFrames = size(bfiq,3);
framerate = 800;  % compounded frame rate = prf / angle num
res = 10;

% Define the ULM parameters structure
ULM = struct( ...                       
    'numberOfParticles', 100,...         % Number of particles extracted per frame (30-100)
    'res', res,...                       % Resolution factor (typically 10)
    'SVD_cutoff', [50, NbFrames - 50],...% SVD cutoff limits for tissue clutter filtering
    'max_linking_distance', 6,...        % Max allowed linking distance between frames (pixels)
    'min_length', 15,...                 % Minimum valid track length (rejects noise)
    'fwhm', [1 1]*3,...                  % Mask size for sub-pixel localization [z, x]
    'max_gap_closing', 0,...             % Max allowed missing frames in a single track
    'size', [size(bfiq(ROI_z,ROI_x,:),1), size(bfiq(ROI_z,ROI_x,:),2), NbFrames],...
    'scale', [mean(diff(z_axis))*1000 mean(diff(x_axis))*1000 1/framerate],...     % Physical scale [z, x, dt]
    'numberOfFramesProcessed', NbFrames,... 
    'interp_factor', 1/res,...           % Interpolation factor for track smoothing
    'butter_CuttofFreq', [50 300],...    % Butterworth filter cutoff frequencies
    'framerate', framerate...            % Frame rate (Hz)
    );

%% ========================================================================
% LOAD DLL & INITIALIZE GPU
% =========================================================================
% Check if the C++ CUDA DLL is loaded; if not, load it
if ~libisloaded('ULM')
    loadlibrary('ULM.dll', 'ULM.h');
end

% Initialize the GPU processing handle with ULM parameters
disp("------ Initializing ULM GPU ------")
gpu_handle = calllib('ULM', 'initializeULMGPU', ...
    ULM.numberOfParticles, ... 
    ULM.res, ...
    ULM.SVD_cutoff(1), ULM.SVD_cutoff(2),...
    ULM.max_linking_distance,...
    ULM.min_length,...
    ULM.fwhm(1), ULM.fwhm(2),...
    ULM.max_gap_closing,...
    ULM.size(1), ULM.size(2), ULM.size(3),...
    ULM.scale(1), ULM.scale(2), ULM.scale(3),...
    ULM.butter_CuttofFreq(1), ULM.butter_CuttofFreq(2),...
    ULM.framerate,...
    flag_vel_norm, flag_vel_z);
disp("------ Initialization Finished ------")

%% ========================================================================
% PROCESS EACH DATA BLOCK AND ACCUMULATE RESULTS
% =========================================================================
% Initialize accumulation matrices for the high-resolution grid
Density_result = zeros(length(ROI_z)*ULM.res, length(ROI_x)*ULM.res);
if(flag_vel_norm)
    Velocity_result_norm = zeros(length(ROI_z)*ULM.res, length(ROI_x)*ULM.res);
end
if(flag_vel_z)
    Velocity_result_z = zeros(length(ROI_z)*ULM.res, length(ROI_x)*ULM.res);
end

disp("------ Processing Start ------")

for idx = 1:IQnum

    load(IQdataList{idx});
    disp("Processing " + IQdataList{idx})
    bfiq = IQ(ROI_z, ROI_x, :); % Crop to ROI

    % Extract real part and create a C-style single-precision pointer
    single_data_real = single(real(bfiq));
    single_data_real = reshape(single_data_real, [], 1);
    single_data_real_ptr = libpointer('singlePtr', single_data_real);

    % Extract imaginary part and create a C-style pointer
    single_data_imag = single(imag(bfiq));
    single_data_imag = reshape(single_data_imag, [], 1);
    single_data_imag_ptr = libpointer('singlePtr', single_data_imag);

    % Prepare a pointer with enough memory to receive the output from GPU
    single_rev = zeros(3 * size(single_data_imag, 1), 1);
    single_data_rev_ptr = libpointer('singlePtr', single_rev);

    % Call GPU DLL (Step 5: Render Maps)
    ret = calllib('ULM', 'processULMGPU', gpu_handle, ...
        single_data_real_ptr, single_data_imag_ptr, single_data_rev_ptr, 5);
    
    aa = single_data_rev_ptr.Value;
    
    % 1. Parse matrix dimensions
    Nz_out = double(aa(1));
    Nx_out = double(aa(2));
    total_pixels = Nz_out * Nx_out;
    ptr = double(3);
    
    % 2. Extract and accumulate Density Map (always exists)
    MatOut_1D = aa(ptr : ptr + total_pixels - 1);
    MatOut = reshape(MatOut_1D, Nz_out, Nx_out); 
    ptr = ptr + total_pixels;

    Density_result = Density_result + MatOut;
    
    % 3. Extract and accumulate Absolute Velocity Map (if enabled)
    if flag_vel_norm
        MatOut_vel_norm_1D = aa(ptr : ptr + total_pixels - 1);
        MatOut_vel_norm = reshape(MatOut_vel_norm_1D, Nz_out, Nx_out);
        ptr = ptr + total_pixels;

        Velocity_result_norm = Velocity_result_norm + MatOut_vel_norm;
    end
    
    % 4. Extract and accumulate Axial Velocity Map (if enabled)
    if flag_vel_z
        MatOut_vel_z_1D = aa(ptr : ptr + total_pixels - 1);
        MatOut_vel_z = reshape(MatOut_vel_z_1D, Nz_out, Nx_out);
        ptr = ptr + total_pixels;

        Velocity_result_z = Velocity_result_z + MatOut_vel_z;
    end
end

disp("------ Processing Finished ------")


%% ========================================================================
% VISUALIZATION: FIGURE 1 - DENSITY MAP
% =========================================================================
figure(1); clf, set(gcf, 'Position', [652 393 941 585]);

% Apply non-linear power compression to enhance the visibility of small capillaries
IntPower = 0.4; 
SigmaGauss = 0;
im = imagesc(x_axis(ROI_x)*1000, z_axis(ROI_z)*1000, Density_result.^IntPower); 
axis image;

% Apply Gaussian smoothing if specified
if SigmaGauss > 0, im.CData = imgaussfilt(im.CData, SigmaGauss); end

title('ULM Intensity Display (Density Map)')
colormap(gca, hot(128))

% Add a colorbar and saturate the top 20% to make vessels pop out
clbar = colorbar; caxis(caxis * .8);  
clbar.Label.String = 'Number of counts';

% Convert the colorbar ticks back to linear counts (undoing the power compression)
clbar.TickLabels = round(clbar.Ticks.^(1/IntPower), 1);
xlabel('mm'); ylabel('mm'); axis equal; axis tight;

% --- Save Figure 1 ---
fig1_path = fullfile(resultFolder, 'ULM_DensityMap');
exportgraphics(figure(1), [fig1_path, '.png'], 'Resolution', 300); % High-res PNG
saveas(figure(1), [fig1_path, '.fig']);                            % Raw MATLAB Figure


%% ========================================================================
% VISUALIZATION: FIGURE 2 - AXIAL DIRECTIONAL VELOCITY MAP
% =========================================================================
if (flag_vel_z)
    % Calculate average velocity per pixel and handle division by zero
    Velocity_z = Velocity_result_z ./ Density_result;
    Velocity_z(isnan(Velocity_z)) = 0;

    figure(2); clf, set(gcf, 'Position', [652 393 941 585]);
    
    % Create a custom diverging colormap (Blue -> Black -> Red)
    velColormap = cat(1, flip(flip(hot(128), 1), 2), hot(128)); 
    velColormap = velColormap(15:end-15, :); % Remove extreme white parts
    
    IntPower = 0.25;
    % Multiply the compressed density by the sign of the smoothed velocity 
    % to encode flow direction (positive = upward, negative = downward)
    im = imagesc(x_axis(ROI_x)*1000, z_axis(ROI_z)*1000, ...
        (Density_result).^IntPower .* sign(imgaussfilt(Velocity_z, .8)));
    
    % Shift data slightly towards zero to filter out near-zero noise
    im.CData = im.CData - sign(im.CData)/2; 
    axis equal; axis tight;
    
    title('ULM Intensity Display with Axial Flow Direction')
    colormap(gca, velColormap)
    caxis([-1 1] * max(caxis) * .7) % Add saturation
    
    clbar = colorbar; clbar.Label.String = 'Flow Direction & Intensity';
    clbar.TickLabels = round(clbar.Ticks.^(1/IntPower), 1);
    ca = gca; ca.Position = [.05 .05 .8 .9];

    % --- Save Figure 2 ---
    fig2_path = fullfile(resultFolder, 'ULM_Velocity_Z_Direction');
    exportgraphics(figure(2), [fig2_path, '.png'], 'Resolution', 300);
    saveas(figure(2), [fig2_path, '.fig']);
end


%% ========================================================================
% VISUALIZATION: FIGURE 3 - ABSOLUTE VELOCITY MAGNITUDE MAP
% =========================================================================
if (flag_vel_norm)
    % Calculate average velocity per pixel
    Velocity_n = Velocity_result_norm ./ Density_result;
    Velocity_n(isnan(Velocity_n)) = 0;
    
    % Determine the max display velocity based on the 98th percentile (ignores outliers)
    vmax_disp = ceil(quantile(Velocity_n(abs(Velocity_n)>0), .98) / 10) * 10;
    
    figure(3); clf, set(gcf, 'Position', [652 393 941 585]);
    clbsize = [180, 50]; % Dimensions for the 2D velocity-density colorbar
    
    % Normalize velocity to [0, 1] range
    Mvel_rgb = Velocity_n / vmax_disp; 
    
    % Manually draw a velocity colorbar at the top-left corner
    Mvel_rgb(1:clbsize(1), 1:clbsize(2)) = repmat(linspace(1,0,clbsize(1))', 1, clbsize(2)); 
    
    % Apply non-linear compression to the velocity values
    Mvel_rgb = Mvel_rgb.^(1/1.5); 
    Mvel_rgb(Mvel_rgb > 1) = 1;
    Mvel_rgb = imgaussfilt(Mvel_rgb, .5);
    
    % Convert the 2D normalized velocity matrix into a 3D RGB image using 'jet'
    Mvel_rgb = ind2rgb(round(Mvel_rgb * 256), jet(256)); 
    
    % Create a density shadow mask (modulates brightness based on microbubble count)
    MatShadow = Density_result; 
    MatShadow = MatShadow ./ max(MatShadow(:) * .3); % Overexpose dense regions
    MatShadow(MatShadow > 1) = 1;
    
    % Manually draw a density gradient colorbar for the X-axis of the top-left legend
    MatShadow(1:clbsize(1), 1:clbsize(2)) = repmat(linspace(0,1,clbsize(2)), clbsize(1), 1);
    
    IntPower = 1/3;
    % Multiply the RGB velocity map by the compressed density shadow mask
    % Background becomes black, while dense vessels become bright and colored
    Mvel_rgb = Mvel_rgb .* (MatShadow.^IntPower);
    Mvel_rgb = brighten(Mvel_rgb, .4); % Brighten the entire image
    
    imagesc(x_axis(ROI_x)*1000, z_axis(ROI_z)*1000, Mvel_rgb);
    axis on; axis equal; axis tight
    title(['Velocity Magnitude (0 - ' num2str(vmax_disp) ' mm/s)'])
    ca = gca; ca.Position = [.05 .05 .8 .9];

    % --- Save Figure 3 ---
    fig3_path = fullfile(resultFolder, 'ULM_Velocity_Magnitude');
    exportgraphics(figure(3), [fig3_path, '.png'], 'Resolution', 300);
    saveas(figure(3), [fig3_path, '.fig']);
end


%% ========================================================================
% SAVE MATRIX DATA (.mat)
% =========================================================================
disp("------ Saving Matrix Data ------")
mat_save_path = fullfile(resultFolder, 'ULM_Final_Data.mat');

% Dynamically build the list of variables to save based on the active flags
varsToSave = {'Density_result', 'x_axis', 'z_axis', 'ROI_x', 'ROI_z', 'framerate', 'res'};
if flag_vel_z
    varsToSave = [varsToSave, {'Velocity_result_z', 'Velocity_z'}];
end
if flag_vel_norm
    varsToSave = [varsToSave, {'Velocity_result_norm', 'Velocity_n', 'vmax_disp'}];
end

% Save to a .mat file for future fast loading
save(mat_save_path, varsToSave{:});
disp(['All results saved successfully to: ', resultFolder]);


%% ========================================================================
% UNLOAD DLL AND CLEANUP
% =========================================================================
% Release GPU memory handles and unload the dynamic library
calllib('ULM', 'deleteULMGPUHandle', gpu_handle);
unloadlibrary('ULM')
disp("------ ULM Processing Completely Finished ------")