% This script is the fourth step of Ultrafast Power Doppler
% Step 1 is data acquisition
% Step 2 is Beamforming Imaging PlaneWave_ADC_BF_IQ_Analysis
% Step 3 is Matrix Concatenation Mat2Mat

clear;clc
% Change this path to the path of the data saved in step 3
load('D:\software_matlab\exampledata\doppler\20251020190356\bfiq\bfiq_com.mat')

[Nz,Nx,Nt] = size(bfiq_com);
% Parameter settings: ord1 and ord2
% Use SVD for filtering, ord1 and ord2 are SVD filter parameters,
% ord1 controls how much tissue signal is filtered out, ord2 controls how much noise signal is filtered out

% check beamforming data
figure(1); 
for i = 1:size(bfiq_com,3)
    img_envelope = abs(bfiq_com(:,:,i));
    img_log = log_compressed(img_envelope);
    imagesc(x_axis,z_axis,img_log,[-60 0]);
    colormap(gray);axis equal;axis tight
    title("frame"+i)
    pause(0.01);
end


%% svd filter
img_c = reshape(bfiq_com,[],size(bfiq_com,3));
IQ_mat = img_c(:,:);
casorati_mat = IQ_mat;
[U,S,V] = svd(casorati_mat,'econ');
clear casorati;

% SVD orders
ord1 = 55;
ord2 = size(IQ_mat,2)-20;
S_filt = zeros(size(S));
for i = ord1:ord2
    S_filt(i,i) = S(i,i);
end
filtered_mat = U*S_filt*V';


%% power doppler
power_doppler = mean(abs(filtered_mat).^2,2);
power_doppler = reshape(power_doppler,[Nz Nx]);
figure(2);imagesc(x_axis,z_axis,log_compressed(power_doppler),[-25 0]);
colormap('hot');axis equal;axis tight;colorbar


