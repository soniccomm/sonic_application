% 此脚本为超快能量多普勒第四步
% 第一步为采集数据
% 第二步为波束合成成像PlaneWave_ADC_BF_IQ_Analysis
% 第三步为矩阵拼接Mat2Mat

clear;clc
% 将此路径改为第三步保存的数据的路径
load('D:\software_matlab\exampledata\doppler\20251020190356\bfiq\bfiq_com.mat')

[Nz,Nx,Nt] = size(bfiq_com);
% 参数设置：ord1和ord2
% 使用svd进行滤波，ord1和ord2为svd滤波的参数，
% ord1用于控制滤掉多少组织信号，ord2用于控制滤掉多少噪声信号

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

% svd阶数
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


