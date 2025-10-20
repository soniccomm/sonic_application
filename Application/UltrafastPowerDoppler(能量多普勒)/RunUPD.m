% 此脚本为超快能量多普勒第三步（第一步为采集数据，第二步为波束合成成像PlaneWave_bf_linear/PlaneWave_bf_convex）

% 将此路径改为第二步波束合成保存的数据的路径
load('D:\software_matlab\exampledata\dog_kidney_ADC_data\20250624164838\bfiq\bfiq_1.mat')

% 参数设置：ord1和ord2
% 使用svd进行滤波，ord1和ord2为svd滤波的参数，
% ord1用于控制滤掉多少组织信号，ord2用于控制滤掉多少噪声信号

% check beamforming data
figure(1); 
for i = 1:size(bfiq,3)
    img_envelope = abs(bfiq(:,:,i));
    img_log = log_compressed(img_envelope);
    imshow(img_log,[-60 0]);
    pause(0.02);
end


%% svd filter
img_c = reshape(bfiq,[],size(bfiq,3));
IQ_mat = img_c(:,:);
casorati_mat = IQ_mat;
[U,S,V] = svd(casorati_mat,'econ');
clear casorati;


%%
ord1 = 20;
ord2 = size(IQ_mat,2)-20;
S_filt = zeros(size(S));
for i = ord1:ord2
    S_filt(i,i) = S(i,i);
end
filtered_mat = U*S_filt*V';
% power doppler
power_doppler = mean(abs(filtered_mat).^2,2);
power_doppler = reshape(power_doppler,[Nz Nx]);
figure;imshow(log_compressed(power_doppler),[-60 0]);colormap('hot');


