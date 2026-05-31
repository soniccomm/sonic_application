% 此脚本为超快能量多普勒第四步
% 第一步为采集数据
% 第二步为波束合成成像DivergingWave_ADC_BF_IQ_Analysis
% 第三步为矩阵拼接Mat2Mat

clear;clc
% 将此路径改为第三步保存的数据的路径
load('E:\20260418154739\bfiq\bfiq_com.mat')

[Nz,Nx,Nt] = size(bfiq_com);
% 参数设置：ord1和ord2
% 使用svd进行滤波，ord1和ord2为svd滤波的参数，
% ord1用于控制滤掉多少组织信号，ord2用于控制滤掉多少噪声信号

% check beamforming data
figure(1); 
for i = 1:size(bfiq_com,3)
    pcolor(xx_grid*100,zz_grid*100,log_compressed(abs(bfiq_com(:,:,i))));
    colormap(gray);clim([-60, 0]);xlabel("cm");ylabel("cm")
    axis equal;axis tight
    set(gca, 'YDir', 'reverse');
    ylim([0,max(zz_grid(:))*100])
    shading interp;
    title("frame"+i)
    colorbar
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
figure(2);
pcolor(xx_grid*100,zz_grid*100,log_compressed(power_doppler));
clim([-60, 0]);xlabel("cm");ylabel("cm")
colormap('hot');axis equal;axis tight;colorbar
set(gca, 'YDir', 'reverse');
ylim([0,max(zz_grid(:))*100])
shading interp;





