% 函数功能：
% 读取波束合成后的IQ数据，拼成指定帧数的矩阵

clear all
clc
close all
% 加载当前环境变量
currentPath = pwd;
parentDir = fileparts(fileparts(fileparts(currentPath)));
addpath(genpath(parentDir));


%% 波束合成数据路径（该路径已有波束合成数据）

data_filepath = 'D:\software_matlab\exampledata\大鼠脑\bfiq';

%% 获取数据文件列表
[load_file_start_idx,min_num,max_num,sorted_files] = getfiles_mat(data_filepath);

% 获取尺寸
load(fullfile(sorted_files(1).folder, sorted_files(1).name))
[H,W,frameperfile] = size(bfdata_iq);

%% 读取拼接
% 需要多少帧计算多普勒
framenum = 400;

% 计算需要读取多少mat
need_filenum = 1;
while (need_filenum*frameperfile<framenum)
    need_filenum = need_filenum + 1;
end

% 预分配
bfiq_com = single(zeros(H, W, need_filenum*frameperfile));

% 判断
if (load_file_start_idx+need_filenum-1)>max_num
    error("数据不足")
end

disp("读取开始...")
for file_i = load_file_start_idx+1-min_num:load_file_start_idx+1-min_num+need_filenum-1
    disp(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name))
    load(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name))
    idx = file_i - load_file_start_idx + min_num;
    bfiq_com(:,:,(idx-1)*frameperfile+1:(idx)*frameperfile) = bfdata_iq;
end
disp("读取完毕...")

% 裁剪
bfiq_com = bfiq_com(:,:,1:framenum);

% 保存
save(fullfile(data_filepath, "bfiq_com.mat"),"bfiq_com","x_axis","z_axis")
disp("保存完毕...")

