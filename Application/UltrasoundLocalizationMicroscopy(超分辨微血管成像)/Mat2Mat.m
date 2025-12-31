% 函数功能：
% 本脚本旨在将原始分段保存的波束合成后IQ数据（通常每文件帧数较少或不固定），
% 按照后续算法（如多普勒计算、ULM成像）所需的指定帧长度（framenum），
% 重新拼接、裁切并打包成新的数据文件序列。

clear all
clc
close all
% 加载当前环境变量
currentPath = pwd;
parentDir = fileparts(fileparts(fileparts(currentPath)));
addpath(genpath(parentDir));


% 波束合成数据路径（读取数据的路径）
data_filepath = 'D:\software_matlab\exampledata\ULM\20250703165301\bfiq';

% 波束合成数据路径（保存数据的路径）
data_save_filepath = 'D:\software_matlab\exampledata\ULM\20250703165301\bfiq_com';
if ~exist(data_save_filepath,'dir')
        mkdir(data_save_filepath);
end


%% 获取数据文件列表
[load_file_start_idx,min_num,max_num,sorted_files] = getfiles_mat(data_filepath);

% 获取尺寸
load(fullfile(sorted_files(1).folder, sorted_files(1).name))
[H,W,frameperfile] = size(bfdata_iq);

%% 读取拼接
% 需要多少帧计算多普勒
framenum = 100;

% 计算从开始索引到结束一共有多少个文件
total_available_files = max_num - load_file_start_idx + 1;

% 计算总帧数
total_frames_all = total_available_files * frameperfile;

% 计算可以完整拼凑的包数 (向下取整，不足一包的丢弃)
num_packages = floor(total_frames_all / framenum);

fprintf('可用总帧数: %d\n', total_frames_all);
fprintf('每包帧数: %d\n', framenum);
fprintf('预计生成包数: %d\n', num_packages);

if num_packages == 0
    error("数据总量不足以拼凑一个完整的包 (%d 帧)", framenum);
end

%% 循环打包处理
disp("---------------------------------")
disp("开始分包处理...")

% 外层循环：对应每一个要生成的包 (Package)
for pkg_idx = 1 : num_packages
    
    % 1. 预分配当前包的内存
    IQ = (zeros(H, W, framenum));
    
    % 2. 计算当前包在全局时间轴上的 起始帧 和 结束帧 (绝对索引)
    global_req_start = (pkg_idx - 1) * framenum + 1;
    global_req_end   = pkg_idx * framenum;
    
    % 3. 计算这一包数据跨越了哪些文件 (相对于 sorted_files 的索引)
    % 索引从0开始算，方便取余和整除
    file_rel_idx_start = floor((global_req_start - 1) / frameperfile);
    file_rel_idx_end   = floor((global_req_end   - 1) / frameperfile);
    
    % 内层循环：遍历覆盖当前包的所有源文件
    for f_rel = file_rel_idx_start : file_rel_idx_end
        
        % 定位实际文件
        % load_file_start_idx 是 getfiles_mat 返回的起始文件在列表中的位置
        % 注意：这里假设 sorted_files 是按顺序排列的
        current_file_list_idx = load_file_start_idx + f_rel; 
        
        % 文件地址
        file_path = fullfile(sorted_files(current_file_list_idx-min_num+1).folder, sorted_files(current_file_list_idx-min_num+1).name);
        
        % 加载源文件
        % disp(['  -> 读取源文件片段: ', sorted_files(current_file_list_idx).name]);
        load(file_path, 'bfdata_iq');
        
        % --- 核心逻辑：计算剪切和拼接的索引 ---
        
        % 当前文件包含的全局帧范围
        file_global_start = f_rel * frameperfile + 1;
        file_global_end   = (f_rel + 1) * frameperfile;
        
        % 计算当前文件与当前包需求的 交集范围
        overlap_start = max(global_req_start, file_global_start);
        overlap_end   = min(global_req_end, file_global_end);
        
        % 计算交集长度
        len = overlap_end - overlap_start + 1;
        
        if len > 0
            % 源数据索引 (在当前读取的 bfdata_iq 中的位置)
            src_idx_start = overlap_start - file_global_start + 1;
            %disp("src_idx_start "+src_idx_start)
            src_idx_end   = src_idx_start + len - 1;
            %disp("src_idx_end "+src_idx_end)
            
            % 目标数据索引 (在当前包 bfiq_com 中的位置)
            dst_idx_start = overlap_start - global_req_start + 1;
            %disp("dst_idx_start "+dst_idx_start)
            dst_idx_end   = dst_idx_start + len - 1;
            %disp("dst_idx_end "+dst_idx_end)
            
            % 填入数据
            IQ(:, :, dst_idx_start:dst_idx_end) = bfdata_iq(:, :, src_idx_start:src_idx_end);
        end
    end
    
    % 4. 保存当前包
    % 命名格式：bfiq_com_1.mat, bfiq_com_2.mat ...
    save_filename = sprintf('bfiq_com_%03d.mat', pkg_idx);
    save_fullpath = fullfile(data_save_filepath, save_filename);
    
    % 保存变量，注意这里保存的变量名是 bfiq_com
    % 如果有 x_axis 和 z_axis，也可以一起保存，假设它们是不变的
    if exist('x_axis', 'var') && exist('z_axis', 'var')
        save(save_fullpath, "IQ", "x_axis", "z_axis");
    else
        save(save_fullpath, "IQ");
    end
    
    fprintf('已保存第 %d 包: %s (包含 %d 帧)\n', pkg_idx, save_filename, framenum);
    
end

disp("---------------------------------")
disp("所有完整数据包处理完毕。");



