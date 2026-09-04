%% PALA_ULM_CustomData.m
% =========================================================================
% 超声超分辨定位显微成像 (ULM) 处理脚本
% adapt from PALA_InVivoULM_example
% 
% 核心流程：
%   1. 读取自采 Beamformed IQ 数据块 (Data Blocks)
%   2. 杂波抑制滤波 (SVD 奇异值滤波 + 时域带通滤波)
%   3. 微泡亚像素中心定位 (默认径向对称法 Radial)
%   4. 帧间微泡运动追踪关联 (Tracking & Pairing)
%   5. 超分辨率网格积累与流速图解算 (Track to Super-Resolution Matrix)
%   6. 多模态可视化渲染 (微血管密度图、轴向流向图、绝对流速图)
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------
%% 0. 环境与工具箱路径初始化
% -------------------------------------------------------------
currentPath = pwd;
parentDir   = fileparts(currentPath);
addpath(genpath(parentDir)); % 递归添加 PALA 工具箱的所有子函数库

% -------------------------------------------------------------
%% 1. 数据路径与系统物理参数配置 (用户核心适配区)
% -------------------------------------------------------------
data_dir   = 'D:\software_matlab\exampledata\R256-ULM\bfiq\';           % 存放 IQ 数据块 (.mat) 的文件夹
save_dir   = 'D:\software_matlab\exampledata\R256-ULM\bfiq\ULM_Results\'; % 处理结果及图像的输出路径
file_list  = dir(fullfile(data_dir, '*.mat'));                         % 检索所有数据文件
assert(~isempty(file_list), '未在数据目录中检索到任何 .mat 文件，请检查路径！');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

% --- 超声系统物理参数 ---
c0       = 1540;            % 组织中平均声速 (m/s)
f0       = 15e6;            % 探头中心发射频率 (Hz)
lambda   = (c0 / f0) * 1e3; % 声学波长 (mm)
PRF      = 1000;            % 脉冲重复频率/超快成像等效帧率 (Hz)

% --- 自动探测数据尺寸与空间网格 ---
sample_data = load(fullfile(file_list(1).folder, file_list(1).name));

% 提取 IQ 矩阵及尺寸: [Nz (轴向/深度点数), Nx (横向点数), Nt (帧数)]
SizeOfBloc = size(sample_data.IQ);

% 读取随数据保存的空间坐标轴 (统一换算为 mm)
z_axis   = sample_data.z_axis * 1000;
x_axis   = sample_data.x_axis * 1000;

% 计算原始图像每个像素对应的物理步长 (mm)
dz       = abs(mean(diff(z_axis))); % 轴向像元分辨率
dx       = abs(mean(diff(x_axis))); % 横向像元分辨率

% -------------------------------------------------------------
%% 2. ULM 算法超参数配置 (Algorithm Parameters)
% -------------------------------------------------------------
% 【分辨率提升倍率】
% res = 10 表示将原始超声网格在轴向和横向上均细分为 10 份
res = 10; 

ULM = struct();
ULM.size     = SizeOfBloc;                % 输入矩阵三维尺寸 [Nz, Nx, Nt]
ULM.scale    = [1, 1, 1/PRF];             % 追踪阶段以像元为基准单位，时间步长为 s
ULM.res      = res;                       % 超分辨率倍率因子
ULM.lambda   = lambda;                    % 波长 (mm)

% 【微泡密度控制】
% numberOfParticles: 单帧图像中提取的最强微泡上限数量
% 调优经验: 浓度高且重叠严重时适当调低 (30~60)；稀疏时可放宽到 (80~120)
ULM.numberOfParticles = 80;

% 【SVD 杂波滤波阶数 [低阶截断, 高阶截断]】
% 组织杂波集中在低阶奇异值分量，高频热噪声集中在高阶分量
ULM.SVD_cutoff = [50, SizeOfBloc(3)];

% 【微泡定位窗口大小 (Point Spread Function 拟合尺寸)】
% 单位: 原始像元，取奇数
pix_fwhm_z = 3; 
pix_fwhm_x = max(3, round(pix_fwhm_z * (dz / dx))); 
if mod(pix_fwhm_x, 2) == 0, pix_fwhm_x = pix_fwhm_x + 1; end 
ULM.fwhm  = [pix_fwhm_x, pix_fwhm_z]; % [横向x像素数, 轴向z像素数]

% 【定位算法选择】
% 'Radial': 径向对称法 (兼顾速度与精度，推荐)
% 'WA': 加权质心法 (Weighted Average)
% 'CurveFitting': 高斯曲面拟合
ULM.LocMethod  = 'Radial';

% 【帧间微泡追踪 (Tracking) 约束】
% max_linking_distance: 连续两帧之间，微泡允许移动的最大位移 (单位: 原始像元)
ULM.max_linking_distance = 3; 

% min_length: 构成一条有效血管轨迹所需的最短持续帧数 (滤除闪烁噪点)
ULM.min_length           = 10; 

% max_gap_closing: 允许微泡中途丢失/未检测到的最大间隔帧数 (0 表示不允许断帧)
ULM.max_gap_closing      = 0;  

% interp_factor: 轨迹内插步长，保证内插轨迹平滑穿过每个超分辨网格
ULM.interp_factor        = 1 / (ULM.max_linking_distance * ULM.res) * 0.8;

% 【时域带通滤波器】
% 进一步滤除低频组织慢漂移和极高频噪声 (30 Hz ~ 300 Hz)
[but_b, but_a] = butter(2, [30, 300] / (PRF / 2), 'bandpass');

% -------------------------------------------------------------
%% 3. 批处理循环：滤波、定位与轨迹追踪
% -------------------------------------------------------------
Nbuffers  = numel(file_list);
Track_tot = cell(Nbuffers, 1); 

fprintf('--- 开始 ULM 批处理: 共 %d 个 Blocks ---\n', Nbuffers);
t_start = tic;

for hhh = 1:Nbuffers
    fprintf('正在处理 Block %d/%d ...\n', hhh, Nbuffers);
    
    % 1. 加载单 Block IQ 数据 (Nz x Nx x Nt 复数矩阵)
    mat_obj = load(fullfile(file_list(hhh).folder, file_list(hhh).name));
    raw_IQ  = mat_obj.IQ;
    
    % 2. 空间 SVD 奇异值分解滤波 (去除静止与慢速组织回波)
    IQ_filt = SVDfilter(raw_IQ, ULM.SVD_cutoff);
    
    % 3. 时域滤波与异常值清理
    IQ_filt = filter(but_b, but_a, IQ_filt, [], 3);
    IQ_filt(~isfinite(IQ_filt)) = 0; 
    
    % 4. 微泡亚像素检测与定位 (输入为包络检波后的强度图 abs(IQ))
    % 输出矩阵 MatTracking: [intensity, z_pix, x_pix, frame_idx]
    MatTracking = ULM_localization2D(abs(IQ_filt), ULM);
    
    if isempty(MatTracking)
        Track_tot{hhh} = {};
        continue;
    end
    
    % 5. 像元尺度微泡配对追踪与轨迹平滑插值
    tracks_pix = ULM_tracking2D(MatTracking, ULM);
    
    % 排除无有效轨迹的情况
    if isempty(tracks_pix) || (numel(tracks_pix) == 1 && size(tracks_pix{1}, 2) < 4)
        Track_tot{hhh} = {};
        continue;
    end
    
    % 6. 将轨迹由像元坐标转换为物理空间单位 (mm 与 mm/s)
    % 存储格式: [z(mm), x(mm), vz(mm/s), vx(mm/s), time(s)]
    Track_tot{hhh} = cellfun(@(t) [ ...
        (t(:, 1) - 1) * dz, ...
        (t(:, 2) - 1) * dx, ...
        t(:, 3) * dz, ...
        t(:, 4) * dx, ...
        t(:, 5) ], tracks_pix, 'UniformOutput', false);
end

% 整合所有数据块的轨迹
Track_tot = Track_tot(~cellfun('isempty', Track_tot));
Track_tot = cat(1, Track_tot{:});
fprintf('批处理完毕！耗时: %.2f 分钟, 共构建 %d 条血管轨迹。\n', toc(t_start)/60, numel(Track_tot));

% 保存轨迹矩阵供离线分析使用
save(fullfile(save_dir, 'ULM_Tracks.mat'), 'Track_tot', 'ULM', '-v7.3');

% -------------------------------------------------------------
%% 4. 超分辨率网格投影 (轨迹转栅格图像)
% -------------------------------------------------------------
fprintf('--- 正在将轨迹栅格化投影至超分辨率像素网格 ---\n');

% 计算超分辨网格像元物理大小 (SRscale) 及矩阵尺寸 (SRsize)
ULM.SRscale = dz / ULM.res; 
ULM.SRsize  = round([SizeOfBloc(1)*dz, SizeOfBloc(2)*dx] / ULM.SRscale);
grid_size   = ULM.SRsize + 1; 

% 将物理连续坐标 (mm) 量化映射为超分辨率网格的整数索引
Track_matout = cellfun(@(x) [ ...
    x(:, 1) / ULM.SRscale + 1, ...
    x(:, 2) / ULM.SRscale + 1, ...
    x(:, 3:end) ], Track_tot, 'UniformOutput', false);

% 通过轨迹穿透累加计算多模态生理学特征图
MatOut      = ULM_Track2MatOut(Track_matout, grid_size);                         % 微泡通量密度 (灌注量)
MatOut_zdir = ULM_Track2MatOut(Track_matout, grid_size, 'mode', '2D_vel_z');   % 轴向速度矢量 (含正负流向)
MatOut_vel  = ULM_Track2MatOut(Track_matout, grid_size, 'mode', '2D_velnorm'); % 标量流动速度大小 (mm/s)

% 物理空间坐标轴向量 (mm)
z_axis_sr = (0:ULM.SRsize(1)) * ULM.SRscale;
x_axis_sr = (0:ULM.SRsize(2)) * ULM.SRscale;

save(fullfile(save_dir, 'ULM_MatOuts.mat'), ...
     'MatOut', 'MatOut_zdir', 'MatOut_vel', 'x_axis_sr', 'z_axis_sr', 'ULM');

% -------------------------------------------------------------
%% 5. 超分辨率图像渲染与成果图导出
% -------------------------------------------------------------
fprintf('--- 正在渲染并导出图像 ---\n');

% 标尺参数 (右下角 1 mm 标尺)
scale_len_mm = 1.0; 
scale_z = z_axis_sr(end) - 0.5; 
scale_x = [x_axis_sr(end) - 0.5 - scale_len_mm, x_axis_sr(end) - 0.5];

% ----------------- 图 1: 血管灌注密度图 (Density Map) -----------------
f1 = figure('Color', 'k', 'Position', [100, 100, 850, 650]);
im_density = MatOut .^ (1/3); % 幂律压缩压制大血管高光，拉伸毛细血管对比度
imagesc(x_axis_sr, z_axis_sr, im_density); 
axis image; colormap(gca, hot(256));
clim([0, quantile(im_density(:), 1)]);

hold on; plot(scale_x, [scale_z, scale_z], 'w-', 'LineWidth', 3);
text(scale_x(1), scale_z - 0.2, sprintf('%.1f mm', scale_len_mm), 'Color', 'w', 'FontSize', 10);

title('ULM Vascular Density Map (Power Law compressed)', 'Color', 'w', 'FontSize', 12);
xlabel('Lateral Distance (mm)', 'Color', 'w'); 
ylabel('Depth (mm)', 'Color', 'w');
set(gca, 'XColor', 'w', 'YColor', 'w', 'Color', 'k');
exportgraphics(f1, fullfile(save_dir, 'Density_Map.png'), 'Resolution', 300);

% ----------------- 图 2: 轴向血流方向图 (Axial Flow Map) -----------------
f2 = figure('Color', 'k', 'Position', [150, 150, 850, 650]);
dir_color = cat(1,flip(flip(hot(128),1),2),hot(128));

% 空间平滑以压制微泡方向符号抖动
im_dir = (MatOut .^ (1/4)) .* sign(imgaussfilt(MatOut_zdir, 0.8));
imagesc(x_axis_sr, z_axis_sr, im_dir); 
axis image; colormap(gca, dir_color);
c_limit = quantile(abs(im_dir(:)), 1);
if c_limit == 0, c_limit = 1; end
clim([-c_limit, c_limit]);

hold on; plot(scale_x, [scale_z, scale_z], 'w-', 'LineWidth', 3);
text(scale_x(1), scale_z - 0.2, sprintf('%.1f mm', scale_len_mm), 'Color', 'w', 'FontSize', 10);

title('ULM Axial Flow Direction (Red/Yellow: Downward, Black/Red: Upward)', 'Color', 'w', 'FontSize', 12);
xlabel('Lateral Distance (mm)', 'Color', 'w'); 
ylabel('Depth (mm)', 'Color', 'w');
set(gca, 'XColor', 'w', 'YColor', 'w', 'Color', 'k');
exportgraphics(f2, fullfile(save_dir, 'Axial_Direction_Map.png'), 'Resolution', 300);

% ----------------- 图 3: 绝对流速分布图 (Velocity Map) -----------------
f3 = figure('Color', 'k', 'Position', [200, 200, 850, 650]);
if any(MatOut_vel(:) > 0)
    v_max = quantile(MatOut_vel(MatOut_vel > 0), 0.98);
else
    v_max = 1;
end
if v_max <= 0, v_max = 1; end

imagesc(x_axis_sr, z_axis_sr, MatOut_vel); 
axis image; 

% 背景无微泡区域设为纯黑
cmap_turbo = turbo(256);
cmap_turbo(1, :) = [0, 0, 0];
colormap(gca, cmap_turbo);
clim([0, v_max]);

cb = colorbar; 
cb.Color = 'w'; 
cb.Label.String = 'Mean Velocity (mm/s)'; 
cb.Label.Color = 'w';

hold on; plot(scale_x, [scale_z, scale_z], 'w-', 'LineWidth', 3);
text(scale_x(1), scale_z - 0.2, sprintf('%.1f mm', scale_len_mm), 'Color', 'w', 'FontSize', 10);

title(sprintf('ULM Absolute Velocity Magnitude (0 ~ %.1f mm/s)', v_max), 'Color', 'w', 'FontSize', 12);
xlabel('Lateral Distance (mm)', 'Color', 'w'); 
ylabel('Depth (mm)', 'Color', 'w');
set(gca, 'XColor', 'w', 'YColor', 'w', 'Color', 'k');
exportgraphics(f3, fullfile(save_dir, 'Velocity_Map.png'), 'Resolution', 300);

fprintf('所有流程执行完毕，结果已保存至: %s\n', save_dir);