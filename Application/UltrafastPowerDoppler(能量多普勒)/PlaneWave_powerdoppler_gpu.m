% 此为能量多普勒（功率多普勒）第二步，第一步为采集数据

% 数据处理流程：
% 带解调的波束合成 -> svd滤波 -> 计算能量

% 可保存波束合成后IQ数据和能量图

clear all
clc
close all

currentPath = pwd;
parentDir = fileparts(fileparts(currentPath));
addpath(genpath(parentDir));

%% 必要参数
filefolder = 'D:\software_matlab\exampledata\colordoppler\20250809H\20250809144618';  % 数据所在文件夹
% filefolder = 'D:\software_matlab\exampledata\colordoppler\20250809Z\20250809144808';

[fs,fc,BF_SampleN,ImageDepth,~,cstartoffset,frame_nums,numperfile,prf,steering_deg] = read_adc_para_pw(strcat(filefolder,'\Param.txt'),'plane_wave');

                            % fs               采样率
                            % BF_SampleN       采样点数
                            % numperfile       保存的数据每一包所含帧数
                            % steering_deg     平面波发射角度
                            % ImageDepth       图像深度
                            % prf              采集数据所用的脉冲重复频率 
                            % fc               发射频率

load_file_num = 10;         % 读取并处理的文件数量

Imagestart = 0.002;         % B模式图像起始深度
ImageDepth = 0.040;         % B模式图像深度
BeamN = 256;                % 线束（波束合成x方向线束数量）       

t_axis_span = 2;            % 秒 定义循环刷新的图的时间跨度
fft_period = 50;            % 计算一次血流所用帧数（Window Size）
lag = 5;                    % 相邻两次血流计算间隔的帧数（Update Interval/Hop Size）

x1_loc_real = -0.00;        % 取样框左上角物理坐标x1
z1_loc_real = 0.0056;       % 取样框左上角物理坐标z1
x2_loc_real = 0.0121;       % 取样框右下角物理坐标x2
z2_loc_real = 0.0194;       % 取样框右下角物理坐标z2

svd_auto = 1;               % 是否使用自适应阈值进行svd滤波，若为1，则svd_ord1和svd_ord2不起作用
svd_ord1 = 20;              % SVD滤波起始阶数
svd_ord2 = 45;              % SVD滤波结束阶数

is_save_BF = 1;             % 是否保存波束合成后IQ数据 1为保存 0为不保存
Bmode_save_index = 1:17;       % 每一包数据保存的帧号
                            % 仅在is_save_BF = 1时生效 例：1表示每一包仅保存第1帧，1:numperfile表示保存第1帧到最后一帧
is_save_power = 1;          % 是否保存血流数据 1为保存 0为不保存

drange_B = 60;              % bmode动态范围
drange_power = 30;          % 能量多普勒动态范围

% 默认参数
probe_name = 'L5-10';       % 探头名称
TxChannel = 128;            % 发射通道数量
RxChannel = 128;            % 接收通道数量
sos = 1540;                 % 声速



%% 波束合成参数计算

probe = Probe_para(probe_name);

% 通道对应关系（不用改）
ch_map =[82 ,114 ,86 ,118 ,90 ,122 ,94 ,126 ,80 ,112 ,84 ,116 ,88 ,120 ,92 ,124 ,83 ,115 ,87 ,119 ,91 ,123 ,95 ,127 ,81 ,113 ,85 ,117 ,89 ,121 ,93 ,125 ,34 ,98 ,38 ,102 ,42 ,106 ,46 ,110 ,32 ,96 ,36 ,100 ,40 ,104 ,44 ,108 ,35 ,99 ,39 ,103 ,43 ,107 ,47 ,111 ,33 ,97 ,37 ,101 ,41 ,105 ,45 ,109 ,18 ,66 ,22 ,70 ,26 ,74 ,30 ,78 ,16 ,64 ,20 ,68 ,24 ,72 ,28 ,76 ,19 ,67 ,23 ,71 ,27 ,75 ,31 ,79 ,17 ,65 ,21 ,69 ,25 ,73 ,29 ,77 ,2 ,50 ,6 ,54 ,10 ,58 ,14 ,62 ,0 ,48 ,4 ,52 ,8 ,56 ,12 ,60 ,3 ,51 ,7 ,55 ,11 ,59 ,15 ,63 ,1 ,49 ,5 ,53 ,9 ,57 ,13 ,61 ];
elex = probe.element_pos.x;
elez = probe.element_pos.z;

x_axis = linspace(elex(1),elex(end),BeamN);
step_z = sos/fs/2;
z_axis = (0:step_z:(BF_SampleN-1) * step_z);
[x_grid, z_grid] = meshgrid(x_axis,z_axis);
scan.xx  =  x_grid;
scan.zz  =  z_grid;

allbeamx = reshape(scan.xx,[BF_SampleN*BeamN,1]);
allbeamz = reshape(scan.zz,[BF_SampleN*BeamN,1]);

probe_type_ptr = libpointer('cstring', "linear");

single_Steering = single(steering_deg);
single_Steering_ptr = libpointer('singlePtr', single_Steering);
single_Steering_len = length(single_Steering);

single_allbeamx = single(allbeamx);
single_allbeamx_ptr = libpointer('singlePtr', single_allbeamx);
single_allbeamx_len = length(single_allbeamx);

single_allbeamz = single(allbeamz);
single_allbeamz_ptr = libpointer('singlePtr', single_allbeamz);
single_allbeamz_len = length(single_allbeamz);

single_cstartoffset = single(0);
single_cstartoffset_ptr = libpointer('singlePtr', single_cstartoffset);
single_cstartoffset_len = length(single_cstartoffset);

single_ch_map = single(ch_map);
single_ch_map_ptr = libpointer('singlePtr', single_ch_map);
single_ch_map_len = length(single_ch_map);

elexz = cat(2, elex, elez)';
single_elexz = single(elexz);
single_elexz_ptr = libpointer('singlePtr', single_elexz);
single_elexz_len = length(single_elexz);

lambda = sos/fc;
probe.element_width = probe.element_pitch/2;
f_number = est_fNumber(probe.element_width,lambda,0.71);

f_mask = zeros(BF_SampleN*BeamN,probe.element_num,'single');
for i = 1:probe.element_num
    f_mask(:,i) = allbeamz./abs(allbeamx -  probe.element_pos.x(i))/2 > f_number;
end
rx_apod = f_mask;
xm = bsxfun(@minus, probe.element_pos.x,allbeamx);
zm = bsxfun(@minus,probe.element_pos.z,allbeamz);
rx_delay = sqrt(xm.^2+zm.^2)/sos*fs;
SteeringNum = numel(steering_deg);
tx_delay = zeros(BF_SampleN*BeamN,SteeringNum);
steer = deg2rad(steering_deg);
for i = 1:SteeringNum
    if steer(i) >= 0
        tx_delay(:,i) = ((allbeamx - min(probe.element_pos.x))*sin(steer(i)) + allbeamz*cos(steer(i)))/sos*fs;
    else
        tx_delay(:,i) = ((max(probe.element_pos.x) - allbeamx)*sin(-steer(i)) + allbeamz*cos(steer(i)))/sos*fs;
    end
end

single_rx_apod = reshape(rx_apod,[],1);
single_rx_apod_ptr = libpointer('singlePtr', single_rx_apod);
single_rx_apod_len = length(single_rx_apod);

single_rx_delay = reshape(rx_delay,[],1);
single_rx_delay_ptr = libpointer('singlePtr', single_rx_delay);
single_rx_delay_len = length(single_rx_delay);

single_tx_delay = reshape(tx_delay,[],1);
single_tx_delay_ptr = libpointer('singlePtr', single_tx_delay);
single_tx_delay_len = length(single_tx_delay);


%% 参数计算

% 解调滤波器
bandwidth = 80; 
Wn = (fc*bandwidth/100)/(fs/2);
% 确定滤波器长度 (经验公式)
M = ceil(6.64 * fs/2 / (fc*bandwidth/100));   % Hamming 窗专用公式
% 确保奇数长度以获得线性相位
if mod(M, 2) == 0
    M = M + 1;
end
% 使用 Hamming 窗设计低通 FIR 滤波器
b_fir = fir1(M - 1, Wn, 'low', hamming(M));
%给gpu传
single_filter = single(b_fir);
single_filter_ptr = libpointer('singlePtr', single_filter);
single_filter_len = length(single_filter);

% 
buffer_num = floor(fft_period/(numperfile/SteeringNum))+2; %缓冲区个数
t_idx = fft_period; % 这个idx是当前在图的idx 初始为fft_period 
t_buffer_idx = fft_period; % 这个idx是在缓冲区的idx 初始为fft_period 


%% 用于求血流的像素的坐标


x_dif = abs(x_axis - x1_loc_real);
z_dif = abs(z_axis - z1_loc_real);
[~, x_left] = min(x_dif);
[~, z_up] = min(z_dif);

x_dif = abs(x_axis - x2_loc_real);
z_dif = abs(z_axis - z2_loc_real);
[~, x_right] = min(x_dif);
[~, z_down] = min(z_dif);

rec_x_num = x_right - x_left + 1;
rec_z_num = z_down - z_up + 1;

x_loc_all = [];
z_loc_all = [];
for i = x_left:x_right
    for j = z_up:z_down
        z_loc_all = [z_loc_all,j];
        x_loc_all = [x_loc_all,i];
    end
end

loc_num = numel(x_loc_all);

single_x_loc_all = single(x_loc_all);
single_x_loc_all_ptr = libpointer('singlePtr', single_x_loc_all);
single_z_loc_all = single(z_loc_all);
single_z_loc_all_ptr = libpointer('singlePtr', single_z_loc_all);
single_points_len = loc_num;

%% 加载

%加载dll
if ~libisloaded('US_APP')
    loadlibrary('US_APP.dll', 'ApplicationMatlabInterface.h');
end

bprocesspara.Demod_AFE_Dynamic = -9;
gpu_handle = calllib('US_APP', 'initializepowerGPU', ...
    prf, ...
    numperfile, ...
    buffer_num, ...
    fft_period, ...
    lag, ...
    t_axis_span, ...
    t_idx, ...
    t_buffer_idx, ...
    probe_type_ptr, ...
    bprocesspara.Demod_AFE_Dynamic, ...
    1, ... %AcqConfig.Tx.FsNum
    TxChannel, ... %AcqConfig.Tx.Channel
    RxChannel, ... %RxChannel
    128, ... %AcqConfig.Probe.element_num
    probe.element_pitch, ... %AcqConfig.Probe.element_pitch
    32, ...
    BF_SampleN, ...
    BeamN, ...
    sos, ...
    fs, ...
    fc, ...
    svd_auto,... svd
    svd_ord1,... svd 
    svd_ord2,... svd
    single_Steering_ptr,single_Steering_len,...
    single_allbeamx_ptr,single_allbeamx_len,...
    single_allbeamz_ptr,single_allbeamz_len,...
    single_cstartoffset_ptr,single_cstartoffset_len,...
    single_ch_map_ptr,single_ch_map_len,...
    single_elexz_ptr,single_elexz_len,...
    single_filter_ptr,single_filter_len,...
    single_rx_apod_ptr,single_rx_apod_len,...
    single_rx_delay_ptr,single_rx_delay_len,...
    single_tx_delay_ptr,single_tx_delay_len, ...
    single_x_loc_all_ptr,single_z_loc_all_ptr,single_points_len, ...
    rec_x_num, rec_z_num);


%% 获取所有有效数据文件
files =dir (fullfile (filefolder ,'**' ,'*bin' ));
num_files =length (files);
file_numbers =zeros (num_files, 1);

for i =1 :num_files
    [~,file_name ]=fileparts (files(i).name );
    try
        file_numbers(i)=str2double(file_name);
    catch
        file_numbers(i)=inf ;
    end
end

valid_indices =~isinf(file_numbers);
valid_files =files (valid_indices);
valid_numbers =file_numbers (valid_indices);
[~,sort_indices ]=sort (valid_numbers);
sorted_files =valid_files (sort_indices);


%% 定义图像

% 要显示的B模式图像（波束合成图）
bfdata = zeros(1, 2* BF_SampleN * BeamN * numperfile / SteeringNum);
bfdata = single(bfdata);
bfdata_ptr = libpointer('singlePtr', bfdata);

valid_indices = find(z_axis >= Imagestart & z_axis <= ImageDepth);
zz_cut = z_axis(valid_indices);
Nz_cut = numel(valid_indices);

% 获取屏幕尺寸
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);
            
real_width = x_grid(end)- x_grid(1);
real_height = z_grid(end) - z_grid(1);

% 规定窗口占据屏幕比例
fig_width_ratio = 0.3;
fig_height_ratio = real_height/real_width*fig_width_ratio*1.5;
% 设置窗口大小
fig_width = screen_width * fig_width_ratio;
fig_height = screen_height * fig_height_ratio;
% 设置窗口位置
fig_left = screen_width * (0.5-fig_width_ratio/2);
fig_bottom = screen_height * (1-fig_height_ratio)/2;

hFig1 = figure('Name',"Bmode",'Position', [fig_left, fig_bottom, fig_width, fig_height]);
hIm1 = imagesc(x_axis, zz_cut, zeros(Nz_cut,BeamN)-100,[-60 0]);
colormap(gray);title("能量多普勒");
axis equal;axis tight

power_min = -drange_power;  % 能量最低 dB
power_max = 0;   % 能量最高 dB
colormap(hot); 
cb = colorbar; % 添加 colorbar
cb.Label.String = 'power (dB)'; % 设置 colorbar 标签
% 设置 colorbar 的刻度范围，与血流数据范围对应
caxis([power_min, power_max]);

hold on
rectangle('Position', [x_axis(x_left), z_axis(z_up), x_axis(x_right) - x_axis(x_left), z_axis(z_down) - z_axis(z_up)], 'EdgeColor', 'r', 'LineWidth', 2);

% 要显示的blood
power_matrix = single(zeros(rec_z_num, rec_x_num));
power_matrix_ptr = libpointer('singlePtr', power_matrix);

% 如果要保存则创建文件夹
if is_save_BF
    if ~exist(fullfile(filefolder, 'bfdata'), 'dir')
        mkdir(fullfile(filefolder, 'bfdata'));
    end
end
if is_save_power
    if ~exist(fullfile(filefolder, 'power'), 'dir')
        mkdir(fullfile(filefolder, 'power'));
    end
end


%% 读取数据并处理
for file_i = 1:load_file_num
    disp(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name))
    fileID = fopen(fullfile(sorted_files(file_i).folder, sorted_files(file_i).name), 'rb'); 

    fseek(fileID,0,1);
    nFileLen = ftell(fileID);
    fseek(fileID,0,-1);
    %注意这里类型控制
    alldata = fread(fileID, nFileLen,'int8=>int8');
    
    sid = 0;
    alldata_len = SteeringNum*(2*BF_SampleN*TxChannel+32);
    cur_data = alldata(1:(numperfile/SteeringNum)*alldata_len);
    
    alldata_ptr = libpointer('int8Ptr', cur_data);
    bag_idx = file_i-1;  %这个bag_idx要从0开始
    ret = calllib('US_APP', 'processPowerDataBeamformingandPostGPU', gpu_handle, alldata_ptr, ...
        alldata_len, bfdata_ptr, power_matrix_ptr, bag_idx); 
    
    rev = bfdata_ptr.Value;
    bfdata_real = reshape(rev(1:2:end),BF_SampleN,BeamN,[]);
    bfdata_imag = reshape(rev(2:2:end),BF_SampleN,BeamN,[]);
    bfdata_iq = bfdata_real + 1i* bfdata_imag;
    B_image = log_compressed(abs(bfdata_iq(valid_indices,:,end)));
    B_image(B_image<-drange_B) = -drange_B;
    B_image(B_image>0) = 0;
    B_image = B_image + drange_B;
    B_image = B_image/drange_B; % 现在 B_image 的范围是 [0, 1]
    B_image_rgb = repmat(B_image, [1, 1, 3]); % 复制三份，创建 [rows, cols, 3] 的 RGB 矩阵

    if is_save_BF
        bfdata_save = bfdata_iq(:,:,Bmode_save_index);
        save(fullfile(filefolder,'bfdata',num2str(bag_idx)+".mat"), 'bfdata_save',"x_axis","z_axis");
    end

    % 获取power套回原图位置 
    power_matrix_update = reshape(power_matrix_ptr.Value, rec_z_num, rec_x_num);
    power_matrix_update = log_compressed(power_matrix_update);
    power_full = zeros(BF_SampleN,BeamN);
    power_full(z_up:z_down,x_left:x_right) = power_matrix_update;
    power_full = power_full(valid_indices,:);

    % 获取区域mask 套回原图位置
    flow_mask_full = zeros(BF_SampleN,BeamN);
    flow_mask_full(z_up:z_down,x_left:x_right) = 1;
    flow_mask_full = logical(flow_mask_full); 
    flow_mask_full = flow_mask_full(valid_indices,:);
    % 将掩码扩展到3个颜色通道
    mask_3D = repmat(flow_mask_full, [1, 1, 3]);

    % 将数据归一化到 [1, 256] 的整数索引
    power_indices = round( (power_full - power_min) / (power_max - power_min) * 255 + 1 );
    % 处理超出范围的值
    power_indices(power_indices < 1) = 1;
    power_indices(power_indices > 256) = 256;
    % 使用 ind2rgb 将索引矩阵和 colormap 转换为 RGB 图像
    power_image_rgb = ind2rgb(power_indices, hot(256));

    % 在掩码区域，用血流图数据替换掉原来的灰度数据
    B_image_rgb(mask_3D) = power_image_rgb(mask_3D);

    hIm1.CData = B_image_rgb;
    pause(0.001)

    if is_save_power
        save(fullfile(filefolder,'power',num2str(bag_idx)+".mat"), 'power_matrix_update');
    end

    fclose(fileID);
end


%% 卸载

%释放显存和内存
calllib('US_APP', 'deleteBeamformingGPUHandle', gpu_handle);

%卸载dll
unloadlibrary('US_APP');







