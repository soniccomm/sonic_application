% 此脚本为超分辨微血管成像第二步（第一步为采集数据，第三步为ULM_process）
% 脚本功能：
% 适用于凸阵探头
% 读取软件采集的ADC数据，做波束合成获取IQ解调数据

% 用法:
% 将data_filepath的路径改为保存数据的（以时间命名的）完整路径
% 将bfdata_filepath的路径改为波束合成结果要保存到的路径
% 程序运行时，将从data_filepath中读取保存的数据，分批进行波束合成以及必要的后处理，分批保存到bfdata_filepath中。
% 总共要进行波束合成的数据大小和每一批进行波束合成的数据大小由以下两个参数控制：
% frame_nums_want和frame_per_bf_file
% frame_nums_want为波束合成的帧数
% frame_per_bf_file为每一批波束合成的帧数
% 程序运行时，每合成好指定的帧数后会将数据保存到bfdata_filepath中，如bfiq_1.mat,bfiq_2.mat,...
% 程序运行时会出现进度提示，"bfiq_m"提示当前处于第m批波束合成， “beamforming frame #n”提示当前处于第m批波束合成的第n帧。

% 波束合成的参数在成像参数中控制

clear all
clc
close all
%% 加载当前环境变量
currentPath = pwd;
parentDir = fileparts(fileparts(currentPath));
addpath(genpath(parentDir));

data_filepath = 'D:\software_matlab\exampledata\convex\20250804173238';   %用于解析的ADC数据文件路径
adc_folder = data_filepath;
bfdata_filepath = 'D:\software_matlab\exampledata\convex\20250804173238\bfiq';   %用于保存波束合成结果的路径
bf_folder = bfdata_filepath;
if ~exist(bfdata_filepath,'dir')
    mkdir(bfdata_filepath);
end

disp("波束合成结果将保存到"+bfdata_filepath)

is_show = 1;

%% 成像参数
frame_nums_want = 1; % 波束合成帧数
frame_per_bf_file = 1; % 保存的一包波束合成数据的帧数

Nx = 512;               % 波束合成线束数量x轴方向
x_start = -0.06;     % 波束合成区域左端x坐标(m)
x_end = 0.06;        % 波束合成区域右端x坐标(m)

Nz = 1024;               % 波束合成线束数量z轴方向
z_start = -0.0074;        % 波束合成区域上端z坐标(m)
z_end = 0.1;           % 波束合成区域下端z坐标(m)

wavetype ='plane_wave';
chmap =[82 ,114 ,86 ,118 ,90 ,122 ,94 ,126 ,80 ,112 ,84 ,116 ,88 ,120 ,92 ,124 ,83 ,115 ,87 ,119 ,91 ,123 ,95 ,127 ,81 ,113 ,85 ,117 ,89 ,121 ,93 ,125 ,34 ,98 ,38 ,102 ,42 ,106 ,46 ,110 ,32 ,96 ,36 ,100 ,40 ,104 ,44 ,108 ,35 ,99 ,39 ,103 ,43 ,107 ,47 ,111 ,33 ,97 ,37 ,101 ,41 ,105 ,45 ,109 ,18 ,66 ,22 ,70 ,26 ,74 ,30 ,78 ,16 ,64 ,20 ,68 ,24 ,72 ,28 ,76 ,19 ,67 ,23 ,71 ,27 ,75 ,31 ,79 ,17 ,65 ,21 ,69 ,25 ,73 ,29 ,77 ,2 ,50 ,6 ,54 ,10 ,58 ,14 ,62 ,0 ,48 ,4 ,52 ,8 ,56 ,12 ,60 ,3 ,51 ,7 ,55 ,11 ,59 ,15 ,63 ,1 ,49 ,5 ,53 ,9 ,57 ,13 ,61 ];
probe_name = 'C2-5';  %探头名称
sos = 1540; % 声速
channel = 128; % 通道数
bandwidth = 70; % 解调时低通的带宽控制



% frame_nums 帧数
[fs,fc,sampleNum,imagedepth,focus_depth,cstartoffset,frame_nums,NumsPerFile,prf,steer] = read_adc_para_pw(strcat(data_filepath,'\Param.txt'),wavetype);
% [fs,sampleNum,imagedepth,focus_depth,cstartoffset,frame_nums,NumsPerFile,steer] = read_adc_para(strcat(data_filepath,'\Param.txt'),wavetype);
% fc = 3.2e6;


frame_pack_num_1 = floor(frame_nums_want/frame_per_bf_file);
frame_pack_num_2 = floor(frame_nums/frame_per_bf_file);
if (frame_pack_num_1>frame_pack_num_2)
    frame_pack_num = frame_pack_num_2;
    frame_nums = frame_pack_num*frame_per_bf_file;
    disp("要做波束合成的帧数超过采集帧数")
    disp("波束合成的帧数为"+frame_nums)
    disp("保存波束合成数据每个文件包含"+frame_per_bf_file+"帧")
    disp("总共"+frame_pack_num+"个")
else
    frame_pack_num = frame_pack_num_1;
    frame_nums = frame_pack_num*frame_per_bf_file;
    disp("波束合成的帧数为"+frame_nums)
    disp("保存波束合成数据每个文件包含"+frame_per_bf_file+"帧")
    disp("总共"+frame_pack_num+"个")
end

% 配置参数
probe = Probe_para(probe_name);
AcqConfig.Probe = probe;
AcqConfig.Tx.channel = channel;
AcqConfig.Tx.fs = fs;
AcqConfig.Tx.sos = sos;
AcqConfig.Tx.focus_depth = focus_depth;
AcqConfig.Rx.fs = fs;
AcqConfig.Rx.sos = sos;
AcqConfig.Rx.sample_num = sampleNum;
lambda = AcqConfig.Rx.sos/fc;


steer = deg2rad(steer);
Filenums =  frame_nums*length(steer);



%% 接收波束合成

z_start = z_start + AcqConfig.Probe.R;
z_end = z_end + AcqConfig.Probe.R;

x_axis = linspace(x_start,x_end,Nx);
z_axis = linspace(z_start,z_end,Nz);

AcqConfig.Rx.rev_line.x = x_axis;
AcqConfig.Rx.rev_line.z = z_axis;

[x_grid, z_grid] = meshgrid(x_axis,z_axis);
rf = zeros(AcqConfig.Rx.sample_num ,Nx);
x_grid = reshape(x_grid,[Nz*Nx,1]);
z_grid = reshape(z_grid,[Nz*Nx,1]);

% 计算发射和接收延时
xm = bsxfun(@minus, AcqConfig.Probe.element_pos.x,x_grid);
zm = bsxfun(@minus,AcqConfig.Probe.element_pos.z,z_grid);
rx_delay = sqrt(xm.^2+zm.^2)/sos*fs ;    

tx_delay = zeros(Nx*Nz,numel(steer));
for i = 1:numel(steer)
    angle_thres = 25/180*pi;
    active_ele = abs(probe.element_pos.theta - steer(i))<angle_thres;
    thetas = steer(i) - probe.element_pos.theta;
    dist_all = cos(thetas)*probe.R;
    dist_comp = min(dist_all(active_ele));
    
    tx_delay(:,i) = (x_grid*sin(steer(i)) + z_grid*cos(steer(i)) - dist_comp)/sos*fs;
end


% 夹角
xz_grid_rad = atan(x_grid./z_grid);



% 高通滤波，去除原始信号中的直流
fstop = 0.1e6;
[b1,a1] = butter(5,fstop/(fs)*2,"high");
% 低通滤波，去除载波信号
Wn = (fc*bandwidth/100)/(AcqConfig.Tx.fs/2);
wn = 0.5;
[b2,a2] = butter(5,Wn);
% [b2,a2] = fir1(63,0.5,'low',hamming(64));


rx_delay = gpuArray(single(rx_delay));
tx_delay = gpuArray(single(tx_delay));


% 根据阵元指向性开窗
f_number = est_fNumber(AcqConfig.Probe.element_pitch,lambda,0.8);
f_mask = zeros(Nx*Nz,AcqConfig.Probe.element_num,'single');
for i = 1:AcqConfig.Probe.element_num
    f_mask(:,i) = z_grid./abs(x_grid -  AcqConfig.Probe.element_pos.x(i))/2 > f_number;
end
rx_apod = f_mask;
rx_apod = gpuArray(rx_apod);

if is_show
    hFig = figure('Name',"planewave reconstruction",'NumberTitle','off');
    hIm = imagesc(x_axis*100,(z_axis-AcqConfig.Probe.R)*100,zeros(Nz,Nx));
    colormap(gray);title("平面波成像");clim([-60, 0]);xlabel("cm");ylabel("cm")
    axis equal;axis tight
end

% dsc
if(AcqConfig.Probe.type == "convex")
    valid_rad = atan(x_grid./(z_grid));
    dsc_table = (valid_rad>AcqConfig.Probe.element_pos.theta(1)) .* (valid_rad<AcqConfig.Probe.element_pos.theta(end));

    valid_range = (sqrt(x_grid.^2 + (z_grid).^2)>=(AcqConfig.Probe.R+0.001));
    dsc_table = dsc_table.*valid_range;

    dsc_table = reshape(dsc_table, Nz, Nx);
end

% weight
for i = 1:numel(steer)
    a = 10/180*3.1415926;
    b = 13/180*3.1415926;
    pixel_rad = atan(x_grid./z_grid);
    theta_rad = abs(pixel_rad - steer(i));

    dsc_weight = (theta_rad<=a).*1 + (theta_rad>a & theta_rad<b).*(1 - (theta_rad - a)/(b - a));
    dsc_weight(theta_rad>b) = 0;

    dsc_weight_table(:,i) = dsc_weight;
end
dsc_weight_table = gpuArray(dsc_weight_table);



for pack_num_idx = 1:frame_pack_num
    
disp("bfiq_"+pack_num_idx)

s_idx = frame_per_bf_file*(pack_num_idx-1)*length(steer)+1;
e_idx = frame_per_bf_file*(pack_num_idx)*length(steer);
    
[temp,adc_scan_length] = ReadADC_period (adc_folder,channel ,NumsPerFile,Filenums, s_idx, e_idx);
temp = reshape(temp,channel ,adc_scan_length(1),length(steer),frame_per_bf_file);
temp = permute(temp,[2,1,3,4]);
adc_data = temp(:,chmap+1,:,:);

t = (0:size(adc_data,1)-1)/AcqConfig.Tx.fs;
img_c = zeros(Nz*Nx,frame_per_bf_file,'single');

img_c = gpuArray(img_c);

for i = 1:frame_per_bf_file
    disp("beamforming frame #"+num2str(i));
    for j = 1:numel(steer)
        idx = (i-1)*numel(steer)+j;
        rf_data = adc_data(:,:,j,i);
        rf_data = gpuArray(double(rf_data));
        rf_data = filtfilt(b1,a1,rf_data);
        % demodulation
        IQdata = double(rf_data).*exp(-1i*2*pi*fc*t');
        IQdata = filtfilt(b2,a2,IQdata)*2;
        % beamforming and compounding
        for k = 1:AcqConfig.Probe.element_num
            temp_delay = (tx_delay(:,j) + rx_delay(:,k))/AcqConfig.Tx.fs;
            if k == 1
                bfiq = rx_apod(:,k).*interp1(t,IQdata(:,k),temp_delay,"linear",0).*exp(2*1j*pi*fc*temp_delay);  % 相位补偿
            else
                bfiq = bfiq + rx_apod(:,k).*interp1(t,IQdata(:,k),temp_delay,"linear",0).*exp(2*1j*pi*fc*temp_delay);  % 相位补偿
            end
        end
        img_c(:,i) = img_c(:,i) + bfiq.*dsc_weight_table(:,j);
    end
end
img_c = reshape(img_c,Nz,Nx,[])/numel(steer);
img_c = gather(img_c);

bfiq = img_c;
bfiq = bfiq.*dsc_table;

if is_show
    hIm.CData = log_compressed(abs(bfiq(:,:,1)));
    pause(0.0001)
end

save(fullfile(bf_folder,"bfiq_"+pack_num_idx+".mat"),"bfiq")


end






















