% Directly invoke PALA toolbox
clear;clc;close all;
addpath(genpath('.\PALA_SRUS'));
currentPath = pwd;
parentDir = fileparts(fileparts(currentPath));
addpath(genpath(parentDir));

%% load data
% if micro bubbles are too much, choose ROI to get better tracking
ROI_z = 1:400;
ROI_x = 1:256;

dataFolder = 'D:\software_matlab\exampledata\S256-ULM\bfiq';
resultFolder = fullfile(dataFolder,'results');

if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

IQdataList = {};
IQnum = 1;
fileList = dir([dataFolder,filesep,'*.mat']);
for i = 1:numel(fileList)
    IQdataList{IQnum} = [fileList(i).folder,filesep,fileList(i).name];
    IQnum = IQnum + 1;
end
IQnum = IQnum - 1;



%% Adapt parameters
load(IQdataList{1});
bfiq = IQ;
PData.PDelta = [1 1 1];
PData.Origin = [0 size(bfiq,2)/2 0];
PData.Size = [size(bfiq(ROI_z,ROI_x,:),1) size(bfiq(ROI_z,ROI_x,:),2) 1];
% PData.Size = [244 244 1];
framerate = 1200;  % compounded frame rate
NbFrames = size(bfiq,3);   % frame num each IQdata contains
lambda = 0.1;
%% ULM parameters
res = 10;
ULM = struct( ...                       
    'numberOfParticles', 100,...         % Number of particles per frame. (30-100)
    'res',10,...                        % Resolution factor. Typically 10 for images at lambda/10.
    'SVD_cutoff',[50 NbFrames - 50],...  % svd filtering
    'max_linking_distance',2,...        % Maximum linking distance between two frames to reject pairing, in pixels units (UF.scale(1)). (2-4 pixel).
    'min_length', 15,...                 % Minimum length of the tracks. (5-20)
    'fwhm',[1 1]*3,...                  % Size of the mask for localization. (3x3 for pixel at lambda, 5x5 at lambda/2). [fmwhz fmwhx]
    'max_gap_closing', 0,...            % Allowed gap in microbubbles pairing. (0)
    'size',[size(bfiq(ROI_z,ROI_x,:),1),size(bfiq(ROI_z,ROI_x,:),2),NbFrames],...
    'scale',[1 1 1/framerate],...       % Scale [z x t]
    'numberOfFramesProcessed',NbFrames,... % Number of processed frames
    'interp_factor',1/res...            % interpfactor
    );
ULM.butter.CuttofFreq = [50 300];       % Cut off frequency (Hz) for additional filter. Typically [20 300] at 1kHz.
ULM.butter.samplingFreq = framerate;    % Sampling frequency (Hz)
[but_b,but_a] = butter(2,ULM.butter.CuttofFreq/(ULM.butter.samplingFreq/2),'bandpass');
ULM.parameters.NLocalMax = 3;           % Safeguard on the number of maxLocal in the fwhm*fwhm grid (3 for fwhm=3, 7 for fwhm=5)
res = ULM.res;
% listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
% listAlgo = {'wa','gaussian_fit','radial'};
listAlgo = {'gaussian_fit'};
Nalgo = numel(listAlgo);

%% select SVD filtering Noise
% load one of the bfiq to test SVD thresholds
bag_idx = 1;
load(IQdataList{1});
bfiq = IQ(ROI_z,ROI_x,:);
bulles = SVDfilter(bfiq,ULM.SVD_cutoff);
bulles = filter(but_b,but_a,bulles,[],3);
bulles(~isfinite(bulles))=0;

figure("Name","Filter visualization");
for i = 1:size(bfiq,3)
    img_envelope = abs(bulles(:,:,i));
    img_log = log_compressed(img_envelope);
    imagesc(x_axis(ROI_x),z_axis(ROI_z),img_log,[-30 0]);
    colormap("gray");axis equal;axis tight
    title("Bag "+bag_idx+" Frame "+i)
    pause(0.01);
end
%% Load and localize data     
fprintf('--- ULM PROCESSING --- \n\n')
clear Track_tot Track_tot_interp ProcessingTime bulles IQ dB
t1=tic;
for iqIdx = 1:IQnum
    fprintf('Processing bloc %d/%d\n',iqIdx,IQnum);
    tmp = load(IQdataList{iqIdx});
    bfiq = tmp.IQ;
    bfiq = bfiq(ROI_z,ROI_x,:);
%     bfiq = imresize(bfiq,[244,244]);
    IQ_filt = SVDfilter(bfiq,ULM.SVD_cutoff);
    IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0;
    [~,~] = PALA_multiULM(IQ_filt,listAlgo,ULM,PData,'savingfilename',[resultFolder filesep 'Tracks' num2str(iqIdx,'%.3d') '.mat']);
%     [Track_raw,Track_interp,ProTime] = PALA_multiULM(IQ_filt,listAlgo,ULM,PData);
%     save([dataFolder filesep 'Tracks' num2str(iqIdx,'%.3d') '.mat'],'Track_raw','Track_interp','ProTime','ULM','UF','PData','-v6')
end
t2=toc(t1);
fprintf('ULM done in %d hours %.1f minutes. (for all localization algorithms) \n', floor(t2/60/60), rem(t2/60,60));

%% Create MatOuts     
% for each algorithm, create the MatOut density with interpolated tracks for visual analysis, and with non interpolated tracks for aliasing index calculation.
fprintf('--- CREATING MATOUTS --- \n\n')
MatOutSat = [];
NbrOfLoc = zeros(Nalgo,1);
MatOut = cell(Nalgo,1);
MatOut(:)={0};
MatOutNoInterp = MatOut;
MatOut_vel = MatOut;

ULM.SRscale = ULM.scale(1)/ULM.res;
ULM.SRsize = round(ULM.size(1:2).*[1 2].*ULM.scale(1:2)/ULM.SRscale);
ULM.lambda = lambda;

for iqIdx = 1:IQnum % Generate MatOut density matrix
    disp("load "+[resultFolder filesep 'Tracks' num2str(iqIdx,'%.3d') '.mat'])
    load([resultFolder filesep 'Tracks' num2str(iqIdx,'%.3d') '.mat'],'Track_raw','Track_interp','ProTime')
    aa = -PData(1).Origin([3 1])+[1 1]*1;  % get origin
    bb = [1./PData(1).PDelta([3 1])*ULM.res];  % fix the size of pixel
    aa(3) = 0;bb(3) = 1; % for velocity
    for ialgo = 1:Nalgo
        % MatOut and MatOutVel rendering with interpolated tracks
        Track_matout = Track_interp{ialgo};
        Track_matout = cellfun(@(x) (x(:,[1 2 3])+aa).*bb,Track_matout,'UniformOutput',0);
        [MatOut_i,MatOut_vel_i] = ULM_Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1,'mode','2D_velmean'); %pos in superpix [z x]
        clear Track_matout
        MatOut_vel{ialgo} = MatOut_vel{ialgo}.*MatOut{ialgo}+MatOut_vel_i.*MatOut_i; % weighted summation
        MatOut{ialgo} = MatOut{ialgo}+MatOut_i;
        MatOut_vel{ialgo}(MatOut{ialgo}>0) = MatOut_vel{ialgo}(MatOut{ialgo}>0)./MatOut{ialgo}(MatOut{ialgo}>0); % average velocity
        MatOutSat(iqIdx,ialgo) = nnz(MatOut{ialgo}>0); % compute saturation curve

        % MatOut without interpolation, for gridding index
        Track_matout = Track_raw{ialgo};
        Track_matout = cellfun(@(x) (x(:,[1 2])+aa(1:2)).*bb(1:2),Track_matout,'UniformOutput',0);
        MatOut_i = ULM_Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1); %pos in superpix [z x]
        MatOutNoInterp{ialgo} = MatOutNoInterp{ialgo}+MatOut_i;
        Track_count = cat(1,Track_matout{:});clear Track_matout
        NbrOfLoc(ialgo) = NbrOfLoc(ialgo)+size(Track_count,1);
    end
end
% clear Track_raw Track_interp Track_count Track_matout MatOut_i
save([resultFolder filesep 'MatOut_multi'],'MatOut','MatOut_vel','MatOutSat','NbrOfLoc','ULM','listAlgo','Nalgo','PData');
save([resultFolder filesep 'MatOut_multi_nointerp'],'MatOutNoInterp','MatOutSat','NbrOfLoc','ULM','listAlgo','PData','Nalgo');

%% Display MatOut Intensity

load([resultFolder filesep 'MatOut_multi'],'MatOut','MatOut_vel','MatOutSat','NbrOfLoc','ULM','listAlgo','Nalgo','PData');
load([resultFolder filesep 'MatOut_multi_nointerp'],'MatOutNoInterp','MatOutSat','NbrOfLoc','ULM','listAlgo','PData','Nalgo');

figure;clf,set(gcf,'Position',[652 393 941 585]);
% imagesc(log_compressed(MatOut{1}),[-72 0]); colormap hot% ialgo-th ULM result
wv = 1540/15*1e-3*1e-3;
imagesc(x_axis(ROI_x),z_axis(ROI_z),MatOut_vel{1}.*wv);colormap jet  % ialgo-th ULM result   
title("mean velocity[mm/s]")
axis equal;axis tight;colorbar
%% MatOut intensity rendering, with compression factor
algo_idx = 1;

fprintf('--- GENERATING IMAGE RENDERINGS --- \n\n')
figure;clf,set(gcf,'Position',[652 393 941 585]);
IntPower = 1/2;SigmaGauss=0;
im=imagesc(x_axis(ROI_x)*1000,z_axis(ROI_z)*1000,MatOut{algo_idx}.^IntPower);axis image
if SigmaGauss>0,im.CData = imgaussfilt(im.CData,SigmaGauss);end

title('ULM intensity display')
colormap(gca,hot(128))
clbar = colorbar;caxis(caxis*.8)  % add saturation in image
clbar.Label.String = 'number of counts';
clbar.TickLabels = round(clbar.Ticks.^(1/IntPower),1);
xlabel('mm');ylabel('mm');axis equal;axis tight;
% ca = gca;
% ca.Position = [.05 .05 .8 .9];
% BarWidth = round(1./(ULM.SRscale*lambda)); % 1 mm
% im.CData(size(MatOut{algo_idx},1)-50+[0:3],60+[0:BarWidth])=max(caxis);
