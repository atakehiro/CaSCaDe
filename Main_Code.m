%% tifファイルの読み取り
tic
[file, file_path] = uigetfile('*.tif');
file_info = imfinfo([file_path, file]);
d1 = file_info(1).Height;
d2 = file_info(1).Width;
T = numel(file_info);
bit = file_info(1).BitDepth;
   
im = zeros(d1,d2,T);
for t = 1:T
    im(:,:,t) = imread([file_path, file], t);
end
disp('データ読み取り完了')
toc

%% パラメータ設定
 p.foffset=10; % how many initial frames to exclude in analysis
 p.norm_signal='std'; % ('std','bkg','sub') % different way to normalize intenisty
 p.spf=1 ; % frame rate at acquisition
 % event detection
 p.min_int_ed=0.5; % minimum intenisty value for start-end of a event;
 p.peak_int_ed=2.0; % minimum peak intesnity value for being considered as signal
 p.min_peak_dist_ed=2;
 p.min_peak_length=2;
 % background trending correction
 p.int_correct= 0; % if 1, correct bkg, if 0, no correction.
 % 3D band pass convolution (bpass3d_v1.m)
 p.lb=1; % low bound size for in-x,y dim
 p.hb=11; % high bound size for in-x,y dim
 p.zlb=1;% low bound size for in-z(t) dim
 p.zhb=21; % high bound size for in-z(t) dim
 
 %% 関数に渡す
res=Cal_anl_main2sa_forreview_x(im, p);

%% 図示
figure
imagesc(res.L)
colorbar
title(['Segmentation Result: ', num2str(res.obnum), 'ROIs'])
figure
imagesc(res.intout0')
title('trace')
