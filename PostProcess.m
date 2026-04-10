%% Initial
clc; clear; close all;

%% Loading data
% CF
[filename1, pathname1] = uigetfile('*.tif', 'Confocal Data:'); 
if isequal(filename1, 0), return; end
CF_raw = imread(fullfile(pathname1, filename1));

% % WF
% [filename2, pathname2] = uigetfile('*.tif', 'Widefield Data:');
% if isequal(filename2, 0), return; end
% WF = tiffreadVolume([pathname2 filename2]);

%% Step 1: CF 预处理 —— 保留前 20% 高强度像素
CF_2D = rgb2gray(CF_raw);

% 计算前 5% 强度阈值（即第 95 百分位数）
threshold_CF = prctile(CF_2D(:), 97);

% 保留强度 >= 阈值的像素，其余置零
CF_processed = CF_2D;
CF_processed(CF_2D < threshold_CF) = 0;

% 预览处理结果
figure('Name', 'CF 预处理结果预览', 'NumberTitle', 'off');
imshow(CF_processed,[]);

%% Step 2: 提取 CF 骨架（非零像素坐标）
skeleton_mask = CF_processed > 0;  % logical 矩阵，true = 骨架点

% 获取骨架点的行列坐标
[skel_rows, skel_cols] = find(skeleton_mask);
fprintf('CF 骨架点数量: %d\n', numel(skel_rows));

if isempty(skel_rows)
    error('CF 骨架为空，请检查阈值设置或输入数据。');
end

figure;
imshow(skeleton_mask);
title('骨架 mask（白色=骨架点）');
%% Step 3: 利用距离变换预计算骨架距离图
% bwdist 计算每个像素到最近非零点（骨架）的欧氏距离（单位：像素）
dist_map = bwdist(skeleton_mask);  % 大小与 CF 相同
%% Step 4: 逐帧计算纳米粒子到 CF 骨架的距离
ROI_X = 362;   
ROI_Y = 60;    
ROI_W = 85;   
ROI_H = 168;   
data = load('D3f4_TrackResult.mat');
spotsAll = data.trackitBatch.results.spotsAll;

num_frames = size(spotsAll, 1);

pixel_size = 0.2654;  % μm/pixel

frames_kept    = [];
dist_um        = [];

for f = 1:num_frames
    if isempty(spotsAll{f}) || numel(spotsAll{f}) < 2
        continue;
    end
    x_full = spotsAll{f}(1);  % 第1列：X
    y_full = spotsAll{f}(2);  % 第2列：Y
    
    % 转换为 ROI 内坐标
    col = x_full - ROI_X + 1;
    row = y_full - ROI_Y + 1;
    
    % 检查是否落在 ROI 范围内
    if col < 1 || col > ROI_W || row < 1 || row > ROI_H
        warning('第 %d 帧坐标 (%.1f, %.1f) 超出 ROI 范围，已跳过。', f, x_full, y_full);
        continue;
    end
    
    % 四舍五入取整用于查表
    row_idx = round(row);
    col_idx = round(col);
    
    % 查距离图
    d_px = dist_map(row_idx, col_idx);
    d_um = d_px;
    
    frames_kept(end+1) = f;    %#ok<AGROW>
    dist_um(end+1)     = d_um; %#ok<AGROW>
end

fprintf('总帧数: %d，有效帧数: %d\n', num_frames, numel(frames_kept));
%% Step 6: 绘制 距离-帧数 图
figure('Name', 'Minimum distance', 'NumberTitle', 'off');

bar(frames_kept, dist_um);

xlabel('Frame Index', 'FontSize', 13);
ylabel('Minimum distance（px）', 'FontSize', 13);
title('The distance between Nanoparticles and skeleton', 'FontSize', 14);

grid on;
xlim([1 num_frames]);

%% Step 6: Distance distribution histogram
figure('Name', 'Distance distribution histogram', 'NumberTitle', 'off');

histogram(dist_um);

xlabel('Minimum distance（px）', 'FontSize', 13);
ylabel('Times', 'FontSize', 13);
title('Distance distribution histogram', 'FontSize', 14);

grid on;
xlim([0, max(dist_um, [], 'omitnan') * 1.05]);