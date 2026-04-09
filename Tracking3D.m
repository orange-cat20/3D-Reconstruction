clear; clc; close all;
 
%% ========== 1. 加载数据 ==========
fprintf('正在加载 trackitBatch.mat ...\n');
data = load('D6F2_TrackResult.mat');
 
% 提取 spotsAll（cell数组，每个cell对应一帧）
spotsAll = data.trackitBatch.results.spotsAll;
nFrames  = length(spotsAll);
fprintf('共加载 %d 帧数据\n', nFrames);
 
%% ========== 2. 提取有效帧的粒子数据 ==========
allX     = [];
allY     = [];
allIntensity = [];
allWidth = [];
allFrame = [];
 
for i = 1:nFrames
    frameData = spotsAll{i};
    
    % 跳过空帧（该帧未追踪到粒子）
    if isempty(frameData)
        continue;
    end
    
    % 提取 x, y（第1、2列）和像素宽度（第5列）
    x     = frameData(:, 1);
    y     = frameData(:, 2);
    intensity = frameData(:, 3);
    width = frameData(:, 5);
    
    % 过滤无效值（NaN / 0 / 负值）
    valid = ~isnan(x) & ~isnan(y) & ~isnan(intensity) & ~isnan(width) & (width > 0);
    
    allX     = [allX;     x(valid)];
    allY     = [allY;     y(valid)];
    allIntensity = [allIntensity; intensity(valid)];
    allWidth = [allWidth; width(valid)];
    allFrame = [allFrame; repmat(i, sum(valid), 1)];
end
 
nParticles = length(allX);
fprintf('有效粒子检测总数: %d\n', nParticles);
fprintf('像素宽度范围: %.4f ~ %.4f pixels\n', min(allWidth), max(allWidth));
%% ========== 3. 绘制双Y轴直方图 ==========
fprintf('正在绘制直方图...\n');

% 定义帧区间（bin边界）
nBins = 50;  % 直方图分组数，可根据需要调整
binEdges = linspace(1, nFrames, nBins + 1);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

% 计算每个帧区间内第3列（像素强度）和第5列（像素宽度）的均值
intensity_mean = zeros(1, nBins);
width_mean     = zeros(1, nBins);

for b = 1:nBins
    % 找到属于当前帧区间的粒子
    inBin = (allFrame >= binEdges(b)) & (allFrame < binEdges(b+1));
    
    if any(inBin)
        intensity_mean(b) = mean(allIntensity(inBin));
        width_mean(b)     = mean(allWidth(inBin));
    else
        intensity_mean(b) = 0;
        width_mean(b)     = 0;
    end
end

% ---- 创建图窗 ----
figure('Name', 'Histogram：Intensity & Width', 'NumberTitle', 'off', ...
       'Position', [100, 100, 900, 500]);

% ---- 左Y轴：像素强度（第3列） ----
yyaxis left
bar(binCenters, intensity_mean, 1, ...
    'FaceColor', [0.2 0.5 0.8], ...
    'FaceAlpha', 0.7, ...
    'EdgeColor', 'none');
ylabel('Intensity', 'FontSize', 12);
ylim([0, max(intensity_mean) * 1.2]);

% ---- 右Y轴：像素宽度（第5列） ----
yyaxis right
plot(binCenters, width_mean, '-o', ...
    'Color', [0.85 0.35 0.1], ...
    'LineWidth', 1.8, ...
    'MarkerSize', 4, ...
    'MarkerFaceColor', [0.85 0.35 0.1]);
ylabel('Width', 'FontSize', 12);
ylim([0, max(width_mean) * 1.2]);

% ---- 坐标轴修饰 ----
xlabel('Frames', 'FontSize', 12);
title('Intensity and width distribution', 'FontSize', 14, 'FontWeight', 'bold');
xlim([1, nFrames]);
legend({'Intensity（left）', 'width（right）'}, ...
       'Location', 'northwest', 'FontSize', 10);
grid on;
box on;

fprintf('直方图绘制完成。\n');
 
% %% ========== 3. Z轴映射 ==========
% % 物理依据：
% %   - PSF宽度（sigma）与z轴离焦量的关系近似抛物线（高斯光学近似）：
% %       sigma(z) = sigma0 * sqrt(1 + (z/d)^2)
% %   - 在有效追踪范围内，假设像素宽度与z轴位置单调对应
% %
% % 映射策略：
% %   最小宽度 → 焦平面（z范围中心）
% %   宽度增大 → 向z范围两端扩展 
% %   此处采用抛物线模型（可切换为线性模型）
% 
% z_min = 5508;
% z_max = 5523;
% z_center = (z_min + z_max) / 2;   % 焦平面对应z
% z_half_range = (z_max - z_min) / 2;
% 
% width_min = min(allWidth);
% width_max = max(allWidth);
% 
% % --- 抛物线映射（推荐）---
% % sigma(z) = sigma_min * sqrt(1 + ((z - z_center)/depth)^2)
% % 反推：z - z_center = ±depth * sqrt((sigma/sigma_min)^2 - 1)
% % depth参数由宽度范围标定
% depth_param = z_half_range;  % 可根据标定曲线调整
% 
% % 归一化宽度（相对于最小宽度）
% width_norm = allWidth / width_min;
% 
% % 抛物线反映射（取绝对值再根据帧序列判断上下方向）
% delta_z_abs = depth_param * sqrt(max(width_norm.^2 - 1, 0));
% 
% % 利用帧编号奇偶模拟轴向扫描（若数据是轴向扫描序列请修改此部分）
% % 如果数据不是轴向扫描，默认将宽度大的粒子映射到z两端，中心为焦平面
% % 此处简化为：随帧序号线性分配上/下半空间的符号
% frame_sign = sign(allFrame - median(allFrame));
% frame_sign(frame_sign == 0) = 1;
% 
% allZ = z_center + frame_sign .* delta_z_abs;
% 
% % 钳位到z范围内
% allZ = max(min(allZ, z_max), z_min);
% 
% fprintf('Z轴映射完成，Z范围: %.2f ~ %.2f\n', min(allZ), max(allZ));
% 
% %% ========== 4. 构建三维体积数据（用于volumeviewer）==========
% % 将离散粒子坐标离散化到体素网格
% 
% % 定义体素网格分辨率
% voxel_xy = 0.207;          % xy方向：1 pixel/voxel（与原始图像一致）
% voxel_z  = 1;        % z方向：0.1单位/voxel（可调整）
% 
% x_range = [floor(min(allX)), ceil(max(allX))];
% y_range = [floor(min(allY)), ceil(max(allY))];
% z_range_vox = [z_min, z_max];
% 
% % 体素坐标轴
% x_axis = x_range(1) : voxel_xy : x_range(2);
% y_axis = y_range(1) : voxel_xy : y_range(2);
% z_axis = z_range_vox(1) : voxel_z : z_range_vox(2);
% 
% nx = length(x_axis);
% ny = length(y_axis);
% nz = length(z_axis);
% 
% fprintf('体积网格尺寸: %d × %d × %d (X×Y×Z)\n', nx, ny, nz);
% 
% % 将粒子坐标索引化
% xi = round((allX - x_axis(1)) / voxel_xy) + 1;
% yi = round((allY - y_axis(1)) / voxel_xy) + 1;
% zi = round((allZ - z_axis(1)) / voxel_z)  + 1;
% 
% % 钳位到有效范围
% xi = max(1, min(xi, nx));
% yi = max(1, min(yi, ny));
% zi = max(1, min(zi, nz));
% 
% % 构建三维体积数组（粒子计数/密度图）
% vol = zeros(ny, nx, nz, 'uint16');   % MATLAB数组顺序：(row=y, col=x, slice=z)
% 
% for k = 1:nParticles
%     vol(yi(k), xi(k), zi(k)) = vol(yi(k), xi(k), zi(k)) + 1;
% end
% 
% % 高斯平滑（使粒子点在体积中可见）
% sigma_blur = 1.5;  % 体素单位
% vol_smooth = imgaussfilt3(single(vol), sigma_blur);
% 
% fprintf('体积数据构建完成，最大密度值: %.2f\n', max(vol_smooth(:)));
% 
% %% ========== 5. 保存结果 ==========
% % 保存粒子坐标表
% T = table(allX, allY, allZ, allWidth, allFrame, ...
%     'VariableNames', {'X_pixel','Y_pixel','Z_pos','PSF_Width','Frame'});
% save('particle_tracks_3D.mat', 'T', 'vol_smooth', 'x_axis', 'y_axis', 'z_axis');
% fprintf('三维坐标已保存至 particle_tracks_3D.mat\n');
% 
% %% ========== 6. 二维预览图 ==========
% figure('Name', 'Particle Distribution Overview', 'Position', [100, 100, 1200, 400]);
% 
% subplot(1,3,1);
% scatter(allX, allY, 5, allZ, 'filled', 'MarkerFaceAlpha', 0.4);
% colormap(gca, jet); colorbar;
% xlabel('X (pixels)'); ylabel('Y (pixels)');
% title('XY Projection（Color=Z position）');
% axis equal tight;
% 
% subplot(1,3,2);
% scatter(allX, allZ, 5, allZ, 'filled', 'MarkerFaceAlpha', 0.4);
% colormap(gca, jet); colorbar;
% xlabel('X (pixels)'); ylabel('Z position');
% title('XZ Projection');
% 
% subplot(1,3,3);
% histogram(allWidth, 50, 'FaceColor', '#0072BD');
% xlabel('PSF Width (pixels)'); ylabel('Count');
% title('PSF Width Distribution');
% 
% sgtitle('Single Particle Tracking - 3D Distribution Preview');
% 
% %% ========== 7. Volume Viewer 可视化 ==========
% fprintf('\n正在启动 volumeviewer...\n');
% fprintf('使用提示：\n');
% fprintf('  - 在 volumeviewer 中调整 "Rendering Style" 为 "MaximumIntensityProjection"\n');
% fprintf('    或 "VolumeRendering" 获得最佳三维效果\n');
% fprintf('  - 调整 "Colormap" 和亮度滑块以突出粒子\n');
% fprintf('  - 可在 "Viewer" 菜单中导出截图\n\n');
% 
% % 创建空间参考对象（指定真实物理坐标）
% R = imref3d(size(vol_smooth), ...
%     [x_axis(1), x_axis(end)], ...
%     [y_axis(1), y_axis(end)], ...
%     [z_axis(1), z_axis(end)]);
% 
% volumeViewer(vol_smooth);
% 
% fprintf('volumeviewer 已启动，共显示 %d 个有效粒子的三维分布\n', nParticles);
% fprintf('Z轴物理范围：%.1f ~ %.1f\n', z_min, z_max);