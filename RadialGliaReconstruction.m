    %% Initial
    clc; clear; close all;
    %% Loading data
    
    [file, path] = uigetfile({'*.tif;*.tiff;*.ome.tif'}, 'choose 3D Tiff data');
    if isequal(file, 0), return; end
    fullPath = fullfile(path, file);
    info = imfinfo(fullPath);
    num_slices = numel(info);
    fprintf('Loading %d slices image...\n', num_slices);
    data = zeros(info(1).Height, info(1).Width, num_slices, 'single');
        for i = 1:num_slices
            data(:,:,i) = single(imread(fullPath, i));
        end
    data_full_disp = data;
    % volumeViewer(data);
    % %% Data Cropped
    % % ref_mip = max(data, [], 3);                         % Z轴最大强度投影
    % % ref_disp = double(ref_mip);
    % % ref_disp = ref_disp / max(ref_disp(:) + eps);       % 归一化到[0,1]便于显示
    % % 
    % % fig_roi = figure('Name', 'ROI Selection — Draw rectangle on MIP', ...
    % %                  'Position', [200 200 900 700]);
    % % imshow(ref_disp, [], 'InitialMagnification', 'fit');
    % % title({'在 MIP 图上拖拽矩形框选 XY ROI 区域', ...
    % %        '框选完成后双击矩形确认（或按 Enter）'}, 'FontSize', 10);
    % % colormap(gca, gray);
    % % 
    % % % --- 2. 交互框选 XY ROI ---
    % % h_rect = drawrectangle('Color','cyan', 'LineWidth', 1.5, ...
    % %                         'Label','ROI','LabelTextColor','cyan');
    % % % 等待用户双击或按 Enter 确认
    % % input_msg = 'ROI 框选后按 Enter 确认...';
    % % fprintf('%s\n', input_msg);
    % % wait(h_rect);                   % 阻塞直到用户完成编辑
    % 
    % % roi_pos = round(h_rect.Position);   % [x, y, width, height]，x/y 为左上角
    % % close(fig_roi);
    roi_pos = [477, 410, 104, 237];

    % 防止 ROI 超出图像边界
    x1 = max(1, roi_pos(1));
    y1 = max(1, roi_pos(2));
    x2 = min(size(data,2), x1 + roi_pos(3) - 1);
    y2 = min(size(data,1), y1 + roi_pos(4) - 1);
    roi_w = x2 - x1 + 1;
    roi_h = y2 - y1 + 1;

    z_start = 1;
    z_end   = 26;
    if z_start > z_end
        error('Z 层范围无效：z_start (%d) > z_end (%d)', z_start, z_end);
    end

    data = data(y1:y2, x1:x2, z_start:z_end);
    % 
    %% Isotropic modification
    pixel_size_xy = 0.207;
    pixel_size_z = 1.0;
    z_scale = pixel_size_z / pixel_size_xy;
    V_iso = imresize3(data, [size(data,1), size(data,2), round(size(data,3)*z_scale)], 'linear');
    % volumeViewer(V_iso);
    %% preprocessing and denoising
    V_median = medfilt3(V_iso, [5 5 5]);
    V_median = medfilt3(V_median, [3 3 3]);
    se = strel('sphere', 4); 
    V_tophat = imtophat(V_median, se);
    volumeViewer(V_tophat);
    %% normalisation
    V_max = max(V_tophat(:));
    V_norm = double(V_tophat) / double(V_max);
    % volumeViewer(V_norm);
    %% Binarization
    thresh_manual = 0.08;
    thresh_otsu   = thresh_manual;
     
    V_bin = V_norm > thresh_otsu;
     
    min_vol_voxels = 50;                       
    V_bin = bwareaopen(V_bin, min_vol_voxels);  
    V_bin = imfill(V_bin, 'holes');             
     
    volumeViewer(V_bin);
    %% skeletonisation
    V_skel = bwskel(V_bin, 'MinBranchLength', 30);
    volumeViewer(V_skel);
    %% overlay
    V_iso_disp = double(V_iso);
    iso_max = max(V_iso_disp(:));
    if iso_max > 0, V_iso_disp = V_iso_disp / iso_max; end

    % -----------------------------------------------------------------
    % 7.3  3D 骨架叠加渲染（isosurface）
    % -----------------------------------------------------------------
    figure('Name','3D骨架叠加渲染','Position',[820 100 700 600]);
    scale = 1;   % 降采样比例（加速渲染）
    iso_small  = imresize3(V_iso_disp,  scale, 'linear');
    bin_small  = imresize3(single(V_bin),  scale, 'nearest');
    skel_small = imresize3(single(V_skel), scale, 'nearest');

    % 原始数据半透明等值面
    p_iso = patch(isosurface(iso_small, 0.3));
    set(p_iso, 'FaceColor',[0.75 0.85 1.0], 'EdgeColor','none', 'FaceAlpha',0.10);
    hold on;

    % 二值化前景轮廓（淡蓝，极低透明度提示范围）
    p_bin = patch(isosurface(bin_small, 0.5));
    set(p_bin, 'FaceColor',[0.3 0.6 1.0], 'EdgeColor','none', 'FaceAlpha',0.08);

    % 骨架（亮橙，不透明）
    p_skel = patch(isosurface(skel_small, 0.5));
    set(p_skel, 'FaceColor',[1.0 0.45 0.1], 'EdgeColor','none', 'FaceAlpha',1.0);

    axis equal tight; lighting gouraud;
    camlight('headlight'); camlight('right');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Rendering：gray=Initial，blue=foreground，orange=skeleton');
    legend([p_iso p_bin p_skel], {'V\_iso','Foreground','Skeleton'}, 'Location','best');
    view(3); rotate3d on;
%% 
    V_viewer = V_iso_disp;        % 复制原始体，保留真实强度
    V_viewer(V_skel) = 1.0;       % 骨架位置强制写为最大值
    volumeViewer(V_viewer);
    mip_skeleton = max(V_skel, [], 3);
    figure
    imshow(mip_skeleton, []);
%% Skeleton MIP
    mip_z = max(V_viewer, [], 3); 
    figure
    imshow(mip_z, []);
   
    %% 骨架叠加回原始全图高亮显示
    orig_H = size(data_full_disp, 1);  % 还不存在，改用已知量

    % 用 imfinfo 重新获取尺寸，避免 info 变量冲突
    img_info  = imfinfo(fullPath);
    orig_H    = img_info(1).Height;
    orig_W    = img_info(1).Width;
    orig_Z    = numel(img_info);
    % --- Step 1：重新加载原始全图（crop前）用作背景 ---
    data_full = zeros(orig_H, orig_W, orig_Z, 'single');
    for i = 1:orig_Z
        data_full(:,:,i) = single(imread(fullPath, i));
    end
    data_full_disp = double(data_full) / (max(double(data_full(:))) + eps);
    
    % --- Step 2：将 V_skel 从各向同性空间缩回 crop 空间 ---
    % crop 空间尺寸：[roi_h, roi_w, z_end-z_start+1]
    % num_z_crop = z_end - z_start + 1;
    % V_skel_crop = imresize3(single(V_skel), ...
    %     [roi_h, roi_w, num_z_crop], 'nearest') > 0;   % logical mask，crop尺寸
    num_z_iso  = size(V_skel, 3);
    num_z_crop = z_end - z_start + 1;
    
    V_skel_crop = false(roi_h, roi_w, num_z_crop);
    for z = 1:num_z_crop
        z_iso_start = round((z-1) * z_scale) + 1;
        z_iso_end   = min(round(z * z_scale), num_z_iso);
        V_skel_crop(:,:,z) = any(V_skel(:,:,z_iso_start:z_iso_end), 3);
    end
    
    % --- Step 3：将 crop 空间骨架贴回全图坐标 ---
    skel_full = false(orig_H, orig_W, orig_Z);
    skel_full(y1:y2, x1:x2, z_start:z_end) = V_skel_crop;
    
    % --- Step 4：骨架写入全图体数据 → volumeViewer ---
    V_full_viewer = data_full_disp;
    V_full_viewer(skel_full) = 1.0;
    volumeViewer(V_full_viewer);
    
    % --- Step 5：MIP 叠加图（全图坐标）---
    mip_data_full = max(data_full_disp, [], 3);
    mip_skel_full = max(skel_full,      [], 3);
    
    rgb_overlay = repmat(mip_data_full, [1, 1, 3]);
    % 骨架高亮为白色
    rgb_overlay(:,:,1) = max(rgb_overlay(:,:,1), double(mip_skel_full));
    rgb_overlay(:,:,2) = max(rgb_overlay(:,:,2), double(mip_skel_full));
    rgb_overlay(:,:,3) = max(rgb_overlay(:,:,3), double(mip_skel_full));
    
    figure('Name', '骨架 MIP 叠加至原始全图', 'Position', [100 100 900 700]);
    imshow(rgb_overlay, 'InitialMagnification', 'fit');
    
    % 标注 ROI 边框位置
    hold on;
    rectangle('Position', [x1, y1, roi_w, roi_h], ...
              'EdgeColor', 'cyan', 'LineWidth', 1.5, 'LineStyle', '--');
    title('骨架 MIP 高亮叠加至原始全图（青色框=ROI区域）', 'FontSize', 11);
    out_tif = fullfile(path, [file(1:end-4) '_skeleton_overlay3.tif']);
    V_save  = uint16(V_full_viewer / max(V_full_viewer(:)) * 65535);  % 转为16bit
    for i = 1:size(V_save, 3)
        if i == 1
            imwrite(V_save(:,:,i), out_tif, 'tif', ...
                'WriteMode', 'overwrite', 'Compression', 'none');
        else
            imwrite(V_save(:,:,i), out_tif, 'tif', ...
                'WriteMode', 'append',    'Compression', 'none');
        end
    end
    fprintf('体数据已保存: %s\n', out_tif);

%%以WF时间帧为x轴，算纳米粒子和骨架的最短距离，距离做histogram，看距离分布。距离分布再做nominal distance
