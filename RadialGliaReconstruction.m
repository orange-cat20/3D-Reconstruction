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
    % ref_mip = max(data, [], 3);                         % Z轴最大强度投影
    % ref_disp = double(ref_mip);
    % ref_disp = ref_disp / max(ref_disp(:) + eps);       % 归一化到[0,1]便于显示
    % 
    % fig_roi = figure('Name', 'ROI Selection — Draw rectangle on MIP', ...
    %                  'Position', [200 200 900 700]);
    % imshow(ref_disp, [], 'InitialMagnification', 'fit');
    % title({'在 MIP 图上拖拽矩形框选 XY ROI 区域', ...
    %        '框选完成后双击矩形确认（或按 Enter）'}, 'FontSize', 10);
    % colormap(gca, gray);
    % 
    % % --- 2. 交互框选 XY ROI ---
    % h_rect = drawrectangle('Color','cyan', 'LineWidth', 1.5, ...
    %                         'Label','ROI','LabelTextColor','cyan');
    % % 等待用户双击或按 Enter 确认
    % input_msg = 'ROI 框选后按 Enter 确认...';
    % fprintf('%s\n', input_msg);
    % wait(h_rect);                   % 阻塞直到用户完成编辑
    % 
    % roi_pos = round(h_rect.Position);   % [x, y, width, height]，x/y 为左上角
    % close(fig_roi);
    roi_pos = [445, 227, 139, 316];

    % 防止 ROI 超出图像边界
    x1 = max(1, roi_pos(1));
    y1 = max(1, roi_pos(2));
    x2 = min(size(data,2), x1 + roi_pos(3) - 1);
    y2 = min(size(data,1), y1 + roi_pos(4) - 1);
    roi_w = x2 - x1 + 1;
    roi_h = y2 - y1 + 1;

    z_start = 1;
    z_end   = num_slices;

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
    thresh_manual = 0.20;
    thresh_otsu   = thresh_manual;
     
    V_bin = V_norm > thresh_otsu;
     
    min_vol_voxels = 50;                       
    V_bin = bwareaopen(V_bin, min_vol_voxels);  
    V_bin = imfill(V_bin, 'holes');             
     
    volumeViewer(V_bin);
    %%  Fill and closing operations
    se_close = strel('sphere', 6);    
    solid_vessels = imclose(V_bin, se_close);    
    solid_vessels = imfill(solid_vessels, 'holes');    
    % 在进行 bwskel 之前    
    solid_vessels = imopen(solid_vessels, strel('sphere', 1)); % 抹平表面    
    solid_vessels = imclose(solid_vessels, strel('sphere', 1)); % 再次连接       
    volumeViewer(solid_vessels); 

    %% skeletonisation
    V_skel = bwskel(solid_vessels, 'MinBranchLength', 15);
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
    set(p_iso, 'FaceColor',[0.75 0.85 1.0], 'EdgeColor','none', 'FaceAlpha',0);
    % hold on;

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
    %% Skeleton ROI Selection — keep only skeleton within polygon
    
    % --- 在骨架 MIP 上手动绘制多边形，保留多边形内的骨架 ---
    mip_skel_select = max(single(V_skel), [], 3);   % 骨架MIP用于选取参考
    
    % 同时显示原始数据MIP作为背景，便于对照解剖结构定位
    mip_iso_bg = max(V_iso, [], 3);
    
    % 构建显示图：背景灰度 + 骨架叠加为亮色
    fig_poly = figure('Name', 'Draw polygon to select skeleton region', ...
                      'Position', [100 100 1000 800]);
    imshow(mip_iso_bg, [], 'InitialMagnification', 'fit');
    hold on;
    % 骨架用红色叠加显示
    skel_rgb_overlay = cat(3, mip_skel_select, zeros(size(mip_skel_select)), ...
                               zeros(size(mip_skel_select)));
    h_skel_layer = imshow(skel_rgb_overlay);
    set(h_skel_layer, 'AlphaData', mip_skel_select * 0.8);  % 半透明叠加
    title({'在图上绘制多边形框选需要保留的骨架区域', ...
           '双击多边形内部或按 Enter 确认'}, 'FontSize', 11);
    
    % --- 交互绘制多边形 ---
    h_poly = drawpolygon('Color', 'cyan', 'LineWidth', 1.5);
    fprintf('多边形绘制完成后按 Enter 确认...\n');
    wait(h_poly);
    
    % 获取多边形顶点坐标
    poly_vertices = h_poly.Position;   % [N×2]，每行为 [x, y]
    close(fig_poly);
    
    % --- 生成多边形 mask（XY平面，与 V_skel XY 尺寸一致）---
    poly_mask = poly2mask(poly_vertices(:,1), poly_vertices(:,2), ...
                          size(V_skel,1), size(V_skel,2));   % [H×W] logical
    
    % --- 将 2D mask 扩展到3D，逐层与骨架做 AND ---
    poly_mask_3d = repmat(poly_mask, [1, 1, size(V_skel,3)]);
    V_skel = V_skel & poly_mask_3d;   % 只保留多边形内的骨架体素
    
    fprintf('多边形选取完成，保留骨架体素数: %d\n', sum(V_skel(:)));
    
    % --- 预览选取结果 ---
    mip_skel_after = max(single(V_skel), [], 3);
    figure('Name', '多边形选取后的骨架 MIP', 'Position', [100 100 900 700]);
    rgb_check = repmat(mip_iso_bg, [1,1,3]);
    rgb_check(:,:,1) = max(rgb_check(:,:,1), mip_skel_after);  % 红色高亮保留部分
    rgb_check(:,:,2) = rgb_check(:,:,2) .* ~logical(mip_skel_after);
    rgb_check(:,:,3) = rgb_check(:,:,3) .* ~logical(mip_skel_after);
    % 多边形边界叠加
    imshow(rgb_check, 'InitialMagnification', 'fit');
    hold on;
    pgon_closed = [poly_vertices; poly_vertices(1,:)];  % 闭合多边形
    plot(pgon_closed(:,1), pgon_closed(:,2), 'c-', 'LineWidth', 1.5);
    title('选取后骨架（红色）+ 多边形边界（青色）', 'FontSize', 11);

%% 
    V_viewer = V_iso_disp;        % 复制原始体，保留真实强度
    V_viewer(V_skel) = 1.0;       % 骨架位置强制写为最大值
    % volumeViewer(V_viewer);
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
    % volumeViewer(V_full_viewer);
    
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
    out_tif = fullfile(path, [file(1:end-4) '_skeleton_overlay7.tif']);
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

