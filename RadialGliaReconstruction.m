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
    % volumeViewer(data);
    %% Isotropic modification
    pixel_size_xy = 0.207;
    pixel_size_z = 1.0;
    z_scale = pixel_size_z / pixel_size_xy;
    V_iso = imresize3(data, [size(data,1), size(data,2), round(size(data,3)*z_scale)], 'linear');
    % volumeViewer(V_iso);
    %% preprocessing and denoising
    V_median = medfilt3(V_iso, [3 3 3]);
    se = strel('sphere', 5); 
    V_tophat = imtophat(V_median, se);
    % volumeViewer(V_tophat);
    %% normalisation
    V_max = max(V_tophat(:));
    V_norm = double(V_tophat) / double(V_max);
    % volumeViewer(V_norm);
    %% Binarization
    thresh_manual = 0.28;
    thresh_otsu   = thresh_manual;
     
    V_bin = V_norm > thresh_otsu;
     
    min_vol_voxels = 50;                        
    V_bin = bwareaopen(V_bin, min_vol_voxels);  
    V_bin = imfill(V_bin, 'holes');             
     
    % volumeViewer(V_bin);
    %% skeletonisation
    V_skel = bwskel(V_bin, 'MinBranchLength', 10);
    % volumeViewer(V_skel);
    %% overlay
    V_iso_disp = double(V_iso);
    iso_max = max(V_iso_disp(:));
    if iso_max > 0, V_iso_disp = V_iso_disp / iso_max; end

    % % -----------------------------------------------------------------
    % % 7.3  3D 骨架叠加渲染（isosurface）
    % % -----------------------------------------------------------------
    % figure('Name','3D骨架叠加渲染','Position',[820 100 700 600]);
    % scale = 1;   % 降采样比例（加速渲染）
    % iso_small  = imresize3(V_iso_disp,  scale, 'linear');
    % bin_small  = imresize3(single(V_bin),  scale, 'nearest');
    % skel_small = imresize3(single(V_skel), scale, 'nearest');
    % 
    % % 原始数据半透明等值面
    % p_iso = patch(isosurface(iso_small, 0.3));
    % set(p_iso, 'FaceColor',[0.75 0.85 1.0], 'EdgeColor','none', 'FaceAlpha',0.10);
    % hold on;
    % 
    % % 二值化前景轮廓（淡蓝，极低透明度提示范围）
    % p_bin = patch(isosurface(bin_small, 0.5));
    % set(p_bin, 'FaceColor',[0.3 0.6 1.0], 'EdgeColor','none', 'FaceAlpha',0.08);
    % 
    % % 骨架（亮橙，不透明）
    % p_skel = patch(isosurface(skel_small, 0.5));
    % set(p_skel, 'FaceColor',[1.0 0.45 0.1], 'EdgeColor','none', 'FaceAlpha',1.0);
    % 
    % axis equal tight; lighting gouraud;
    % camlight('headlight'); camlight('right');
    % xlabel('X'); ylabel('Y'); zlabel('Z');
    % title('3D Rendering：gray=Initial，blue=foreground，orange=skeleton');
    % legend([p_iso p_bin p_skel], {'V\_iso','Foreground','Skeleton'}, 'Location','best');
    % view(3); rotate3d on;
%% 
    V_viewer = V_iso_disp;        % 复制原始体，保留真实强度
    V_viewer(V_skel) = 1.0;       % 骨架位置强制写为最大值
    % volumeViewer(V_viewer);
%% Skeleton MIP
    mip_z = max(V_viewer, [], 3); 
    figure
    imshow(mip_z, []);
