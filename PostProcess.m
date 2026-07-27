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

%% Step 1: Confocal preprocessing - keep high-intensity pixels
CF_2D = rgb2gray(CF_raw);

% Calculate the high-intensity threshold.
threshold_CF = prctile(CF_2D(:), 96);

% Keep pixels above threshold and set the rest to zero.
CF_processed = CF_2D;
CF_processed(CF_2D < threshold_CF) = 0;

% Preview the processed confocal image.
figure('Name', 'CF preprocessing preview', 'NumberTitle', 'off');
imshow(CF_processed,[]);

%% Step 2: Extract confocal skeleton candidate pixels
skeleton_mask = CF_processed > 0;  % true = candidate skeleton pixel

% Get row/column coordinates of candidate skeleton pixels.
[skel_rows, skel_cols] = find(skeleton_mask);
fprintf('Number of CF skeleton candidate pixels: %d\n', numel(skel_rows));

if isempty(skel_rows)
    error('CF skeleton is empty. Check the threshold or input image.');
end

figure;
imshow(skeleton_mask);
title('Skeleton mask (white = candidate skeleton pixels)');
%% Step 2-New: Manually select valid skeleton ROI
figure('Name', 'Select valid skeleton ROI', 'NumberTitle', 'off');
imshow(skeleton_mask, []);
title('Draw the valid skeleton ROI, then double-click to finish');

% Draw a polygon ROI around the valid structure/skeleton region.
hROI = drawpolygon('Color', 'y', 'LineWidth', 1.5);

% Create an ROI mask from the polygon.
validROI_mask = createMask(hROI);

% Keep only skeleton pixels inside the selected ROI.
skeleton_mask_roi = skeleton_mask & validROI_mask;

% Remove small noisy regions.
skeleton_mask_roi = bwareaopen(skeleton_mask_roi, 5);

% Check whether any skeleton pixels remain inside the ROI.
if ~any(skeleton_mask_roi(:))
    error('No valid skeleton points inside the ROI. Reselect ROI or lower the threshold.');
end

% Replace the original skeleton with the ROI-filtered skeleton.
skeleton_mask = skeleton_mask_roi;

% Recalculate skeleton coordinates after ROI filtering.
[skel_rows, skel_cols] = find(skeleton_mask);

fprintf('Number of valid skeleton points inside ROI: %d\n', numel(skel_rows));

figure('Name', 'Valid skeleton ROI result', 'NumberTitle', 'off');
imshow(skeleton_mask, []);
title('Valid skeleton retained inside ROI');
%% Step 3: Precompute distance map to the skeleton
% bwdist computes the Euclidean distance to the nearest skeleton pixel.
dist_map = bwdist(skeleton_mask);  % Same size as the confocal image
%% ===== Step 3-New: Estimate local skeleton direction vectors =====
% Parameters for local skeleton direction and trajectory alignment analysis.
SKEL_WINDOW    = 10;   % Skeleton direction estimation window (px), suggested 5-15
TRAJ_SMOOTH    = 5;    % Trajectory smoothing window (frames), suggested 3-7
DIST_THRESHOLD = 20;   % Distance threshold for near-skeleton analysis (px)
ANGLE_THRESHOLD = 30;  % Parallel-alignment threshold (deg), suggested 20-45

fprintf('Estimating local skeleton direction vectors...\n');
skeletonPts = [skel_cols, skel_rows];   % [x, y] format
N_skel = size(skeletonPts, 1);
skeletonDirs = zeros(N_skel, 3);        % [x, y, angle_rad]

for i = 1:N_skel
    cx = skeletonPts(i, 1);
    cy = skeletonPts(i, 2);
    dists_skel = sqrt((skeletonPts(:,1)-cx).^2 + (skeletonPts(:,2)-cy).^2);
    neighbors = skeletonPts(dists_skel <= SKEL_WINDOW & dists_skel > 0, :);

    if size(neighbors, 1) >= 2
        centered = neighbors - mean(neighbors);
        [~, ~, V] = svd(centered, 'econ');
        tangent = V(:, 1);
        angle = atan2(tangent(2), tangent(1));
    else
        angle = NaN;
    end
    skeletonDirs(i, :) = [cx, cy, angle];
end

% Remove skeleton points without a valid local direction estimate.
skeletonDirs = skeletonDirs(~isnan(skeletonDirs(:, 3)), :);
fprintf('Number of valid skeleton direction points: %d\n', size(skeletonDirs, 1));
%% Step 4: Calculate UCNP-to-skeleton distance for each frame
ROI_X = 258;   
ROI_Y = 168;    
ROI_W = 28;   
ROI_H = 93;   
TRACKIT_COORDINATES_ARE_ROI_LOCAL = true;
data = load('TrackResult.mat');
spotsAll = data.trackitBatch.results.spotsAll;

num_frames = size(spotsAll, 1);

pixel_size = 0.2654;  % um/pixel

% Extract trajectory coordinates for distance and direction analysis.
traj_x = NaN(num_frames, 1);
traj_y = NaN(num_frames, 1);

frames_kept = [];
dist_px     = [];

for f = 1:num_frames
    if isempty(spotsAll{f}) || numel(spotsAll{f}) < 2
        continue;
    end
    x_full = spotsAll{f}(1);
    y_full = spotsAll{f}(2);

    if TRACKIT_COORDINATES_ARE_ROI_LOCAL
        col = x_full;
        row = y_full;
    else
        col = x_full - ROI_X + 1;
        row = y_full - ROI_Y + 1;
    end

    if col < 1 || col > ROI_W || row < 1 || row > ROI_H
        warning('Frame %d coordinate (%.1f, %.1f) is outside ROI and was skipped.', f, x_full, y_full);
        continue;
    end

    row_idx = round(row);
    col_idx = round(col);
    d_px    = dist_map(row_idx, col_idx);

    frames_kept(end+1) = f;       %#ok<AGROW>
    dist_px(end+1)     = d_px;    %#ok<AGROW>
    traj_x(f)          = col;
    traj_y(f)          = row;
end

fprintf('Total frames: %d, valid frames: %d\n', num_frames, numel(frames_kept));
%% ===== Step 4-New: Visualize local skeleton directions =====
% Map local skeleton direction to 0-180 deg for color display.
skel_angle_deg = mod(rad2deg(skeletonDirs(:, 3)), 180);

figure('Name', 'Skeleton Direction Map', 'NumberTitle', 'off');

% Panel 1: color-code skeleton points by local direction.
subplot(1, 2, 1);
imshow(skeleton_mask, []); hold on;

% Scatter color indicates local direction angle.
scatter(skeletonDirs(:,1), skeletonDirs(:,2), ...
    8, skel_angle_deg, 'filled', 'MarkerEdgeColor', 'none');

colormap(hsv);
cb = colorbar;
clim([0 180]);
cb.Label.String = 'Local direction (deg)';
cb.Label.FontSize = 11;
cb.Ticks = [0 45 90 135 180];
cb.TickLabels = {'0', '45', '90', '135', '180'};

title('Skeleton: local direction (color-coded)', 'FontSize', 12);
xlabel('X (px)'); ylabel('Y (px)');
axis image;

% Panel 2: show subsampled bidirectional skeleton tangents.
subplot(1, 2, 2);
imshow(skeleton_mask, []); hold on;

ARROW_STEP   = 8;    % Plot one tangent for every 8 skeleton direction points.
ARROW_LEN    = 6;    % Half-length of each displayed tangent (px).

idx_show = 1:ARROW_STEP:size(skeletonDirs, 1);
sx  = skeletonDirs(idx_show, 1);
sy  = skeletonDirs(idx_show, 2);
ang = skeletonDirs(idx_show, 3);

% Convert direction angle to x/y vector components.
dx_arr = ARROW_LEN * cos(ang);
dy_arr = ARROW_LEN * sin(ang);

% Draw bidirectional local tangent segments.
quiver(sx, sy,  dx_arr,  dy_arr, 0, ...
    'Color', [1 0.6 0.1], 'LineWidth', 0.8, ...
    'MaxHeadSize', 0.6, 'AutoScale', 'off');
quiver(sx, sy, -dx_arr, -dy_arr, 0, ...
    'Color', [1 0.6 0.1], 'LineWidth', 0.8, ...
    'MaxHeadSize', 0, 'AutoScale', 'off');

title('Skeleton: direction arrows (subsampled)', 'FontSize', 12);
xlabel('X (px)'); ylabel('Y (px)');
axis image;

%% ===== Step 5-New: Calculate trajectory-skeleton angle differences =====
% Extract valid trajectory frames and coordinates.
valid_frames = frames_kept;
tx_raw = traj_x(valid_frames);
ty_raw = traj_y(valid_frames);

% Smooth trajectory coordinates to reduce localization noise.
if TRAJ_SMOOTH > 1 && numel(tx_raw) > TRAJ_SMOOTH
    tx_sm = movmean(tx_raw, TRAJ_SMOOTH, 'omitnan');
    ty_sm = movmean(ty_raw, TRAJ_SMOOTH, 'omitnan');
else
    tx_sm = tx_raw;
    ty_sm = ty_raw;
end

% Calculate frame-to-frame trajectory directions and compare with local skeleton direction.
N_valid = numel(valid_frames);
angleDiff_all  = NaN(N_valid-1, 1);  % Angle difference folded to 0-90 deg.
nearestDist_ad = NaN(N_valid-1, 1);  % Distance to nearest skeleton direction point.
isParallel_all = false(N_valid-1, 1);
midX_all = NaN(N_valid-1, 1);
midY_all = NaN(N_valid-1, 1);

for t = 1:N_valid-1
    dx = tx_sm(t+1) - tx_sm(t);
    dy = ty_sm(t+1) - ty_sm(t);
    stepLen = sqrt(dx^2 + dy^2);

    % Skip nearly stationary steps.
    if stepLen < 0.1 || isnan(dx) || isnan(dy)
        continue;
    end

    % Midpoint of the trajectory step.
    mx = (tx_sm(t) + tx_sm(t+1)) / 2;
    my = (ty_sm(t) + ty_sm(t+1)) / 2;
    midX_all(t) = mx;
    midY_all(t) = my;

    % Trajectory direction angle.
    traj_angle = atan2(dy, dx);

    % Find the nearest skeleton point with a valid local direction.
    d_to_skel = sqrt((skeletonDirs(:,1) - mx).^2 + (skeletonDirs(:,2) - my).^2);
    [minDist, idx] = min(d_to_skel);
    nearestDist_ad(t) = minDist;

    % Only analyze steps close to the skeleton.
    if minDist <= DIST_THRESHOLD
        skel_angle = skeletonDirs(idx, 3);

        % Fold the undirected angle difference into 0-90 deg.
        raw_diff  = abs(traj_angle - skel_angle);
        raw_diff  = mod(raw_diff, pi);
        angleDiff = min(raw_diff, pi - raw_diff);
        angleDiff_all(t)  = rad2deg(angleDiff);
        isParallel_all(t) = angleDiff_all(t) <= ANGLE_THRESHOLD;
    end
end

% Summary statistics for the frame-to-frame angle analysis.
validAngle_idx = ~isnan(angleDiff_all);
validAngles    = angleDiff_all(validAngle_idx);
n_parallel     = sum(isParallel_all(validAngle_idx));
n_total_angle  = sum(validAngle_idx);

fprintf('\n===== Frame-to-frame direction analysis =====\n');
fprintf('Number of analyzed steps: %d\n', n_total_angle);
if n_total_angle > 0
    fprintf('Median angle difference: %.2f deg\n', median(validAngles, 'omitnan'));
    fprintf('Mean angle difference: %.2f deg\n', mean(validAngles, 'omitnan'));
    fprintf('Parallel-step fraction (<= %d deg): %.1f%%\n', ...
        ANGLE_THRESHOLD, 100*n_parallel/n_total_angle);
end
fprintf('========================\n\n');

%% ===== Step 6: Visualize distance and angle analyses =====

% Figure 1: distance to skeleton over time.
figure('Name', 'Minimum distance', 'NumberTitle', 'off');
bar(frames_kept, dist_px);
xlabel('Frame Index', 'FontSize', 13);
ylabel('Minimum distance (px)', 'FontSize', 13);
title('Distance between Nanoparticles and Skeleton', 'FontSize', 14);
grid on;
xlim([1 num_frames]);

% Figure 2: distribution of UCNP-to-skeleton distance.
figure('Name', 'Distance distribution histogram', 'NumberTitle', 'off');
histogram(dist_px);
xlabel('Minimum distance (px)', 'FontSize', 13);
ylabel('Count', 'FontSize', 13);
title('Distance Distribution Histogram', 'FontSize', 14);
grid on;
maxDistForXlim = max(dist_px, [], 'omitnan');
if isempty(maxDistForXlim) || isnan(maxDistForXlim) || maxDistForXlim <= 0
    maxDistForXlim = 1;
end
xlim([0, maxDistForXlim * 1.05]);

% Figure 3: spatial map of angle difference.
figure('Name', 'Angle difference spatial map', 'NumberTitle', 'off');
imshow(skeleton_mask); hold on;
ax = gca;
scatter(ax, midX_all(validAngle_idx), midY_all(validAngle_idx), ...
    30, validAngles, 'filled');
colormap(ax, jet);
cb = colorbar(ax);
clim([0 90]);
cb.Label.String = 'Angle difference (deg)';
cb.Label.FontSize = 12;
title('Spatial Map: Trajectory-Skeleton Angle Difference', 'FontSize', 13);
xlabel('X (px)'); ylabel('Y (px)');

% Figure 4: angle difference histogram.
figure('Name', 'Angle difference histogram', 'NumberTitle', 'off');
histogram(validAngles, 18, 'FaceColor', [0.25 0.55 0.85], ...
    'EdgeColor', 'white', 'Normalization', 'probability');
hold on;
xline(ANGLE_THRESHOLD, 'r--', 'LineWidth', 2, ...
    'Label', sprintf('Parallel threshold %d deg', ANGLE_THRESHOLD), ...
    'LabelVerticalAlignment', 'bottom');
xlabel('Angle difference (deg)', 'FontSize', 13);
ylabel('Probability', 'FontSize', 13);
title('Trajectory-Skeleton Angle Difference Distribution', 'FontSize', 14);
xlim([0 90]);
grid on;

% Figure 5: parallel vs non-parallel fraction.
figure('Name', 'Parallelism assessment', 'NumberTitle', 'off');
n_nonparallel = n_total_angle - n_parallel;
pie([n_parallel, n_nonparallel], ...
    {sprintf('Parallel (<=%d deg)\n%.1f%%', ANGLE_THRESHOLD, 100*n_parallel/n_total_angle), ...
     sprintf('Non-parallel (>%d deg)\n%.1f%%', ANGLE_THRESHOLD, 100*n_nonparallel/n_total_angle)});
colormap([0.20 0.70 0.40; 0.90 0.30 0.30]);
title('Parallelism Assessment', 'FontSize', 14);

% Figure 6: joint distance-angle analysis.
figure('Name', 'Distance vs Angle', 'NumberTitle', 'off');
scatter(nearestDist_ad(validAngle_idx), validAngles, ...
    30, validAngles, 'filled', 'MarkerEdgeColor', 'none');
colormap(jet); colorbar;
clim([0 90]);
hold on;
xline(DIST_THRESHOLD, 'k--', 'LineWidth', 1.5, ...
    'Label', sprintf('Distance threshold %d px', DIST_THRESHOLD));
yline(ANGLE_THRESHOLD, 'r--', 'LineWidth', 1.5, ...
    'Label', sprintf('Parallel threshold %d deg', ANGLE_THRESHOLD));
xlabel('Distance to skeleton (px)', 'FontSize', 13);
ylabel('Angle difference (deg)', 'FontSize', 13);
title('Distance vs Angle Difference (Joint Analysis)', 'FontSize', 14);
grid on;
%% ===== Step 6-New: PCA/window summary for UCNP-structure alignment =====
% This block preserves the previous figures and adds a 2x3 summary figure.
% Goals:
% 1) UCNP points stay close to the structure skeleton.
% 2) Windowed UCNP motion is directionally aligned with local skeleton direction.
% 3) The global skeleton PCA axis and UCNP PCA motion axis are aligned.

TRAJ_WINDOW = 10;  % non-empty UCNP trajectory points per displacement window

if numel(tx_raw) < 2 || numel(ty_raw) < 2 || size(skeletonDirs, 1) < 2
    warning(['Not enough valid UCNP trajectory points or skeleton direction points ', ...
             'for PCA/window summary. Check ROI_X, ROI_Y, ROI_W, ROI_H and coordinate mapping.']);
else

% --- Global PCA axis of the skeleton ---
skelXY_pca = [skel_cols(:), skel_rows(:)];
skelCenter_pca = mean(skelXY_pca, 1);
skelCentered_pca = skelXY_pca - skelCenter_pca;
[skelEigVec, skelEigVal] = eig(cov(skelCentered_pca));
[~, skelMainIdx] = max(diag(skelEigVal));
skelMainAxis = skelEigVec(:, skelMainIdx);

% --- Global PCA axis of the UCNP trajectory ---
trajXY_pca = [tx_raw(:), ty_raw(:)];
validTrajPca = all(~isnan(trajXY_pca), 2);
trajXY_pca = trajXY_pca(validTrajPca, :);
trajCenter_pca = mean(trajXY_pca, 1);
trajCentered_pca = trajXY_pca - trajCenter_pca;
[trajEigVec, trajEigVal] = eig(cov(trajCentered_pca));
[~, trajMainIdx] = max(diag(trajEigVal));
trajMainAxis = trajEigVec(:, trajMainIdx);

firstToLast = [tx_raw(end) - tx_raw(1); ty_raw(end) - ty_raw(1)];
if dot(trajMainAxis, firstToLast) < 0
    trajMainAxis = -trajMainAxis;
end
if dot(skelMainAxis, trajMainAxis) < 0
    skelMainAxis = -skelMainAxis;
end

trajMainAngleDeg = atan2d(trajMainAxis(2), trajMainAxis(1));
skelMainAngleDeg = atan2d(skelMainAxis(2), skelMainAxis(1));
mainAxisDiffDeg = abs(mod(trajMainAngleDeg - skelMainAngleDeg + 90, 180) - 90);
trajProjectedPosition = ([tx_raw(:), ty_raw(:)] - trajXY_pca(1, :)) * trajMainAxis;

% --- 10-point windowed UCNP displacement vs nearest skeleton tangent ---
window_mid_frame = [];
window_mid_x = [];
window_mid_y = [];
window_projected_position = [];
window_dist_to_skel_px = [];
window_angle_diff_deg = [];
window_speed_px_per_frame = [];
window_is_near = [];
window_is_parallel = [];

for t = 1:(N_valid - TRAJ_WINDOW)
    t2 = t + TRAJ_WINDOW;
    dxw = tx_sm(t2) - tx_sm(t);
    dyw = ty_sm(t2) - ty_sm(t);
    frameGapW = valid_frames(t2) - valid_frames(t);
    stepLenW = hypot(dxw, dyw);

    if stepLenW < 0.1 || frameGapW <= 0 || isnan(dxw) || isnan(dyw)
        continue;
    end

    mxw = (tx_sm(t) + tx_sm(t2)) / 2;
    myw = (ty_sm(t) + ty_sm(t2)) / 2;
    trajAngleW = atan2(dyw, dxw);

    d_to_skel_w = hypot(skeletonDirs(:,1) - mxw, skeletonDirs(:,2) - myw);
    [minDistW, idxW] = min(d_to_skel_w);
    skelAngleW = skeletonDirs(idxW, 3);

    rawDiffW = abs(trajAngleW - skelAngleW);
    rawDiffW = mod(rawDiffW, pi);
    angleDiffW = min(rawDiffW, pi - rawDiffW);
    angleDiffDegW = rad2deg(angleDiffW);

    window_mid_frame(end+1, 1) = mean([valid_frames(t), valid_frames(t2)]); %#ok<AGROW>
    window_mid_x(end+1, 1) = mxw; %#ok<AGROW>
    window_mid_y(end+1, 1) = myw; %#ok<AGROW>
    window_projected_position(end+1, 1) = ([mxw, myw] - trajXY_pca(1, :)) * trajMainAxis; %#ok<AGROW>
    window_dist_to_skel_px(end+1, 1) = minDistW; %#ok<AGROW>
    window_angle_diff_deg(end+1, 1) = angleDiffDegW; %#ok<AGROW>
    window_speed_px_per_frame(end+1, 1) = stepLenW / frameGapW; %#ok<AGROW>
    window_is_near(end+1, 1) = minDistW <= DIST_THRESHOLD; %#ok<AGROW>
    window_is_parallel(end+1, 1) = angleDiffDegW <= ANGLE_THRESHOLD; %#ok<AGROW>
end

nearWindowIdx = window_is_near == 1;
nearAndParallelIdx = nearWindowIdx & window_is_parallel == 1;

if any(nearWindowIdx)
    nearWindowMedianAngle = median(window_angle_diff_deg(nearWindowIdx), 'omitnan');
    nearParallelFraction = 100 * sum(nearAndParallelIdx) / sum(nearWindowIdx);
else
    nearWindowMedianAngle = NaN;
    nearParallelFraction = NaN;
end

fprintf('\n===== PCA main-axis and 10-point window summary =====\n');
fprintf('UCNP PCA axis angle: %.2f deg\n', trajMainAngleDeg);
fprintf('Skeleton PCA axis angle: %.2f deg\n', skelMainAngleDeg);
fprintf('UCNP-skeleton PCA axis difference: %.2f deg\n', mainAxisDiffDeg);
fprintf('Window count: %d\n', numel(window_angle_diff_deg));
fprintf('Near-window fraction (<= %d px): %.1f%%\n', DIST_THRESHOLD, 100 * mean(nearWindowIdx));
if any(nearWindowIdx)
    fprintf('Near-window median angle difference: %.2f deg\n', nearWindowMedianAngle);
    fprintf('Near and parallel window fraction: %.1f%%\n', nearParallelFraction);
end
fprintf('===============================================\n\n');

% --- Export numerical summary ---
summaryTable = table( ...
    mainAxisDiffDeg, ...
    100 * mean(dist_px <= DIST_THRESHOLD), ...
    median(dist_px, 'omitnan'), ...
    100 * mean(nearWindowIdx), ...
    nearWindowMedianAngle, ...
    nearParallelFraction, ...
    'VariableNames', { ...
        'pca_axis_difference_deg', ...
        'frame_fraction_within_distance_threshold_percent', ...
        'median_frame_distance_px', ...
        'window_fraction_within_distance_threshold_percent', ...
        'near_window_median_angle_difference_deg', ...
        'near_and_parallel_window_fraction_percent' ...
    });
writetable(summaryTable, 'UCNP_structure_alignment_summary.xlsx', 'Sheet', 'summary');

windowTable = table( ...
    window_mid_frame, ...
    window_mid_x, ...
    window_mid_y, ...
    window_projected_position, ...
    window_dist_to_skel_px, ...
    window_angle_diff_deg, ...
    window_speed_px_per_frame, ...
    window_is_near, ...
    window_is_parallel, ...
    'VariableNames', { ...
        'window_mid_frame', ...
        'window_mid_x_px', ...
        'window_mid_y_px', ...
        'position_on_ucnp_pca_axis_px', ...
        'distance_to_skeleton_px', ...
        'angle_difference_deg', ...
        'speed_px_per_frame', ...
        'within_distance_threshold', ...
        'within_angle_threshold' ...
    });
writetable(windowTable, 'UCNP_structure_alignment_summary.xlsx', 'Sheet', 'window_analysis');

N_POSITION_BINS = 20;
if numel(window_projected_position) >= 2
    posEdgesForTable = linspace(min(window_projected_position), max(window_projected_position), N_POSITION_BINS + 1);
    if min(window_projected_position) == max(window_projected_position)
        posEdgesForTable = [min(window_projected_position) - 0.5, max(window_projected_position) + 0.5];
    end

    posCentersForTable = posEdgesForTable(1:end-1) + diff(posEdgesForTable) / 2;
    nWindowsByPosition = zeros(numel(posCentersForTable), 1);
    medianAngleByPositionTable = nan(numel(posCentersForTable), 1);
    q25AngleByPositionTable = nan(numel(posCentersForTable), 1);
    q75AngleByPositionTable = nan(numel(posCentersForTable), 1);

    for b = 1:numel(posCentersForTable)
        if b < numel(posCentersForTable)
            inBin = window_projected_position >= posEdgesForTable(b) & window_projected_position < posEdgesForTable(b + 1);
        else
            inBin = window_projected_position >= posEdgesForTable(b) & window_projected_position <= posEdgesForTable(b + 1);
        end

        binAngles = window_angle_diff_deg(inBin);
        nWindowsByPosition(b) = nnz(inBin);
        if ~isempty(binAngles)
            medianAngleByPositionTable(b) = median(binAngles, 'omitnan');
            q25AngleByPositionTable(b) = prctile(binAngles, 25);
            q75AngleByPositionTable(b) = prctile(binAngles, 75);
        end
    end

    positionBinTable = table( ...
        posCentersForTable(:), ...
        nWindowsByPosition(:), ...
        medianAngleByPositionTable(:), ...
        q25AngleByPositionTable(:), ...
        q75AngleByPositionTable(:), ...
        'VariableNames', { ...
            'position_on_ucnp_pca_axis_bin_center_px', ...
            'n_windows', ...
            'median_angle_difference_deg', ...
            'q25_angle_difference_deg', ...
            'q75_angle_difference_deg' ...
        });
    writetable(positionBinTable, 'UCNP_structure_alignment_summary.xlsx', 'Sheet', 'position_angle_bins');
end

% --- 2x3 summary figure ---
figure('Name', 'UCNP-Structure Alignment Summary', 'NumberTitle', 'off', 'Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: skeleton, UCNP trajectory, and PCA axes
nexttile;
imshow(skeleton_mask, []); hold on;
scatter(tx_raw, ty_raw, 14, valid_frames, 'filled');
plot(tx_raw, ty_raw, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
axisLen = max([range(skel_cols), range(skel_rows), range(tx_raw), range(ty_raw)]) * 0.35;
skelStart = skelCenter_pca(:) - skelMainAxis * axisLen;
skelEnd = skelCenter_pca(:) + skelMainAxis * axisLen;
trajStart = trajCenter_pca(:) - trajMainAxis * axisLen;
trajEnd = trajCenter_pca(:) + trajMainAxis * axisLen;
plot([skelStart(1), skelEnd(1)], [skelStart(2), skelEnd(2)], 'r-', 'LineWidth', 2);
plot([trajStart(1), trajEnd(1)], [trajStart(2), trajEnd(2)], 'c-', 'LineWidth', 2);
cb = colorbar; cb.Label.String = 'Frame';
title(sprintf('PCA axes, diff = %.1f deg', mainAxisDiffDeg));
xlabel('X (px)'); ylabel('Y (px)'); axis image;

% Panel 2: distance over time
nexttile;
plot(trajProjectedPosition, dist_px, '-o', 'LineWidth', 1.1, 'MarkerSize', 3); hold on;
yline(DIST_THRESHOLD, 'k--', sprintf('%d px', DIST_THRESHOLD));
xlabel('Position on UCNP PCA axis (px)'); ylabel('Distance to skeleton (px)');
title('UCNP-structure distance'); grid on;

% Panel 3: distance distribution
nexttile;
histogram(dist_px, 'FaceColor', [0.35 0.55 0.85], 'EdgeColor', 'white'); hold on;
xline(DIST_THRESHOLD, 'k--', sprintf('%d px', DIST_THRESHOLD));
xlabel('Distance to skeleton (px)'); ylabel('Count');
title('Distance distribution histogram'); grid on;

% Panel 4: window angle difference along the UCNP PCA axis
nexttile;
scatter(window_projected_position, window_angle_diff_deg, 18, ...
    [0.65 0.65 0.65], 'filled', ...
    'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45); hold on;

N_POSITION_BINS = 20;
if numel(window_projected_position) >= 2
    posEdges = linspace(min(window_projected_position), max(window_projected_position), N_POSITION_BINS + 1);
    if min(window_projected_position) == max(window_projected_position)
        posEdges = [min(window_projected_position) - 0.5, max(window_projected_position) + 0.5];
    end

    posCenters = posEdges(1:end-1) + diff(posEdges) / 2;
    medianAngleByPosition = nan(numel(posCenters), 1);
    q25AngleByPosition = nan(numel(posCenters), 1);
    q75AngleByPosition = nan(numel(posCenters), 1);

    for b = 1:numel(posCenters)
        if b < numel(posCenters)
            inBin = window_projected_position >= posEdges(b) & window_projected_position < posEdges(b + 1);
        else
            inBin = window_projected_position >= posEdges(b) & window_projected_position <= posEdges(b + 1);
        end

        binAngles = window_angle_diff_deg(inBin);
        if ~isempty(binAngles)
            medianAngleByPosition(b) = median(binAngles, 'omitnan');
            q25AngleByPosition(b) = prctile(binAngles, 25);
            q75AngleByPosition(b) = prctile(binAngles, 75);
        end
    end

    validPositionBins = ~isnan(medianAngleByPosition);
    if any(validPositionBins)
        errorbar(posCenters(validPositionBins), medianAngleByPosition(validPositionBins), ...
            medianAngleByPosition(validPositionBins) - q25AngleByPosition(validPositionBins), ...
            q75AngleByPosition(validPositionBins) - medianAngleByPosition(validPositionBins), ...
            '-o', ...
            'Color', [0.00 0.35 0.75], ...
            'MarkerFaceColor', [0.00 0.35 0.75], ...
            'MarkerEdgeColor', [0.00 0.35 0.75], ...
            'LineWidth', 1.5, ...
            'MarkerSize', 4, ...
            'CapSize', 5);
    end
end
yline(ANGLE_THRESHOLD, 'r--', sprintf('%d deg', ANGLE_THRESHOLD));
xlabel('Position on UCNP PCA axis (px)'); ylabel('Angle difference (deg)');
title(sprintf('%d-point window angle difference', TRAJ_WINDOW));
ylim([0 90]); grid on;

% Panel 5: angle difference distribution for near windows
nexttile;
if any(nearWindowIdx)
    histogram(window_angle_diff_deg(nearWindowIdx), 0:5:90, ...
        'FaceColor', [0.25 0.65 0.45], 'EdgeColor', 'white', 'Normalization', 'probability'); hold on;
    xline(ANGLE_THRESHOLD, 'r--', sprintf('%d deg', ANGLE_THRESHOLD));
    title(sprintf('Near-window median = %.1f deg', median(window_angle_diff_deg(nearWindowIdx), 'omitnan')));
else
    text(0.5, 0.5, 'No near windows', 'HorizontalAlignment', 'center');
    title('Near-window angle distribution');
end
xlabel('Angle difference (deg)'); ylabel('Probability');
xlim([0 90]); grid on;

% Panel 6: joint distance-angle analysis
nexttile;
scatter(window_dist_to_skel_px, window_angle_diff_deg, 28, window_speed_px_per_frame, ...
    'filled', 'MarkerEdgeColor', 'none'); hold on;
xline(DIST_THRESHOLD, 'k--', sprintf('%d px', DIST_THRESHOLD));
yline(ANGLE_THRESHOLD, 'r--', sprintf('%d deg', ANGLE_THRESHOLD));
colormap(parula); cb = colorbar; cb.Label.String = 'Window speed (px/frame)';
xlabel('Distance to skeleton (px)'); ylabel('Angle difference (deg)');
title('Distance-angle joint analysis');
ylim([0 90]); grid on;

sgtitle('UCNP proximity and directional alignment with structure skeleton');

savefig('UCNP_structure_alignment_summary.fig');
exportgraphics(gcf, 'UCNP_structure_alignment_summary.png', 'Resolution', 300);
end
