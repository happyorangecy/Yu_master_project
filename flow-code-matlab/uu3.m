% Define Paths
image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
d_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\after_tiling';
outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRtiling\';
ground_truth_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1';

% Initialize
block_size = [120, 120];
num_iterations = 5;

% Initialize error accumulators
tiling_error = 0;
updated_error = 0;
total_pixels = 0;

% Iteratively modify blocks
for img_idx = 1:49
    % Read the image, its corresponding flow field and ground truth flow
    img1 = imread(fullfile(image_path, sprintf('frame_%04d.png', img_idx)));
    img1 = img1(:,:,1);  % only take one channel
    img2 = imread(fullfile(image_path, sprintf('frame_%04d.png', img_idx+1)));
    img2 = img2(:,:,1);  % only take one channel

    file_d = sprintf('%s\\flow_000%02d_000%02d.flo', d_flow_path, img_idx, img_idx+1);
    tiling_flow = readFlowFile(file_d);
    
    file_gt = sprintf('%s\\frame_%04d.flo', ground_truth_path, img_idx);
    gt_flow = readFlowFile(file_gt);

     %MAE for flow_tiling and flow_gt
    abs_error_tiling = abs(tiling_flow - gt_flow);
    mae_tiling = mean(abs_error_tiling(:));
    mae_values_tiling(img_idx) = mae_tiling;


    updated_flow = tiling_flow;
    
    % Break image into 3x4 blocks and apply flow fields
    for row_start = 1:block_size(1):size(img1, 1)
        for col_start = 1:block_size(2):size(img1, 2)
            % Check if we have a full block; if not, continue
            if (row_start+block_size(1)-1 > size(img1, 1)) || (col_start+block_size(2)-1 > size(img1, 2))
                continue;
            end
            
            % Extract blocks
            block_prev = double(img1(row_start:row_start+block_size(1)-1, col_start:col_start+block_size(2)-1, 1)); % if you are taking the first channel
            
            % Compute block motion using flow field
            blockflow = updated_flow(row_start:row_start+block_size(1)-1, col_start:col_start+block_size(2)-1, :);
            
            % Fetch the block from the next frame using the block motion
            block_next = fetch_block_from_frame_using_interp2(img2, 1:block_size(2), 1:block_size(1), blockflow);

            % Perform gradient-based block matching to get new motion vector
            new_u = calcgradest(block_prev, block_next);
            
            % Update the flow field
            updated_flow(row_start:row_start+block_size(1)-1, col_start:col_start+block_size(2)-1, :) = repmat(reshape(new_u, [1, 1, 2]), [block_size(1), block_size(2), 1]);
        end
    end
    
    % Calculate errors
    tiling_error = tiling_error + sum(abs(tiling_flow(:) - gt_flow(:)));
    updated_error = updated_error + sum(abs(updated_flow(:) - gt_flow(:)));
    total_pixels = total_pixels + numel(gt_flow);
end

% Compute MAE
tiling_MAE = tiling_error / total_pixels;
updated_MAE = updated_error / total_pixels;

fprintf('MAE between ground truth and tiling flow: %f\n', tiling_MAE);
fprintf('MAE between ground truth and updated flow: %f\n', updated_MAE);

function [update] = calcgradest(block_prev, block_curr)
    [gx, gy] = gradient(block_prev);
    E = block_curr - block_prev;
    G = [gx(:) gy(:)]';
    GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];
    update = -pinv(GtG) * G * E(:);
    [x,y] = size(block_prev);
    update = update /(x*y);
end

function mc_block = fetch_block_from_frame_using_interp2(imgI_n, bx, by, blockflow)
    % imgI_n: Image frame (2D or 3D array)
    % bx, by: Block coordinates (1D arrays)
    % blockflow: Block motion vectors (2D array)

    % Know image dimensions
    [height, width, channels] = size(imgI_n);

    % Compute the block coordinates using the flow vectors
    brows = length(by);
    bcols = length(bx);

    new_bx_coords = ones(brows, 1) * bx + blockflow(:, :, 1);
    new_by_coords = by' * ones(1, bcols)  + blockflow(:, :, 2);

    % Initialize the warped image
    mc_block = zeros(brows, bcols, channels);

    % Warp each channel using interp2
    for ch = 1:channels
        mc_block(:, :, ch) = interp2(1:width, 1:height, double(imgI_n(:, :, ch)), new_bx_coords, new_by_coords, 'linear', 0);
    end
end

