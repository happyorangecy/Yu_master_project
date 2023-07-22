close all;
clear;

image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
gt_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1';
d_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\after_tiling';
outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRtiling\';

startFrame = 1;
endFrame = 5;
blockSize = [128, 128];
max_iters = 5;  % Increased iteration count
convergence_threshold = 1e-6;  % Convergence threshold for the motion field difference

num_files = endFrame - startFrame;

mae_values_gt = zeros(1, num_files);
mae_values_d_original = zeros(1, num_files);
mae_values_d_updated = zeros(1, num_files);

% Iterate over frames
for i = startFrame:endFrame-1
    fprintf('Processing frame %d\n', i);
    
    % Read the frames
    image_file1 = sprintf('frame_%04d.png', i);
    image_file2 = sprintf('frame_%04d.png', i+1);
    In_pre = double(imread([image_path image_file1])) / 255;
    I_n_curr = double(imread([image_path image_file2])) / 255;
    
    % Read the ground truth and estimated motion flow fields
    file_gt = sprintf('%s\\frame_%04d.flo', gt_flow_path, i);
    flow_gt = readFlowFile(file_gt);
    
    file_d = sprintf('%s\\flow_000%02d_000%02d.flo', d_flow_path, i, i+1);
    flow_d = readFlowFile(file_d);
    
    % Get the size of the frames
    [rows, cols, ~] = size(In_pre);

    
    % Convert frames to grayscale
    I_n_prev_gray =  255 * double(im2gray(In_pre));
    I_n_curr_gray = 255 * double(im2gray(I_n_curr));
    
    % Initialize the updated motion field
    flow_updated = flow_d;
    
    % Perform gradient-based block matching
    for iter = 1:max_iters
        fprintf('Iteration %d\n', iter);
        
        % Initialize the previous motion field
        flow_prev = flow_updated;
        
        % Iterate over blocks
        for y = 1+blockSize(1):blockSize(1):rows-2*blockSize(1)
            for x = 1+blockSize(2):blockSize(2):cols-2*blockSize(2)

                % Get the block flow vectors
                block_flow_gt = flow_gt(y:y+blockSize(1)-1, x:x+blockSize(2)-1, :);
                block_flow_d = flow_d(y:y+blockSize(1)-1, x:x+blockSize(2)-1, :);
                block_flow_updated = flow_updated(y:y+blockSize(1)-1, x:x+blockSize(2)-1, :);
                
                % Compute the predicted frame based on the previous frame and motion vectors
                predicted_frame = fetch_block_from_frame_using_interp2(I_n_prev_gray, x:x+blockSize(2)-1, y:y+blockSize(1)-1, block_flow_updated);
                
                % Compute the gradient of the current frame
                [grad_x, grad_y] = gradient(I_n_curr_gray(y:y+blockSize(1)-1, x:x+blockSize(2)-1));
                
                % Compute the difference between the current frame and the predicted frame
                Z = I_n_curr_gray(y:y+blockSize(1)-1, x:x+blockSize(2)-1) - predicted_frame;
                
                % Compute the updated motion vector using the formula
                G = [grad_x(:), grad_y(:)];
                u = (G' * G) \ (G' * Z(:))  ; 
             %  u= u /((rows*cols)/((blockSize(1)*blockSize(2))));

               
                % Update the motion vectors
                block_flow_updated(:, :, 1) = block_flow_d(:, :, 1) + u(1) ;
                block_flow_updated(:, :, 2) = block_flow_d(:, :, 2) + u(2)  ;
                disp(u);
                 
                % Store the updated motion vectors
                flow_updated(y:y+blockSize(1)-1, x:x+blockSize(2)-1, :) = block_flow_updated ;
            end
        end
        
%         % Calculate the difference between the updated and previous motion fields
%         diff = sum(abs(flow_updated - flow_prev), 'all');
%         
%         % Check the convergence criterion
%         if diff < convergence_threshold
%             fprintf('Convergence reached at iteration %d\n', iter);
%             break;
%         end
    end
    
    % Store the updated motion field
    file_updated = sprintf('%s\\updated_motion_flow_field_%04d.flo', outrefine_flow_path, i);
    writeFlowFile(flow_updated, file_updated);
    
    % Calculate the Mean Absolute Error (MAE) for ground truth and estimated motion fields within specified range
    flow_gt_cropped = flow_gt(blockSize(1):end-blockSize(1), blockSize(2):end-blockSize(2), :);
    flow_d_cropped = flow_d(blockSize(1):end-blockSize(1), blockSize(2):end-blockSize(2), :);
    flow_updated_cropped = flow_updated(blockSize(1):end-blockSize(1), blockSize(2):end-blockSize(2), :);
    mae_values_gt(i) = mean(abs(flow_gt_cropped(:) - flow_d_cropped(:)));
    mae_values_d_original(i) = mean(abs(flow_gt_cropped(:) - flow_d_cropped(:)));
    mae_values_d_updated(i) = mean(abs(flow_gt_cropped(:) - flow_updated_cropped(:)));
end

% Calculate the average MAE over all frames
average_mae_gt = mean(mae_values_gt);
average_mae_d_original = mean(mae_values_d_original);
average_mae_d_updated = mean(mae_values_d_updated);

fprintf('Average MAE for ground truth motion field: %f\n', average_mae_gt);
fprintf('Average MAE for original estimated motion field: %f\n', average_mae_d_original);
fprintf('Average MAE for updated motion field: %f\n', average_mae_d_updated);

% Plot the average MAE values
figure;
t = startFrame:endFrame-1;
%plot(t, mae_values_gt, 'x-', 'LineWidth', 2);
plot(t, mae_values_d_original, '-o', 'LineWidth', 2);
hold on;

plot(t, mae_values_d_updated, '-+', 'LineWidth', 2);
xlabel('Frame numbers');
ylabel('Mean Absolute Error (MAE)');
title('Average MAE Comparison');
%legend('Ground Truth vs. Estimated', 'Ground Truth vs. Original Estimated', 'Ground Truth vs. Updated');
legend( 'Ground Truth vs. Original Estimated', 'Ground Truth vs. Updated');
grid on;

function mc_block = fetch_block_from_frame_using_interp2(imgI_n, bx, by, blockflow)

    % know image dimensions
    [height, width, channels] = size(imgI_n);

    % Create grids for the x and y coordinates
    [x_grid, y_grid] = meshgrid(1:width, 1:height);

    % Compute the block coordinates using the flow vectors
    brows = length(by);
    bcols = length(bx);
    new_bx_coords = ones(brows, 1) * bx + blockflow(:, :, 1);
    new_by_coords = by' * ones(1, bcols)  + blockflow(:, :, 2);

    % Initialize the warped image
    brows = length(by);
    bcols = length(bx);
    mc_block = zeros(brows, bcols, channels);

    % Warp each channel using interp2
    for ch = 1:channels
        mc_block(:, :, ch) = interp2(x_grid, y_grid, imgI_n(:, :, ch), new_bx_coords, new_by_coords, 'linear', 0);
    end

end
