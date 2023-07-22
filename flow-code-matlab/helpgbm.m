close all;
clear;

image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
gt_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1';
d_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\after_tiling';
%d_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBI\';

outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRtiling\';

startFrame = 1;
endFrame = 50;
blockSize = [16, 16];
max_iters = 5;

num_files = endFrame - startFrame;

mae_values_gt = zeros(1, num_files);
mae_values_d_original = zeros(1, num_files);
mae_values_d_updated = zeros(1, num_files);

for i = startFrame:endFrame-1
    fprintf('Processing frame %d\n', i);
    
    image_file1 = sprintf('frame_%04d.png', i);
    image_file2 = sprintf('frame_%04d.png', i+1);
    In_pre = double(imread([image_path image_file1])) / 255;
    I_n_curr = double(imread([image_path image_file2])) / 255;
    
    file_gt = sprintf('%s\\frame_%04d.flo', gt_flow_path, i);
    flow_gt = readFlowFile(file_gt);
    
    file_d = sprintf('%s\\flow_000%02d_000%02d.flo', d_flow_path, i, i+1);
%   file_d = sprintf('%s\\flow_000%02d_000%02d_BIupsampled.flo', d_flow_path, i, i+1);

    flow_d = readFlowFile(file_d);
    
    [rows, cols, ~] = size(In_pre);
    
    I_n_prev_gray =  255 * double(rgb2gray(In_pre));
    I_n_curr_gray = 255 * double(rgb2gray(I_n_curr));
    
    flow_updated = flow_d ;
    
%    for iter = 1:max_iters
%         fprintf('Iteration %d\n', iter);
%         flow_prev = flow_updated;
         for y = 1+blockSize(1) :blockSize(1):rows-blockSize(1)
             for x = 1+blockSize(2) :blockSize(2):cols-blockSize(2)
                block_flow_d = flow_d(y:y+blockSize(1)-1, x:x+blockSize(2)-1, :);
                block_flow_updated = flow_updated(y:y+blockSize(1)-1, x:x+blockSize(2)-1, :);
                
                SY = (y:y+blockSize(1)-1)' * ones(1, blockSize(2));
                SX = ones(blockSize(1), 1) * (x:x+blockSize(2)-1);

                    for iter = 1:max_iters
                        fprintf('Iteration %d\n', iter);
                        flow_prev = flow_updated;
                        block_prev = interp2(SX, SY, I_n_prev_gray(y:y+blockSize(1)-1, x:x+blockSize(2)-1), SX + block_flow_updated(:,:,1), SY + block_flow_updated(:,:,2), 'linear', 0);

                    [gx, gy] = gradient(block_prev);
                    E = I_n_curr_gray(y:y+blockSize(1)-1, x:x+blockSize(2)-1) - block_prev;
                    G = [gx(:) gy(:)]';
                    GtG = [sum(sum(gx .* gx)) sum(sum(gx .* gy)); sum(sum(gx .* gy)) sum(sum(gy .* gy))];
                   % update = GtG \ (-G * E(:));
                    update = -pinv(GtG) * G * E(:); 

                    end
                  num_blocks = ((rows+1-2*blockSize(1))/blockSize(1)) * ((cols+1-2*blockSize(2))/blockSize(2));
                 update = update /  num_blocks;

                block_flow_updated(:,:,1) = block_flow_d(:,:,1) + update(1)   ;
                block_flow_updated(:,:,2) = block_flow_d(:,:,2) + update(2);
                
                flow_updated(y:y+blockSize(1)-1, x:x+blockSize(2)-1, :) = block_flow_updated;
             end
         end
%    end
    
    if any(isnan(flow_updated(:)))
        fprintf('Warning: flow_updated contains NaN values\n');
    end
    
    file_out = sprintf('%s\\flow_000%02d_000%02d.flo', outrefine_flow_path, i, i+1);
    writeFlowFile(flow_updated, file_out);
    
    flow_gt_cropped = flow_gt(blockSize(1):end-blockSize(1), blockSize(2):end-blockSize(2), :);
    flow_d_cropped = flow_d(blockSize(1):end-blockSize(1), blockSize(2):end-blockSize(2), :);
    flow_updated_cropped = flow_updated(blockSize(1):end-blockSize(1), blockSize(2):end-blockSize(2), :);
    
    mae_values_gt(i) = mean(abs(flow_gt_cropped(:) - flow_d_cropped(:)));
    mae_values_d_original(i) = mean(abs(flow_gt_cropped(:) - flow_d_cropped(:)));
    mae_values_d_updated(i) = mean(abs(flow_gt_cropped(:) - flow_updated_cropped(:)));
end

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