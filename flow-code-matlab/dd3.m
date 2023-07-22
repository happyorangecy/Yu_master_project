close all
clear all

folder_raft = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\estimation_flo5\';
folder_tiling = 'C:\\Users\\23839\\Desktop\\RAFT-master\\RAFT-master\\after_tiling5';
folder_gt = 'C:\\Users\\23839\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\market_2';
frame_numbers = 1:49;
num_files = 49;
mae_values_raft = zeros(num_files, 1);
mae_values_tiling = zeros(num_files, 1);
R_MCMAE = zeros(num_files, 1);
for i = 1:num_files
%% MAE between the estimated raft motion and the actual ground truth motion
    % Read flow data from folder_raft (no divided sequences)
    file_raft = sprintf('%s\\flow_000%02d_000%02d.flo', folder_raft, i, i+1);
    flow_raft = readFlowFile(file_raft);

    % Read flow data from folder_tiling
    file_tiling = sprintf('%s\\flow_000%02d_000%02d.flo', folder_tiling, i, i+1);
    flow_tiling = readFlowFile(file_tiling);
    
    % Read flow data from the second folder_gt
    file_gt = sprintf('%s\\frame_%04d.flo', folder_gt, i);
    flow_gt = readFlowFile(file_gt);
    
    % MAE for flow_raft and flow_gt
    abs_error_raft = abs(flow_raft - flow_gt);
    mae_raft = mean(abs_error_raft(:));
    mae_values_raft(i) = mae_raft;

    %MAE for flow_tiling and flow_gt
    abs_error_tiling = abs(flow_tiling - flow_gt);
    mae_tiling = mean(abs_error_tiling(:));
    mae_values_tiling(i) = mae_tiling;

%% MAE between the motion compensated frames using estimated RAFT motion OR GT
    image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\market_2\';
    image_file1 = sprintf('frame_%04d.png', i);
    image_file2 = sprintf('frame_%04d.png', i+1);
    frame_n = double(imread([image_path image_file1])) / 255;
    frame_np1 = double(imread([image_path image_file2])) / 255;


   % RAFT ONLY
    I_np1_raft = warp_image_using_interp2(frame_np1, flow_raft); 
    R_abs_diff = abs(frame_n - I_np1_raft); 
%   fprintf('MC Mean Absolute Error F_%d to F_%d : %f\n', i,i+1,  mean(F_abs_diff(:)));
%   image(R_abs_diff);
%   drawnow
    N = size(R_abs_diff, 1); M = size(R_abs_diff, 2);
    R_MCMAE(i) = sum(R_abs_diff(:)) / (N * M);
%     fprintf('With MC MAE raft only %d to %d: %f\n', i,i+1, R_MCMAE(i));

    % RAFT with tiling
     I_np1_tiling = warp_image_using_interp2(frame_np1, flow_tiling); 
     T_abs_diff = abs(frame_n - I_np1_tiling); 
     N = size(T_abs_diff, 1); M = size(T_abs_diff, 2);
     T_MCMAE(i) = sum(T_abs_diff(:)) / (N * M);

     
     % ground truth
     I_np1_gt = warp_image_using_interp2(frame_np1, flow_gt); 
     gt_abs_diff = abs(frame_n - I_np1_gt); 
     N = size(gt_abs_diff, 1); M = size(gt_abs_diff, 2);
     gt_MCMAE(i) = sum(gt_abs_diff(:)) / (N * M);


     % non motion compensation
     nmc_abs_diff = abs(frame_n - frame_np1); 
     nmc_MAE(i) = sum(nmc_abs_diff(:)) / (N * M);


end

figure(1);
plot(1:num_files, mae_values_raft, 'o-', 'LineWidth', 2);
hold on
plot(1:num_files, mae_values_tiling, 'x-', 'LineWidth', 2);
xlabel('frame numbers');
ylabel('Mean Absolute Error (MAE)');
title('MAE for Optical Flow File Pairs');
legend('RAFT vs. Ground Truth', 'Tiling vs. Ground Truth')
grid on;

figure(2);
plot(frame_numbers, R_MCMAE, 'LineWidth', 1);
xlabel('Frame Number');
ylabel('Motion Compensation MAE ');
grid on;
hold on;
plot(frame_numbers, T_MCMAE, 'LineWidth', 1);
plot(frame_numbers, gt_MCMAE, 'LineWidth', 1); 
plot(frame_numbers, nmc_MAE, 'LineWidth', 1); 
title('MAE with MC using raft, tiling and gt & without MC')
legend('raft only', 'raft with tiling', 'gt', 'no MC')
hold off;

