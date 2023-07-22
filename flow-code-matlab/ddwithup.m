close all
clear all

folder_raft = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\estimation_flo\';
folder_tiling = 'C:\Users\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\after_tiling';
folder_gt = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\alley_1';
folder_up = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwice\';
folder_BIup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBI\';
folder_BCup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBC\';
outrefine_BI = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRBI\';
%outrefine_AT = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRtiling\';

outrefine_AT = 'C:\Users\Chengyu\Desktop\R64\';
outrefine_UP = 'C:\Users\Chengyu\Desktop\newRUP16';


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
    
    % Read flow data from the folder_upsample_twice
    file_up = sprintf('%s\\flow_000%02d_000%02d_upsampled.flo', folder_up, i, i+1);
    flow_up = readFlowFile(file_up);

    % Read flow data from the folder_upsample_twice_Bilinear_interpolate
    file_BIup = sprintf('%s\\flow_000%02d_000%02d_BIupsampled.flo', folder_BIup, i, i+1);
    flow_BIup = readFlowFile(file_BIup);

    % Read flow data from the folder_upsample_twice_BiCUBIC_interpolate
    file_BCup = sprintf('%s\\flow_000%02d_000%02d_BCupsampled.flo', folder_BCup, i, i+1);
    flow_BCup = readFlowFile(file_BCup);

    % Read flow data from the output refinement BI
    file_outBI = sprintf('%s\\refined_flow_%04d.flo', outrefine_BI, i);
    flow_outBI = readFlowFile(file_outBI);

    % Read flow data from the output refinement UP
    file_outUP = sprintf('%s\\refined_flow_%04d.flo', outrefine_UP, i);
    flow_outUP = readFlowFile(file_outUP);


    % Read flow data from the output refinement AT
    file_outAT = sprintf('%s\\refined_flow_%04d.flo', outrefine_AT, i);
    flow_outAT = readFlowFile(file_outAT);

    % MAE for flow_raft and flow_gt
    abs_error_raft = abs(flow_raft - flow_gt);
    mae_raft = mean(abs_error_raft(:));
    mae_values_raft(i) = mae_raft;

    %MAE for flow_tiling and flow_gt
    abs_error_tiling = abs(flow_tiling - flow_gt);
 %   mae_tiling = mean(abs_error_tiling(:));
      mae_tiling = mean(mean(abs_error_tiling(190:230,:)));
    mae_values_tiling(i) = mae_tiling;

    %MAE for flow_upsampling_duplication and flow_gt
    abs_error_up = abs(flow_up - flow_gt);
   % mae_up = mean(abs_error_up(:));
mae_up = mean(mean(abs_error_up(64:200,:)));

    mae_values_up(i) = mae_up;

    %MAE for flow_upsampling_Bilinear and flow_gt
    abs_error_BIup = abs(flow_BIup - flow_gt);
    mae_BIup = mean(abs_error_BIup(:));
    mae_values_BIup(i) = mae_BIup;

    %MAE for flow_upsampling_Bicubic and flow_gt
    abs_error_BCup = abs(flow_BCup - flow_gt);
    mae_BCup = mean(abs_error_BCup(:));
    mae_values_BCup(i) = mae_BCup;

    %MAE for OUTREFINE_BI and flow_gt
    abs_error_outBI = abs(flow_outBI - flow_gt);
    mae_outBI = mean(abs_error_outBI(:));
    mae_values_outBI(i) = mae_outBI;

    %MAE for OUTREFINE_UP and flow_gt
    abs_error_outUP = abs(flow_outUP - flow_gt);
 %   mae_outUP = mean(abs_error_outUP(:));
   mae_outUP = mean(mean(abs_error_outUP(190:230,:)));

    mae_values_outUP(i) = mae_outUP;


    %MAE for OUTREFINE_BI and flow_gt
    abs_error_outAT = abs(flow_outAT - flow_gt);
 %   mae_outAT = mean(abs_error_outAT(:));
    mae_outAT = mean(mean(abs_error_outAT(64:200,:)));
    mae_values_outAT(i) = mae_outAT;






%% MAE between the motion compensated frames using estimated RAFT motion OR GT
    image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
    image_file1 = sprintf('frame_%04d.png', i);
    image_file2 = sprintf('frame_%04d.png', i+1);
    frame_n = double(imread([image_path image_file1])) ;
    frame_np1 = double(imread([image_path image_file2])) ;


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

    % RAFT with upsample-twice
     I_np1_up = warp_image_using_interp2(frame_np1, flow_up); 
     up_abs_diff = abs(frame_n - I_np1_up); 
     N = size(up_abs_diff, 1); M = size(up_abs_diff, 2);
     up_MCMAE(i) = sum(up_abs_diff(:)) / (N * M);

    % RAFT with upsample-twice-Bilinear-interpolate
     I_np1_BIup = warp_image_using_interp2(frame_np1, flow_BIup); 
     BIup_abs_diff = abs(frame_n - I_np1_BIup); 
     N = size(BIup_abs_diff, 1); M = size(BIup_abs_diff, 2);
     BIup_MCMAE(i) = sum(BIup_abs_diff(:)) / (N * M); 

     % RAFT with upsample-twice-Bicubic-interpolate
     I_np1_BCup = warp_image_using_interp2(frame_np1, flow_BCup); 
     BCup_abs_diff = abs(frame_n - I_np1_BCup); 
     N = size(BCup_abs_diff, 1); M = size(BCup_abs_diff, 2);
     BCup_MCMAE(i) = sum(BCup_abs_diff(:)) / (N * M);

     % RAFT with re
     I_np1_outAT = warp_image_using_interp2(frame_np1, flow_outAT); 
     outAT_abs_diff = abs(frame_n - I_np1_outAT); 
     N = size(outAT_abs_diff, 1); M = size(outAT_abs_diff, 2);
 %    outAT_MCMAE(i) = sum(sum(outAT_abs_diff(64:200,:))) / (N * M);
    outAT_MCMAE(i) = sum(outAT_abs_diff(:)) / (N * M);

         % RAFT with re
     I_np1_outup = warp_image_using_interp2(frame_np1, flow_outUP); 
     outup_abs_diff = abs(frame_n - I_np1_outup); 
     N = size(outup_abs_diff, 1); M = size(outup_abs_diff, 2);
 %    outup_MCMAE(i) = sum(sum(outup_abs_diff(64:200,:))) / (N * M);
    outup_MCMAE(i) = sum(outup_abs_diff(:)) / (N * M);


end

%frame difference  
figure(1);
plot(1:num_files, mae_values_raft, 'o-', 'LineWidth', 2);
hold on
plot(1:num_files, mae_values_tiling, 'x-', 'LineWidth', 2);
xlabel('frame numbers');
ylabel('Mean Absolute Error (MAE)');
grid on;
hold on;
plot(1:num_files, mae_values_up, '+-', 'LineWidth', 2);
plot(1:num_files, mae_values_BIup, '-*', 'LineWidth', 2);
plot(1:num_files, mae_values_BCup, '-o', 'LineWidth', 2);
%plot(1:num_files, mae_values_outBI, '-o', 'LineWidth', 2);
%plot(1:num_files, mae_values_outAT, '-o', 'LineWidth', 2);
plot(1:num_files, mae_values_outUP, '-o', 'LineWidth', 2);

title('frame difference error');
%legend('RAFT vs. Ground Truth', 'Tiling vs. GT', 'UP vs. GT','BIUP vs. GT', 'BCUP vs. GT','outBI vs. GT', 'outAT vs. GT')
legend('RAFT vs. Ground Truth', 'Tiling vs. GT', 'UP vs. GT','BIUP vs. GT', 'BCUP vs. GT', 'outAT vs. GT')
grid on;




%pixel motion compensated MAE
figure(2);
plot(frame_numbers, R_MCMAE, 'LineWidth', 2);
xlabel('Frame Number');
ylabel('Motion Compensation MAE ');
grid on;
hold on;
plot(frame_numbers, T_MCMAE, 'LineWidth', 2);
plot(frame_numbers, gt_MCMAE, 'LineWidth', 2); 
%plot(frame_numbers, nmc_MAE, 'LineWidth', 1); 
plot(frame_numbers, up_MCMAE, 'LineWidth', 2);
plot(frame_numbers, BIup_MCMAE, 'LineWidth', 2);
plot(frame_numbers, BCup_MCMAE, 'LineWidth', 2);


plot(frame_numbers, outup_MCMAE, 'r', 'LineWidth', 2);
%plot(frame_numbers, outAT_MCMAE, 'r', 'LineWidth', 2);
%title('MAE with MC using raft, tiling , gt , without MC , Upsampled , BIupsampled & BCupsampled ')
title('MAE with MC using raft, tiling , gt , Upsampled , BIupsampled & BCupsampled ')
legend('raft only', 'raft with tiling', 'gt' ,'upsampled','BIupsampled', 'BCupsampled' ,'outup')
%legend('raft only', 'raft with tiling', 'gt' ,'upsampled','BIupsampled', 'BCupsampled' ,'outAT')
%legend('raft only', 'raft with tiling', 'gt', 'no MC','upsampled','BIupsampled', 'BCupsampled' )
hold off;


