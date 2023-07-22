close all
clear all

folder_raft = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\estimation_flo\';
folder_tiling = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\after_tiling';
folder_gt = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\alley_1';
folder_up = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwice\';
folder_BIup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBI\';
folder_BCup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBC\';
outrefine_BI = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRBI\';
outrefine_BC = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRBC\';
outrefine_AT = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRtiling\';
% outrefine_AT = 'C:\Users\Chengyu\Desktop\uuR\';

%frame_numbers = 1:10;
num_files = 5;
mae_values_raft = zeros(num_files, 1);
mae_values_tiling = zeros(num_files, 1);
R_MCMAE = zeros(num_files, 1);

% Read flow data from the second folder_gt
file_gt = sprintf('%s\\frame_%04d.flo', folder_gt, 1);
flow_gt = readFlowFile(file_gt);

for i = 1:num_files
  % for row = 1 +blockSize(1) :blockSize(1) : 218     %rows - blockSize(1)
  %    for col = 1 + blockSize(2) : blockSize(2) : 218    % cols  - blockSize(2)

%     [rows, cols, ~] = size(flow_gt);
%     for row = 1 + blockSize(1) :blockSize(1) :  rows - 1*blockSize(1)
%         for col = 1 + blockSize(2) :blockSize(2) :  cols  - 1*blockSize(2)     
%             rowEnd = min(rows - 1*blockSize(1), rows);
%             colEnd = min(cols  -  1*blockSize(2) , cols);


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
        
            % Read flow data from the output refinement BC
            file_outBC = sprintf('%s\\refined_flow_%04d.flo', outrefine_BC, i);
            flow_outBC = readFlowFile(file_outBC);
        
            % Read flow data from the output refinement BC
            file_outAT = sprintf('%s\\refined_flow_%04d.flo', outrefine_AT, i);
            flow_outAT = readFlowFile(file_outAT);
        
%             block_flow_gt = flow_gt(row:rowEnd, col:colEnd, :);
%             block_flow_tiling = flow_tiling(row:rowEnd, col:colEnd, :);
%             block_flow_outAT = flow_outAT(row:rowEnd, col:colEnd, :);
    
            
            %% MAE for flow_raft and flow_gt

        
            abs_error_raft = abs(flow_raft - flow_gt);
            mae_raft = mean(abs_error_raft(:));
            mae_values_raft(i) = mae_raft;
        
            %MAE for flow_tiling and flow_gt
            abs_error_tiling = abs(flow_tiling - flow_gt);
            mae_tiling = mean(mean(abs_error_tiling(16:100,:)));
            mae_values_tiling(i) = mae_tiling;
        
            %MAE for flow_upsampling_duplication and flow_gt
            abs_error_up = abs(flow_up - flow_gt);
            mae_up = mean(abs_error_up(:));
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
        
            %MAE for OUTREFINE_BC and flow_gt
            abs_error_outBC = abs(flow_outBC - flow_gt);
            mae_outBC = mean(abs_error_outBC(:));
            mae_values_outBC(i) = mae_outBC;
        
            %MAE for OUTREFINE_BC and flow_gt
            abs_error_outAT = abs(flow_outAT - flow_gt);
            mae_outAT = mean(mean(abs_error_outAT(16:100,:)));
            mae_values_outAT(i) = mae_outAT;
 
%             abs_error_tiling = abs(block_flow_tiling - block_flow_gt);
%             mae_tiling = mean(abs_error_tiling(:));
%     
%             abs_error_outAT = abs(block_flow_outAT - block_flow_gt);
%             mae_outAT = mean(abs_error_outAT(:));
%     
%             mae_values_tiling(i) = mae_tiling;
%             mae_values_outAT(i) = mae_outAT;



 %    end
 % end

end
average_mae_d_original = mean(mae_values_tiling);
average_mae_d_updated = mean(mae_values_outAT);

fprintf('Average MAE for original estimated motion field: %f\n', average_mae_d_original);
fprintf('Average MAE for updated motion field: %f\n', average_mae_d_updated);



%frame difference  
figure(1);
%plot(1:num_files, mae_values_raft, 'o-', 'LineWidth', 2);
hold on
plot(1:num_files, mae_values_tiling, 'x-', 'LineWidth', 2);
xlabel('frame numbers');
ylabel('Mean Absolute Error (MAE)');
grid on;
hold on;

plot(1:num_files, mae_values_outAT, '-o', 'LineWidth', 2);
title('vector difference error');

legend( 'Tiling vs. GT',  'outAT vs. GT' )
grid on;






