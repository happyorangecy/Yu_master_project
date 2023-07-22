close all;
clear;

% image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
% folder_gt = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\alley_1';
folder_up = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwice\';
folder_BIup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBI\';
folder_BCup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBC\';
folder_tiling = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\after_tiling';

f2= 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\after_tiling5';
u2 = 'C:\Users\Chengyu\Desktop\uptwice5';
%  image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\ambush_5\';
%  folder_gt = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\ambush_5';

%  image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\bamboo_1\';
%  folder_gt = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\bamboo_1';

%  image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\cave_2\';
%  folder_gt = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\cave_2';

 image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\market_2\';
 folder_gt = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\market_2';

outrefine_flow_path = 'C:\Users\Chengyu\Desktop\q6';
% Define output path for refined flow field
%outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRUP\';
%outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRBI\';
%outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRBC\';
%outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRtiling\';




% Define the start and end frame indices
startFrame = 1;
endFrame = 49;
mae_values_AT = zeros(1, endFrame);
mae_values_AT0 = zeros(1, endFrame);

blockSize = [128 , 128];


for i = startFrame:endFrame
    fprintf('I am processing the previous frame %d\n', i);
    % Load the previous image
    image_file1 = sprintf('frame_%04d.png', i);
    In_pre = double(imread([image_path image_file1])) / 255;
  %  In_pre = double(imread([image_path image_file1])) ;

    % Load the current image
    image_file2 = sprintf('frame_%04d.png', i+1);
    I_n_curr = double(imread([image_path image_file2])) / 255;
  %   I_n_curr = double(imread([image_path image_file2])) ;   

    [rows, cols, ~] = size(In_pre);
    x = 1:cols;
    y = 1:rows;
    %SSS = (blockSize(1)*blockSize(2))/rows *cols;



    figure(1);
    image(x, y, In_pre);
    axis image;
    title('previous frame');
    hold on;
    figure(2);
    image(x, y, I_n_curr);
    axis image;
    title('current frame');
    hold on;
    shg;

    % Read flow data from folder_bilinear_interpolate
    file_BIup = sprintf('%s\\flow_000%02d_000%02d_BIupsampled.flo', folder_BIup, i, i+1);
    %d = readFlowFile(file_BIup);


    file_BCup = sprintf('%s\\flow_000%02d_000%02d_BCupsampled.flo', folder_BCup, i, i+1);
    %d = readFlowFile(file_BCup);
%     d0 = readFlowFile(file_BCup);
%     d = readFlowFile(file_BCup);

    % Read flow data from the folder_upsample_twice
 %   file_up = sprintf('%s\\flow_000%02d_000%02d_upsampled.flo', folder_up, i, i+1);
    file_upu = sprintf('%s\\flow_000%02d_000%02d_upsampled.flo', u2, i, i+1);
   % flow_up = readFlowFile(file_up);
      d0 = readFlowFile(file_upu);
      d = readFlowFile(file_upu);


    % Read flow data from folder_tiling
    file_tiling = sprintf('%s\\flow_000%02d_000%02d.flo', folder_tiling, i, i+1);
   
    file_tiling2 = sprintf('%s\\flow_000%02d_000%02d.flo', f2, i, i+1);
%      d0 = readFlowFile(file_tiling2);
%      d = readFlowFile(file_tiling2);





    % Read flow data from the second folder_gt
    file_gt = sprintf('%s\\frame_%04d.flo', folder_gt, i);
    flow_gt = readFlowFile(file_gt);


        I_n_prev_gray = 255 * double(im2gray(In_pre));
        I_n_curr_gray = 255 * double(im2gray(I_n_curr));
%        I_n_prev_gray = double(rgb2gray(In_pre));
%        I_n_curr_gray = double(rgb2gray(I_n_curr));

   % I_n_pred = warp_image_using_interp2(I_n_prev_gray, d);
    %E = I_n_curr_gray - I_n_pred;

    % Iterate through each block in the image
    max_iters = 5;
    par_hist = zeros(max_iters, 4);

% for iteration = 1:max_iters
     for row = 1 +blockSize(1) :blockSize(1) :  rows - blockSize(1)
         for col = 1 + blockSize(2) : blockSize(2) :      cols  - blockSize(2)

            % Define block limits
            rowEnd = min(row+blockSize(1)-1, rows);
            colEnd = min(col+blockSize(2)-1, cols);

            block_prev = I_n_prev_gray(row:rowEnd, col:colEnd);
           % block_prev =   I_n_pred(row:rowEnd, col:colEnd);
            block_curr = I_n_curr_gray(row:rowEnd, col:colEnd);
            bx = col : colEnd;
            by = row :rowEnd;

            %gt = flow_gt(row:rowEnd, col:colEnd,:);
            d_block =  d(row:rowEnd, col:colEnd,:);
            d0_block = d0(row:rowEnd, col:colEnd,:);
            gt_block = flow_gt(row:rowEnd, col:colEnd,:);


            figure(3);subplot(1, 2, 1); imshow(I_n_prev_gray);
            title('Previous Image');
            hold on;
            rectangle('Position', [col, row, blockSize(2), blockSize(1)], 'EdgeColor', 'r', 'LineWidth', 2);
            hold off;
            subplot(1, 2, 2);imshow(I_n_curr_gray); title('Current Frame');   
            hold on;
            rectangle('Position', [col, row, blockSize(2), blockSize(1)], 'EdgeColor', 'r', 'LineWidth', 2);
            hold off;
            sgtitle('watching the process');

           % d_block(:, :, 1)=0; d_block(:, :, 2)=0;
            % Iteration
            for iteration = 1:max_iters
                block_prev = fetch_block_from_frame_using_interp2(I_n_prev_gray, bx, by, d_block);
                [G, E, update] = calcgradest(block_prev, block_curr);
                %d_block = zeros(size(block_prev));
                d_block(:, :, 1) = d_block(:, :, 1) - update(1) ;
                d_block(:, :, 2) = d_block(:, :, 2) - update(2);

                par_hist(iteration, :) = [d_block(blockSize(1)/2, blockSize(2)/2), update(1), update(2),mean(abs(E(:)))];
           %     block_prev = fetch_block_from_frame_using_interp2(I_n_prev_gray, bx, by, d_block);
               % block_prev = fetch_block_from_frame_using_interp2(I_n_pred, bx, by, d_block)

             % Update the block motion vector to the whole image motion vector
%                 d(row:rowEnd, col:colEnd, 1) = d(row:rowEnd, col:colEnd, 1) + d_block(1);
%                 d(row:rowEnd, col:colEnd, 2) = d(row:rowEnd, col:colEnd, 2) + d_block(2);
               d(row:rowEnd, col:colEnd, 1) =  d_block(:, : ,1);
               d(row:rowEnd, col:colEnd, 2) =  d_block(:, : ,2);


           end 

        end

    end
    % Save the refined flow field
    outrefine_flow_file = sprintf('%s\\refined_flow_%04d.flo', outrefine_flow_path, i);
    writeFlowFile(d , outrefine_flow_file);
    disp(i);        % should display a frame number
    
    % MAE for flow_raft and flow_gt
    abs_error_AT = abs(d - flow_gt);
    mae_AT = mean(abs_error_AT(:));
    mae_values_AT(i) = mae_AT;

    abs_error_AT0 = abs(d0 - flow_gt);
    mae_AT0 = mean(abs_error_AT0(:));
    mae_values_AT0(i) = mae_AT0;


end

average_mae_d_original = mean(mae_values_AT0);
average_mae_d_updated = mean(mae_values_AT);

fprintf('Average MAE for original estimated motion field: %f\n', average_mae_d_original);
fprintf('Average MAE for updated motion field: %f\n', average_mae_d_updated);



figure(4);
plot((1 : max_iters), par_hist(:, 1), 'r-x', 'linewidth' , 1.5);
title('Covergence of motion');
legend('displacement');
xlabel('Iteration');
ylabel('Motion');
shg;

figure(5);
plot((1 : max_iters), par_hist(:, 4), 'r-x', 'linewidth' , 1.5);
title('Covergence of Motion compensated MAE');
xlabel('Iteration');
ylabel('vector difference MAE');
shg;

figure(6);
plot(1:endFrame, mae_values_AT, 'o-',  'LineWidth', 2);
hold on;
plot(1:endFrame, mae_values_AT0, '+-', 'LineWidth', 2);
xlabel('frame numbers');
ylabel('Mean Absolute Error (MAE)');
legend('MAE with refinement strategy', 'Original MAE')
%legend('UP', 'UP0')
%legend('BIUP', 'BIUP0')

title('vector difference error');
shg;



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


 function [G, E, update] = calcgradest(block_prev, block_curr)
     [gx, gy] = gradient(block_prev);
     E = block_curr - block_prev;
     G = [gx(:) gy(:)]';
     GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];
     update = -pinv(GtG) * G * E(:);

 end
