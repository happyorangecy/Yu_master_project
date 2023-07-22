close all;
clear;

image_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
folder_BIup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBI\';
folder_BCup = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\uptwiceBC\';
folder_tiling = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\after_tiling';


% Define output path for refined flow field
%outrefine_flow_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\newRBI\';
%outrefine_flow_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\newRBC\';
outrefine_flow_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\newRtiling\';




% Define the start and end frame indices
startFrame = 1;
endFrame = 5;
blockSize = [64, 64];
max_iters = 5;
par_hist = zeros(max_iters, 5);
for iteration = 1:max_iters
    for i = startFrame:endFrame-1
    fprintf('I am processing the previous frame %d\n', i);
    % Load the previous image
    image_file1 = sprintf('frame_%04d.png', i);
    In_pre = double(imread([image_path image_file1])) / 255;
    % Load the current image
    image_file2 = sprintf('frame_%04d.png', i+1);
    I_n_curr = double(imread([image_path image_file2])) / 255;

    [rows, cols, ~] = size(In_pre);
    x = 1:cols;
    y = 1:rows;
    SSS = (blockSize(1)*blockSize(2))/rows *cols;



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

    % Read flow data from folder_tiling
    file_tiling = sprintf('%s\\flow_000%02d_000%02d.flo', folder_tiling, i, i+1);
    d = readFlowFile(file_tiling);

     I_n_prev_gray = 255 * double(im2gray(In_pre));
     I_n_curr_gray = 255 * double(im2gray(I_n_curr));
%     I_n_prev_gray = double(rgb2gray(In_pre));
%     I_n_curr_gray = double(rgb2gray(I_n_curr));

    I_n_pred = warp_image_using_interp2(I_n_prev_gray, d);
    E = I_n_curr_gray - I_n_pred;

    % Iterate through each block in the image


%         for row = 1 + blockSize(1) : 20 : rows - blockSize(1)
%             for col = 1 + blockSize(2) : 20 : cols  - blockSize(2)

        for row = 1 + blockSize(1) : 5 : 10 + blockSize(1)
            for col = 1 + blockSize(2) : 5 : 10  + blockSize(2)
            % Define block limits
%             rowEnd = min(row+blockSize(1)-1, rows);
%             colEnd = min(col+blockSize(2)-1, cols);
            rowEnd = row-blockSize(1);
            colEnd = col-blockSize(2);

           % block_prev = I_n_prev_gray(row:rowEnd, col:colEnd);
            block_prev = I_n_pred(row:rowEnd, col:colEnd);
            block_curr = I_n_curr_gray(row:rowEnd, col:colEnd);
            bx = col : colEnd;
            by = row :rowEnd;

            %d_block = zeros(size(block_prev));
            d_block =  d(row:rowEnd, col:colEnd,:);

%             % Draw and display block_prev and block_curr
%             figure(3);
%             subplot(1,2,1);
%             imshow(block_prev);
%             title('Block Previous');
%             subplot(1,2,2);
%             imshow(block_curr);
%             title('Block Current');
%             %pause(0.5);

                figure(3);
                subplot(1, 2, 1);
                imshow(In_pre);
                title('Previous Image');
                hold on;
                rectangle('Position', [col, row, blockSize(2), blockSize(1)], 'EdgeColor', 'r', 'LineWidth', 2);
                hold off;
                subplot(1, 2, 2);
                imshow(I_n_curr);
                title('Current Frame');   
                hold on;
                rectangle('Position', [col, row, blockSize(2), blockSize(1)], 'EdgeColor', 'r', 'LineWidth', 2);
                hold off;
                sgtitle('watching the process');


            % Iteration
        %    for iteration = 1:max_iters
                % block_prev is updated with the warped version from previous iteration
                block_prev = fetch_block_from_frame_using_interp2(I_n_prev_gray, bx, by, d_block);
                
                [update] = calcgradest(block_prev, block_curr);
                
                % d_block is updated with the sum of current and previous values
                d_block(:, :, 1) = d_block(:, :, 1) - update(1);
                d_block(:, :, 2) = d_block(:, :, 2) - update(2);
                
            %    par_hist(iteration, :) = [d_block(blockSize(1)/2, blockSize(2)/2), d_block(blockSize(1)/2, blockSize(2)/2), update(1), update(2),mean(abs(E(:)))];
            end
            
            d(row:rowEnd, col:colEnd, 1) = d_block(:, : ,1) ;
            d(row:rowEnd, col:colEnd, 2) = d_block(:, : ,2) ;

        end
    end
    % Save the refined flow field
    outrefine_flow_file = sprintf('%s\\refined_flow_%04d.flo', outrefine_flow_path, i);
    writeFlowFile(d , outrefine_flow_file);
    disp(i);        % should display a frame number



end

figure(4);
plot((1 : max_iters), par_hist(:, 1), 'r-x', (1 : max_iters), par_hist(:, 2), 'b-x', 'linewidth' , 1.5);
title('Covergence of motion');
legend('DX', 'DY');
xlabel('Iteration');
ylabel('Motion');
shg;

figure(5);
plot((1 : max_iters), par_hist(:, 5), 'r-x', 'linewidth' , 1.5);
title('Covergence of Motion compensated MAE');
xlabel('Iteration');
ylabel('frame difference MAE');
shg;




% function [G, E, update] = calcgradest(block_prev, block_curr)
%     [gx, gy] = gradient(block_prev);
%     E = block_curr - block_prev;
%     G = [gx(:) gy(:)]';
%     GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];
%     update = -pinv(GtG) * G * E(:);
%     [x,y] = size(block_prev);
%     update = update /(x*y);
% end

function [update] = calcgradest(block_prev, block_curr)
    [gx, gy] = gradient(block_prev);
    E = block_curr - block_prev;
    G = [gx(:) gy(:)]';
    GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];
    update = -pinv(GtG) * G * E(:);
    [x,y] = size(block_prev);
    update = update /(x*y);
end

