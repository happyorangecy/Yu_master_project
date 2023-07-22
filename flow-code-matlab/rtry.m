close all;
clear;

image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
folder_BIup = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBI\';
outrefine_flow_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\newRBI\';

startFrame = 1;
endFrame = 5;
blockSize = [128 128];

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

    % Visualization of the previous and current frames
    figure(1);
    image(x, y, In_pre);
    axis image;
    title('Previous Frame');
    hold on;
    figure(2);
    image(x, y, I_n_curr);
    axis image;
    title('Current Frame');
    hold on;
    shg;

    % Read flow data from folder_BIup
    file_BIup = sprintf('%s\\flow_000%02d_000%02d_BIupsampled.flo', folder_BIup, i, i+1);
    d = readFlowFile(file_BIup);

    % Convert frames to grayscale
    I_n_prev_gray = double(rgb2gray(In_pre));
    I_n_curr_gray = double(rgb2gray(I_n_curr));

    % Initialize motion compensation
    I_n_pred = warp_image_using_interp2(I_n_prev_gray, d);
    E = I_n_curr_gray - I_n_pred;

    % Iterate through each block in the image
    max_iters = 5;
    par_hist = zeros(max_iters, 5);

    for row = 1:blockSize(1):rows
        for col = 1:blockSize(2):cols
            % Define block limits
            rowEnd = min(row+blockSize(1)-1, rows);
            colEnd = min(col+blockSize(2)-1, cols);

            block_prev = I_n_prev_gray(row:rowEnd, col:colEnd);
            block_curr = I_n_curr_gray(row:rowEnd, col:colEnd);

            d_block = d(row:rowEnd, col:colEnd, :);

            % Visualization of the block
            figure(3);
            subplot(1, 2, 1);
            imshow(block_prev);
            title('Block Previous');
            subplot(1, 2, 2);
            imshow(block_curr);
            title('Block Current');
            sgtitle('Block Visualization');
            %pause(0.5);

            % Iteration
            for iteration = 1:max_iters
                block_prev_warped = warp_image_using_interp2(block_prev, d_block);
                [G, E, update] = calcgradest(block_prev_warped, block_curr);

                d_block(1) = d_block(1) - update(1);
                d_block(2) = d_block(2) - update(2);

                par_hist(iteration, :) = [d_block(1), d_block(2), update(1), update(2), mean(abs(E(:)))];

                block_prev = block_prev_warped;

                % Update the block motion vectors in the flow field
                d(row:rowEnd, col:colEnd, 1) = d_block(1);
                d(row:rowEnd, col:colEnd, 2) = d_block(2);
            end
        end
    end

    % Save the refined flow field
    outrefine_flow_file = sprintf('%s\\refined_flow_%04d.flo', outrefine_flow_path, i);
    writeFlowFile(d, outrefine_flow_file);
    disp(outrefine_flow_path);
    disp(i);
    disp(outrefine_flow_file);
end

figure(4);
plot((1:max_iters), par_hist(:, 1), 'r-x', (1:max_iters), par_hist(:, 2), 'b-x', 'linewidth', 1.5);
title('Convergence of Motion');
legend('DX', 'DY');
xlabel('Iteration');
ylabel('Motion');
shg;

figure(5);
plot((1:max_iters), par_hist(:, 5), 'r-x', 'linewidth', 1.5);
title('Convergence of Motion Compensated MAE');
xlabel('Iteration');
ylabel('Frame Difference MAE');
shg;

function [G, E, update] = calcgradest(block_prev, block_curr)
    [gx, gy] = gradient(block_prev);
    E = block_curr - block_prev;
    G = [gx(:) gy(:)]';
    GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];
    update = -inv(GtG) * G * E(:);
end
