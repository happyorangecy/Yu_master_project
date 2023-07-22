close all;
clear;

image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
folder_BIup = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBI\';
% Define output path for refined flow field
outrefine_flow_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\newR\';

% Define the start and end frame indices
startFrame = 1;
endFrame = 50;
blockSize = [32 32];

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
    d = readFlowFile(file_BIup);

    I_n_prev_gray = double(im2gray(In_pre));
    I_n_curr_gray = double(im2gray(I_n_curr));
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

            % Iteration
            for iteration = 1:max_iters
                %block_prev = interp2(x, y, double(I_n_prev_gray), x + dx, y + dy);
                block_prev = warp_image_using_interp2(block_prev , d(row:rowEnd, col:colEnd, :));
                [G, E, update] = calcgradest(block_prev, block_curr);
                d(1) = d(1) - update(1);
                d(2) = d(2) - update(2);
                par_hist(iteration, :) = [d(1, 1, 1), d(1, 1, 2), update(1), update(2), mean(mean(abs(E(:))))];

                % Update the block motion vector to the whole image motion vector
                d(row:rowEnd, col:colEnd, 1) = d(row:rowEnd, col:colEnd, 1) + d(1);
                d(row:rowEnd, col:colEnd, 2) = d(row:rowEnd, col:colEnd, 2) + d(2);
            end
        end
    end
    % Save the refined flow field
    outrefine_flow_file = sprintf('%s\\refined_flow_%04d.flo', outrefine_flow_path, i);
    writeFlowFile(d, outrefine_flow_file);
    disp(outrefine_flow_path); % should display a valid directory path as a one-row string
    disp(i);        % should display a frame number
    disp(outrefine_flow_file); % should display the full path to the output .flo file as a one-row string


end

figure(3);
plot((1 : max_iters), par_hist(:, 1), 'r-x', (1 : max_iters), par_hist(:, 2), 'b-x', 'linewidth' , 1.5);
title('Covergence of motion');
legend('DX', 'DY');
xlabel('Iteration');
ylabel('Motion');
shg;

figure(4);
plot((1 : max_iters), par_hist(:, 5), 'r-x', 'linewidth' , 1.5);
title('Covergence of Motion compensated MAE');
xlabel('Iteration');
ylabel('frame difference MAE');
shg;




function [G, E, update] = calcgradest(block_prev, block_curr)
    [gx, gy] = gradient(block_prev);
    E = block_curr - block_prev;
    G = [gx(:) gy(:)]';
    GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];
    update = -pinv(GtG) * G * E(:);
end


