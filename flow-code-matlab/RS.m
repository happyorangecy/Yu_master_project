close all;
clear;

image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
folder_gt = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1';
folder_BIup = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBI2\';

folder_raft = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\estimation_flo\';
folder_tiling = 'C:\\Users\\23839\\Desktop\\RAFT-master\\RAFT-master\\after_tiling';
folder_up = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwice2\';
folder_BCup = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBC2\';

folders = {folder_BIup, folder_raft, folder_tiling, folder_up, folder_BCup};
labels = {'BIup', 'raft', 'tiling', 'up', 'BCup'};

% Define the start and end frame indices
startFrame = 1;
endFrame = 5;

for f = 1:numel(folders)
    folder = folders{f};
    label = labels{f};
    for i = startFrame:endFrame-1
        fprintf('I am processing the previous frame %d\n', i);
        %disp(i);
        % Load the previous image
        image_file1 = sprintf('frame_%04d.png', i);
        In_pre = double(imread([image_path image_file1])) / 255;
        % Load the current image
        image_file2 = sprintf('frame_%04d.png', i+1);
        I_n_curr = double(imread([image_path image_file2])) / 255;
    
        [rows, cols, ~] = size(In_pre);
        x = 1:cols;
        y = 1:rows;
    %     X = ones(rows, 1) * x;
    %     Y = (1:rows)' * ones(1, cols);
    
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
    
    
        % Read flow data from current folder
        file_flow = sprintf('%s\\flow_000%02d_000%02d.flo', folder, i, i+1);
        d = readFlowFile(file_flow);
    %     % Read flow data from the second folder_gt
    %     file_gt = sprintf('%s\\frame_%04d.flo', folder_gt, i);
    %     flow_gt = readFlowFile(file_gt);
        
        I_n_prev_gray = double(im2gray(In_pre));
        I_n_curr_gray = double(im2gray(I_n_curr));

        I_n_pred = warp_image_using_interp2(I_n_prev_gray, d);
        E = I_n_curr_gray - I_n_pred;
    
        % Iteration
        max_iters = 5;
        par_hist = zeros(max_iters, 5);
        correlation = d;  % Initialize correlation with initial motion vectors
        for iteration = 1:max_iters
            [G, E, update] = calcgradest(I_n_pred, I_n_curr_gray);
            correlation(1) = correlation(1) - update(1); % Update motion vectors
            correlation(2) = correlation(2) - update(2);
            par_hist(iteration, :) = [d(1, 1, 1), d(1, 1, 2), correlation(1), correlation(2), mean(E(:).^2)];
            I_n_preu = warp_image_using_interp2(I_n_curr_gray, correlation);
            I_n_pred = I_n_preu;
        end
    
        % Update the estimated motion vectors
        d = correlation;
        par_hist_all(:,:,f) = par_hist;
    end
end 

figure(3);
hold on;
for f = 1:numel(folders)
    plot((1 : max_iters), par_hist_all(:, 1, f), 'x-', 'linewidth' , 1);
end
title('Covergence of motion');
legend(labels);
xlabel('Iteration');
ylabel('Motion');
shg;

figure(4);
hold on;
for f = 1:numel(folders)
    plot((1 : max_iters), par_hist_all(:, 2, f), 'x-', 'linewidth' , 1);
end
title('Covergence of motion');
legend(labels);
xlabel('Iteration');
ylabel('Motion');
shg;


figure(5);
hold on;
for f = 1:numel(folders)
    plot((1 : max_iters), par_hist_all(:, 5, f), 'x-', 'linewidth' , 1.5);
end
title('Covergence of Motion compensated MSE');
legend(labels);
xlabel('Iteration');
ylabel('frame difference MSE');
shg;





function [G, E, update] = calcgradest(block_prev, block_curr)
    [gx, gy] = gradient(block_prev);
    E = block_curr - block_prev;
    G = [gx(:) gy(:)]';
    GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];
    update = -pinv(GtG) * G * E(:);
end
