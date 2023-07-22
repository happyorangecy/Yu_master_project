close all
clear all

%% render the vectors to arrows
im = imread('frame_0002.png');
% downsample the image
%im_downsampled = imresize(im, [height, width]); 
imb = imread('frame_0001.png');
imf = imread('frame_0003.png');
% read forward optical flow
img_flo_f = readFlowFile('flow_00002_00003.flo');
out_f = flowToColor(img_flo_f);

% read backward optical flow
img_flo_b = readFlowFile('flow_00002_00001.flo');
out_b = flowToColor(img_flo_b);

% visualize forward optical flow
% figure(1)
% imshow(out_f)
% title('Forward Optical Flow')

% visualize backward optical flow
% figure(2)
% imshow(out_b);
% title('Backward Optical Flow')

% forward flow vectors
u_f = img_flo_f(:, :, 1);
v_f = img_flo_f(:, :, 2);

% backward flow vectors
u_b = img_flo_b(:, :, 1);
v_b = img_flo_b(:, :, 2);

% Define a meshgrid for the image
[xx, yy] = meshgrid(1:size(im, 2), 1:size(im, 1)); % col and row

% downsample the meshgrid to plot arrows
xq = xx(1:10:end, 1:10:end);
yq = yy(1:10:end, 1:10:end);
uq_f = u_f(1:10:end, 1:10:end);
vq_f = v_f(1:10:end, 1:10:end);
uq_b = u_b(1:10:end, 1:10:end);
vq_b = v_b(1:10:end, 1:10:end);

% Visualize the optical flow vectors as arrows on the image frame
figure(1); imshow(im);
hold on; 
quiver(xq, yq, uq_f, vq_f,0, 'color', 'g', 'LineWidth', 1.5);
quiver(xq, yq, uq_b, vq_b,0, 'color', 'r', 'LineWidth', 1.5);
title('Optical Flow foreward and backward Vectors')

figure (5); imshow(imb);
title('frame B');

figure(6); imshow(imf);
title('frame F');

%% MAE
% Initialize variables
estflopath = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\estimation_flo\';
num_files = 49;
mae_values = zeros(num_files, 1);
mse_values = zeros(num_files, 1);
frame_numbers = 1:49;


figure(100)
% Loop through the desired frames
for i = 1:num_files
    % Read the estimated forward .flo file
    festflo_file = sprintf('flow_%05d_%05d.flo', i, i+1);
    festflo_forward = readFlowFile([estflopath festflo_file]);

    % Read the estimated backward .flo file
    bestflo_file = sprintf('flow_%05d_%05d.flo', i+1, i);
    bestflo_backward = readFlowFile([estflopath bestflo_file]);

    % Read the corresponding ground truth .flo file
    truth_file = sprintf('frame_%04d.flo', i); %  get the corresponding frame number
    truth_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1\';
    truth = readFlowFile([truth_path truth_file]);

    %% computation MAE between F and gt
    %%
    mean_abs_diff = mean(abs(festflo_forward(:) - truth(:)));
    mae_values1(i) = mean_abs_diff;
    fprintf('Mean Absolute Error for F_%d and gt_%d: %f\n', i,i, mean_abs_diff);

    %% Calculate MAE with motion compensation
    %%  F
    image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
    F_image_file1 = sprintf('frame_%04d.png', i);
    F_image_file2 = sprintf('frame_%04d.png', i+1);
    F_frame_n = double(imread([image_path F_image_file1])) / 255;
    frame_np1 = double(imread([image_path F_image_file2])) / 255;
    % Warp frame n+1 using forward flow vectors
    I_np1_xd = warp_image_using_interp2(frame_np1, festflo_forward); 
    % Calculate the motion compensation MAE with 1/(NM) factor
    F_abs_diff = abs(F_frame_n - I_np1_xd); %plot vis
    %fprintf('MC Mean Absolute Error F_%d to F_%d : %f\n', i,i+1,  mean(F_abs_diff(:)));
     image(F_abs_diff);
     drawnow
    N = size(F_abs_diff, 1);
    M = size(F_abs_diff, 2);
    F_MCMAE(i) = sum(F_abs_diff(:)) / (N * M);
    fprintf('With MC Mean Absolute Error F %d to %d: %f\n', i,i+1, F_MCMAE(i));


    % B
    B_image_file1 = sprintf('frame_%04d.png', i+1);
    B_image_file2 = sprintf('frame_%04d.png', i);
    B_frame_n = double(imread([image_path B_image_file1])) / 255;
    frame_nm1 = double(imread([image_path B_image_file2])) / 255;
    % Warp frame n-1 using forward flow vectors
    I_nm1_xd = warp_image_using_interp2(frame_nm1, bestflo_backward); 
    % Calculate the motion compensation MAE with 1/(NM) factor
    B_abs_diff = abs(B_frame_n - I_nm1_xd);
    N = size(B_abs_diff, 1);
    M = size(B_abs_diff, 2);
    B_MCMAE(i) = sum(B_abs_diff(:)) / (N * M);
    %fprintf('With MC Mean Absolute Error B %d to %d: %f\n', i,i-1, B_MCMAE(i));


    % no compensated
    no_abs_diff = abs(F_frame_n -frame_np1 );
    no_abs_diff_mae(i) = sum(no_abs_diff(:)) / (N * M);

%gt with compensation
I_np1_gt = warp_image_using_interp2(frame_np1, truth); 
    % Calculate the motion compensation MAE with 1/(NM) factor
    gt_abs_diff = abs(F_frame_n - I_np1_gt); %plot vis
    N = size(gt_abs_diff, 1);
    M = size(gt_abs_diff, 2);
    gt_MCMAE(i) = sum(gt_abs_diff(:)) / (N * M);
    fprintf('With MC Mean Absolute Error F %d to %d: %f\n', i,i+1, gt_MCMAE(i));


end

% Calculate the MAE between F and gt
mae_avg1 = mean(mae_values1);
disp(['Average MAE between F and gt over frames ' num2str(frame_numbers(1)) '-' num2str(frame_numbers(end)) ': ' num2str(mae_avg1)]);
 figure(2);    
 plot(1:num_files, mae_values1, 'o-');
 xlabel('frame number');
 ylabel('Mean Absolute Error');
 title('Mean Absolute Error for Foreground and Ground Truth .flo Files');
 grid on;


 % Display the motion compensation MAE values
 figure(3);
 % Plot the F_MCMAE values
 yyaxis left; % Use the left y-axis for F_MCMAE
 plot(frame_numbers, F_MCMAE, 'LineWidth', 1);
 xlabel('Frame Number');
 ylabel('Motion Compensation MAE F');
 grid on;
 hold on;
 % Plot the B_MCMAE values
 yyaxis right; % Use the right y-axis for B_MCMAE
 plot(frame_numbers, B_MCMAE, 'LineWidth', 1);
 ylabel('Motion Compensation MAE B'); 
 title('Motion Compensation MAE F and B vs. Frame Number');
 grid on;




% Calculate the average MAE and MSE over all the frames
FMCmae_avg = mean(F_MCMAE);
BMCmae_avg = mean(B_MCMAE);
% Display the result
disp(['Average F_MCMAE ' num2str(frame_numbers(1)) '-' num2str(frame_numbers(end)) ': ' num2str(FMCmae_avg)]);
disp(['Average B_MCMAE ' num2str(frame_numbers(1)) '-' num2str(frame_numbers(end)) ': ' num2str(BMCmae_avg)]);


figure(4);
%plot(1:num_files, mae_values1, 'o-', 'LineWidth', 2);
plot(frame_numbers, F_MCMAE, 'LineWidth', 1);

xlabel('Frame Number');
ylabel('Mean Absolute Error');
grid on;
hold on;
%plot(frame_numbers, F_MCMAE, 'LineWidth', 1);
plot(frame_numbers, B_MCMAE, 'LineWidth', 1);
plot(frame_numbers, no_abs_diff_mae, 'LineWidth', 1);
plot(frame_numbers, gt_MCMAE, 'LineWidth', 1);

title('Mean Absolute Error and Motion Compensation MAE');
%legend('Ground Truth', 'Motion Compensation F', 'Motion Compensation B', 'non compensation','gt MC');
legend( 'Motion Compensation F', 'Motion Compensation B', 'non compensation','gt MC');
hold off;









