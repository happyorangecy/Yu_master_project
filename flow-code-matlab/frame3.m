close all
clear all

%% render the vectors to arrows
%frame n-1
im_nm1= imread('frame_0001.png');
figure(1)
subplot(3,1,1)
imshow(im_nm1);
%frame n
im_n = imread('frame_0002.png');
figure(1)
subplot(3,1,2)
imshow(im_n);
%frame n+1
im_np1 = imread('frame_0003.png');
figure(1)
subplot(3,1,3)
imshow(im_np1);

img_flo_f = readFlowFile('flow_00002_00003.flo');
out_f = flowToColor(img_flo_f);

img_flo_b = readFlowFile('flow_00002_00003.flo');
out_b = flowToColor(img_flo_b);

% %flow vectors
% u = img_flo(:, :, 1);
% v = img_flo(:, :, 2);
% 
% % Define a meshgrid for the image
% [xx, yy] = meshgrid(1:size(im_n, 2), 1:size(im_n, 1)); % col and row
% 
% % downsample the meshgrid to plot arrows
% xq = xx(1:10:end, 1:10:end);
% yq = yy(1:10:end, 1:10:end);
% uq = u(1:10:end, 1:10:end);
% vq = v(1:10:end, 1:10:end);
% 
% % Visualize the optical flow vectors as arrows on the image frame
% figure(2)
% %figure(2);subplot(3,1,1); 
% imshow(im_nm1);title('quiver on frame n-1');
% hold on; quiver(xq, yq, uq, vq, 'color', 'g', 'LineWidth', 1.5);
% 
% figure(3)
% %figure(2);subplot(3,1,2); 
% imshow(im_n);title('quiver on frame n');
% hold on; quiver(xq, yq, uq, vq, 'color', 'g', 'LineWidth', 1.5);
% 
% figure(4)
% %figure(2);subplot(3,1,3); 
% imshow(im_np1);title('quiver on frame n+1');
% hold on; quiver(xq, yq, uq, vq, 'color', 'g', 'LineWidth', 1.5);

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
