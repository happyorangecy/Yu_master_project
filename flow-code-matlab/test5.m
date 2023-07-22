close all
clear all

% test if it has been tiled from 4 parts into a whole one
 flo_filepath = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\after_tiling\flow_00025_00026.flo';
 img_flo_f = readFlowFile(flo_filepath);
 out_f = flowToColor(img_flo_f);
 figure(1)
 imshow(out_f);
%%
raft_filepath = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\estimation_flo\flow_00025_00026.flo';
img_flo_r = readFlowFile(raft_filepath);
out_r = flowToColor(img_flo_r);
figure(2)
imshow(out_r);
 %% gt
truth_path = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1\frame_0025.flo';
               % C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1
img_flo_gt = readFlowFile(truth_path);
out_gt = flowToColor(img_flo_gt);
 figure(3)
 imshow(out_gt)

 %%
close all
clear all

im_n = imread('frame_0025.png');
im_np1 = imread('frame_0026.png');
[xx, yy] = meshgrid(1:size(im_n, 2), 1:size(im_n, 1)); % col and row

 flo_filepath = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\after_tiling\flow_00025_00026.flo';
 img_flo_f = readFlowFile(flo_filepath);
 out_f = flowToColor(img_flo_f);
 u_f = img_flo_f(:, :, 1);
 v_f = img_flo_f(:, :, 2);

xq = xx(1:10:end, 1:10:end);
yq = yy(1:10:end, 1:10:end);
uq_f = u_f(1:10:end, 1:10:end);
vq_f = v_f(1:10:end, 1:10:end);
figure(4); imshow(im_n);title('frame n_raft-tiling');
hold on; 
quiver(xq, yq, uq_f, vq_f,0, 'color', 'g', 'LineWidth', 1.5);
quiver(xq, yq, uq_f, zeros(size(vq_f)),0, 'color', 'r', 'LineWidth', 1.5);
quiver(xq, yq, zeros(size(uq_f)), vq_f, 0, 'color', 'b', 'LineWidth', 1.5);




figure(5); imshow(im_np1);title('frame n+1');


raft_filepath = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\estimation_flo\flow_00025_00026.flo';
img_flo_r = readFlowFile(raft_filepath);
out_r = flowToColor(img_flo_r);
u_r = img_flo_r(:, :, 1); v_r = img_flo_r(:, :, 2);
xq = xx(1:10:end, 1:10:end); yq = yy(1:10:end, 1:10:end);
uq_r = u_r(1:10:end, 1:10:end); vq_r = v_r(1:10:end, 1:10:end);
figure(6); imshow(im_n);title('frame n_raft-only');
hold on; 
quiver(xq, yq, uq_r, vq_r,0, 'color', 'm', 'LineWidth', 1.5);
quiver(xq, yq, uq_r, zeros(size(vq_r)),0, 'color', 'r', 'LineWidth', 1.5);
quiver(xq, yq, zeros(size(uq_r)), vq_r, 0, 'color', 'b', 'LineWidth', 1.5);


gt_filepath = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\alley_1\\frame_0025.flo';
img_flo_gt = readFlowFile(gt_filepath);
out_gt = flowToColor(img_flo_gt);
u_gt = img_flo_gt(:, :, 1); v_gt = img_flo_gt(:, :, 2);
xq = xx(1:10:end, 1:10:end); yq = yy(1:10:end, 1:10:end);
uq_gt = u_gt(1:10:end, 1:10:end); vq_gt = v_gt(1:10:end, 1:10:end);
figure(7); imshow(im_n);title('frame gt');
hold on; 
quiver(xq, yq, uq_gt, vq_gt,0, 'color', 'y', 'LineWidth', 1.5);
quiver(xq, yq, uq_gt, zeros(size(vq_gt)),0, 'color', 'r', 'LineWidth', 1.5);
quiver(xq, yq, zeros(size(uq_gt)), vq_gt, 0, 'color', 'b', 'LineWidth', 1.5);




% Calculate the new coordinates after applying the flow
%xx_new = xx + u_f;
%yy_new = yy + v_f;
% Calculate the displacement vectors
%dxx = xx_new - xx;
%dyy = yy_new - yy;
 dxx =u_f ;
 dyy =v_f;
% Select every 10th element for visualization purposes
xq = xx(1:10:end, 1:10:end);
yq = yy(1:10:end, 1:10:end);
dxx_q = dxx(1:10:end, 1:10:end);
dyy_q = dyy(1:10:end, 1:10:end);
% Plot the image and quiver only the dx displacement
figure(8); imshow(im_n); title('frame n');
hold on;
quiver(xq, yq, dxx_q, zeros(size(dxx_q)), 0, 'color', 'g', 'LineWidth', 1.5);

%dy
figure(9); imshow(im_n); title('frame n');
hold on;
quiver(xq, yq,  zeros(size(dyy_q)),dyy_q, 0, 'color', 'r', 'LineWidth', 1.5);