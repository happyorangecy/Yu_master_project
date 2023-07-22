close all
clear;

im_n = imread('frame_0025.png');
im_np1 = imread('frame_0026.png');
[xx, yy] = meshgrid(1:size(im_n, 2), 1:size(im_n, 1)); % col and row

flo_filepath = 'C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\after_tiling\flow_00025_00026.flo';
img_flo_f = readFlowFile(flo_filepath);
out_f = flowToColor(img_flo_f);
u_f = img_flo_f(:, :, 1);
v_f = img_flo_f(:, :, 2);

xq = xx(380:396, 516:532);
yq = yy(380:396, 516:532);
uq_f = u_f(380:396, 516:532);
vq_f = v_f(380:396, 516:532);
figure(4); imshow(im_n);title('after-tiling frame 25');
hold on; 
quiver(xq, yq, uq_f, vq_f,0, 'color', 'g', 'LineWidth', 1.5);
axis([516 532 380 396]);

raft_filepath = 'C:\Users\Chengyu\Desktop\R64\flow_00025_00026.flo';

img_flo_r = readFlowFile(raft_filepath);
out_r = flowToColor(img_flo_r);
u_r = img_flo_r(:, :, 1); v_r = img_flo_r(:, :, 2);
xq = xx(380:396, 516:532); yq = yy(380:396, 516:532);
uq_r = u_r(380:396, 516:532); vq_r = v_r(380:396, 516:532);
figure(6); imshow(im_n);title('after-tiling frame 25 with refinement');
hold on; 
quiver(xq, yq, uq_r, vq_r,0, 'color', 'm', 'LineWidth', 1.5);
axis([516 532 380 396]);

gt_filepath = 'C:\\Users\\Chengyu\\Desktop\\RAFT-master\\RAFT-master\\datasets\\Sintel\\training\\flow\\alley_1\\frame_0025.flo';
img_flo_gt = readFlowFile(gt_filepath);
out_gt = flowToColor(img_flo_gt);
u_gt = img_flo_gt(:, :, 1); v_gt = img_flo_gt(:, :, 2);
xq = xx(380:396, 516:532); yq = yy(380:396, 516:532);
uq_gt = u_gt(380:396,516:532); vq_gt = v_gt(380:396, 516:532);
figure(7); imshow(im_n);title('ground truth frame 25');
hold on; 
quiver(xq, yq, uq_gt, vq_gt,0, 'color', 'b', 'LineWidth', 1.5);
axis([516 532 380 396]);
