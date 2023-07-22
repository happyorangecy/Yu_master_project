close all
clear all

% Groundtruth
 truth_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1\frame_0025.flo';
               % C:\Users\Chengyu\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\flow\alley_1
 img_flo_gt = readFlowFile(truth_path);
 out_gt = flowToColor(img_flo_gt);
 figure(1)
 imshow(out_gt)

 % Tiling: test if it has been tiled from 4 parts into a whole one
 flo_filepath = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\after_tiling\flow_00025_00026.flo';
 img_flo_f = readFlowFile(flo_filepath);
 out_f = flowToColor(img_flo_f);
 figure(2)
 imshow(out_f);

% Downsampled
 draft_filepath = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\DhalfR\flow_00025_00026.flo';
 img_flo_r = readFlowFile(draft_filepath);
 out_r = flowToColor(img_flo_r);
 figure(3)
 imshow(out_r);

 % Upsampled
 up_filepath = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwice\flow_00025_00026_upsampled.flo';
 img_flo_up = readFlowFile(up_filepath);
 out_up = flowToColor(img_flo_up);
 figure(4)
 imshow(out_up);