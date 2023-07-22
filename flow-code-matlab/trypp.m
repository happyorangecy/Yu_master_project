clear
close all
clear all

image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
hres        = 1024;
vres        = 436;
B           = 16;
w           = 16;
mae_t       = 10;
start_frame = 15;
nframes     = 30;

fprintf('processing the sequence\n')

for frame = start_frame+1:start_frame + nframes-1
    F_image_file1 = sprintf('frame_%04d.png', frame);
    F_image_file2 = sprintf('frame_%04d.png', frame+1);
    F_frame_n = imread([image_path F_image_file1]) ;
    F_frame_np1 = imread([image_path F_image_file2]) ;
    non_mc_dfd = abs(F_frame_np1-F_frame_n); 
    non_mc_all = mae(F_frame_n,F_frame_np1);
    non_mc_mae(frame) = mean(non_mc_all(:));


    [bdx, bdy, dfd] = blockmatching(F_frame_np1, F_frame_n, B, w, mae_t);
    [x, y] = meshgrid((B/2):B:hres-(B/2), (B/2):B:vres-(B/2)); 
    ly = length((B/2):B:vres-(B/2));
    lx = length((B/2):B:hres-(B/2));
    x = x(1:ly,1:lx);
    y = y(1:ly,1:lx);
%     figure; image((1:hres),(1:vres),F_frame_np1);colormap(gray(256)); axis image; 
%     hold on; title('Motion vectors for each block superimposed on current frame');
%     h = quiver(x(:), y(:), bdx(:), bdy(:), 0, 'b-');
%     set(h,'linewidth',1.5); 
%     xlabel('Columns'); ylabel('Rows'); hold off;drawnow;
    mae_mc_dfd = mae(dfd);
    mc_mae(frame) = mean(mae_mc_dfd(:));
end

figure(3);plot(non_mc_mae, 'o-');
hold on;
plot(mc_mae,'+-');
title('motion  compensated MAE ');


function [motion_x, motion_y, dfd] = blockmatching(curr_frame, other_frame, B, w, mae_t)
% [motion_x, motion_y, dfd] = blockmatching(curr_frame, other_frame, B, w, mae_t)
% 
% This function implements a simple Block Matching algorithm
% The block size is B (eg. 16), Edges of the image are
% ignored. The search width is w (eg. 4).
% The algorithm searches for a match between blocks only when the MAE for a
% block with its co-located block in the previous frame, exceeds mae_t the
% motion threshold.

[vres, hres] = size(curr_frame);

% lx, ly are the number of blocks across the picture
lx = length((B/2):B:hres-(B/2));
ly = length((B/2):B:vres-(B/2));

motion_x = zeros(ly,lx); %the horizontal component of motion
motion_y = zeros(ly,lx); %the vertical component of motion
dfd = zeros(vres, hres);

search_range_x = -w:w; % horizontal search range
search_range_y = -w:w; % vertical search range

non_mc_dfd = abs(curr_frame-other_frame);

% leave out a border of BxB pels so you don't have to bother about borders
ny = 2;
for j = B:B:vres-B+1-B+1
    nx = 2;
    for i = B:B:hres-B+1-B+1
        bx = i:i+B-1; by = j:j+B-1;
        ref_block = curr_frame(by,bx);
        non_mc_dfd_block = non_mc_dfd(by,bx);
        block_mae = mean(non_mc_dfd_block(:));
        if ( block_mae > mae_t )
            
            % searching for each possible block in the search window
            % in the past_frame and measure the mean abs dfd
            % for every offset block.
            
            mc_block = ref_block;
            min_error_ = +inf;
            for jj = search_range_y
                for ii = search_range_x
                    other_block = fetch_block(other_frame, by+jj, bx+ii);
                    % we use fetch_block (see implementation at end of file) 
                    % to deal with when outside of the image boundaries.

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Find the block with minimum DFD, save it to mc_block
                    % and assign its offset to 'motion_x(ny,nx)' 
                    % and 'motion_y(ny,nx)'
                                        
                    %write your code here
                    current_diff = sum(abs(ref_block-other_block));
                    current_diff = mean(current_diff(:));
              
                    if current_diff < min_error_
                        mc_block = other_block;
                        motion_y(ny,nx) = jj;
                        motion_x(ny,nx) = ii;
                        min_error_ = current_diff;
                    end
                end
            end
          
            dfd(by,bx) = ref_block - mc_block;
                
        else
            motion_x(ny,nx) = 0;
            motion_y(ny,nx) = 0;
        end
        nx = nx+1;
    end % end of horizontal scan
    ny = ny+1;
end % end of vertical scan

% fetch block with index range bx, by in frame. handle boundaries by
% repeating the boundary values.
function block = fetch_block(frame, by, bx)
    block = frame(max(min(by,end),1),max(min(bx,end),1));
end

end