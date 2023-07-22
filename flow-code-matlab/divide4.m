close all
clear all

image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\market_2\';
K = 50; % Number of frames to process
output_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\DecomposedFrames5\';

% for i = 1:4
%     if ~exist([output_path, 'sequence_', num2str(i)], 'dir')
%         mkdir([output_path, 'sequence_', num2str(i)]);
%     end
% Create output directories if they do not exist
dir_names = {'topleft', 'topright', 'bottomleft', 'bottomright'};
for i = 1:4
    if ~exist([output_path, dir_names{i}], 'dir')
        mkdir([output_path, dir_names{i}]);
    end
end

% Loop through K frames
for k = 1:K
    % Read the frame
    frame = imread([image_path, sprintf('frame_%04d.png', k)]);
    
    % Get frame dimensions
    [M, N, ~] = size(frame);
    
    % Decompose frame into 4 smaller frames
    topleft = frame(1:M/2, 1:N/2, :);
    topright = frame(1:M/2, N/2+1:N, :);
    bottomleft = frame(M/2+1:M, 1:N/2, :);
    bottomright = frame(M/2+1:M, N/2+1:N, :);
    
    % Save the smaller frames into their respective directories
    imwrite(topleft, [output_path, 'topleft\', sprintf('frame_%04d.png', k)]);
    imwrite(topright, [output_path, 'topright\', sprintf('frame_%04d.png', k)]);
    imwrite(bottomleft, [output_path, 'bottomleft\', sprintf('frame_%04d.png', k)]);
    imwrite(bottomright, [output_path, 'bottomright\', sprintf('frame_%04d.png', k)]);
end

fprintf('Decomposition of %d frames into 4 sequences of smaller frames is completed.\n', K);
