close all
clear all

image_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\datasets\Sintel\training\clean\alley_1\';
output_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\Downhalf\';

k = 50;
for i = 1:k
    image_name = strcat(image_path, 'frame_', sprintf('%04d', i), '.png');
    image_data = imread(image_name);
    
    % Check if the image has correct size
    if size(image_data, 1) == 436 && size(image_data, 2) == 1024
        % Downsample the image
        downsampled_image = imresize(image_data, [218 512]);

        % Save the downsampled image
        output_image_name = strcat(output_path, 'frame_', sprintf('%04d', i), '_downsampled.png');
        imwrite(downsampled_image, output_image_name);
    else
        fprintf('Image %d is not of the expected size. Skipping...\n', i);
    end
end
