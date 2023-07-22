
function upsampled = Bcubicupsample(flow_data)
    % Get the size of the input data
    [h, w, c] = size(flow_data);

    % Twice both horizontally and vertically
    upsampled = zeros(h*2, w*2, c);

    % Apply bicubic interpolation for each channel
    for z = 1:c
        upsampled(:,:,z) = 2*imresize(flow_data(:,:,z), [h*2, w*2], 'bicubic');
    end
end
