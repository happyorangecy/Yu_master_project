function upsampled = upsample_twice(flow_data)
    % Get the size of the input data
    [h, w, c] = size(flow_data);

    % Twice both horizontally and vertically
    upsampled = zeros(h*2, w*2, c);

    for i = 1:h
        for j = 1:w
            for z = 1:c
                upsampled(2*i-1:2*i, 2*j-1:2*j, z) = 2*flow_data(i, j, z);
            end
        end
    end
end
