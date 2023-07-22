function I_np1_xd = warp_image_using_interp2(imgI_n, flow)

    % know image dimensions
    [height, width, channels] = size(imgI_n);

    % Create grids for the x and y coordinates
    [x_grid, y_grid] = meshgrid(1:width, 1:height);

    % Compute the new coordinates using the flow vectors
    new_x_coords = x_grid + flow(:, :, 1);
    new_y_coords = y_grid + flow(:, :, 2);

    % Initialize the warped image
    I_np1_xd = zeros(height, width, channels);

    % Warp each channel using interp2
    for ch = 1:channels
        I_np1_xd(:, :, ch) = interp2(x_grid, y_grid, imgI_n(:, :, ch), new_x_coords, new_y_coords, 'linear', 0);
    end

end
