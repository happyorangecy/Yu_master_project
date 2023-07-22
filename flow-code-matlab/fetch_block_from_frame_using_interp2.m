function mc_block = fetch_block_from_frame_using_interp2(imgI_n, bx, by, blockflow)

    % know image dimensions
    [height, width, channels] = size(imgI_n);

    % Create grids for the x and y coordinates
    [x_grid, y_grid] = meshgrid(1:width, 1:height);

    % Compute the block coordinates using the flow vectors
    brows = length(by);
    bcols = length(bx);
    new_bx_coords = ones(brows, 1) * bx + blockflow(:, :, 1);
    new_by_coords = by' * ones(1, bcols)  + blockflow(:, :, 2);

    % Initialize the warped image
    brows = length(by);
    bcols = length(bx);
    mc_block = zeros(brows, bcols, channels);

    % Warp each channel using interp2
    for ch = 1:channels
        mc_block(:, :, ch) = interp2(x_grid, y_grid, imgI_n(:, :, ch), new_bx_coords, new_by_coords, 'linear', 0);
    end

end
