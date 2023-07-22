function BIupsampled = bilineartwice(flow_data)
    % Get the size of the input data
    [h, w, c] = size(flow_data);

    % Initialize a new matrix with twice the size in both dimensions
    BIupsampled = zeros(h*2, w*2, c);

    % Get new matrix with interpolated values
    for z = 1:c
        for i = 1:(h*2)
            for j = 1:(w*2)
                % Find corresponding location in original data
                orig_i = ceil(i/2);
                orig_j = ceil(j/2);

                A = flow_data(orig_i, orig_j, z);

                % If we are not on a boundary, interpolate, else keep the value
                if mod(i,2) == 0 && mod(j,2) == 0 % Center point
                    if orig_j < w
                        B = flow_data(orig_i, orig_j+1, z);
                    else
                        B = A;
                    end

                    if orig_i < h
                        C = flow_data(orig_i+1, orig_j, z);
                    else
                        C = A;
                    end

                    if orig_i < h && orig_j < w
                        D = flow_data(orig_i+1, orig_j+1, z);
                    else
                        D = A;
                    end

                    BIupsampled(i, j, z) = 2 * (A + B + C + D) / 4;

                elseif mod(i,2) == 0 % Vertical edge point
                    if orig_i < h
                        C = flow_data(orig_i+1, orig_j, z);
                    else
                        C = A;
                    end
                    BIupsampled(i, j, z) = 2 * (A + C) / 2;
                    
                elseif mod(j,2) == 0 % Horizontal edge point
                    if orig_j < w
                        B = flow_data(orig_i, orig_j+1, z);
                    else
                        B = A;
                    end
                    BIupsampled(i, j, z) = 2 * (A + B) / 2;
                    
                else % Corner point
                    BIupsampled(i, j, z) = 2 * A;
                end
            end
        end
    end
end
