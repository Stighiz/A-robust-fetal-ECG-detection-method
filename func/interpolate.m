function interpolated_data = interpolate(data)
    interpolated_data = data;
    
    for ch=1:size(data, 2)
        channel = data(:, ch);

        if any(isnan(channel))

            valid_idx = find(~isnan(channel));
            nan_idx = isnan(channel);
            valid_vals = channel(valid_idx);
            
            interpolated_ch = interp1(valid_idx, valid_vals, find(nan_idx), 'linear', 'extrap');
            channel(nan_idx) = interpolated_ch;
            interpolated_data(:, ch) = channel;
        end 

    end
end 