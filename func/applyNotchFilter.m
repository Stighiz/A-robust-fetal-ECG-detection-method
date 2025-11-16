function filtered_data = applyNotchFilter(data, b, a, num_ch)
    for j = 1:num_ch
        data(:, j) = filtfilt(b, a, data(:, j));
    end
    filtered_data = data;
end