function filtered_data = filterData(ecg, hpf, num_ch)
    filtered_data = zeros(size(ecg));
    for i = 1:num_ch
        filtered_data(:, i) = filter(hpf, 1, ecg(:, i));
    end
end