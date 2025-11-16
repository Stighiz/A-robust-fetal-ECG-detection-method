function QRS_segments = extractQRSComplexes(singleChannelECG, Rpeaks_locs, target_fs)
    mecg_duration = 0.25 + 0.45;
    num_samples_mecg = round(mecg_duration * target_fs);
    QRS_segments = [];
    for i = 1:length(Rpeaks_locs)
        if Rpeaks_locs(i) > 0.25 * target_fs && Rpeaks_locs(i) + 0.45 * target_fs <= length(singleChannelECG)
            start_idx = Rpeaks_locs(i) - round(0.25 * target_fs);
            end_idx = Rpeaks_locs(i) + round(0.45 * target_fs);
            QRS_segments = [QRS_segments; singleChannelECG(start_idx:end_idx - 1)'];
        end
    end
end