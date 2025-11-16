clearvars;  clc;    close all; 
addpath '.\func'

fs = 1000;                     
target_fs = 2000;  

count_success = 0;

files = string(compose('%02d', 1:1:25));

results = struct('Found_MHR', {}, 'Found_FHR', {}, 'Target_FHR', {}, 'Success', {});
result_index=1;


for file_id = 1:length(files)
                                                                                           
    % Load solutions for later testing
    Rpeaks_file = './dataset/a' + files(file_id) + '.fqrs.txt'; 
    fid = fopen(Rpeaks_file, 'r');
    if fid == -1
        error('Unable to open the file');
    end
    target_Rpeaks_locs = fscanf(fid, '%d');
    target_FHR = length(target_Rpeaks_locs);
    fclose(fid);

    % Load data
    data_file = './dataset/a' + files(file_id) + '.csv';
    data = readmatrix(data_file);
    
    data = interpolate(data(:, 2:end));     % Replace NaN values using interpolation
   
    [N, num_ch] = size(data);    
    t = (0:N-1) / fs;
    f = (-N/2:N/2-1)*(fs/N);
    
    % figure('NumberTitle', 'off', 'Name', ('Signal no. ' + string(file_id)));
    % az(1) = subplot(3, 2, 1); 
    % plotData(t, data, 'Original signal', 'Time (s)', 'Amplitude');
    % ax(1) = subplot(3, 2, 2); 
    % data_spec = abs(fftshift(fft(data)));
    % plotData(f, data_spec, 'Original signal', 'Frequency', 'Magnitude', 1);  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BASELINE WANDER REMOVAL
    
    fc = 3;                                             % Cutoff frequency in Hz
    filter_order = 1000;                                % Filter order
    hpf = fir1(filter_order, fc/(fs/2), 'high');        % High-pass filter design
    filtered_data = filterData(data, hpf, num_ch);      % High-pass filter application
    
    % filtered_data_spec = abs(fftshift(fft(filtered_data)));
    % az(2) = subplot(3, 2, 3);
    % plotData(t, filtered_data, 'Highpass filtered signal at 3Hz', 'Time (s)', 'Amplitude');
    % ax(2) = subplot(3, 2, 4);
    % plotData(f, filtered_data_spec, 'Highpass filtered signal at 3Hz', 'Frequency', 'Magnitude', 1);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % POWER LINE INTERFERENCE CANCELLATION
   
    f0 = 50;                      % Power line frequency                                    
    bw = 0.06;                    % Filter bandwidth
    
    wo = f0 / (fs/2);            % Normalized frequency
    [b, a] = iirnotch(wo, bw);                                         % Notch filter design
    filtered_data = applyNotchFilter(filtered_data, b, a, num_ch);     % Notch filter application

    % filtered_data_spec = abs(fftshift(fft(filtered_data)));
    % az(3) = subplot(3, 2, 5);
    % plotData(t, filtered_data, 'Notch filtered signal', 'Time (s)', 'Amplitude');
    % ax(3) = subplot(3, 2, 6); 
    % linkaxes(ax,'x');
    % linkaxes(az,'x');
    % plotData(f, filtered_data_spec, 'Notch filtered signal',  'Frequency', 'Magnitude', 1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATERNAL QRS ENHANCED DETECTION USING PCA
    
    data_mean = mean(filtered_data);    
    data_std = std(filtered_data);
    data_norm = (filtered_data - data_mean) ./ data_std;    

    upsampled_data = resample(data_norm, target_fs, fs);
    [N, num_ch] = size(upsampled_data);    
    t = (0:N-1) / target_fs;

    data_score = performPCA(upsampled_data);
     
    [~, Rpeaks_locs, ~] = panTompkin(data_score, target_fs, 0);
    MHR = length(Rpeaks_locs); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MECG TEMPLATE ESTIMATION FOR EACH CHANNEL
    
    mecg_template = cell(1, 4);
    for ch=1:num_ch
        QRS_segments = extractQRSComplexes(upsampled_data(:,ch), Rpeaks_locs, target_fs);   
        mecg_template{ch} =  mean(QRS_segments, 1);
    end
    
   % figure('NumberTitle', 'off', 'Name', ('MECG TEMPLATE, Signal no. ' + string(file_id)));
   %  for ch = 1:num_ch
   %      ax_t = subplot(2, 2, ch);
   %      plot(mecg_template{ch});
   %      title(['MECG template - channel ', num2str(ch)]);
   %      xlabel('Samples');
   %      ylabel('Amplitude');
   %      grid on; 
   %  end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MECG CANCELLATION FOR EACH CHANNEL

    mecg_cancelled_data = zeros(size(upsampled_data));

    for ch = 1:num_ch
        mecg_cancelled_ch = upsampled_data(:, ch);

        for i = 1:MHR
            start_idx = Rpeaks_locs(i) - round(0.25 * target_fs);
            end_idx = Rpeaks_locs(i) + round(0.45 * target_fs) - 1;

            if start_idx > 0 && end_idx <= length(mecg_cancelled_ch)    
                mecg_cancelled_ch(start_idx:end_idx) = mecg_cancelled_ch(start_idx:end_idx) - mecg_template{ch}';
            end
        end
        mecg_cancelled_data(:, ch) = mecg_cancelled_ch;
    end
    
    % figure('NumberTitle', 'off', 'Name', ('MECG CENCELLATION, Signal no. ' + string(file_id)));
    % for ch = 1:num_ch
    %     ax_canc(1) = subplot(num_ch, 2, 2 * ch - 1);
    %     plotData(t, upsampled_data(:, ch), ['MECG NON Cancelled - Channel ', num2str(ch)], 'Time (s)', 'Amplitude');
    %     ax_canc(2)= subplot(num_ch, 2, 2 * ch);
    %     plotData(t, mecg_cancelled_data(:, ch), ['MECG Cancelled - Channel ', num2str(ch)], 'Time (s)', 'Amplitude');
    %     linkaxes(ax_canc,'x');
    % end
   


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FETAL QRS ENHANCED DETECTION USING PCA

    data_mean = mean(mecg_cancelled_data);
    data_std = std(mecg_cancelled_data);
    mecg_cancelled_data_norm = (mecg_cancelled_data - data_mean) ./ data_std;

    mecg_cancelled_data_score = performPCA(mecg_cancelled_data_norm);
    mecg_cancelled_data_score = filloutliers(mecg_cancelled_data_score, 'linear', 'ThresholdFactor', 6);
    
    [~, found_Rpeaks_locs, ~] = panTompkin(mecg_cancelled_data_score, target_fs, 0);
    found_FHR = length(found_Rpeaks_locs);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESULTS

    results(file_id).Found_FHR = found_FHR;
    results(file_id).Found_MHR = MHR;
    results(file_id).Target_FHR = target_FHR;
    
    error = abs(target_FHR - found_FHR);
    
    if error < 11
        count_success = count_success + 1;
        results(file_id).Success = true;
        mark = "✔";
    else 
        results(file_id).Success = false;
        mark = "✗";
    end

    disp("Signal no." + string(compose('%02d', file_id)) + ": " + char(9) + ...
        "Target FHR: " + string(target_FHR) + "," + char(9)+  'Found FHR: ' + string(found_FHR) ...
        + ',' + char(9) + char(9) + 'Error: ' + error + " bpm" + char(9) + mark);
    
    
end

count_failed = length(files) - count_success;
disp(newline + "No. of signals with correctly identified FHR: " + char(9) + string(count_success))
disp("No. of signals with wrongly identified FHR: " + char(9) + string(count_failed))
disp("SUCCESS RATE = " + string(count_success/length(files))+ "%")
disp("FAILURE RATE = " + string(count_failed/length(files))+ "%")

        
                       














                



