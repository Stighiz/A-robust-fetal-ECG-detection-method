function wave_sig=wavelet(sig, gr)
    wt = modwt(sig, 4, 'sym4');
    wtwrec= zeros(size(wt));
    wtwrec(3:4, :) = wt(3:4, :); 
    y = imodwt(wtwrec, 'sym4');
    wave_sig = abs(y').^2;
        wave_sig = filloutliers(wave_sig,"linear","percentiles",[0,95]);

    if gr
        figure; subplot(2,1,1);
        plot(abs(sig));
        title("no wavelet")
        subplot(2,1,2);
        plot(wave_sig);
        title("wavelet")
    end
end 