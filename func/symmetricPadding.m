function padded_sig = symmetricPadding(sig, pad_len)
    left_pad = fliplr(sig(1:pad_len));   
    right_pad = fliplr(sig(end-pad_len+1:end)); 
    padded_sig = vertcat(left_pad, sig, right_pad); 
end