function [qam] = ofdm_demod_stereo(out_aligned,symbol_count,CP,H12)

fftsize = (symbol_count+1)*2;

%% Cut off CP;
out_aligned = reshape(out_aligned, fftsize+CP, []);
out_aligned = out_aligned((CP+1):(fftsize+CP), :);

%% fft
Qam_data = fft(out_aligned,fftsize);

%% equalization
Qam_data = Qam_data./H12;

%% output
qam = Qam_data(2:fftsize/2,:);
qam = reshape(qam,[],1);

end

