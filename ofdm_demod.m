function [H_est] = ofdm_demod(Signal,N,CP,ofdmDummy)
    % How many Complex numbers are in each frame is known
    fftsize = (N+1)*2;
    
    %% cut the signal
    Signal = Signal(1:length(ofdmDummy));
    
    P = length(Signal)/(fftsize+CP);
    
    %Let's reshape since we know fftsize*P = Signal length
    ofdm = reshape(Signal,(fftsize+CP),P);
    ofdmx = reshape(ofdmDummy, fftsize+CP, P);
    
    %Remove Cyclic prefix before we take FFT
    ofdm = ofdm(((CP+1):(fftsize+CP)),:);
    ofdmx = ofdmx(((CP+1):(fftsize+CP)),:);
    
    %Take the FFT
    fft_ofdm = fft(ofdm, fftsize);
    fft_ofdmx = fft(ofdmx, fftsize);
    
    %estimation
    for i = 1:fftsize
        if (i == 1 || i == fftsize/2+1)
            H_est(i,1) = 0;
        else
            H_est(i,1) = mean(fft_ofdm(i,:) ./ fft_ofdmx(i,:));
        end
    end
    
    figure(21)
    plot(abs(H_est));
    xlim([1 fftsize/2]);
    
 end