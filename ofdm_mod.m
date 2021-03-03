function [ofdm] = ofdm_mod(qam, N, CP)% H, Threshold)
    %% Original without on-off bit
    % Find how many frames we need with respect to N
    L = length(qam);
    P = L/N;
    
    ofdm = reshape(qam,N,P);
    Null = zeros(1,P);
    ofdm = [Null;ofdm;Null;flipud(conj(ofdm))];
    
    fftsize = (N+1)*2;
    ifft_ofdm = ifft(ofdm,fftsize);
    
    if not(CP == 0)
        ifft_ofdm = [ifft_ofdm(((fftsize-CP+1):fftsize),:);ifft_ofdm];
    end
    
    % Go from matrix (N,P) to series data (P*(2N+2),1) (column vector)
    ofdm = reshape(ifft_ofdm,P*(fftsize+CP),1);
 end