function [ofdmStream1, ofdmStream2] = ofdm_mod_stereo(a, b, qam1, qam2, N, CP)
%% Original without on-off bit
% Find how many frames we need with respect to N
L1 = length(qam1);
P1 = L1/N;
P2 = P1;

ofdm1 = reshape(qam1,N,P1);
ofdm2 = reshape(qam2,N,P2);

Null_1 = zeros(1,P1);
Null_2 = zeros(1,P2);

ofdm1 = [Null_1;ofdm1;Null_1;flipud(conj(ofdm1))];
ofdm2 = [Null_2;ofdm2;Null_2;flipud(conj(ofdm2))];

fftsize = (N+1)*2;
ofdm1 = a.*ofdm1;
ofdm2 = b.*ofdm2;
ifft_ofdm1 = ifft(ofdm1, fftsize);
ifft_ofdm2 = ifft(ofdm2, fftsize);

if not(CP == 0)
    ifft_ofdm1 = [ifft_ofdm1(((fftsize-CP+1):fftsize),:);ifft_ofdm1];
    ifft_ofdm2 = [ifft_ofdm2(((fftsize-CP+1):fftsize),:);ifft_ofdm2];
end

% Go from matrix (N,P) to series data (P*(2N+2),1) (column vector)
ofdmStream1 = reshape(ifft_ofdm1,P1*(fftsize+CP),1);
ofdmStream2 = reshape(ifft_ofdm2,P2*(fftsize+CP),1);

end