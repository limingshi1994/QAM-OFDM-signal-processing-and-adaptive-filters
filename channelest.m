function [H_est1,H_est2] = channelest(received, CP, symbol_count, frameCount, ofdmStream1, ofdmStream2)
% How many Complex numbers are in each frame is known
fftsize = (symbol_count+1)*2;

%% cut the signal
startIndex = 1;
received = received(startIndex : startIndex+(fftsize+CP)*frameCount-1);

%% h1
ofdm1 = received(1:length(received)/2);
ofdmStream1_half = ofdmStream1(1:length(ofdmStream1)/2);

%% h2
ofdm2 = received(length(received)/2+1:end);
ofdmStream2_half = ofdmStream2(1:length(ofdmStream2)/2);

%Let's reshape since we know fftsize*P = Signal length
ofdm1y = reshape(ofdm1,(fftsize+CP),[]);
ofdm2y = reshape(ofdm2,(fftsize+CP),[]);
ofdm1x = reshape(ofdmStream1_half,(fftsize+CP),[]);
ofdm2x = reshape(ofdmStream2_half,(fftsize+CP),[]);

%Remove Cyclic prefix before we take FFT
ofdm1y = ofdm1y(((CP+1):(fftsize+CP)),:);
ofdm2y = ofdm2y(((CP+1):(fftsize+CP)),:);
ofdm1x = ofdm1x(((CP+1):(fftsize+CP)),:);
ofdm2x = ofdm2x(((CP+1):(fftsize+CP)),:);

%Take the FFT
fft_ofdm1y = fft(ofdm1y, fftsize);
fft_ofdm2y = fft(ofdm2y, fftsize);
fft_ofdm1x = fft(ofdm1x, fftsize);
fft_ofdm2x = fft(ofdm2x, fftsize);

%estimation
for i = 1:fftsize
    if (i == 1 || i == fftsize/2+1)
        H_est1(i,1) = 1e-6;
        H_est2(i,1) = 1e-6;
    else
        H_est1(i,1) = mean(fft_ofdm1y(i,:) ./ fft_ofdm1x(i,:));
        H_est2(i,1) = mean(fft_ofdm2y(i,:) ./ fft_ofdm2x(i,:));
    end
end

figure(21)
plot(1:1:length(H_est1), 20*log(abs(H_est1)));
xlim([1 length(H_est1)/2]);
figure(22)
plot(1:1:length(H_est2), 20*log(abs(H_est2)));
xlim([1 length(H_est2)/2]);

end

