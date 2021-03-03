clear all;
close all;

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
save('imageData.mat', 'imageData');
save('colorMap.mat', 'colorMap');

%% Parameters
N_P = 4;
M = 2^N_P;
symbol_count = 1024;
fftsize = (symbol_count+1)*2;
CP = 2047;
fs = 16e3;
Ld = 5;    %each packet of data contains Ld frames;
Lt = 5;    %each packet of train contains Lt frames;

%% BitStream padding;
bitBlock_size = symbol_count*N_P;
if mod(length(bitStream), bitBlock_size) ~= 0
    bitFrames = fix(length(bitStream) / bitBlock_size)+1;
    residue = bitStream((bitFrames-1)*bitBlock_size+1:end);
    topad = zeros(bitBlock_size-length(residue), 1);
    bitStream_padded = [bitStream; topad];
else
    bitFrames = length(bitStream) / bitBlock_size;
    bitStream_padded = bitStream;
end

%% Calculation of packet count;
bitPerPacket = Ld*bitBlock_size;
packetCount = ceil(length(bitStream_padded)/bitPerPacket);

%% Trainblock
trainBlock = randi([0 1], bitBlock_size, 1);

%% Dummy trainpacket
qamTrainBlock = qam_mod(trainBlock, M);
qamDummy = repmat(qamTrainBlock, Lt, 1);

%% Ofdm for dummy
ofdmDummy = ofdm_mod(qamDummy, symbol_count, CP);

%% Play signal and record, align
pulse = [0 1 1 1 1 1 1 1 1 0];
IRlength = 4500;
ofdmDummy = transpose(ofdmDummy);

[simin,nbsecs,fs] = initparams(ofdmDummy,fs,pulse,IRlength);
figure(1);
plot(simin);
title('simin');

sim('recplay');
received = simout.signals.values;
received = transpose(received);

%% Alignment;
out_aligned = alignIO(received,pulse,IRlength);

%% Estimation of channel
H_est = ofdm_demod(out_aligned, symbol_count, CP, ofdmDummy);
h_estDum = ifft(H_est, fftsize);
H_estDum = H_est;
save('H_estDum.mat', 'H_estDum');
save('h_estDum.mat', 'h_estDum');

%% Adaptive bit loading
% Calculate M size for adaptive QAM
k = 1:1:fftsize;
k = k.';
g = 1e-2; %gamma
noise = wgn(fftsize,1,0,1,0);
P = fft(noise, fftsize);
P = abs(P).^2/fftsize;
figure(22)
plot(k, P);
bk = floor(log2(1+(abs(H_est).^2)./(g*P)));
figure(23)
plot(k, bk);

%% Adaptive indices
% find indices of useful bins
ik_on = find(bk ~= 0);
ik_off = find(bk == 0);
save('ik_off.mat', 'ik_off');

% QAM value
bk_on = bk(ik_on);
bk_on(bk_on>6) = 6;

% bits per frame
bitsPerFrame_adapt = sum(bk_on(1:length(bk_on)/2));

% adaptive frame count
frameCount_adapt = ceil(length(bitStream)/bitsPerFrame_adapt);

% padding zeros to bitstream
topad_adapt = zeros(frameCount_adapt*bitsPerFrame_adapt-length(bitStream), 1); 
bitStreamPadded_adapt = [bitStream; topad_adapt];

% bitStream to parallel, each bin with different M=2^N
bitStreamPar_adapt = reshape(bitStreamPadded_adapt, bitsPerFrame_adapt, frameCount_adapt);

%% Qam modulation for bitStream in parallel
% split matrix into matrices with N rows for each bin
bk_on = bk_on(1:length(bk_on)/2);
bins_adapt = mat2cell(bitStreamPar_adapt, bk_on);
qamBitStreamPar_adapt = [];
for i = 1:length(ik_on)/2 %bin count
    M = 2^(bk_on(i));
    bin_adapt_i = bins_adapt{i};
    bin_adapt_i_Ser = reshape(bin_adapt_i, [], 1);
    qambin_adapt_i = qam_mod(bin_adapt_i_Ser, M);
    qambin_adapt_i = reshape(qambin_adapt_i, 1, frameCount_adapt);
    qamBitStreamPar_adapt = [qamBitStreamPar_adapt; qambin_adapt_i];
end

%% Insert useful qam bins into empty parallel matrix
ik_on = ik_on(1:length(ik_on)/2);
qamAll = zeros(symbol_count, frameCount_adapt);
for i = 1:length(ik_on)
    qamAll(ik_on(i), :) = qamBitStreamPar_adapt(i, :); 
end

qamAll_adapt = reshape(qamAll, [], 1);

%% OFDM modulation
ofdmStream = ofdm_mod_adapt(qamAll_adapt, symbol_count, CP);

%% Play signal and record, align
pulse = [0 1 1 1 1 1 1 1 1 0];
IRlength = 4500;
ofdmStream = transpose(ofdmStream);

[simin,nbsecs,fs] = initparams(ofdmStream,fs,pulse,IRlength);
figure(1);
plot(simin);
title('simin');

sim('recplay');
received = simout.signals.values;
received = transpose(received);

%% Alignment;
out_aligned = alignIO(received,pulse,IRlength);

%% OFDM demodulation
[rxQamStream] = ofdm_demod_eq_adapt(out_aligned, symbol_count, CP, H_est, ik_on, frameCount_adapt, N_P);

%% QAM demodulation
rxQamPar = reshape(rxQamStream, length(ik_on), frameCount_adapt);
rxBitPar = [];
for i = 1:length(ik_on)
    M = 2^bk_on(i);
    rxBitPar_i = qam_demod(rxQamPar(i,:), M);
    rxBitPar = [rxBitPar; rxBitPar_i];
end
rxBitStream = reshape(rxBitPar, [], 1);
rxBitStream = rxBitStream(1:length(bitStream));

%% Compute BER and show
berTransmission = ber(bitStream,rxBitStream)

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);
save('imageRx.mat', 'imageRx');

%% Plot images
figure(12)
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); 
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); 

%% visualize
close all;
visualize_demod_adapt(fs, fftsize, ik_on, colorMap, imageData, frameCount_adapt, bk_on, bitStream, imageSize, bitsPerPixel);

