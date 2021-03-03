clear all;
close all;

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
save('imageData.mat', 'imageData');
save('colorMap.mat', 'colorMap');

%% Parameters
N_P = 4;
M = 2^N_P;
symbol_count = 512;
fftsize = (symbol_count+1)*2;
CP = 1023;
fs = 16e3;
Lt = 15;    %each packet of train contains Lt frames;

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

%% Trainblock
trainBlock = randi([0,1], bitBlock_size, 1);

%% Qam modulation for trainblock
qamTrainBlock = qam_mod(trainBlock, M);
qamTrainStream = repmat(qamTrainBlock, Lt, 1);

%% Qam modulation for bitStream and reshape
qamBitStream = qam_mod(bitStream_padded, M);
qamBitParallel = reshape(qamBitStream, symbol_count, []);
Frame_count = size(qamBitParallel,2)+Lt;

%% OFDM modulation
ofdmTrain = ofdm_mod(qamTrainStream,symbol_count,CP);
ofdmData  = ofdm_mod(qamBitStream,symbol_count,CP);

ofdmStream = [ofdmTrain;ofdmData];


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
[rxQamStream, channel_est] = ofdm_demod_eq(out_aligned, symbol_count, CP, Lt, Frame_count, ofdmTrain,N_P);
rxQamStream = rxQamStream(1 : length(qamBitStream));

%% QAM demodulation
rxBitStream = qam_demod(rxQamStream, M);
rxBitStream = rxBitStream(1:length(bitStream));

%% Compute BER and show
berTransmission = ber(bitStream,rxBitStream)

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure(12)
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); 
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); 

%% Visualize
close all;
visualize_demod(fs, fftsize, Lt, colorMap, imageData, Frame_count, M, bitStream, imageSize, bitsPerPixel);

