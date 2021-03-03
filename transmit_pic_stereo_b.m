clear;
close all;

%% Parameters
N_P = 4;
M = 2^N_P;
symbol_count = 1024;
bitPerFrame = symbol_count * N_P;
fftsize = (symbol_count+1)*2;
CP = 2047;
fs = 16e3;
Lt = 10;    %each packet of train contains Lt frames;
Ld = 10;

%% train signals
trainBlock = randi([0,1], bitPerFrame, 1);
trainPacket = repmat(trainBlock, 1, Lt);
empty = zeros(bitPerFrame, Lt);

%% toplay packets(in bit)
toplay1 = [trainPacket, empty];
toplay2 = [empty, trainPacket];

%% toplay streams(in bit)
toplay1 = reshape(toplay1, [], 1);
toplay2 = reshape(toplay2, [], 1);

%% qam_mod(in symbol)
qamTrain1 = qam_mod(toplay1, M);
qamTrain2 = qam_mod(toplay2, M);

%% ofdm_mod
ofdmStream1 = ofdm_mod(qamTrain1, symbol_count, CP);
ofdmStream2 = ofdm_mod(qamTrain2, symbol_count, CP);

%% impulse generation
pulse = [0 1 1 1 1 1 1 1 1 0];
IRlength = 4500;
ofdmStream1_t = transpose(ofdmStream1);
ofdmStream2_t = transpose(ofdmStream2);

[simin,nbsecs,fs] = initparams_stereo(ofdmStream1_t, ofdmStream2_t, fs, pulse, IRlength);
figure(3);
plot(simin);
title('simin');

sim('recplay');
received = simout.signals.values;
received = transpose(received);

%% Alignment;
out_aligned = alignIO(received,pulse,IRlength);

%% ofdm_demod and channel estimation
frameCount_est = Lt*2;
[H1, H2] = channelest(out_aligned, CP, symbol_count, frameCount_est, ofdmStream1, ofdmStream2);
H1 = transpose(H1);
H2 = transpose(H2);

%% beamformer
[a, b, H12] = fixed_transmitter_side_beamformer(H1, H2);

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
save('imageData.mat', 'imageData');
save('colorMap.mat', 'colorMap');

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
bitPerPacket = Ld*bitPerFrame;
packetCount = ceil(length(bitStream_padded)/bitPerPacket);

%% Qam modulation for trainblock
qamTrainBlock = qam_mod(trainBlock, M);
qamTrainStream = repmat(qamTrainBlock, Lt*packetCount, 1);

%% Qam modulation for bitStream and reshape
qamBitStream = qam_mod(bitStream_padded, M);
qamBitParallel = reshape(qamBitStream, symbol_count, bitFrames);

%% Pair Lt Qamtrain blocks with Ld Qambit blocks to form one packet;
qam = [];
for i = 1:packetCount
    qam_train_ith = repmat(qamTrainBlock, 1, Lt);
    if i == packetCount
        qam_bit_ith = qamBitParallel(:, (i-1)*Ld+1:end);
    else
        qam_bit_ith = qamBitParallel(:, (i-1)*Ld+1:i*Ld);
    end
    qam_ith = [qam_train_ith qam_bit_ith];
    qam = [qam qam_ith];
end

%% Serialize stream pairs;
qam_pairs_serial = reshape(qam, [], 1);
qam1 = qam_pairs_serial;
qam2 = qam_pairs_serial;

%% OFDM modulation
L = size(qam, 2);

% In order to avoid unexpected confliction, 
% I renamed ofdm_mod_stereo to ofdm_mod_stereo_b
%  a = 1;
%  b = 0;
[ofdmStream1_trans, ofdmStream2_trans] = ofdm_mod_stereo_b(a, b, qam1, qam2, symbol_count, CP);
ofdmTrainStream = ofdm_mod(qamTrainStream, symbol_count, CP);

%% impulse generation
pulse = [0 1 1 1 1 1 1 1 1 0];
IRlength = 4500;
ofdmStream1_trans = transpose(ofdmStream1_trans);
ofdmStream2_trans = transpose(ofdmStream2_trans);

[simin,nbsecs,fs] = initparams_stereo(ofdmStream1_trans, ofdmStream2_trans, fs, pulse, IRlength);
figure(3);
plot(simin);
title('simin');

sim('recplay');
received = simout.signals.values;
received = transpose(received);

%% Alignment;
out_aligned = alignIO(received,pulse,IRlength);

%% ofdm_demod
frameCount = (Lt+Ld)*packetCount;
[rxQamStream] = ofdm_demod_stereo_b(out_aligned,symbol_count,CP,Ld,Lt,frameCount,packetCount,ofdmTrainStream,H12);
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

