clear;
close all;

%% Parameters
N_P = 4;
M = 2^N_P;
symbol_count = 512;
bitPerFrame = symbol_count * N_P;
fftsize = (symbol_count+1)*2;
CP = 1023;
fs = 16e3;
Lt = 15;    %each packet of train contains Lt frames;

%% Impulse response 1 & 2
h_length = 400;
hmin=1e-4;
hmax=2e-4;
hn=fftsize-h_length;
h0=hmin+rand(1,hn)*(hmax-hmin);
h1 = load('h1.mat', 'h').h;
h1 = h1(101:500);
h1 = [h1, h0];
h2 = load('h2.mat', 'h').h;
h2 = h2(101:500);
h2 = [h2, h0];
figure(1)
subplot(2,1,1)
plot(1:1:length(h1), h1);
title('h1');
subplot(2,1,2)
plot(1:1:length(h2), h2);
title('h2');

H1 = fft(h1, fftsize);
H2 = fft(h2, fftsize);

%% beamformer
[a, b, H12] = fixed_transmitter_side_beamformer(H1, H2);

% %% mono-transmission
%   a = 0; %a = 1;
%   b = 1; %b = 0;

%% dummy signals
sig1 = randi([0,1], bitPerFrame, 1);
sig2 = sig1; %randi([0,1], bitPerFrame, 1);

%% Qam
qam1 = qam_mod(sig1, M);
qam2 = qam_mod(sig2, M);

%% ofdm
[ofdmStream1, ofdmStream2] = ofdm_mod_stereo(a, b, qam1, qam2, symbol_count, CP);

%% convolve(vitual transmission)
receive1 = conv(h1, ofdmStream1);
receive1 = receive1(1:length(ofdmStream1));
receive2 = conv(h2, ofdmStream2);
receive2 = receive2(1:length(ofdmStream2));
receive = receive1+receive2;
figure(2)
plot(1:1:length(receive), receive);

%% noise
SNR = 100;
receive = awgn(receive, SNR);

%% ofdm_demod
[qam] = ofdm_demod_stereo(receive,symbol_count,CP,H12);

%% qam_demod
rxBitStream = qam_demod(qam, M);

%% BER
ber = ber(rxBitStream, sig1)
