function visualize_demod_adapt(fs, fftsize, ik_on, colorMap, imageData, Frame_count, bk_on, bitStream, imageSize, bitsPerPixel)

% fs = 16e3;
% fftsize = 2050;
% Lt = 5;
% Ld = 5;
t_perframe = fftsize*(1/fs);
t = t_perframe;

figure('Name','ResultDisplay','NumberTitle','off');
load('h_i.mat', 'h');
load('HLS.mat', 'H');
load('data_adapt.mat', 'qamp');

for i = 1:size(h,2)
    %% IR
    subplot(2,2,1)
    plot(1:1:size(h,1), h(:,i));
%    ylim([min(min(h)) max(max(h))]);
    
    %% Professor Moonen
    subplot(2,2,2)
    colormap(colorMap);
    image(imageData);
    axis image;
    
    %% CFR
    subplot(2,2,3)
    H = abs(H);
%   H = H(1:fftsize/2, :);
    x = normalize(2:fftsize/2,'range',[0 fs/2]);
    plot(x, 20*log10(H(:,i)));
    
    %% Transmitted
    subplot(2,2,4)
    data = qamp(:,1:i);
    empty = zeros(length(ik_on), Frame_count-size(data,2));
    data = [data empty];
    rx = [];
    for ik = 1:length(ik_on)
        M = 2^bk_on(ik);
        rx_ik = qam_demod(data(ik,:), M);
        rx = [rx; rx_ik];
    end
    bit = reshape(rx,[],1);
    bit = bit(1:length(bitStream));
    imageRx = bitstreamtoimage(bit, imageSize, bitsPerPixel);
    colormap(colorMap);
    image(imageRx);
    axis image;
    
    pause(t);
end

