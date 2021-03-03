function visualize_demod(fs, fftsize, Lt, colorMap, imageData, Frame_count, M, bitStream, imageSize, bitsPerPixel)

% fs = 16e3;
% fftsize = 2050;
% Lt = 5;
% Ld = 5;
t_perframe = fftsize*(1/fs);
t = t_perframe;

figure('Name','ResultDisplay','NumberTitle','off');

load('h_adaptfilter_1.mat', 'h');
load('H_adaptfilter_2.mat', 'H');
load('qam_adaptfilter.mat', 'qam');

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
    data = qam(:,1:i);
    empty = zeros(fftsize/2-1, Frame_count-Lt-size(data,2));
    data = [data empty];
    qam_toprint = reshape(data,[],1);
    bit = qam_demod(qam_toprint, M);
    bit = bit(1:length(bitStream));
    imageRx = bitstreamtoimage(bit, imageSize, bitsPerPixel);
    colormap(colorMap);
    image(imageRx);
    axis image;
    
    pause(t);
end

