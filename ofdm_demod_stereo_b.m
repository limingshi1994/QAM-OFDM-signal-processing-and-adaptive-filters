function [qam] = ofdm_demod_stereo_b(out_aligned,symbol_count,CP,Ld,Lt,frameCount,packetCount,ofdmTrainStream,H12)
fftsize = (symbol_count+1)*2;

%% pickout training and data packets;
startIndex = 1;
out_aligned = out_aligned(startIndex : startIndex+(fftsize+CP)*frameCount-1);


%% Cut off CP;
out_aligned = reshape(out_aligned, fftsize+CP, []);
out_aligned = out_aligned((CP+1):(fftsize+CP), :);
ofdmTrainStream = reshape(ofdmTrainStream, fftsize+CP, []);
ofdmTrainStream = ofdmTrainStream((CP+1):(fftsize+CP), :);

%% Pick out train and data packets;
ofdm_train_y = [];
ofdm_train_x = [];
ofdm_data = [];
ofdm_train_x_ref = ofdmTrainStream;
h = [];
HLS = [];
fft_ofdm_data_eq = [];
data_cell = cell(packetCount, 1);

for i = 1:packetCount
    iTrainStart = (i-1)*(Lt+Ld)+1;
    iTrainEnd = iTrainStart+Lt-1;
    iDataStart = iTrainEnd+1;
    if i == packetCount
        iDataEnd = frameCount;
    else
        iDataEnd = iDataStart+Ld-1;
    end
    ofdm_train_y_i = out_aligned(:, iTrainStart:iTrainEnd);
    ofdm_train_x_i = ofdm_train_x_ref(:, ((i-1)*Lt+1:i*Lt));
    ofdm_data_i = out_aligned(:, iDataStart:iDataEnd);
    ofdm_train_y = [ofdm_train_y ofdm_train_y_i];
    ofdm_train_x = [ofdm_train_x ofdm_train_x_i];
    ofdm_data = [ofdm_data ofdm_data_i];
    
    %% take fft;
    fft_ofdm_trainy_i = fft(ofdm_train_y_i,fftsize);
    fft_ofdm_trainx_i = fft(ofdm_train_x_i,fftsize);
    fft_ofdm_data_i = fft(ofdm_data_i, fftsize);
    
    %% (IF USE TRAINING PACKETS)
    %% H for each training packet;
    HLS_i = zeros(fftsize, 1);
    for iH = 1 : fftsize
        if (iH == 1 || iH == fftsize/2+1)
            HLS_i(iH, 1) = 0;
        else
            HLS_i(iH, 1) = mean(fft_ofdm_trainy_i(iH, :) ./ fft_ofdm_trainx_i(iH, :));
        end
    end
    HLS = [HLS HLS_i];
    
    %% IR;
    h_i = ifft(HLS_i, fftsize);
    h(:, i) = h_i;
    
    %% Plot;
    figure(13);
    plot(1:1:fftsize, abs(HLS_i(:,1)));
    xlim([0 fftsize/2])
    figure(14);
    plot(1:1:length(h_i), h_i);
    
    %% Equalizer;
    %HLS_i = repmat(HLS_i, 1, size(fft_ofdm_data_i, 2));
    fft_ofdm_data_eq_i = fft_ofdm_data_i ./ HLS_i;
    fft_ofdm_data_eq = [fft_ofdm_data_eq fft_ofdm_data_eq_i];

%     %% (USE H12 estimation previously)
%     fft_ofdm_data_eq_i = fft_ofdm_data_i ./ H12;
%     fft_ofdm_data_eq = [fft_ofdm_data_eq fft_ofdm_data_eq_i];
    
    %% cell for each loop
    data_cell{i} = fft_ofdm_data_eq;
    save('fft_ofdm_data.mat', 'data_cell');
    
end

%% save impulse response
save('h_i.mat', 'h');
save('HLS.mat', 'HLS');

ofdm = fft_ofdm_data_eq(2:fftsize/2,:);
%Get qam in a series vector
%FFTSIZE/2 = N+1 so we do -1 to get N.
qam = reshape(ofdm,[],1);
%     qam = qam(1 : end-pad_size);
%calculate
end