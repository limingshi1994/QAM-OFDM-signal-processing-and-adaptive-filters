function [qam, HLS] = ofdm_demod_eq(out_aligned, symbol_count, CP, Lt, Frame_count, ofdmTrainStream,N_P)

fftsize = (symbol_count+1)*2;

%% pickout training and data packets;
startIndex = 1;
out_aligned = out_aligned(startIndex : startIndex+(fftsize+CP)*Frame_count-1);

%% Cut off CP;
out_aligned = reshape(out_aligned, fftsize+CP, []);
out_aligned = out_aligned((CP+1):(fftsize+CP), :);
ofdmTrainStream = reshape(ofdmTrainStream, fftsize+CP, []);
ofdmTrainStream = ofdmTrainStream((CP+1):(fftsize+CP), :);

%% Training input x,Training output y and data
ofdm_train_y = out_aligned(:,(1:Lt));
ofdm_train_x = ofdmTrainStream;
ofdm_data = out_aligned(:,(Lt+1:end));

%% 
Qam_train_y = fft(ofdm_train_y,fftsize);
Qam_train_x = fft(ofdm_train_x,fftsize); 
Qam_data = fft(ofdm_data,fftsize);
M = 2^N_P;

%% Yay
H = Qam_train_y./Qam_train_x;

for k = 1:fftsize
    Hmean(k,1) = mean(H(k,:));
end 

Hmean = Hmean(2:fftsize/2,:);

%%
Qam_data = Qam_data(2:fftsize/2,:);
stepsize = 0.1; % mu/alpha
Q_T = conj(Qam_data);
Qam_dataQam_data = Q_T.*Qam_data;
Mult = stepsize./Qam_dataQam_data;

iter = 1;
W = zeros(symbol_count, iter+1); %See what happens in filter coeffecients
W(:,1) = 1./conj(Hmean); % Inititialization of coeffecients
W_initial = W;
error = [];
W_final = [];

for i = 1:(Frame_count-Lt) % symbol_count = 1000
    
    W = W_initial;
    Qam_data_i = Qam_data(:,i).*conj(W(:,1));
    Qam_bit_Desired = qam_demod(Qam_data_i,M);
    Qam_Desired_data_i = qam_mod(Qam_bit_Desired,M);
    for j = 1:iter
        Qam_data_i = Qam_data(:,i).*conj(W(:,j));
        error(:,j) = conj(Qam_Desired_data_i - Qam_data_i);
        W(:,j+1) = W(:,j) + Mult(:,i).*Qam_data(:,i).*error(:,j);
    end
    W_final(:,i) = W(:,j+1);
end
%ploterr(error, Frame_count, Lt);
HLS = W_final(:,1:(Frame_count-Lt));

qam = Qam_data.*conj(HLS);

H = conj(1./HLS);
h = ifft(H,fftsize);

%% save impulse response
save('h_adaptfilter_1.mat', 'h');
save('H_adaptfilter_2.mat', 'H');
save('qam_adaptfilter.mat', 'qam');

%ofdm = fft_ofdm_data_eq(2:fftsize/2,:);
%Get qam in a series vector
%FFTSIZE/2 = N+1 so we do -1 to get N.
qam = reshape(qam,[],1);
%     qam = qam(1 : end-pad_size);
%calculate

end