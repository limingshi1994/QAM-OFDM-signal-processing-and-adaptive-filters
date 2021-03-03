clear all;
close all;

%% Parameters
N_P = 6;
M = 2^N_P;
symbol_count = N_P*1000;
Hk = 3+2i;

%% Qam 
bitBlock = randi([0,1], symbol_count, 1);
Xk = qam_mod(bitBlock, M);
Xk = Xk';
%% Output and filter updates with NLMS
Yk = Xk*Hk;

stepsize = 0.5; % mu/alpha
yy = Yk.*conj(Yk);
Mult = stepsize./yy;

iteration = 20; % Number of loops in the filter
W = zeros(iteration,1000); %See what happens in filter coeffecients
Werror = W; % See what's happening in the error
W(1,:) = (ones(1,1000)./conj(Hk)) + 0.001*rand(1,1000); % Inititialization of coeffecients


for i=1:iteration-1
    error = conj((Xk - Yk.*conj(W(i,:))));
    W(i+1,:) = W(i,:) + Mult.*Yk.*error;
    Werror(i,:) = W(i) - 1./conj(Hk);
end
%Final error value
Werror(iteration,:) = W(iteration) - 1/conj(Hk);


subplot(2,1,1);
y = 1:1:iteration;
plot(y,W(:,1));
title('filter coffecients with respect to iterations')
xlabel('Iterations')
ylabel('Filter coeffecients')

subplot(2,1,2);
plot(y,max(Werror'));
title('Error with respect to iterations')
xlabel('Iterations')
ylabel('Error')