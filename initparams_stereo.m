function [simin, nbsecs, fs] = initparams_stereo(toplay1, toplay2, fs, pulse, IRlength)
%Period and frequency;
dt = 1/fs;

%silence times;
t_1 = 2/dt;
t_2 = 1/dt;

%silence vectors;
s1 = zeros([1, t_1]);
s2 = zeros([1, t_2]);

%Synchronization pulse
zeroIR = zeros([1, IRlength]);
                   
%concatenate;
simin_a = [s1 pulse zeroIR toplay1 s2];
simin_b = [s1 pulse zeroIR toplay2 s2];

%transverse;
simin = [simin_a; simin_b];
simin = simin';

%total time;
t = size(simin, 1);
nbsecs = dt * t;  

end
