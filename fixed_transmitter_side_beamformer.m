function [a,b,H] = fixed_transmitter_side_beamformer(H1,H2)

H = sqrt(H1.*conj(H1) + H2.*conj(H2));
a = conj(H1)./H;
b = conj(H2)./H;
a = conj(a)';
b = conj(b)';
H = H';

end
