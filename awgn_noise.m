function [y , Noise_Power] = awgn_noise(hx,SNR)


[p , q] = size(hx);

Noise_Power = 10^(-SNR/10);

n = sqrt(Noise_Power/2) * (randn(p,q) + 1j*randn(p,q));

y = hx + n;