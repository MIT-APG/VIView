function [freq,Y] = FFT_complex(t,P)
%FFT_PRESSURE Summary of this function goes here
%   Detailed explanation goes here
period=t(2)-t(1);
L=length(t);
N=ceil(L/2);
freq=(0:N-1)/period/L;
Y=fft(P)/L;

end
