function [freq,amp,phase] = FFT_processed(t,P)
%FFT_PRESSURE Summary of this function goes here
%   Detailed explanation goes here
period=t(2)-t(1);
L=length(t);
N=ceil(L/2);
freq=(0:N-1)/period/L;
Y=fft(P)/L;
amp=abs(Y(1:N));
amp(1)=Y(1); % make sure get the sign of this one right
amp(2:end)=amp(2:end)*2;
phase=atan2(imag(Y(1:N)),real(Y(1:N)));


end

