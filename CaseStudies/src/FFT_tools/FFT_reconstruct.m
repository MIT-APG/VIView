function [Y] = FFT_reconstruct(t,freq,amp,phase, N)
%FFT_RECONSTRUCT Summary of this function goes here
%   Detailed explanation goes here

Y=0*t+amp(1);
for i=2:N
    Y=Y+amp(i)*cos(freq(i)*2*pi*t + phase(i));
end


end

