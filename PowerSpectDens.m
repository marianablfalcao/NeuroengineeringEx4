function PowerSpectDens(data, Fs)

rng default

N = length(data);
xdft = fft(data);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(data):Fs/2;

figure
% plot(freq,10*log10(psdx))
plot(freq,psdx)
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')






