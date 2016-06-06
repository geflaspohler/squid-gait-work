% PSD tst
Fs = 25;
x = OBDA(1:(floor(length(OBDA)/2))*2);
figure; plot(x)

N = length(x);
xdft = fft(x);

%Take first half of the spectrum
xdft = xdft(1:N/2+1);

psdx = (1/(Fs*N))*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

figure; subplot(311);
plot(freq, 10*log10(psdx));
grid on;
title('Overall PSD')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
axis([0 4 -100 0])

Fs = 25;
x = falling(1:(floor(length(falling)/2))*2);
figure; plot(x)

N = length(x);
xdft = fft(x);

%Take first half of the spectrum
xdft = xdft(1:N/2+1);

psdx = (1/(Fs*N))*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

subplot(312);
plot(freq, 10*log10(psdx));
grid on;
title('Rising PSD')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
axis([0 4 -100 0])

x = falling(1:(floor(length(falling)/2))*2);
figure; plot(x)

N = length(x);
xdft = fft(x);

%Take first half of the spectrum
xdft = xdft(1:N/2+1);

psdx = (1/(Fs*N))*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

subplot(313);
plot(freq, 10*log10(psdx));
grid on;
title('Falling PSD')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
axis([0 4 -100 0])
