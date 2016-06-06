function [ fft_matrix_x, fft_matrix_y ] = FFT_data_set( data_matrix, sampling_freq )
%Outputs Single Sided FFT of input data_matrix  sampled at specified
%sampling_freq and outputs as fft_matrix

starta = 1; % Define the start time for the fft
enda = length(data_matrix); % Define the end time for the fft

Fs = sampling_freq;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = length(data_matrix(starta:enda,1));                     % Length of signal
t = (1:length(data_matrix(starta:enda,1)))/Fs;                % Time vector
y = data_matrix(starta:enda,1)';

NFFT = (2^nextpow2(L)); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f_trans = Fs/2*linspace(0,1,NFFT/2+1);

fft_matrix_x = f_trans;
fft_matrix_y = 2*abs(Y(1:NFFT/2+1));

end

