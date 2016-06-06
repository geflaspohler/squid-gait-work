load lf086a_decimated A_dec pressure P_dec R_dec file_prefix %M_dec

conv_min_to_samp = 60*25; % Converst from min to samples


%Uncomment for lf086a
a = 1; %Tag started
b = 40*conv_min_to_samp; % Animal released (40*60*25)
c = 162*conv_min_to_samp; % Sunset (162*60*25)
d = 858*conv_min_to_samp; % Sunrise
e = 61291*25; % Reversal3
f = 1532821; % coasting1
g = 61319*25; % jetting1
h = 61761*25; % finning3
i = 61774*25; % turning3
j = 62149*25; % jetting2
k = 62368*25; % jetting3
l = 1652694;  % Tag off
m = a+1073124; % Finning

lf086a_b = struct('sunset', c, 'sunrise', d, 'reversal3', e, 'coasting1', f, 'jetting1', g, 'finning3', h, 'turning3', i, 'jetting2', j, 'jetting3', k, 'tagoff', l, 'finning', m);


if strcmp(file_prefix,'lf086a')
%     pressure = pressure(b:lf086a_b.tagoff);
%     A_dec = A_dec(b:lf086a_b.tagoff, :);
%     R_dec = R_dec(b:lf086a_b.tagoff);
%     P_dec = P_dec(b:lf086a_b.tagoff);
end

%Uncomment for lf082a
lf082a_b = struct('sunset', c, 'sunrise', d, 'jetting1', 82668*25, 'jetting2', 62668.288*25, 'forwards_finning',  82775.729*25, 'reversal', 82949.336*25);

if strcmp(file_prefix,'lf082a')
%     pressure = pressure(b:length(pressure));
%     A_dec = A_dec(b:length(A_dec), :);
%     R_dec = R_dec(b:length(R_dec));
%     P_dec = P_dec(b:length(P_dec));
end

%Uncomment for lf080a
lf080a_b = struct('sunset', c, 'sunrise', d, 'jetting1', 62284.472*25);
if strcmp(file_prefix,'lf080a')
%     pressure = pressure(b:1.59e6);
%     A_dec = A_dec(b:1.59e6, :);
%     R_dec = R_dec(b:1.59e6);
%     P_dec = P_dec(b:1.59e6);
end

figure; plot(A_dec)

%A_corrected = Rotation_correct(A_dec, P_dec, R_dec);

%%
A_x = A_dec(:, 1);
A_y = A_dec(:, 2);
A_z = A_dec(:, 3);
time = 1:length(A_x);
time = time./25; time = time';
figure; plot(time, A_x)

FIND_FLAG = 0;
if FIND_FLAG == 1
    running_average = OBDA_window_finder( A_x, 1, 250 );
    windows = 1:length(running_average);
    figure; plot(running_average, 'LineWidth',2);
    hold on
    scatter(windows, running_average);
end

%% Calculate OBDA
pressure_smooth = moving(pressure, 50);
pressure_normalized = pressure_smooth*-1;
pressure_diff = moving(diff(pressure_normalized)*50, 25);

STATIC_ACCEL_x = moving(A_x, 50);
STATIC_ACCEL_y = moving(A_y, 50);
STATIC_ACCEL_z = moving(A_z, 50);
STATIC_TOTAL = [STATIC_ACCEL_x STATIC_ACCEL_y STATIC_ACCEL_z];

DYNAMIC_ACCEL_x = abs(A_x-STATIC_ACCEL_x);
DYNAMIC_ACCEL_y = abs(A_y-STATIC_ACCEL_y);
DYNAMIC_ACCEL_z = abs(A_z-STATIC_ACCEL_z);
DYNAMIC_TOTAL = [DYNAMIC_ACCEL_x DYNAMIC_ACCEL_y DYNAMIC_ACCEL_z];

figure; plot(time, A_dec, time, DYNAMIC_TOTAL, time, STATIC_TOTAL)

OBDA = DYNAMIC_ACCEL_x+DYNAMIC_ACCEL_y+DYNAMIC_ACCEL_z;
dynamics = struct('obda', OBDA, 's_x', STATIC_ACCEL_x, 's_y', STATIC_ACCEL_y, 's_z', STATIC_ACCEL_z, ...
       'd_x', DYNAMIC_ACCEL_x, 'd_y', DYNAMIC_ACCEL_y, 'd_z', DYNAMIC_ACCEL_z, 'a_x', A_x, 'a_y', A_y, 'a_z', A_z);

figure; hold on
plot(time, A_dec, time, OBDA)

n = 1691*25; %finning1
finning1_end = 1695;
plot(n/25, A_x(n), '.r')
plot(finning1_end, A_x(finning1_end), '.r')
text(n/25, A_x(n), 'start')
text(finning1_end, A_x(finning1_end), 'end')

%% Plot Movements
figure; hold on;
zero = zeros(1, length(OBDA));
plot(time, STATIC_ACCEL_x, time, pressure_normalized, time, OBDA, time, zero)
%plot(time, A_x, time, A_y, time, A_z)
Combined = A_dec;

%offset = b;

% %Uncomment for lf080a
% lf080a_b.jetting1 = lf080a_b.jetting1-offset;
% plot(floor(lf080a_b.jetting1/25),  Combined(floor(lf080a_b.jetting1/25)), '.r')
% text(floor(lf080a_b.jetting1/25), Combined(floor(lf080a_b.jetting1/25)), 'jetting1')

% Uncomment for lf086a
% lf086a_b.sunset = c - offset;
% lf086a_b.sunrise = d - offset;
% lf086a_b.reversal3 = e - offset;
% lf086a_b.coasting1 = f - offset;
% lf086a_b.jetting1 = g - offset;
% lf086a_b.finning3 = h - offset;
% lf086a_b.turning3 = i - offset;
% lf086a_b.jetting2 = j - offset;
% lf086a_b.jetting3 = k - offset;
% lf086a_b.tagoff = l - offset;

plot(floor(lf086a_b.sunset/25), Combined(floor(lf086a_b.sunset/25)), '.r')
plot(floor(lf086a_b.sunrise/25), Combined(floor(lf086a_b.sunrise/25)), '.r')
plot(floor(lf086a_b.reversal3/25), Combined(floor(lf086a_b.reversal3/25)), '.r')
plot(floor(lf086a_b.jetting1/25), Combined(floor(lf086a_b.jetting1/25)), '.r')
plot(floor(lf086a_b.finning3/25), Combined(floor(lf086a_b.finning3/25)), '.r')
plot(floor(lf086a_b.turning3/25), Combined(floor(lf086a_b.turning3/25)), '.r')
plot(floor(lf086a_b.jetting2/25), Combined(floor(lf086a_b.jetting2/25)), '.r')
plot(floor(lf086a_b.jetting3/25), Combined(floor(lf086a_b.jetting3/25)), '.r')
plot(floor(lf086a_b.tagoff/25), Combined(floor(lf086a_b.tagoff/25)), '.r')

text(floor(lf086a_b.sunset/25), Combined(floor(lf086a_b.sunset/25)), 'sunset')
text(floor(lf086a_b.sunrise/25), Combined(floor(lf086a_b.sunrise/25)), 'sunrise')
text(floor(lf086a_b.reversal3/25), Combined(floor(lf086a_b.reversal3/25)), 'reversal2')
text(floor(lf086a_b.jetting1/25), Combined(floor(lf086a_b.jetting1/25)), 'jetting1')
text(floor(lf086a_b.finning3/25), Combined(floor(lf086a_b.finning3/25)), 'finning3')
text(floor(lf086a_b.turning3/25), Combined(floor(lf086a_b.turning3/25)), 'turning3')
text(floor(lf086a_b.jetting2/25), Combined(floor(lf086a_b.jetting2/25)), 'jetting2')
text(floor(lf086a_b.jetting3/25), Combined(floor(lf086a_b.jetting3/25)), 'jetting3')
text(floor(lf086a_b.tagoff/25), Combined(floor(lf086a_b.tagoff/25)), 'tagoff')


%% 
time = 1:length(STATIC_ACCEL_x);

falling = OBDA(pressure_diff < 0);  %falling = falling(falling < 1);
rising = OBDA(pressure_diff > 0); %rising = rising(rising < 1);
falling_pitch = P_dec(pressure_diff < 0);  %falling_pitch = falling_pitch(falling_pitch < 1);
rising_pitch = P_dec(pressure_diff > 0); %rising_pitch = rising_pitch(rising_pitch < 1);

length(falling)/length(rising)
rising_obda = mean(rising)
falling_obda = mean(falling)

%% FFT Filter Data: Day vs Night

[A_x_x, A_x_y] = FFT_data_set(A_x, 25);
[A_y_x, A_y_y] = FFT_data_set(A_y, 25);

figure;
subplot(211);
plot(A_x_x,A_x_y,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of X')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 6 0 5e-3])

subplot(212);
plot(A_y_x,A_y_y,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of Y')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 6 0 5e-3])

%% PSD tst
Fs = 25;
x = OBDA(1:(floor(length(OBDA)/2))*2);
%figure; plot(x)

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
%figure; plot(x)

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
%figure; plot(x)

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


%% Pitch Analysis
figure; subplot(311);plot(P_dec)
subplot(312); plot(falling_pitch);
subplot(313); plot(rising_pitch);
mean(falling_pitch)
mean(rising_pitch)

%% Segmenting data
ZERO = zeros(1,length(STATIC_ACCEL_x));

%Limit, time, pressure, p_time
UPPER_LIMIT = mean(abs(dynamics.s_x))*0.80; %Limit to be considered arms first
LOWER_LIMIT = mean(abs(dynamics.s_x))*0.60;

TIME = 1:length(STATIC_ACCEL_x);

%Set up divisions; up, down, heads, tails, change
up_idx = pressure_diff >= 0;
down_idx = pressure_diff < 0;
h_idx = (STATIC_ACCEL_x > UPPER_LIMIT);
t_idx = (STATIC_ACCEL_x < LOWER_LIMIT);
c_idx = ((STATIC_ACCEL_x <= UPPER_LIMIT) & (STATIC_ACCEL_x >= LOWER_LIMIT)) ;

% heads_up = up_idx & h_idx;
% heads_down = down_idx & h_idx;
% tails_up = up_idx & t_idx;
% tails_down = down_idx & t_idx; 
% total = sum(heads_up)+sum(heads_down)+sum(tails_up)+sum(tails_down);

% heads_up_total = sum(heads_up)/total
% heads_down_total = sum(heads_down)/total
% tails_up_total = sum(tails_up)/total
% tails_down_total = sum(tails_down)/total

%Plot 
% figure; hold on;
%plot(1:length(STATIC_ACCEL_x),ZERO, 'black');
% plot(TIME(t_idx),STATIC_ACCEL_x(t_idx), '.r', TIME(h_idx), STATIC_ACCEL_x(h_idx), '.b', TIME(c_idx), STATIC_ACCEL_x(c_idx), '.g') 
% plot(p_time(up_idx), pressure_normalized(up_idx), '.m', p_time(down_idx), pressure_normalized(down_idx), '.c');
% figure; hold on;
% plot(TIME(heads_up), STATIC_ACCEL_x(heads_up), '.r', TIME(heads_down), STATIC_ACCEL_x(heads_down), '.b')
% plot(TIME(tails_up), STATIC_ACCEL_x(tails_up), '.g', TIME(tails_down), STATIC_ACCEL_x(tails_down), '.m')
% plot(1:length(STATIC_ACCEL_x),ZERO, 'black');

%OBDA divisions
% total = length(STATIC_ACCEL_x);
% time_heads = sum(h_idx)/total
% time_tails = sum(t_idx)/total
% heads_to_tails = time_heads/time_tails
% 
% OBDA_heads_up = mean(OBDA(heads_up))
% OBDA_heads_down = mean(OBDA(heads_down))
% OBDA_tails_up = mean(OBDA(tails_up))
% OBDA_tails_down = mean(OBDA(tails_down))

%% Flap vs Finning
% Uncomment to recalculate FIN/FLAP IDX
% Somewhat time consuming, so avioid re-calculating
%close all;

mean_obda = mean(OBDA);
low_limiter = mean_obda*1.5;
med_limiter = mean_obda*3;
high_limiter = mean_obda*5;
window_size = 25;
obda_process = zeros(1, length(OBDA))';
fqs = zeros(1, length(OBDA));

smooth_obda = moving(OBDA,5);

step = 10;
for i=1:step:(length(obda_process)-window_size)
    [x, y] = FFT_data_set(OBDA(i:i+window_size), 25);
%     diff_y = diff(y);
%     y_pos = diff_y > 0; y_pos = [0 y_pos]; y_pos = (y_pos == 1);
%     x_sub_pos = x(y_pos);
% 
%     if(~isempty(x_sub_pos))
%         dom_fq = x_sub_pos(1);
%     else
%         dom_fq = 0;
%     end
    [pks, locs] = max(y(2:length(y)));
    dom_fq = x(locs+1);
    fqs(i:i+step)=dom_fq;
    fqs(i:i+window_size)=dom_fq;
    peak = max(findpeaks(OBDA(i:i+window_size)));
    
    if((peak > low_limiter) & (peak < med_limiter))% & (dom_fq > y(1)))
        obda_process(i:i+step) = 1;
        %obda_process(i) = 1;
        1;
    elseif((peak > med_limiter) & (peak < high_limiter))
        obda_process(i:i+step) = 2;
        %obda_process(i) = 1;
        1;
    elseif (peak > high_limiter)
        if(sum(h_idx(i+floor(window_size/2):i+window_size)))
            obda_process(i:i+step) = 2;
            %obda_process(i) = 1;
            1;
        else
            obda_process(i:i+step) = 3;
            %obda_process(i) = 2;
            2;

        end
    else
        obda_process(i:i+step) = 0;
        %obda_process(i) = 0;
        0;
    end
end

fin_idx = (obda_process == 0);
fin_flap_idx = (obda_process == 1);
flap_idx = (obda_process == 2);
jet_idx = (obda_process == 3);

idx = struct('head', h_idx, 'tail', t_idx, 'change', c_idx, 'up', up_idx, 'down', down_idx, 'fin', fin_idx, ...
    'flap', flap_idx, 'jet', jet_idx, 'fin_flap', fin_flap_idx);

%% Generage histogram of peak heights for day and night
%   computationally intensive, so preform as little as possible

[window, binnumber, loc_diff] = Histogram(OBDA);
[window_fin, binnumber_fin, loc_diff_fin] = Histogram(OBDA(fin_idx));
[window_fin_flap, binnumber_fin_flap, loc_diff_fin_flap] = Histogram(OBDA(fin_flap_idx));
[window_flap, binnumber_flap, loc_diff_flap] = Histogram(OBDA(flap_idx));
[window_jet, binnumber_jet, loc_diff_jet] = Histogram(OBDA(jet_idx));
save lf086a_hist 
%% Plot the histogram

%load lf086a_hist

h = figure;
subplot(321)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
loc_std = zeros(round(2/window),1);
bar(binnumber,loc_diff,'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Overall OBDA')

subplot(322)
loc_std = zeros(round(2/window_fin),1);
bar(binnumber_fin,loc_diff_fin,'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Fin')

subplot(323)
loc_std = zeros(round(2/window_jet),1);
bar(binnumber_jet,loc_diff_jet,'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Jet')

subplot(324)
loc_std = zeros(round(2/window_flap),1);
bar(binnumber_flap,loc_diff_flap,'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Flap')

subplot(325)
loc_std = zeros(round(2/window_fin_flap),1);
bar(binnumber_fin_flap,loc_diff_fin_flap,'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Fin_flap')
%% GOPRO001.MP4 lf086a
close all;
film_start = 61318-(6*60+30.7-0.159);
M_S = 5;
S_S = 30;
M_E = 7;
S_E = 15;

Segment_plotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, fqs, mean_obda);

%% GOPRO002.MP4 lf086a
close all;
film_start = 12*60+08-4;
M_S = 16;
S_S = 38;
M_E = 16;
S_E = 57;

Segment_plotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, fqs, mean_obda);

%% GP010002.MP4 lf086a
close all;
film_start = 29*60+33-4;
M_S = 0;
S_S = 00;
M_E = 29;
S_E = 25;

Segment_plotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, fqs, mean_obda);

%% CRX-700 0036 lf086a, creat animation
%close all;
%Start time and end time, minutes, second, in 00036.MTS
film_start = 61227-457-6.2;
% M_S = 10;
% S_S = 30;
% M_E = 12;
% S_E = 30;

%For trial
M_S = 7;
S_S = 30;
M_E = 9;
S_E = 30;

Segment_plotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx);

%% Tails vs Heads FFT]
close all;
A_mag = (A_dec(:,1).^2+A_dec(:,2).^2+A_dec(:,3).^2).^0.5-1;
[heads_x, heads_y] = FFT_data_set(STATIC_ACCEL_y, 25);
[tails_x, tails_y] = FFT_data_set(A_mag, 25);
%figure; plot(OBDA)

figure;
subplot(211);
plot(heads_x, heads_y,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of up OBDA')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 5 0 .5e-2])

subplot(212);
plot(tails_x, tails_y,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of down OBDA')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([ 5 0 0.1e-3])


%% CRX-700 0034 lf086a

film_start = 61227-(7*60+36);
M_S = 1;
S_S = 50;
M_E = 2;
S_E = 30;

START = (M_S*60+S_S)+film_start;
END = (M_E*60+S_E)+film_start;
S_SAMP = START*25;
E_SAMP = END*25;
figure; hold on;

ZERO = zeros(1,END-START+1);

%axis([0 END-START -0.5 0.5])

%plot(time(START:END), STATIC_ACCEL_x(START:END));%, time, pressure_smooth, time, DYNAMIC_ACCEL_x);
TIME = (1:(E_SAMP-S_SAMP+1))/25';
STATIC_SHORT = STATIC_ACCEL_x(S_SAMP:E_SAMP);
pressure_short = pressure_diff(S_SAMP:E_SAMP);
OBDA_SHORT = OBDA(S_SAMP:E_SAMP);
h_idx = (STATIC_SHORT < UPPER_LIMIT);
t_idx = (STATIC_SHORT > UPPER_LIMIT);
c_idx = (STATIC_SHORT < UPPER_LIMIT) & (STATIC_SHORT >= 0) ;
plot(TIME(t_idx),STATIC_SHORT(t_idx), '.r', TIME(h_idx), STATIC_SHORT(h_idx), '.b', TIME(c_idx), STATIC_SHORT(c_idx), '.g') 

plot(0:END-START,ZERO, 'black');





