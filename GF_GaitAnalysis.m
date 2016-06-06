%% Load Data Files
file_prefix = '86';
load(strcat('lf0', file_prefix, 'a_decimated'),'A_dec', 'pressure',  'P_dec',  'R_dec'); %M_dec
%load('testprh200.mat','A_dec', 'pressure',  'P_dec',  'R_dec', 'M_dec'); %M_dec

figure;hold on;
plot(A_dec); plot(pressure, 'c')

if strcmp(file_prefix,'77')
    pressure = pressure(4.868e4:length(pressure));
    A_dec = A_dec(4.868e4:length(A_dec), :);
    R_dec = R_dec(4.868e4:length(R_dec));
    P_dec = P_dec(4.868e4:length(P_dec));
end

if strcmp(file_prefix,'79')
    pressure = pressure(6904:1.199e6);
    A_dec = A_dec(6904:1.199e6, :);
    R_dec = R_dec(6904:1.199e6);
    P_dec = P_dec(6904:1.199e6);
end

if strcmp(file_prefix,'80')
    pressure = pressure(3.92e4:1.590e6);
    A_dec = A_dec(3.92e4:1.590e6, :);
    R_dec = R_dec(3.92e4:1.590e6);
    P_dec = P_dec(3.92e4:1.590e6);
end

if strcmp(file_prefix,'82')
    pressure = pressure(5.679e4:1.95e6);
    A_dec = A_dec(5.679e4:1.95e6, :);
    R_dec = R_dec(5.679e4:1.95e6);
    P_dec = P_dec(5.679e4:1.95e6);
end

if strcmp(file_prefix,'86')
%     pressure = pressure(5.224e4:1.859e6);
%     A_dec = A_dec(5.224e4:1.859e6, :);
%     R_dec = R_dec(5.224e4:1.859e6);
%     P_dec = P_dec(5.224e4:1.859e6);
end

if strcmp(file_prefix,'89')
    pressure = pressure(1.632e5:2.077e6);
    A_dec = A_dec(1.632e5:2.077e6, :);
    R_dec = R_dec(1.632e5:2.077e6);
    P_dec = P_dec(1.632e5:2.077e6);
end

if strcmp(file_prefix,'90')
%     pressure = setdiff(pressure(1:length(pressure)), pressure(1.186e6:1.2085e6), 'stable');
%     A_dec = setdiff(A_dec(1:length(A_dec), :), A_dec(1.186e6:1.2085e6, :), 'stable');
%     R_dec = setdiff(R_dec(1:length(R_dec)), R_dec(1.186e6:1.2085e6), 'stable');
%     P_dec = setdiff(P_dec(1:length(P_dec)), P_dec(1.186e6:1.2085e6), 'stable');
%     pressure = pressure(1:1.186e6);
%     A_dec = A_dec(1.632e5:2.077e6, :);
%     R_dec = R_dec(1.632e5:2.077e6);
%     P_dec = P_dec(1.632e5:2.077e6);
end

%A_corrected = Rotation_correct(A_dec, P_dec, R_dec);

%% Split Data; Change FIND_FLAG to 1 to find the optimal OBDA window
A_x = A_dec(:, 1);
A_y = A_dec(:, 2);
A_z = A_dec(:, 3);
time = 1:length(A_x);
time = time./25; time = time';
figure;hold on;
plot(A_dec); plot(pressure, 'c')

figure; hold on;
plot(A_x, 'b'); plot(A_y, 'g'); plot(A_z, 'r'); plot(light, 'y')
FIND_FLAG = 0;
if FIND_FLAG == 1
    running_average = OBDA_window_finder( A_x, 1, 250 );
    windows = 1:length(running_average);
    figure; plot(running_average, 'LineWidth',2);
    hold on;
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

% figure; plot(time, A_dec, time, DYNAMIC_TOTAL, time, STATIC_TOTAL)

OBDA = DYNAMIC_ACCEL_x+DYNAMIC_ACCEL_y+DYNAMIC_ACCEL_z;
dynamics = struct('obda', OBDA, 's_x', STATIC_ACCEL_x, 's_y', STATIC_ACCEL_y, 's_z', STATIC_ACCEL_z, ...
       'd_x', DYNAMIC_ACCEL_x, 'd_y', DYNAMIC_ACCEL_y, 'd_z', DYNAMIC_ACCEL_z, 'a_x', A_x, 'a_y', A_y, 'a_z', A_z);

% figure; hold on
% plot(time, A_dec, time, OBDA)


%% Calculate Falling vs Rising Stats
time = 1:length(STATIC_ACCEL_x);

falling = OBDA(pressure_diff < 0);  %falling = falling(falling < 1);
rising = OBDA(pressure_diff > 0); %rising = rising(rising < 1);
falling_pitch = nanmean(P_dec(pressure_diff < 0));  %falling_pitch = falling_pitch(falling_pitch < 1);
rising_pitch = nanmean(P_dec(pressure_diff > 0)); %rising_pitch = rising_pitch(rising_pitch < 1);

ratio_falling_rising = length(falling)/length(rising);
rising_obda = nanmean(rising);
falling_obda = nanmean(falling);


%% Create Mean vector
window_len = 25*60;
temp_static = dynamics.s_x + abs(min(dynamics.s_x));
mean_vector = zeros(1, length(temp_static));
for i=1:window_len
   mean_vector(i) = mean(temp_static(i:i+window_len)); 
end
for i=window_len+1:length(temp_static)-window_len-1
    mean_vector(i) = mean(temp_static(i-window_len:i+window_len));
    
end
for i=length(temp_static)-window_len:length(temp_static)
   mean_vector(i) = mean(temp_static(i-window_len:i)); 
end

mean_vector = mean_vector';

%% Segmenting data
ZERO = zeros(1,length(STATIC_ACCEL_x));

%Limit, time, pressure, p_time
% UPPER_LIMIT = nanmean(abs(dynamics.s_x))*0.80; %Limit to be considered arms first
% LOWER_LIMIT = nanmean(abs(dynamics.s_x))*0.60;
% temp_static = STATIC_ACCEL_x;

UPPER_LIMIT = mean_vector*1.00;
LOWER_LIMIT = mean_vector*0.95;
TIME = 1:length(STATIC_ACCEL_x);

%Set up divisions; up, down, heads, tails, change
up_idx = pressure_diff >= 0;
down_idx = pressure_diff < 0;

h_idx = (temp_static(1:length(pressure_diff)) > UPPER_LIMIT(1:length(pressure_diff)));
t_idx = (temp_static(1:length(pressure_diff)) < LOWER_LIMIT(1:length(pressure_diff)));
c_idx = ~(h_idx | t_idx);
%c_idx = ((STATIC_ACCEL_x(1:length(pressure_diff)) <= UPPER_LIMIT(1:length(pressure_diff))) ...
%& (STATIC_ACCEL_x(1:length(pressure_diff)) >= LOWER_LIMIT(1:length(pressure_diff)))) ;

heads_up = up_idx & h_idx;
heads_down = down_idx & h_idx;
tails_up = up_idx & t_idx;
tails_down = down_idx & t_idx; 
total = sum(heads_up)+sum(heads_down)+sum(tails_up)+sum(tails_down);

heads_up_ratio = sum(heads_up)/total;
heads_down_ratio = sum(heads_down)/total;
tails_up_ratio = sum(tails_up)/total;
tails_down_ratio = sum(tails_down)/total;

%OBDA divisions
total = length(STATIC_ACCEL_x);
time_heads = sum(h_idx)/total;
time_tails = sum(t_idx)/total;
heads_to_tails = time_heads/time_tails;

OBDA_heads_up = nanmean(OBDA(heads_up));
OBDA_heads_down = nanmean(OBDA(heads_down));
OBDA_tails_up = nanmean(OBDA(tails_up));
OBDA_tails_down = nanmean(OBDA(tails_down));
OBDA_heads = nanmean(OBDA(h_idx));
OBDA_tails = nanmean(OBDA(t_idx));

f = figure('Position',[200 200 400 400]);
%f = figure('units','normalized','outerposition',[0 0 1 1]);
dat = [ratio_falling_rising; rising_obda; falling_obda; rising_pitch; falling_pitch; heads_up_ratio; heads_down_ratio; ...
    tails_up_ratio; tails_down_ratio; heads_to_tails; OBDA_heads_up; OBDA_heads_down; OBDA_tails_up; OBDA_tails_down; OBDA_heads; ...
    OBDA_tails];
cnames = {'Value'};
rnames = {'Ratio Falling/Rising','Mean Rising OBDA','Mean Falling OBDA', 'Mean Rising Pitch', 'Mean Falling Pitch' ...
    'Ratio Heads Up', 'Ratio Heads Down', 'Ratio Tails Up', 'Ratio Tails Down', 'Ratio Heads to Tails', 'Mean OBDA heads up' ...
    'Mean OBDA heads down', 'Mean OBDA tails up', 'Mean OBDA tails down', 'Mean OBDA heads', 'Mean OBDA tails'};
t = uitable('Parent',f,'Data',dat,'ColumnName',cnames,... 
            'RowName',rnames ,'Position',[10 10 360 360]);
set(t,'ColumnWidth',{100})



summary_stats = struct('ratio_falling_rising', ratio_falling_rising, 'rising_obda', rising_obda, 'falling_obda',  falling_obda,...
    'rising_pitch', rising_pitch, 'falling_pitch', falling_pitch, 'heads_up_ratio', heads_up_ratio, 'heads_down_ratio', heads_down_ratio, ...
    'tails_up_ratio', tails_up_ratio, 'tails_down_ratio', tails_down_ratio, 'heads_to_tails',  heads_to_tails, 'OBDA_heads_up',  OBDA_heads_up,...
    'OBDA_heads_down', OBDA_heads_down, 'OBDA_tails_up', OBDA_tails_up, 'OBDA_tails_down', OBDA_tails_down, 'OBDA_heads', OBDA_heads, ...
    'OBDA_tails', OBDA_tails); 
%% FFT Filter Data: Day vs Night

[A_x_x, A_x_y] = GF_FftDataSet(OBDA, 25);
[A_y_x, A_y_y] = GF_FftDataSet(A_x, 25);

figure; subplot(211);
plot(A_x_x,A_x_y,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of OBDA')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 4 0 2e-3])

subplot(212);
plot(A_y_x,A_y_y,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of X')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 4 0 2e-3])

%% Analyze finning frequency vs vertical velocity
window_size = 100;
fqs = zeros(1, length(OBDA)-1);
for i=1:window_size:(length(OBDA)-window_size)
     [x, y] = FFT_data_set(OBDA(i:i+window_size), 25);
     diff_y = diff(y);

    [pks, locs] = max(y(2:length(y)));
    dom_fq = x(locs+1);
    fqs(i)=dom_fq;
end
vert_vel= diff(pressure_normalized); vert_vel = [0; vert_vel];
figure;plot(fqs(fqs ~= 0), vert_vel(fqs ~= 0), '.b'); axis([0 2.5 -2e-2 2e-2]);

 %save lf078a_obda OBDA summary_stats_lf078a
save(strcat('lf0', file_prefix, 'a_obda'),'OBDA', 'summary_stats');
%% Flap vs Finning
% Uncomment to recalculate FIN/FLAP IDX
% Somewhat time consuming, so avioid re-calculating

mean_obda = nanmean(OBDA);
low_limiter = mean_obda*1.5;
med_limiter = mean_obda*3;
high_limiter = mean_obda*5;
window_size = 25;
obda_process = zeros(1, length(OBDA))';
smooth_obda = moving(OBDA,5);

for i=1:window_size:(length(obda_process)-window_size)
    peak = max(findpeaks(OBDA(i:i+window_size)));
    %peak_smooth = max(findpeaks(smooth_obda(i:i+window_size)));
    
    if((peak >= low_limiter) & (peak < med_limiter))% & (dom_fq > y(1)))
        obda_process(i:i+window_size) = 1; %Slow_flap
    elseif((peak >= med_limiter) & (peak < high_limiter))
        obda_process(i:i+window_size) = 2; %Fast flap
    elseif (peak > high_limiter)
        if(sum(h_idx(i+floor(window_size/2):i+window_size)))
            obda_process(i:i+window_size) = 2; % Fast flap
        else
            obda_process(i:i+window_size) = 3; % Jet
        end
    else
        obda_process(i:i+window_size) = 0; %Fin
    end
end

fin_idx = (obda_process == 0);
slow_flap_idx = (obda_process == 1);
fast_flap_idx = (obda_process == 2);
jet_idx = (obda_process == 3);

idx = struct('head', h_idx, 'tail', t_idx, 'change', c_idx, 'up', up_idx, 'down', down_idx, 'fin', fin_idx, ...
    'flap', fast_flap_idx, 'jet', jet_idx, 'fin_flap', slow_flap_idx);
%% Plot Data
figure; hold on;
plot(TIME(fin_idx), OBDA(fin_idx), '.',  'Color',[   0.4000    0.4000    1.0000], 'MarkerFaceColor', ...
    [   0.4000    0.4000    1.0000], 'Markersize', 15);
plot(TIME(slow_flap_idx), OBDA(slow_flap_idx), '.', 'Color',[ 0.2000    0.8000         0], 'MarkerFaceColor', ...
    [ 0.2000    0.8000         0], 'Markersize', 15);
plot(TIME(fast_flap_idx), OBDA(fast_flap_idx), '.', 'Color',[  0.9255    0.9255    0.0745], 'MarkerFaceColor',...
    [  0.9255    0.9255    0.0745], 'Markersize', 15);
plot(TIME(jet_idx), OBDA(jet_idx), '.', 'Color',[ 0.9412         0         0],  'MarkerFaceColor',[ 0.9412         0         0],...
    'Markersize', 15);

pressure_lowered = pressure_normalized - 0.7;
plot(TIME, OBDA, 'black');
plot(TIME(up_idx), pressure_lowered(up_idx), '.m', TIME(down_idx), pressure_lowered(down_idx), '.c');
plot(TIME(t_idx),STATIC_ACCEL_x(t_idx), '.b', TIME(h_idx), STATIC_ACCEL_x(h_idx), '.r', TIME(c_idx), STATIC_ACCEL_x(c_idx), '.g') ;
%plot(TIME, UPPER_LIMIT, 'k');
%plot(TIME, LOWER_LIMIT, 'y');
save(strcat('lf0', file_prefix, 'a_complete'));
