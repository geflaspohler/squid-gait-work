%% Assemble data matrixes from files
% %User file input
file_prefix = 'lf086a';
file_suffix = 'prh.mat';
file_number = 10;

%Load files and concatante into large matrixes A_all, P_all, R_all, H_all
%Concat_files(file_prefix, file_suffix, file_number);
[A_all, P_all, R_all, H_all, pressure, M_all] = Concat_files(file_prefix, file_suffix, file_number);
% 
% %Decimate data to make for faster processing
A_dec = decdc(A_all,10);
P_dec = decdc(P_all,10);
R_dec = decdc(R_all,10);
H_dec = decdc(H_all,10);
M_dec = decdc(M_all, 10);
pressure = decdc(pressure,10);

save lf086a_decimated

%% Plot decimated data
figure;
subplot(221)
plot(A_dec)
subplot(222)
plot(P_dec)
subplot(223)
plot(R_dec)
subplot(224)
plot(H_dec)

%% Enter event time data
%load data_decimated
%load lf086a_decimated

conv_min_to_samp = 60*25; % Converst from min to samples

%For lf080a
% a = 1; %Tag Armed
% b = 38635; % Animal in the water 18:03
% c = 223500; % Sunset
% d = 1267500; % Sunrise
% e = 461927; % Change in Direction
% l = 1590769; % Tag off

%For lf086a
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
lf082a_b = struct('sunset', c, 'sunrise', d, 'jetting1', 82668*25, 'jetting2', 62668.288*25, 'forwards_finning',  82775.729*25, 'reversal', 82949.336*25);
lf080a_b = struct('sunset', c, 'sunrise', d, 'jetting1', 62284.472*25);

%A_corrected = Rotation_correct(A_dec, P_dec, R_dec);


%User input event times in minutes. 
%Format: event_in_miunutes = [animal_released sunset sunrise tag_removed]
event_in_miunutes = [40 162 858 1101];

%Convert minutes to samples
%For lf086a
event_in_samples = event_in_miunutes.*conv_min_to_samp;

%For lf080a
%event_in_samples = [b c d l];

%A_dec = A_corrected;
%Find the +/- magnitude of the three axis of accelerometer data
A_mag = (A_dec(:,1).^2+A_dec(:,2).^2+A_dec(:,3).^2).^0.5-1;
A_x = A_dec(:,1);
A_y = A_dec(:,2);
A_z = A_dec(:,3);

%Indicies used to access event matrix
ANIMAL_RELEASED = 1;
SUNSET = 2;
SUNRISE = 3;
TAG_REMOVED = 4;

%Create day and night matrixes
Afternoon = A_mag(event_in_samples(ANIMAL_RELEASED):event_in_samples(SUNSET));
Morning = A_mag(event_in_samples(SUNRISE):event_in_samples(TAG_REMOVED));
Night = A_mag(event_in_samples(SUNSET):event_in_samples(SUNRISE));

Afternoon_x = A_x(event_in_samples(ANIMAL_RELEASED):event_in_samples(SUNSET));
Morning_x = A_x(event_in_samples(SUNRISE):event_in_samples(TAG_REMOVED));
Night_x = A_x(event_in_samples(SUNSET):event_in_samples(SUNRISE));

Afternoon_y = A_y(event_in_samples(ANIMAL_RELEASED):event_in_samples(SUNSET));
Morning_y = A_y(event_in_samples(SUNRISE):event_in_samples(TAG_REMOVED));
Night_y = A_y(event_in_samples(SUNSET):event_in_samples(SUNRISE));

Afternoon_z = A_z(event_in_samples(ANIMAL_RELEASED):event_in_samples(SUNSET));
Morning_z = A_z(event_in_samples(SUNRISE):event_in_samples(TAG_REMOVED));
Night_z = A_z(event_in_samples(SUNSET):event_in_samples(SUNRISE));

Afternoon_P = P_dec(event_in_samples(ANIMAL_RELEASED):event_in_samples(SUNSET));
Morning_P = P_dec(event_in_samples(SUNRISE):event_in_samples(TAG_REMOVED));
Night_P = P_dec(event_in_samples(SUNSET):event_in_samples(SUNRISE));

Afternoon_R = R_dec(event_in_samples(ANIMAL_RELEASED):event_in_samples(SUNSET));
Morning_R = R_dec(event_in_samples(SUNRISE):event_in_samples(TAG_REMOVED));
Night_R = R_dec(event_in_samples(SUNSET):event_in_samples(SUNRISE));

Afternoon_H = H_dec(event_in_samples(ANIMAL_RELEASED):event_in_samples(SUNSET));
Morning_H = H_dec(event_in_samples(SUNRISE):event_in_samples(TAG_REMOVED));
Night_H = H_dec(event_in_samples(SUNSET):event_in_samples(SUNRISE));

%save lf086a_data_event_times

%%

figure; hold on;
plot(time, A_mag, time, light);
plot(c/25, Combined(c/25), '.r')
plot(d/25, Combined(d/25), '.r')

text(c, Combined(c), 'sunset')
text(d, Combined(d), 'sunrise')
%% Plot Movements
figure; hold on;
zero = zeros(1, length(OBDA));
plot(time, STATIC_ACCEL_x, time, pressure_normalized, time, OBDA, time, zero)
plot(time, A_x, time, A_y, time, A_z)
Combined = A_dec;

offset = b;

%Uncomment for lf080a
lf080a_b.jetting1 = lf080a_b.jetting1-offset;
plot(floor(lf080a_b.jetting1/25),  Combined(floor(lf080a_b.jetting1/25)), '.r')
text(floor(lf080a_b.jetting1/25), Combined(floor(lf080a_b.jetting1/25)), 'jetting1')

%Uncomment for lf086a
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

%% Exlude areas of artifical activity
%load lf086a_data_event_times

%User inputed begenning and end of artificial probing (samples) for lf086a
morning_ex_s = 2.448e5;
morning_ex_e = 2.877e5;
afternoon_ex_s = 1.184e5;
afternoon_ex_e = 1.787e5;

Morning_with = Morning;
Morning = [Morning(1:morning_ex_s);Morning(morning_ex_e:length(Morning))];
Morning_x = [Morning_x(1:morning_ex_s);Morning_x(morning_ex_e:length(Morning_x))];
Morning_y = [Morning_y(1:morning_ex_s);Morning_y(morning_ex_e:length(Morning_y))];
Morning_z = [Morning_z(1:morning_ex_s);Morning_z(morning_ex_e:length(Morning_z))];
Morning_P = [Morning_P(1:morning_ex_s);Morning_P(morning_ex_e:length(Morning_P))];
Morning_R = [Morning_R(1:morning_ex_s);Morning_R(morning_ex_e:length(Morning_R))];
Morning_H = [Morning_H(1:morning_ex_s);Morning_H(morning_ex_e:length(Morning_H))];

Afternoon_with = Afternoon;
Afternoon = [Afternoon(1:afternoon_ex_s);Afternoon(afternoon_ex_e:length(Afternoon))];
Afternoon_x = [Afternoon_x(1:afternoon_ex_s);Afternoon_x(afternoon_ex_e:length(Afternoon_x))];
Afternoon_y = [Afternoon_y(1:afternoon_ex_s);Afternoon_y(afternoon_ex_e:length(Afternoon_y))];
Afternoon_z = [Afternoon_z(1:afternoon_ex_s);Afternoon_z(afternoon_ex_e:length(Afternoon_z))];
Afternoon_P = [Afternoon_P(1:afternoon_ex_s);Afternoon_P(afternoon_ex_e:length(Afternoon_P))];
Afternoon_R = [Afternoon_R(1:afternoon_ex_s);Afternoon_R(afternoon_ex_e:length(Afternoon_R))];
Afternoon_H = [Afternoon_H(1:afternoon_ex_s);Afternoon_H(afternoon_ex_e:length(Afternoon_H))];

Day_f = [Morning; Afternoon];
Day_x = [Morning_x; Afternoon_x];
Day_y = [Morning_y; Afternoon_y];
Day_z = [Morning_z; Afternoon_z];
Day_P = [Morning_P; Afternoon_P];
Day_R = [Morning_R; Afternoon_R];
Day_H = [Morning_H; Afternoon_H];

Night_f = Night;

Day_nf = [Morning_with; Afternoon_with];
Night_nf = Night;

time_day = ((((1:length(Day_f))/25)/60)/60)';
time_night = ((((1:length(Night_f))/25)/60)/60)';

figure; subplot(221), plot(Day_f), axis([0 4e5 0 1]), title('Day Filtered')
subplot(223), plot(Night_f), axis([0 4e5 0 1]), title('Night Filtered')
subplot(222), plot(Day_nf), axis([0 4e5 0 1]), title('Day Not Filtered')
subplot(224), plot(Night_nf), axis([0 4e5 0 1]), title('Night Not Filtered')

%save lf086a_data_artifical_removed

%% %% Find Filtered Day vs Night RMS

Day_rms_f = sqrt((sum(Day_f.^2))/length(Day_f));
Night_rms_f = sqrt((sum(Night_f.^2))/length(Night_f));
Day_rms_nf = sqrt((sum(Day_nf.^2))/length(Day_nf));
Night_rms_nf = sqrt((sum(Night_nf.^2))/length(Night_nf));

h = figure;
set(h,'Color',[1,1,1], 'Position', [1 1 840 330])

subplot(212)
barweb([Night_rms_f Day_rms_f], [0 0], [], [], [], [], [], bone, [], [], 1, [])
axis([.5 1.5 0 0.02])
xlabel('Night vs Day Filtered')
ylabel('A RMS (g)')

subplot(211)
barweb([Night_rms_nf Day_rms_nf], [0 0], [], [], [], [], [], bone, [], [], 1, [])
axis([.5 1.5 0 0.02])
xlabel('Night vs Day Not Filtered')
ylabel('A RMS (g)')

%% FFT Filter Data: Day vs Night

[X_Night, Y_Night] = FFT_data_set(Night_f, 25);

figure;
subplot(211);
plot(X_Night,Y_Night,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of Night')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 10 0 1e-3])

[X_Day, Y_Day] = FFT_data_set(Day_f, 25);

subplot(212)
plot(X_Day,Y_Day,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of Day')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 10 0 1e-3])

%% Generage histogram of peak heights for day and night
%   computationally intensive, so preform as little as possible

[window_day, binnumber_day, loc_diff_day] = Histogram(Day_f);
[window_night, binnumber_night, loc_diff_night] = Histogram(Night_f);
%save lf086a_hist window_day window_night binnumber_day binnumber_night loc_diff_day loc_diff_night

%% Plot the histogram

%load lf086a_hist

h = figure;
subplot(221)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
loc_std = zeros(round(2/window_day),1);
bar(binnumber_day,loc_diff_day,'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Day')

subplot(222)
loc_std = zeros(round(2/window_day),1);
bar(binnumber_day,loc_diff_day*(length(Night_f)/length(Day_f)),'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Day Normalized')

subplot(223)
loc_std = zeros(round(2/window_night),1);
bar(binnumber_night,loc_diff_night,'DisplayName','loc_diff');
figure(gcf)
axis([0 .2 0 7000])
xlabel('Bins')
ylabel('# of Peaks')
title('Night')

%% Begin finning discovery 
% Low pass filter the data
d = fdesign.lowpass('Fp,Fst,Ap,Ast',5,10,0.5,40,100);
Hd = design(d,'equiripple');

Day_f = filter(Hd,Day_f);
Day_x = filter(Hd,Day_x);
Day_y = filter(Hd,Day_y);
Day_z = filter(Hd,Day_z);
Day_P = filter(Hd,Day_P);
Day_R = filter(Hd,Day_R);
Day_H = filter(Hd,Day_H);

Night_f = filter(Hd,Night_f);
Night_x = filter(Hd, Night_x);
Night_y = filter(Hd, Night_y);
Night_z = filter(Hd, Night_z);
Night_P = filter(Hd, Night_P);
Night_R = filter(Hd, Night_R);
Night_H = filter(Hd, Night_H);

info(Hd)                       % View information about filter
fvtool(Hd)

%save lf086a_filtered_data
%% Find night/day areas of high activity
frame_size = 25;
seconds_spent_not_finning_night = 0;
seconds_spent_not_finning_day = 0;
[ Day_totals, Day_peaks ] = Sum_frames( frame_size, Day_f );
[ Night_totals, Night_peaks] = Sum_frames( frame_size, Night_f );


figure; 
subplot(221)
plot(Day_totals), title('Day Totals'), axis([0 5e5 0 1])
subplot(223)
plot(Night_totals), title('Night Totals'), axis([0 10e5 0 1])
subplot(222)
plot(Day_f), title('Day filtered'), axis([0 5e5 0 1])
subplot(224)
plot(Night_f), title('Night filtered'), axis([0 10e5 0 1])

%save lf086a_peaks

%% Remove large events for finning finding
%load lf086a_peaks
%User set threshold for what constitudes an area of high
%activity and should be removed when considering finning
min_threshold = 0.1;

figure;
hold on
plot(Night_f, 'black'), title('Night'),axis([0 length(Night_f) -1 1]);
Night_changed = Night_f;
for i=1:length(Night_totals)
    if abs(Night_totals(i)) > min_threshold
       seconds_spent_not_finning_night = seconds_spent_not_finning_night+1;
       plot(i-frame_size:i, Night_f(i-frame_size:i), '.r', 'markersize', 10), title('Night totals'), axis([0 10e5 -1 1]);
       Night_changed(i-frame_size:i) = 0;
    end
end

figure;
hold on
plot(Day_f, 'black'), title('Day'),axis([0 length(Day_f) -1 1]);
Day_changed = Day_f;
for i=1:length(Day_totals)
    if abs(Day_totals(i)) > min_threshold
       seconds_spent_not_finning_day = seconds_spent_not_finning_day+1;
       plot(i-frame_size:i, Day_f(i-frame_size:i), '.r', 'markersize', 10), title('Day totals')%, axis([0 10e5 -1 1])
       Day_changed(i-frame_size:i) = 0;
    end
end

figure; plot(Day_changed), axis([0 10e5 -1 1]); 
figure; plot(Night_changed), axis([0 10e5 -1 1]);
%save lf086a_large_activity_removed

%% Find Bouts of finning
%close all
start_cut = 0.961e4; 
end_cut = 3.961e4;
%Day_cut = Day_changed(start_cut:end_cut);
%Night_cut = Night_changed(start_cut:end_cut);
Day_cut = Day_f;
Night_cut = Night_f;
day_std = std(Day_cut);
night_std = std(Night_cut);
resting_mag_day = day_std/3;
resting_mag_night = night_std/3;%-0.10*night_std;

[average_wait_night, average_peaks_night, finning_matrix_night, swimming_matrix_night, wait_matrix_night, count_matrix_night] = Gait_parse( Night_cut, resting_mag_night);

prop_swimming_night = length(swimming_matrix_night)/length(Night_cut)*100
prop_finning_night = length(finning_matrix_night)/length(Night_cut)*100
prop_finning_to_swiming_night = length(finning_matrix_night)/length(swimming_matrix_night)
figure; plot(swimming_matrix_night), title('Swimming Matrix Night');
figure; plot(finning_matrix_night), title('Finning Marix Night');

 
 [ average_wait_day, average_peaks_day, finning_matrix_day, swimming_matrix_day, wait_matrix_day, count_matrix_day] = Gait_parse( Day_cut, resting_mag_day);
 
 time_a_swim_day = ((1:length(Day_cut))/25)';
 time_a_fin_day = ((1:length(finning_matrix_day))/25)';
 prop_swimming_day = length(swimming_matrix_day)/length(Day_cut)*100;
 prop_finning_day = length(finning_matrix_day)/length(Day_cut)*100;
 prop_finning_to_swiming_day = length(finning_matrix_day)/length(swimming_matrix_day);
 figure; plot(swimming_matrix_day), title('Swimming Matrix Day');
 figure; plot(time_a_fin_day, finning_matrix_day), title('Finning Marix Day');
 

%% Plot Calculations
figure;
subplot(211)
barweb([prop_finning_night prop_finning_day ], [0 0], [0.75], [], [], [], [], bone, [], [], 1, [])
%axis([.5 1.5 0 0.02])
xlabel('Night vs Day Prop Finning')
ylabel('% of time')

subplot(212)
barweb([prop_swimming_night prop_swimming_day], [0 0], [0.75], [], [], [], [], bone, [], [], 1, [])
axis([.5 1.5 0 10])
xlabel('Night vs Day Prop Swimming')
ylabel('% of time')

%save lf086a_bout_counting

%% Load bout counting data
%load lf086a_bout_counting
count_matrix_sorted_night = sort(count_matrix_night);
count_matrix_sorted_day = sort(count_matrix_day);

wait_matrix_sorted_night = sort(wait_matrix_night);
wait_matrix_sorted_day = sort(wait_matrix_day);

%% Bin Fin bout lens
bin_len = 10;

%Setup for night histogram
max_peaks_night = ceil((max(count_matrix_night)/10))*10;
binnumber_night = linspace(bin_len,max_peaks_night, max_peaks_night/bin_len);
peaks_binned_night = zeros(max_peaks_night/bin_len, 1);

for i=1:(max_peaks_night/bin_len)
  peaks_binned_night(i) = sum((count_matrix_sorted_night <= i*bin_len) & (count_matrix_sorted_night > (i-1)*bin_len)) ;  
end

%Plot night histogram
h = figure;
subplot(232)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
bar(binnumber_night-(bin_len/2), peaks_binned_night, 1);
set(gca,'XTick',[0:max_peaks_night/bin_len]*bin_len*5) 
xlabel('Bins')
axis([0 max_peaks_night+bin_len 0 max(peaks_binned_night)+bin_len*10])
ylabel('# of Peaks/Bout')
title('Night')

%Add for numerical disaply of bins
% for b = 1:length(peaks_binned_night)
%     text(binnumber_night(b)-bin_len/2,peaks_binned_night(b),num2str(peaks_binned_night(b)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
% end

text(max_peaks_night/2,max(peaks_binned_night)/1.5,strcat('Avg pks/bout: ',num2str(average_peaks_night)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
%-------------------------------End night------------------------%

%-------------------------------Begin day------------------------%
max_peaks_day = ceil((max(count_matrix_day)/10))*10;
binnumber_day = linspace(bin_len,max_peaks_day, max_peaks_day/bin_len);
peaks_binned_day = zeros(max_peaks_day/bin_len, 1);

for i=1:(max_peaks_day/bin_len)
  peaks_binned_day(i) = sum((count_matrix_sorted_day <= i*bin_len) & (count_matrix_sorted_day > (i-1)*bin_len)) ;  
end

%Plot day histogram
subplot(231)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
bar(binnumber_day-(bin_len/2), peaks_binned_day, 1);
set(gca,'XTick',[0:max_peaks_day/bin_len]*bin_len*10) 
xlabel('Bins')
axis([0 max_peaks_day+bin_len 0 max(peaks_binned_day)+bin_len*5])
ylabel('# of Peaks/Bout')
title('Day ')

%Add for numerical disaply of bins
% for b = 1:length(peaks_binned_day)
%     text(binnumber_day(b)-bin_len/2,peaks_binned_day(b),num2str(peaks_binned_day(b)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
% end

text(max_peaks_day/2,max(peaks_binned_day)*2,strcat('Avg pks/bout: ',num2str(average_peaks_day)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
%-------------------------------End day------------------------%

%-------------------------------Begin day normalize------------%
peaks_binned_day_norm = floor(peaks_binned_day*(length(Night_cut)/length(Day_cut)));
max_peaks_day_norm = max(peaks_binned_day_norm);

subplot(233)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
bar(binnumber_day-(bin_len/2), peaks_binned_day_norm, 1);
set(gca,'XTick',[0:max_peaks_day/bin_len]*bin_len*10) 
xlabel('Bins')
axis([0 max_peaks_day+bin_len 0 max(peaks_binned_day_norm)+bin_len*5])
ylabel('# of Waits/Bout Norm')
title('Day Normalized')

%Add for numerical disaply of bins
% for b = 1:length(peaks_binned_day)
%     text(binnumber_day(b)-bin_len/2,peaks_binned_day_norm(b),num2str(peaks_binned_day_norm(b)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
% end

text(max_peaks_day/2,max(peaks_binned_day)*2,strcat('Avg pks/bout: ',num2str(average_peaks_day)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);

%--------------------------------------------------------------------%
% Bin Fin wait lens
bin_len = 1;

%Setup night wait times
max_peaks_night = ceil((max(wait_matrix_night)/bin_len))*bin_len;
binnumber_night = linspace(bin_len,max_peaks_night, max_peaks_night/bin_len);
peaks_binned_night = zeros(max_peaks_night/bin_len, 1);

for i=1:(max_peaks_night/bin_len)
  peaks_binned_night(i) = sum((wait_matrix_sorted_night <= i*bin_len) & (wait_matrix_sorted_night > (i-1)*bin_len)) ;  
end

%Plot night wait times
subplot(235);
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
bar(binnumber_night-(bin_len/2), peaks_binned_night, 1);
set(gca,'XTick',[0:max_peaks_night/bin_len]*bin_len*2) 
xlabel('Bins')
axis([0 max_peaks_night+bin_len 0 max(peaks_binned_night)+bin_len*50])
ylabel('# of Waits/Bout')
title('Night')

%Add for numerical disaply of bins
% for b = 1:length(peaks_binned_night)
%     text(binnumber_night(b)-bin_len/2,peaks_binned_night(b),num2str(peaks_binned_night(b)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
% end

text(max_peaks_night/2,max(peaks_binned_night)/1.5,strcat('Avg waits/bout: ',num2str(average_wait_night)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
%End night

%Setup day wait times
max_peaks_day = ceil((max(wait_matrix_day)/bin_len))*bin_len;
binnumber_day = linspace(bin_len,max_peaks_day, max_peaks_day/bin_len);
peaks_binned_day = zeros(max_peaks_day/bin_len, 1);

for i=1:(max_peaks_day/bin_len)
  peaks_binned_day(i) = sum((wait_matrix_sorted_day <= i*bin_len) & (wait_matrix_sorted_day > (i-1)*bin_len)) ;  
end

peaks_binned_day_norm = floor(peaks_binned_day*(length(Night_cut)/length(Day_cut)));
max_peaks_day_norm = max(peaks_binned_day_norm);

%Plot day wait tmes
subplot(234)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
bar(binnumber_day-(bin_len/2), peaks_binned_day, 1);
set(gca,'XTick',[0:max_peaks_day/bin_len]*bin_len*2) 
xlabel('Bins')
axis([0 max_peaks_day+bin_len 0 max(peaks_binned_day)+bin_len*50])
ylabel('# of Waits/Bout')
title('Day')

%Add for numerical disaply of bins
% for b = 1:length(peaks_binned_day)
%     text(binnumber_day(b)-bin_len/2,peaks_binned_day(b),num2str(peaks_binned_day(b)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
% end

text(max_peaks_day/2,max(peaks_binned_day)/1.5,strcat('Avg waits/bout: ',num2str(average_wait_day)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);

subplot(236)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
bar(binnumber_day-(bin_len/2), peaks_binned_day_norm, 1);
set(gca,'XTick',[0:max_peaks_day_norm/bin_len]*bin_len*2) 
xlabel('Bins')
axis([0 max_peaks_day+bin_len 0 max(peaks_binned_day_norm)+bin_len*50])
ylabel('# of Waits/Bout Norm')
title('Day Normalized')

%Add for numerical disaply of bins
% for b = 1:length(peaks_binned_day)
%     text(binnumber_day(b)-bin_len/2,peaks_binned_day_norm(b),num2str(peaks_binned_day_norm(b)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);
% end

text(max_peaks_day/2,max(peaks_binned_day)/1.5,strcat('Avg waits/bout: ',num2str(average_wait_day)),'VerticalAlignment','top', 'BackgroundColor', [.7 .9 .7]);

save lf080a_complete
%%
conv_min_to_samp = 60*25; % Converst from min to samples
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

%Combined = [[Night_x;Day_x] [Night_y-1; Day_y-1] [Night_z; Day_z]];
Combined = A_dec;
Combined_mag = (Combined(:,1).^2+Combined(:,2).^2+Combined(:,3).^2).^0.5-1;
figure;
subplot(211); hold on; plot(Combined)
plot(a, Combined_mag(a), '.r')
plot(b, Combined_mag(b), '.r')
plot(c, Combined_mag(c), '.r')
plot(d, Combined_mag(d), '.r')
plot(e, Combined_mag(e), '.r')
plot(f, Combined_mag(f), '.r')
plot(g, Combined_mag(g), '.r')
plot(h, Combined_mag(h), '.r')
plot(i, Combined_mag(i), '.r')
plot(j, Combined_mag(j), '.r')
plot(k, Combined_mag(k), '.r')
plot(l, Combined_mag(l), '.r')
plot(m, Combined_mag(m), '.r')
text(a, Combined_mag(a), 'a')
text(b, Combined_mag(a), 'b')
text(c, Combined_mag(a), 'c')
text(d, Combined_mag(a), 'd')
text(e, Combined_mag(a), 'e')
text(f, Combined_mag(a), 'f')
text(g, Combined_mag(a), 'g')
text(h, Combined_mag(a), 'h')
text(i, Combined_mag(a), 'i')
text(j, Combined_mag(a), 'g')
text(k, Combined_mag(a), 'h')

subplot(212); hold on; plot(Combined_mag)
plot(a, Combined_mag(a), '.r')
plot(b, Combined_mag(b), '.r')
plot(c, Combined_mag(c), '.r')
plot(d, Combined_mag(d), '.r')
plot(e, Combined_mag(e), '.r')
plot(f, Combined_mag(f), '.r')
plot(g, Combined_mag(g), '.r')
plot(h, Combined_mag(h), '.r')
plot(i, Combined_mag(i), '.r')
plot(j, Combined_mag(j), '.r')
plot(k, Combined_mag(k), '.r')
plot(l, Combined_mag(l), '.r')
plot(m, Combined_mag(m), '.r')
text(a, Combined_mag(a), 'a')
text(b, Combined_mag(a), 'b')
text(c, Combined_mag(a), 'c')
text(d, Combined_mag(a), 'd')
text(e, Combined_mag(a), 'e')
text(f, Combined_mag(a), 'f')
text(g, Combined_mag(a), 'g')
text(h, Combined_mag(a), 'h')
text(i, Combined_mag(a), 'i')
text(j, Combined_mag(a), 'g')
text(k, Combined_mag(a), 'h')


% Finning and mainting state
% Finning and translating
% Finning and jettings
% Depth and pressure; dives
% Gilly: Humbolt Squid tagging