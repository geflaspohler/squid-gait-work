
[x, y] = FFT_data_set(OBDA, 25);

film_start = 61319-(6*60+25);
mins_start = 13;
secs_start = 02;
mins_end = 14;
secs_end = 02;

start = (mins_start*60+secs_start)+film_start;
ender = (mins_end*60+secs_end)+film_start;
start_samp = start*25;
end_samp = ender*25;

OBDA_short_1 = OBDA(S_SAMP:E_SAMP);
figure; hold on;

%figure; plot(OBDA)

figure;
subplot(211)
plot(x, y,'k','linewidth',2) 
title('Single-Sided Amplitude Spectrum of up OBDA')
xlabel('Frequency (Hz)')
ylabel('|Accelleration|')
axis([0 5 0 .5e-2])


%Segement plotter plots a stylized segment of a graph, from a specified start to end
%time, within the data set.
%M_S: minute in video to start, S_S: second in video to start
%M_E: minute in video to end, S_E: second in video to end
%film_start: time (in samples) of the start of the video within the data set.
%dynamics: structure containing the acceleration dynamics data 
%pressure_normalized, pressure_diff: vectors containing pressure data

% START = (M_S*60+S_S)+film_start;
% END = (M_E*60+S_E)+film_start;
% S_SAMP = START*25;
% E_SAMP = END*25;
% figure; hold on;
% 
% ZERO = zeros(1,END-START+1); %Zero vector
% TIME = (1:(E_SAMP-S_SAMP+1))/25'; 
% 
% a_mag = sqrt(dynamics.a_x.^2+dynamics.a_y.^2+dynamics.a_z.^2)-1;
% 
% STATIC_SHORT = dynamics.s_x(S_SAMP:E_SAMP);
% pressure_short = pressure_normalized(S_SAMP:E_SAMP)-0.5;
% OBDA_short = dynamics.obda(S_SAMP:E_SAMP);
% pitch_short = pitch(S_SAMP:E_SAMP);
% a_z_short = a_mag(S_SAMP:E_SAMP);
% 
% up_idx = idx.up(S_SAMP:E_SAMP);
% down_idx = idx.down(S_SAMP:E_SAMP);
% h_idx = idx.head(S_SAMP:E_SAMP);
% t_idx = idx.tail(S_SAMP:E_SAMP);
% c_idx = idx.change(S_SAMP:E_SAMP);
% fin_idx = idx.fin(S_SAMP:E_SAMP);
% flap_idx = idx.flap(S_SAMP:E_SAMP);
% jet_idx = idx.jet(S_SAMP:E_SAMP);
% 
% 
% xlabel('Time(s)');
% subplot(511) 
% hold on;
% pitch_short_rad = (pitch_short./(2*pi))*360;
% plot(TIME, pitch_short_rad, 'black', 'linewidth',2); axis([0 2 -2.5 2.5])
% plot(0:END-START,ZERO, 'black'); ylabel('Angle(degrees)');
% subplot(514)
% hold on;
% plot(0:END-START,ZERO, 'black'); ylabel('Velocity');
% subplot(515)
% hold on;
% plot(TIME, a_z_short, 'black', 'linewidth',2);
% plot(0:END-START,ZERO, 'black'); axis([0 2 -0.025 0.025]); ylabel('Acceleration(cm/s^2)');



