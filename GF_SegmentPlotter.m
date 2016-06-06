function [  ] = GF_SegmentPlotter( M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda, P_dec, t, v)
%Segement plotter plots a stylized segment of a graph, from a specified start to end
%time, within the data set.
%M_S: minute in video to start, S_S: second in video to start
%M_E: minute in video to end, S_E: second in video to end
%film_start: time (in samples) of the start of the video within the data set.
%dynamics: structure containing the acceleration dynamics data 
%pressure_normalized, pressure_diff: vectors containing pressure data

START = (M_S*60+S_S)+film_start;
END = (M_E*60+S_E)+film_start;
S_SAMP = START*25;
E_SAMP = END*25;

ZERO = zeros(1,END-START+1); %Zero vector
TIME = (1:(E_SAMP-S_SAMP+1))/25'; 

a_mag = sqrt(dynamics.a_x.^2+dynamics.a_y.^2+dynamics.a_z.^2)-1;

STATIC_SHORT = dynamics.s_x(S_SAMP:E_SAMP);
pressure_short = pressure_normalized(S_SAMP:E_SAMP)-0.5;
OBDA_short = dynamics.obda(S_SAMP:E_SAMP);
x_short = dynamics.a_x(S_SAMP:E_SAMP);
y_short = dynamics.a_y(S_SAMP:E_SAMP);
z_short = dynamics.a_z(S_SAMP:E_SAMP);
p_short = P_dec(S_SAMP:E_SAMP);


up_idx = idx.up(S_SAMP:E_SAMP);
down_idx = idx.down(S_SAMP:E_SAMP);
h_idx = idx.head(S_SAMP:E_SAMP);
t_idx = idx.tail(S_SAMP:E_SAMP);
c_idx = idx.change(S_SAMP:E_SAMP);
fin_idx = idx.fin(S_SAMP:E_SAMP);
flap_idx = idx.flap(S_SAMP:E_SAMP);
jet_idx = idx.jet(S_SAMP:E_SAMP);


low_limiter = mean_obda*2.5;
med_limiter = mean_obda*3.5;
high_limiter = mean_obda*5;
window_size = 25;
obda_process = zeros(1, length(OBDA_short))';
fqs = zeros(1, length(OBDA_short));

smooth_obda = moving(OBDA_short,5);

step = 10;
for i=1:window_size:(length(obda_process)-window_size)
%     [x, y] = FFT_data_set(OBDA_short(i:i+window_size), 25);
%     diff_y = diff(y);
%     y_pos = diff_y > 0; y_pos = [0 y_pos]; y_pos = (y_pos == 1);
%     x_sub_pos = x(y_pos);
%     if(mod(i, 10) == 0)
%         figure('units','normalized','outerposition',[0 0 1 1]);
%         subplot(411); hold on;
%         plot(TIME, smooth_obda, 'b')
%         plot(TIME((i:i+window_size)), smooth_obda((i:i+window_size)), '.g')
%        
%         subplot(412); hold on;
%         plot(TIME, OBDA_short, 'b')
%         plot(TIME((i:i+window_size)), OBDA_short((i:i+window_size)), '.g')
% 
%         subplot(413); hold on;
%         diff_y = diff(y);
%         plot(x(2:length(diff_y)+1), diff_y, '-*red')
%         plot(x(1:length(y_pos)), y_pos/15, 'magenta')
% 
%         plot(x, y,'-*k','linewidth',2) 
%         title('Single-Sided Amplitude Spectrum of X')
%         xlabel('Frequency (Hz)')
%         ylabel('|Accelleration|')
%         axis([0 5 0 0.2])
% 
%         subplot(414); hold on;
%         plot(TIME((i:i+window_size)), smooth_obda((i:i+window_size)), 'b')
%     end
% 
%     if(~isempty(x_sub_pos))
%         dom_fq = x_sub_pos(1);
%     else
%         dom_fq = 0;
%     end
%     [pks, locs] = max(y(2:length(y)));
%     pks
%     locs
%     dom_fq = x(locs+1)
%     fqs(i:i+step)=dom_fq;
    %fqs(i)=dom_fq;
    peak = max(findpeaks(OBDA_short(i:i+window_size)));
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
    
     if(mod(i, 10) == 0)
        [x,y] = ginput(1); %keep to deley the start of animation until click
        close all;
     end
end

fin_idx = (obda_process == 0);
slow_flap_idx = (obda_process == 1);
fast_flap_idx = (obda_process == 2);
jet_idx = (obda_process == 3);

% for i=1:length(fqs_short)
%     i
%    text(i, abs(sin(i))*0.2, num2str(fqs_short(i)))  
% end
figure; 
subplot(411); 
hold on;
plot(TIME(fin_idx), OBDA_short(fin_idx), '.',  'Color',[   0.4000    0.4000    1.0000], 'MarkerFaceColor', ...
    [   0.4000    0.4000    1.0000], 'Markersize', 15);
plot(TIME(slow_flap_idx), OBDA_short(slow_flap_idx), '.', 'Color',[ 0.4000    0.4000    1.0000], 'MarkerFaceColor', ...
    [ 0.4000    0.4000    1.0000], 'Markersize', 15);
plot(TIME(fast_flap_idx), OBDA_short(fast_flap_idx), '.', 'Color',[  0.9255    0.9255    0.0745], 'MarkerFaceColor',...
    [  0.9255    0.9255    0.0745], 'Markersize', 15);
plot(TIME(jet_idx), OBDA_short(jet_idx), '.', 'Color',[ 0.9412         0         0],  'MarkerFaceColor',[ 0.9412         0         0],...
    'Markersize', 15);
plot(TIME, OBDA_short, 'black');
xlabel('Time (s)');
ylabel('ODBA (g)');
for i=0:5:(END-START)
   hx = graph2d.constantline(i, 'Color','black');
   changedependvar(hx,'x'); 
end
axis([0 END-START 0 max(OBDA_short)+0.1]);

%pressure_short_lowered = pressure_short-max(pressure_short);
pressure_diff_norm = diff(pressure_short)*10;

subplot(412); hold on; 
plot(TIME(up_idx), pressure_short(up_idx), '.m', TIME(down_idx), pressure_short(down_idx), '.c');
xlabel('Time (s)');
ylabel('Depth (m)');
for i=0:5:(END-START)
   hx = graph2d.constantline(i, 'Color','black');
   changedependvar(hx,'x'); 
end
axis([0 END-START min(pressure_short)-0.1 max(pressure_short)+0.1]);

subplot(413);
hold on;
p_short = p_short-mean(p_short(c_idx));
plot(TIME(t_idx),p_short(t_idx), '.b', TIME(h_idx), p_short(h_idx), '.r', TIME(c_idx), p_short(c_idx), '.g') ;
%plot(TIME, p_short, 'k');
xlabel('Time (s)');
ylabel('Pitch (degrees)');

for i=0:5:(END-START)
   hx = graph2d.constantline(i, 'Color','black');
   changedependvar(hx,'x'); 
end
axis([0 END-START min(p_short)-2 max(p_short)+2]);

subplot(414)
for i=0:5:(END-START)
   hx = graph2d.constantline(i, 'Color','black');
   changedependvar(hx,'x'); 
end

if exist('v', 'var')
    hold on;
    
%     plot(t(fin_idx), v(fin_idx), '.',  'Color',[   0.4000    0.4000    1.0000], 'MarkerFaceColor', ...
%     [   0.4000    0.4000    1.0000], 'Markersize', 15);
% plot(t(slow_flap_idx), v(slow_flap_idx), '.', 'Color',[ 0.4000    0.4000    1.0000], 'MarkerFaceColor', ...
%     [ 0.4000    0.4000    1.0000], 'Markersize', 15);
% plot(t(fast_flap_idx), v(fast_flap_idx), '.', 'Color',[  0.9255    0.9255    0.0745], 'MarkerFaceColor',...
%     [  0.9255    0.9255    0.0745], 'Markersize', 15);
% plot(t(jet_idx), v(jet_idx), '.', 'Color',[ 0.9412         0         0],  'MarkerFaceColor',[ 0.9412         0         0],...
%     'Markersize', 15);
%   

    plot(t(t<=6), v(t<=6), '.', 'Color',[  0.9255    0.9255    0.0745], 'MarkerFaceColor', ...
     [  0.9255    0.9255    0.0745], 'Markersize', 15);
    plot(t(t>6 & t<=10), v(t>6 & t<=10), '.', 'Color',[   0.4000    0.4000    1.0000], 'MarkerFaceColor', ...
     [   0.4000    0.4000    1.0000], 'Markersize', 15);
    plot(t(t>10 & t<=11), v(t>10 & t<=11), '.', 'Color',[   0.9412         0         0], 'MarkerFaceColor', ...
     [  0.9412         0         0], 'Markersize', 15);
    plot(t(t>11 & t<=12), v(t>11 & t<=12), '.', 'Color',[  0.9255    0.9255    0.0745], 'MarkerFaceColor', ...
     [  0.9255    0.9255    0.0745], 'Markersize', 15);
    plot(t(t>12), v(t>12), '.', 'Color',[   0.4000    0.4000    1.0000], 'MarkerFaceColor', ...
     [   0.4000    0.4000    1.0000], 'Markersize', 15); 
 
    plot(t, v, 'k');
    
    plot(t(1), v(1), '*k', 'Markersize', 18);
    plot(t(length(t)), v(length(v)), '*k', 'Markersize', 18);
end
axis([0 END-START min(v)-0.1 max(v)+0.1]);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
% for i=0:1:(END-START)
%    hx = graph2d.constantline(i, 'Color',[  0.5020    0.5020    0.5020], 'LineStyle', ':');
%    changedependvar(hx,'x'); 
% end
% for i=0:1:(END-START)
%    hx = graph2d.constantline(i, 'Linestyle', ':', 'Color','black');
%    changedependvar(hx,'x'); 
% end
% Uncomment to simulate data point moving through graph
% 
% [x,y] = ginput(1); %keep to deley the start of animation until click
% 
% control animation speed
% DELAY = 0.058; %wihtout lines
% DELAY = 0.0315; %good for trial
% numPoints = length(STATIC_SHORT);
% 
% Setup animation for STATIC ACCELERATION
% hLine = line('XData',TIME(1), 'YData',STATIC_SHORT(1), 'Color','r', ...
%     'Marker','o', 'MarkerSize',6, 'LineWidth',2);
% hTxt = text(TIME(1),STATIC_SHORT(1), sprintf('(%.3f,%.3f)',TIME(1),STATIC_SHORT(1)), ...
%     'Color',[0.2 0.2 0.2], 'FontSize',8, ...
%     'HorizontalAlignment','left', 'VerticalAlignment','top');
% 
% Setup animation for pressure
% jLine = line('XData',TIME(1), 'YData',pressure_short(1), 'Color','r', ...
%     'Marker','o', 'MarkerSize',6, 'LineWidth',2);
% jTxt = text(TIME(1),pressure_short(1), sprintf('(%.3f,%.3f)',TIME(1),pressure_short(1)), ...
%     'Color',[0.2 0.2 0.2], 'FontSize',8, ...
%     'HorizontalAlignment','left', 'VerticalAlignment','top');
% 
% infinite loop
% i = 1;                                       %# index
% while true      
%     update point & text
%     set(hLine, 'XData',TIME(i), 'YData',STATIC_SHORT(i))   
%     set(hTxt, 'Position',[TIME(i) STATIC_SHORT(i)], ...
%         'String',sprintf('(%.3f,%.3f)',[TIME(i) STATIC_SHORT(i)]))        
%    
%     set(jLine, 'XData',TIME(i), 'YData',pressure_short(i))   
%     set(jTxt, 'Position',[TIME(i) pressure_short(i)], ...
%         'String',sprintf('(%.3f,%.3f)',[TIME(i) pressure_short(i)]))        
%     
%     drawnow                                  %# force refresh
%     pause(DELAY)                           %# slow down animation
% 
%     i = rem(i+1,numPoints)+1;                %# circular increment
%     if ~ishandle(hLine), break; end          %# in case you close the figure
% end


end

