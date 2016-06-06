%% Save OBDA and IDX individaully so that all 6 specimen can be worked on similtaneously 
% Only needs to be run once; very time consuming

% load lf077a_complete
% obda_lf077a = OBDA;
% idx_lf077a = idx;
% summary_stats_lf077a = summary_stats;
% save lf077a_complete
% 
% load lf078a_complete
% obda_lf078a = OBDA;
% idx_lf078a = idx;
% summary_stats_lf078a = summary_stats;
% save lf078a_complete
% 
% load lf079a_complete
% obda_lf079a = OBDA;
% idx_lf079a = idx;
% summary_stats_lf079a = summary_stats;
% save lf079a_complete
% 
% load lf080a_complete
% obda_lf080a = OBDA;
% idx_lf080a = idx;
% summary_stats_lf080a = summary_stats;
% save lf080a_complete
% 
% load lf082a_complete
% obda_lf082a = OBDA;
% idx_lf082a = idx;
% summary_stats_lf082a = summary_stats;
% save lf082a_complete
% 
% load lf086a_complete
% obda_lf086a = OBDA;
% idx_lf086a = idx;
% summary_stats_lf086a = summary_stats;
% save lf086a_complete
% 
% load lf089a_complete
% obda_lf089a = OBDA;
% idx_lf089a = idx;
% summary_stats_lf089a = summary_stats;
% save lf089a_complete
% 
% load lf090a_complete
% obda_lf090a = OBDA;
% idx_lf090a = idx;
% save lf090a_complete

%% Generate Summary Statistics
numSamples = 7;
load lf077a_obda summary_stats_lf077a
load lf078a_obda summary_stats_lf078a
load lf079a_obda summary_stats_lf079a
load lf080a_obda summary_stats_lf080a
load lf082a_obda summary_stats_lf082a
load lf086a_obda summary_stats_lf086a
load lf089a_obda summary_stats_lf089a
load lf077a_complete obda_lf077a idx_lf077a
load lf078a_complete obda_lf078a idx_lf078a
load lf079a_complete obda_lf079a idx_lf079a
load lf080a_complete obda_lf080a idx_lf080a
load lf082a_complete obda_lf082a idx_lf082a
load lf086a_complete obda_lf086a idx_lf086a
load lf089a_complete obda_lf089a idx_lf089a

%% Build Comprehensive Stat Structs
falling_to_rising = [summary_stats_lf077a.ratio_falling_rising summary_stats_lf078a.ratio_falling_rising summary_stats_lf079a.ratio_falling_rising summary_stats_lf080a.ratio_falling_rising ...
    summary_stats_lf082a.ratio_falling_rising summary_stats_lf086a.ratio_falling_rising summary_stats_lf089a.ratio_falling_rising];
rising_obda = [summary_stats_lf077a.rising_obda summary_stats_lf078a.rising_obda summary_stats_lf079a.rising_obda summary_stats_lf080a.rising_obda ...
    summary_stats_lf082a.rising_obda summary_stats_lf086a.rising_obda summary_stats_lf089a.rising_obda];
falling_obda = [summary_stats_lf077a.falling_obda summary_stats_lf078a.falling_obda summary_stats_lf079a.falling_obda summary_stats_lf080a.falling_obda ...
    summary_stats_lf082a.falling_obda summary_stats_lf086a.falling_obda summary_stats_lf089a.falling_obda];
rising_pitch = [summary_stats_lf077a.rising_pitch summary_stats_lf078a.rising_pitch summary_stats_lf079a.rising_pitch summary_stats_lf080a.rising_pitch ...
    summary_stats_lf082a.rising_pitch summary_stats_lf086a.rising_pitch summary_stats_lf089a.rising_pitch];
falling_pitch = [summary_stats_lf077a.falling_pitch summary_stats_lf078a.falling_pitch summary_stats_lf079a.falling_pitch summary_stats_lf080a.falling_pitch ...
    summary_stats_lf082a.falling_pitch summary_stats_lf086a.falling_pitch summary_stats_lf089a.falling_pitch];
heads_up_ratio = [summary_stats_lf077a.heads_up_ratio summary_stats_lf078a.heads_up_ratio summary_stats_lf079a.heads_up_ratio summary_stats_lf080a.heads_up_ratio ...
    summary_stats_lf082a.heads_up_ratio summary_stats_lf086a.heads_up_ratio summary_stats_lf089a.heads_up_ratio];
heads_down_ratio = [summary_stats_lf077a.heads_down_ratio summary_stats_lf078a.heads_down_ratio summary_stats_lf079a.heads_down_ratio summary_stats_lf080a.heads_down_ratio ...
    summary_stats_lf082a.heads_down_ratio summary_stats_lf086a.heads_down_ratio summary_stats_lf089a.heads_down_ratio];
tails_up_ratio = [summary_stats_lf077a.tails_up_ratio summary_stats_lf078a.tails_up_ratio summary_stats_lf079a.tails_up_ratio summary_stats_lf080a.tails_up_ratio ...
    summary_stats_lf082a.tails_up_ratio summary_stats_lf086a.tails_up_ratio summary_stats_lf089a.tails_up_ratio];
tails_down_ratio = [summary_stats_lf077a.tails_down_ratio summary_stats_lf078a.tails_down_ratio summary_stats_lf079a.tails_down_ratio summary_stats_lf080a.tails_down_ratio ...
    summary_stats_lf082a.tails_down_ratio summary_stats_lf086a.tails_down_ratio summary_stats_lf089a.tails_down_ratio];
heads_to_tails = [summary_stats_lf077a.heads_to_tails summary_stats_lf078a.heads_to_tails summary_stats_lf079a.heads_to_tails summary_stats_lf080a.heads_to_tails ...
    summary_stats_lf082a.heads_to_tails summary_stats_lf086a.heads_to_tails summary_stats_lf089a.heads_to_tails];
OBDA_heads_up = [summary_stats_lf077a.OBDA_heads_up summary_stats_lf078a.OBDA_heads_up summary_stats_lf079a.OBDA_heads_up summary_stats_lf080a.OBDA_heads_up ...
    summary_stats_lf082a.OBDA_heads_up summary_stats_lf086a.OBDA_heads_up summary_stats_lf089a.OBDA_heads_up];
OBDA_heads_down = [summary_stats_lf077a.OBDA_heads_down summary_stats_lf078a.OBDA_heads_down summary_stats_lf079a.OBDA_heads_down summary_stats_lf080a.OBDA_heads_down ...
    summary_stats_lf082a.OBDA_heads_down summary_stats_lf086a.OBDA_heads_down summary_stats_lf089a.OBDA_heads_down];
OBDA_tails_up = [summary_stats_lf077a.OBDA_tails_up summary_stats_lf078a.OBDA_tails_up summary_stats_lf079a.OBDA_tails_up summary_stats_lf080a.OBDA_tails_up ...
    summary_stats_lf082a.OBDA_tails_up summary_stats_lf086a.OBDA_tails_up summary_stats_lf089a.OBDA_tails_up];
OBDA_tails_down = [summary_stats_lf077a.OBDA_tails_down summary_stats_lf078a.OBDA_tails_down summary_stats_lf079a.OBDA_tails_down summary_stats_lf080a.OBDA_tails_down ...
    summary_stats_lf082a.OBDA_tails_down summary_stats_lf086a.OBDA_tails_down summary_stats_lf089a.OBDA_tails_down];
OBDA_heads = [summary_stats_lf077a.OBDA_heads summary_stats_lf078a.OBDA_heads summary_stats_lf079a.OBDA_heads summary_stats_lf080a.OBDA_heads ...
    summary_stats_lf082a.OBDA_heads summary_stats_lf086a.OBDA_heads summary_stats_lf089a.OBDA_heads];
OBDA_tails = [summary_stats_lf077a.OBDA_tails summary_stats_lf078a.OBDA_tails summary_stats_lf079a.OBDA_tails summary_stats_lf080a.OBDA_tails ...
    summary_stats_lf082a.OBDA_tails summary_stats_lf086a.OBDA_tails summary_stats_lf089a.OBDA_tails];

figure; 
subplot(211);hold on; title('Ratio of Falling to Fising')
plot([0 1 2 3 4 5 6 7 8], [1 1 1 1 1 1 1 1 1],'k', 'Markersize', 20); axis([0 8 0 1.5])
plot(falling_to_rising,'.b', 'Markersize', 30); axis([0 8 0 1.5])
ax = gca;
set(ax,  'XTickLabel', {'lf077a', 'lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});
set(ax,'XTick',[1 2 3 4 5 6 7]);
avg_string = strcat('  Avg: ', num2str(mean(falling_to_rising), 4));
std_string = strcat('  Std: ', num2str(std(falling_to_rising), 3));
text(0, 0.3, avg_string)
text(0, 0.5, std_string)

subplot(212);hold on; title('Ratio of Heads to Tails')
plot([0 1 2 3 4 5 6 7 8], [1 1 1 1 1 1 1 1 1],'k', 'Markersize', 20); axis([0 8 0 1.5])
plot(heads_to_tails,'.b', 'Markersize', 30); axis([0 8 0 1.5])
ax = gca;
set(ax,  'XTickLabel', {'lf077a', 'lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});
set(ax,'XTick',[1 2 3 4 5 6 7]);
avg_string = strcat('  Avg: ', num2str(mean(heads_to_tails), 4));
std_string = strcat('  Std: ', num2str(std(heads_to_tails), 3));
text(0, 0.3, avg_string)
text(0, 0.5, std_string)


figure;
subplot(221); hold on;title('Rising OBDA vs Falling OBDA')
% plot(rising_obda,'.m', 'Markersize', 30); axis([0 7 0 0.1])
% plot(falling_obda, '.c', 'Markersize', 30); axis([0 7 0 0.1])


C = [obda_lf077a(idx_lf077a.up); obda_lf077a(idx_lf077a.down);...
    obda_lf078a(idx_lf078a.up); obda_lf078a(idx_lf078a.down);...
    obda_lf079a(idx_lf079a.up); obda_lf079a(idx_lf079a.down);...
    obda_lf080a(idx_lf080a.up); obda_lf080a(idx_lf080a.down);...
    obda_lf082a(idx_lf082a.up); obda_lf082a(idx_lf082a.down);...
    obda_lf086a(idx_lf086a.up); obda_lf086a(idx_lf086a.down);...
    obda_lf089a(idx_lf089a.up); obda_lf089a(idx_lf089a.down)];

obda_up_all = [obda_lf077a(idx_lf077a.up);...
    obda_lf078a(idx_lf078a.up);...
    obda_lf079a(idx_lf079a.up);...
    obda_lf080a(idx_lf080a.up);...
    obda_lf082a(idx_lf082a.up);...
    obda_lf086a(idx_lf086a.up);...
    obda_lf089a(idx_lf089a.up)];

obda_down_all = [obda_lf077a(idx_lf077a.down);...
    obda_lf078a(idx_lf078a.down);...
    obda_lf079a(idx_lf079a.down);...
    obda_lf080a(idx_lf080a.down);...
    obda_lf082a(idx_lf082a.down);...
    obda_lf086a(idx_lf086a.down);...
    obda_lf089a(idx_lf089a.down)];

up_78 = zeros(1, sum(idx_lf078a.up)); up_78(1, :) = 0;
down_78 = zeros(1, sum(idx_lf078a.down)); down_78(1, :) = 1;

up_79 = zeros(1, sum(idx_lf079a.up)); up_79(1, :) = 2;
down_79 = zeros(1, sum(idx_lf079a.down)); down_79(1, :) = 3;

up_80 = zeros(1, sum(idx_lf080a.up)); up_80(1, :) = 4;
down_80 = zeros(1, sum(idx_lf080a.down)); down_80(1, :) = 5;

up_82 = zeros(1, sum(idx_lf082a.up)); up_82(1, :) = 6;
down_82 = zeros(1, sum(idx_lf082a.down)); down_82(1, :) = 7;

up_86 = zeros(1, sum(idx_lf086a.up)); up_86(1, :) = 8;
down_86 = zeros(1, sum(idx_lf086a.down)); down_86(1, :) = 9;

up_89 = zeros(1, sum(idx_lf089a.up)); up_89(1, :) = 10;
down_89 = zeros(1, sum(idx_lf089a.down)); down_89(1, :) = 11;

up_77 = zeros(1, sum(idx_lf077a.up)); up_77(1, :) = 12;
down_77 = zeros(1, sum(idx_lf077a.down)); down_77(1, :) = 13;

up_all = [up_77 up_78 up_79 up_80 up_82 up_86 up_89];
down_all = [down_77 down_78 down_79 down_80 down_82 down_86 down_89];

grp = [up_77 down_77 up_78 down_78 up_79 down_79 up_80 down_80 up_82 down_82 ...
    up_86 down_86 up_89 down_89]';

positions = [0.90 1.10 1.90 2.10 2.90 3.10 3.90 4.10 4.90 5.10 5.90 6.10 6.90 7.10];
color = ['g', 'y','g', 'y','g', 'y','g', 'y','g', 'y','g', 'y', 'g', 'y'];

boxplot(C, grp, 'positions', positions, 'symbol', 'r'); axis([0 8 0 0.25])
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
 end

avg_string_up = strcat('Up Avg: ', num2str(mean(obda_up_all), 3));
std_string_up = strcat('Up Std: ', num2str(std(obda_up_all), 3));

avg_string_down = strcat('Down Avg: ', num2str(mean(obda_down_all), 3));
std_string_down = strcat('Down Std: ', num2str(std(obda_down_all), 3));

text(2.5, 0.22, avg_string_up, 'EdgeColor', 'y')
text(2.5, 0.20, std_string_up, 'EdgeColor', 'y') 

text(3.5, 0.22, avg_string_down, 'EdgeColor', 'g')
text(3.5, 0.20, std_string_down, 'EdgeColor', 'g') 
 
 
%legend('Rising', 'Falling', 'Location','south', 'Orientation', 'horizontal'); 
ax = gca;
set(ax,  'XTickLabel', {'lf077a', 'lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});
set(ax,'XTick',[1 2 3 4 5 6 7]);


subplot(222); hold on;title('Heads OBDA vs Tails OBDA')

obda_head_all = [obda_lf077a(idx_lf077a.head);...
    obda_lf078a(idx_lf078a.head);...
    obda_lf079a(idx_lf079a.head);...
    obda_lf080a(idx_lf080a.head);...
    obda_lf082a(idx_lf082a.head);...
    obda_lf086a(idx_lf086a.head);...
    obda_lf089a(idx_lf089a.head)];

obda_tail_all = [obda_lf077a(idx_lf077a.tail);...
    obda_lf078a(idx_lf078a.tail);...
    obda_lf079a(idx_lf079a.tail);...
    obda_lf080a(idx_lf080a.tail);...
    obda_lf082a(idx_lf082a.tail);...
    obda_lf086a(idx_lf086a.tail);...
    obda_lf089a(idx_lf089a.tail)];
 
C = [obda_lf077a(idx_lf077a.head); obda_lf077a(idx_lf077a.tail);...
    obda_lf078a(idx_lf078a.head); obda_lf078a(idx_lf078a.tail);...
    obda_lf079a(idx_lf079a.head); obda_lf079a(idx_lf079a.tail);...
    obda_lf080a(idx_lf080a.head); obda_lf080a(idx_lf080a.tail);...
    obda_lf082a(idx_lf082a.head); obda_lf082a(idx_lf082a.tail);...
    obda_lf086a(idx_lf086a.head); obda_lf086a(idx_lf086a.tail);...
    obda_lf089a(idx_lf089a.head); obda_lf089a(idx_lf089a.tail)];

head_78 = zeros(1, sum(idx_lf078a.head)); head_78(1, :) = 0;
tail_78 = zeros(1, sum(idx_lf078a.tail)); tail_78(1, :) = 1;

head_79 = zeros(1, sum(idx_lf079a.head)); head_79(1, :) = 2;
tail_79 = zeros(1, sum(idx_lf079a.tail)); tail_79(1, :) = 3;

head_80 = zeros(1, sum(idx_lf080a.head)); head_80(1, :) = 4;
tail_80 = zeros(1, sum(idx_lf080a.tail)); tail_80(1, :) = 5;

head_82 = zeros(1, sum(idx_lf082a.head)); head_82(1, :) = 6;
tail_82 = zeros(1, sum(idx_lf082a.tail)); tail_82(1, :) = 7;

head_86 = zeros(1, sum(idx_lf086a.head)); head_86(1, :) = 8;
tail_86 = zeros(1, sum(idx_lf086a.tail)); tail_86(1, :) = 9;

head_89 = zeros(1, sum(idx_lf089a.head)); head_89(1, :) = 10;
tail_89 = zeros(1, sum(idx_lf089a.tail)); tail_89(1, :) = 11;

head_77 = zeros(1, sum(idx_lf077a.head)); head_77(1, :) = 12;
tail_77 = zeros(1, sum(idx_lf077a.tail)); tail_77(1, :) = 13;

head_all = [head_77 head_78 head_79 head_80 head_82 head_86 head_89];
tail_all = [tail_77 tail_78 tail_79 tail_80 tail_82 tail_86 tail_89];

grp = [head_77 tail_77 head_78 tail_78 head_79 tail_79 head_80 tail_80 head_82 tail_82 ...
    head_86 tail_86 head_89 tail_89]';
positions = [0.90 1.10 1.90 2.10 2.90 3.10 3.90 4.10 4.90 5.10 5.90 6.10 6.90 7.10];
color = ['b', 'r','b', 'r','b', 'r','b', 'r','b', 'r','b', 'r' , 'b', 'r' ];
h = boxplot(C, grp, 'positions', positions, 'symbol', 'r'); axis([0 8 0 0.25])
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
 end
 
avg_string_head = strcat('Head Avg: ', num2str(mean(obda_head_all), 3));
std_string_head = strcat('Head Std: ', num2str(std(obda_head_all), 3));

avg_string_tail = strcat('Tail Avg: ', num2str(mean(obda_tail_all), 3));
std_string_tail = strcat('Tail Std: ', num2str(std(obda_tail_all), 3));

text(2.5, 0.22, avg_string_head, 'EdgeColor', 'r')
text(2.5, 0.20, std_string_head, 'EdgeColor', 'r') 

text(3.8, 0.22, avg_string_tail, 'EdgeColor', 'b')
text(3.8, 0.20, std_string_tail, 'EdgeColor', 'b') 
 
%legend('Heads', 'Tails', 'Location','south', 'Orientation', 'horizontal'); 
ax = gca;
set(ax,  'XTickLabel', {'lf077a', 'lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});
set(ax,'XTick',[1 2 3 4 5 6 7]);

subplot(223); hold on;title('Time Spent in Orientations')
time_orientations = [heads_up_ratio; heads_down_ratio; tails_up_ratio; tails_down_ratio]';
bar(time_orientations, 'stacked')
legend('Heads Upwards', 'Heads Downwards', 'Tail Upwards', 'Tails Downwards'); legend('Location', 'south', 'Orientation', 'horizontal');
ax = gca;
set(ax,  'XTickLabel', {'lf077a','lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});
set(ax,'XTick',[1 2 3 4 5 6 7]);

subplot(224); hold on;title('OBDA in Orientations')
orientations = [OBDA_heads_up; OBDA_heads_down; OBDA_tails_up; OBDA_tails_down]';

bar(orientations)
legend_handle = legend('Heads Upwards', 'Heads Downwards', 'Tail Upwards', 'Tails Downwards'); legend('Location', 'south', 'Orientation', 'horizontal');
legend_markers = findobj();
ax = gca;
set(ax,  'XTickLabel', {'lf077a', 'lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});
set(ax,'XTick',[1 2 3 4 5 6 7]);

%% Plot grouped comparison between heads/tails/up/down orientations time and obda
% Hatched graphs are OBDA, solid is ratio of time spent

figure;
subplot(121); hold on;title('Time Spent in Orientations')
time_orientations = [heads_up_ratio; heads_down_ratio; tails_up_ratio; tails_down_ratio]';
bar(time_orientations, 'stacked')
legend('Heads Upwards', 'Heads Downwards', 'Tail Upwards', 'Tails Downwards'); legend('Location', 'south', 'Orientation', 'horizontal');
ax = gca;
set(ax,  'XTickLabel', {'lf077a','lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});
set(ax,'XTick',[1 2 3 4 5 6 7]);


subplot(122); hold on;title('OBDA Ratio in Orientations')
orientations_overall = sum(orientations,2);
orientations_ratio = orientations;
for i=1:numSamples
    orientations_ratio(i, :) = orientations(i, :)./orientations_overall(i);
end

figure;
NumStacksPerGroup = 2;
NumGroupsPerAxis = 7;
NumStackElements = 4;

%Lables to use on tick marks
groupLabels = {'lf077a', 'lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'};
stackData = zeros(NumGroupsPerAxis, NumStacksPerGroup, NumStackElements);
stackData(:, 1, :) = orientations_ratio;
stackData(:, 2, :) = time_orientations;
plotBarStackGroups(stackData, groupLabels);
axis([0 8 0 1]);

%% Plot rising vs Falling Pitch
figure; hold on; title('Rising Pitch vs Falling Pitch')
plot(rising_pitch,'.m', 'Markersize', 30); axis([0 8 0 0.25])
plot(falling_pitch, '.c', 'Markersize', 30); axis([0 8 0 0.25])

legend('Rising', 'Falling', 'Location','south', 'Orientation', 'horizontal');
ax = gca;
set(ax,'XTick',[1 2 3 4 5 6 7]);
set(ax,  'XTickLabel', {'lf077a', 'lf078a', 'lf079a', 'lf080a', 'lf082a', 'lf086a', 'lf089a'});

%% Bout Length Detection and Graph Bout Lengths Histogram
head_lf077a = +idx_lf077a.head;
head_lf078a = +idx_lf078a.head;
head_lf079a = +idx_lf079a.head;
head_lf080a = +idx_lf080a.head;
head_lf082a = +idx_lf082a.head;
head_lf086a = +idx_lf086a.head;
head_lf089a = +idx_lf089a.head;

out_lf077a = findseq(head_lf077a);
out_lf078a = findseq(head_lf078a);
out_lf079a = findseq(head_lf079a);
out_lf080a = findseq(head_lf080a);
out_lf082a = findseq(head_lf082a);
out_lf086a = findseq(head_lf086a);
out_lf089a = findseq(head_lf089a);

out_tails_lf077a = (out_lf077a(:, 1) == 0);
out_tails_lf078a = (out_lf078a(:, 1) == 0);
out_tails_lf079a = (out_lf079a(:, 1) == 0);
out_tails_lf080a = (out_lf080a(:, 1) == 0);
out_tails_lf082a = (out_lf082a(:, 1) == 0);
out_tails_lf086a = (out_lf086a(:, 1) == 0);
out_tails_lf089a = (out_lf089a(:, 1) == 0);

bouts_lf077a = out_lf077a(:, 4);
bouts_lf078a = out_lf078a(:, 4);
bouts_lf079a = out_lf079a(:, 4);
bouts_lf080a = out_lf089a(:, 4);
bouts_lf082a = out_lf082a(:, 4);
bouts_lf086a = out_lf086a(:, 4);
bouts_lf089a = out_lf089a(:, 4);

tail_bouts_lf077a = bouts_lf077a(out_tails_lf077a); tail_bouts_lf077a = tail_bouts_lf077a./25;
head_bouts_lf077a = bouts_lf077a(~out_tails_lf077a); head_bouts_lf077a = head_bouts_lf077a./25;
tail_bouts_lf078a = bouts_lf078a(out_tails_lf078a); tail_bouts_lf078a = tail_bouts_lf078a./25;
head_bouts_lf078a = bouts_lf078a(~out_tails_lf078a); head_bouts_lf078a = head_bouts_lf078a./25;
tail_bouts_lf079a = bouts_lf079a(out_tails_lf079a); tail_bouts_lf079a = tail_bouts_lf079a./25;
head_bouts_lf079a = bouts_lf079a(~out_tails_lf079a); head_bouts_lf079a = head_bouts_lf079a./25;
tail_bouts_lf080a = bouts_lf080a(out_tails_lf080a); tail_bouts_lf080a = tail_bouts_lf080a./25;
head_bouts_lf080a = bouts_lf080a(~out_tails_lf080a); head_bouts_lf080a = head_bouts_lf080a./25;
tail_bouts_lf082a = bouts_lf082a(out_tails_lf082a); tail_bouts_lf082a = tail_bouts_lf082a./25;
head_bouts_lf082a = bouts_lf082a(~out_tails_lf082a); head_bouts_lf082a = head_bouts_lf082a./25;
tail_bouts_lf086a = bouts_lf086a(out_tails_lf086a); tail_bouts_lf086a = tail_bouts_lf086a./25;
head_bouts_lf086a = bouts_lf086a(~out_tails_lf086a); head_bouts_lf086a = head_bouts_lf086a./25;
tail_bouts_lf089a = bouts_lf089a(out_tails_lf089a); tail_bouts_lf089a = tail_bouts_lf089a./25;
head_bouts_lf089a = bouts_lf089a(~out_tails_lf089a); head_bouts_lf089a = head_bouts_lf089a./25;
%% Plot Histograms of Bout Length
figure;
subplot(421)
hist(tail_bouts_lf077a, 50), title('lf077a Tail Bout Length')
avg_string_tail = strcat('Avg: ', num2str(mean(tail_bouts_lf077a), 3));
std_string_tail = strcat('Std: ', num2str(std(tail_bouts_lf077a), 3));
text(10, 200, avg_string_tail, 'EdgeColor', 'b')
text(10, 150, std_string_tail, 'EdgeColor', 'b') 

subplot(422)
hist(head_bouts_lf077a, 50), title('lf077a Head Bout Length')
avg_string_head = strcat('Avg: ', num2str(mean(head_bouts_lf077a), 3));
std_string_head = strcat('Std: ', num2str(std(head_bouts_lf077a), 3));
text(10, 200, avg_string_head, 'EdgeColor', 'r')
text(10, 150, std_string_head, 'EdgeColor', 'r')

subplot(423)
hist(tail_bouts_lf078a, 50), title('lf078a Tail Bout Length')
avg_string_tail = strcat('Avg: ', num2str(mean(tail_bouts_lf078a), 3));
std_string_tail = strcat('Std: ', num2str(std(tail_bouts_lf078a), 3));
text(10, 200, avg_string_tail, 'EdgeColor', 'b')
text(10, 150, std_string_tail, 'EdgeColor', 'b') 

subplot(424)
hist(head_bouts_lf078a, 50), title('lf078a Head Bout Length')
avg_string_head = strcat('Avg: ', num2str(mean(head_bouts_lf078a), 3));
std_string_head = strcat('Std: ', num2str(std(head_bouts_lf078a), 3));
text(10, 200, avg_string_head, 'EdgeColor', 'r')
text(10, 150, std_string_head, 'EdgeColor', 'r')

subplot(425)
hist(tail_bouts_lf079a, 50), title('lf079a Tail Bout Length')
avg_string_tail = strcat('Avg: ', num2str(mean(tail_bouts_lf079a), 3));
std_string_tail = strcat('Std: ', num2str(std(tail_bouts_lf079a), 3));
text(10, 200, avg_string_tail, 'EdgeColor', 'b')
text(10, 150, std_string_tail, 'EdgeColor', 'b') 

subplot(426)
hist(head_bouts_lf079a, 50), title('lf079a Head Bout Length')
avg_string_head = strcat('Avg: ', num2str(mean(head_bouts_lf079a), 3));
std_string_head = strcat('Std: ', num2str(std(head_bouts_lf079a), 3));
text(10, 200, avg_string_head, 'EdgeColor', 'r')
text(10, 150, std_string_head, 'EdgeColor', 'r')

subplot(427)
hist(tail_bouts_lf080a, 50), title('lf080a Tail Bout Length')
avg_string_tail = strcat('Avg: ', num2str(mean(tail_bouts_lf080a), 3));
std_string_tail = strcat('Std: ', num2str(std(tail_bouts_lf080a), 3));
text(10, 200, avg_string_tail, 'EdgeColor', 'b')
text(10, 150, std_string_tail, 'EdgeColor', 'b') 

subplot(428)
hist(head_bouts_lf080a, 50), title('lf080a Head Bout Length')
avg_string_head = strcat('Avg: ', num2str(mean(head_bouts_lf080a), 3));
std_string_head = strcat('Std: ', num2str(std(head_bouts_lf080a), 3));
text(10, 200, avg_string_head, 'EdgeColor', 'r')
text(10, 150, std_string_head, 'EdgeColor', 'r')

figure;
subplot(421)
hist(tail_bouts_lf082a, 50), title('lf082a Tail Bout Length')
avg_string_tail = strcat('Avg: ', num2str(mean(tail_bouts_lf082a), 3));
std_string_tail = strcat('Std: ', num2str(std(tail_bouts_lf082a), 3));
text(10, 200, avg_string_tail, 'EdgeColor', 'b')
text(10, 150, std_string_tail, 'EdgeColor', 'b') 

subplot(422)
hist(head_bouts_lf082a, 50), title('lf082a Head Bout Length')
avg_string_head = strcat('Avg: ', num2str(mean(head_bouts_lf082a), 3));
std_string_head = strcat('Std: ', num2str(std(head_bouts_lf082a), 3));
text(10, 200, avg_string_head, 'EdgeColor', 'r')
text(10, 150, std_string_head, 'EdgeColor', 'r')

subplot(423)
hist(tail_bouts_lf086a, 50), title('lf086a Tail Bout Length')
avg_string_tail = strcat('Avg: ', num2str(mean(tail_bouts_lf086a), 3));
std_string_tail = strcat('Std: ', num2str(std(tail_bouts_lf086a), 3));
text(10, 200, avg_string_tail, 'EdgeColor', 'b')
text(10, 150, std_string_tail, 'EdgeColor', 'b') 

subplot(424)
hist(head_bouts_lf086a, 50), title('lf086a Head Bout Length')
avg_string_head = strcat('Avg: ', num2str(mean(head_bouts_lf086a), 3));
std_string_head = strcat('Std: ', num2str(std(head_bouts_lf086a), 3));
text(10, 200, avg_string_head, 'EdgeColor', 'r')
text(10, 150, std_string_head, 'EdgeColor', 'r')

subplot(425)
hist(tail_bouts_lf089a, 50), title('lf089a Tail Bout Length')
avg_string_tail = strcat('Avg: ', num2str(mean(tail_bouts_lf089a), 3));
std_string_tail = strcat('Std: ', num2str(std(tail_bouts_lf089a), 3));
text(10, 200, avg_string_tail, 'EdgeColor', 'b')
text(10, 150, std_string_tail, 'EdgeColor', 'b') 

subplot(426)
hist(head_bouts_lf089a, 50), title('lf089a Head Bout Length')
avg_string_head = strcat('Avg: ', num2str(mean(head_bouts_lf089a), 3));
std_string_head = strcat('Std: ', num2str(std(head_bouts_lf089a), 3));
text(10, 200, avg_string_head, 'EdgeColor', 'r')
text(10, 150, std_string_head, 'EdgeColor', 'r')


%% Generage histogram of peak heights for day and night
  %%computationally intensive, so preform as little as possible
% load lf077a_obda
% [window_lf077a, binnumber_lf077a, loc_diff_lf077a] = Histogram(OBDA);
% save lf077a_hist
%
% load lf078a_obda
% [window_lf078a, binnumber_lf078a, loc_diff_lf078a] = Histogram(OBDA);
% save lf078a_hist
% 
% load lf079a_obda
% [window_lf079a, binnumber_lf079a, loc_diff_lf079a] = Histogram(OBDA);
% save lf079a_hist
% 
% load lf080a_obda
% [window_lf080a, binnumber_lf080a, loc_diff_lf080a] = Histogram(OBDA);
% save lf080a_hist
% 
% load lf082a_obda
% [window_lf082a, bi,  nnumber_lf082a, loc_diff_lf082a] = Histogram(OBDA);
% save lf082a_hist
% 
% load lf086a_obda
% [window_lf086a, binnumber_lf0786a, loc_diff_lf086a] = Histogram(OBDA);
% save lf086a_hist
% 
% load lf089a_obda
% [window_lf089a, binnumber_lf089a, loc_diff_lf089a] = Histogram(OBDA);
% save lf089a_hist
% 
% [window_fin, binnumber_fin, loc_diff_fin] = Histogram(OBDA(fin_idx));
% [window_fin_flap, binnumber_fin_flap, loc_diff_fin_flap] = Histogram(OBDA(slow_flap_idx));
% [window_flap, binnumber_flap, loc_diff_flap] = Histogram(OBDA(fast_flap_idx));
% [window_jet, binnumber_jet, loc_diff_jet] = Histogram(OBDA(jet_idx));
 
%% Plot the histogram
load lf077a_hist window_lf077a binnumber_lf077a loc_diff_lf077a
load lf078a_hist window_lf078a binnumber_lf078a loc_diff_lf078a
load lf079a_hist window_lf079a binnumber_lf079a loc_diff_lf079a
load lf080a_hist window_lf080a binnumber_lf080a loc_diff_lf080a
load lf082a_hist window_lf082a binnumber_lf082a loc_diff_lf082a
load lf086a_hist window_lf086a binnumber_lf0786a loc_diff_lf086a
load lf089a_hist window_lf089a binnumber_lf089a loc_diff_lf089a

h = figure;
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
loc_std = zeros(round(2/window_lf077a),1);
bar(binnumber_lf077a,loc_diff_lf077a,'DisplayName','loc_diff');
figure(gcf)
axis([0 0.3 0 2500])
xlabel('OBDA Amplitude')
ylabel('# of Peaks')
title('lf077a')

h = figure;
subplot(321)
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
loc_std = zeros(round(2/window_lf078a),1);
bar(binnumber_lf078a,loc_diff_lf078a,'DisplayName','loc_diff');
figure(gcf)
axis([0 0.3 0 2500])
xlabel('OBDA Amplitude')
ylabel('# of Peaks')
title('lf078a')

subplot(322)
loc_std = zeros(round(2/window_lf079a),1);
bar(binnumber_lf079a,loc_diff_lf079a,'DisplayName','loc_diff');
figure(gcf)
axis([0 0.3 0 2500])
xlabel('OBDA Amplitude')
ylabel('# of Peaks')
title('lf079a')

subplot(323)
loc_std = zeros(round(2/window_lf080a),1);
bar(binnumber_lf080a,loc_diff_lf080a,'DisplayName','loc_diff');
figure(gcf)
axis([0 0.3 0 2500])
xlabel('OBDA Amplitude')
ylabel('# of Peaks')
title('lf080a')

subplot(324)
loc_std = zeros(round(2/window_lf082a),1);
bar(binnumber_lf082a,loc_diff_lf082a,'DisplayName','loc_diff');
figure(gcf)
axis([0 0.3 0 2500])
xlabel('OBDA Amplitude')
ylabel('# of Peaks')
title('lf082a')

subplot(325)
loc_std = zeros(round(2/window_lf086a),1);
bar(binnumber_lf0786a,loc_diff_lf086a,'DisplayName','loc_diff');
figure(gcf)
axis([0 0.3 0 2500])
xlabel('OBDA Amplitude')
ylabel('# of Peaks')
title('lf086a')

subplot(326)%figure;
loc_std = zeros(round(2/window_lf089a),1);
bar(binnumber_lf089a,loc_diff_lf089a,'DisplayName','loc_diff');
figure(gcf)
axis([0 0.3 0 2500])
xlabel('OBDA Amplitude')
ylabel('# of Peaks')
title('lf089a')

