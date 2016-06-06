%% Video Plotter

%% GOPRO001.MP4 lf086a
%load lf086a_data;
close all;
film_start = 61318-(6*60+30.7-0.159);

% Finning with direction change
% M_S = 5;
% S_S = 55;
% M_E = 6;
% S_E = 00;

% Jetting mantle first
% M_S = 6;
% S_S = 23;
% M_E = 6;
% S_E = 28;

% Flapping mantle first
% M_S = 6;
% S_S = 00;
% M_E = 6;
% S_E = 05;
%     
% Total thing
% M_S = 6;
% S_S = 15;
% M_E = 6;
% S_E = 30;

% Total thing
M_S = 6;
S_S = 35;
M_E = 6;
S_E = 42;

% t = M(:,1);
% a = M(:,2);
% v = M(:,3);
% 
% zero1 = zeros(
% v = [ zero1 v zero2];


%GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda, p_deg, t, v);
GF_SegmentPlotter_NoV(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda, p_deg);


%%
M = csvread('velcoity_600_630.csv');
t = M(:,1);
a = M(:,2);
v = M(:,3);
t = t+1;
v = v(1:(length(v)-5));
t = t(1:(length(t)-5));
a = a(1:(length(a)-5));
figure; plot(t, v, t, a);

%% GP010001_lf086a_140327.MP4 lf086a 
load lf086a_data;
close all;
film_start = 29*60+33-3;

M_S = 1;
S_S = 40;
M_E = 3;
S_E = 27;

GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda, P_dec, t, v);
%M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda, P_dec, t, v
%% GOPRO002.MP4 lf086a
load lf086a_data;
close all;
film_start = 12*60+08-4;
M_S = 15;
S_S = 04;
M_E = 15;
S_E = 21;

GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda, P_dec, t, v);

%% GP010002.MP4 lf086a
load lf086a_data;
close all;
film_start = 29*60+33-4;

M_S = 3;
S_S = 19;
M_E = 3;
S_E = 23;

GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda);

%% GOPR002_lf080a.MP4 lf080a
close all;
load lf080a_data;
film_start = 16*3600+49*60+37-8;
%film_start = 39648.00;

M_S = 14;
S_S = 28;
M_E = 14;
S_E = 48;

GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda);

%% GO02001_lf082a.MP4 lf082a
close all;
load lf082a_data;

%NOTE: lf082a does not have good pressure data; therefore, up/down cannot
%be determined.  Vido SYNC is bad; ignore

film_start = 22*3600+55*60+57;
%film_start = 39648.00;

M_S = 8;
S_S = 34;
M_E = 8;
S_E = 41;

GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, mean_obda);


%% CRX-700 0036 lf086a, creat animation
%close all;
%Start time and end time, minutes, second, in 00036.MTS
film_start = 61227-457-6.2;

%For trial
M_S = 15;
S_S = 31;
M_E = 15;
S_E = 45;

GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx, fqs, mean_obda);

%% CRX-700 0034 lf086a

film_start = 61227-(7*60+36);
M_S = 1;
S_S = 50;
M_E = 2;
S_E = 30;

GF_SegmentPlotter(M_S, S_S, M_E, S_E, film_start, dynamics, pressure_normalized, idx);