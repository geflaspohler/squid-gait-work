function [window, binnumber, loc_diff] = Histogram( data_matrix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Remove parts of matrix below zero
data_matrix_zero = (data_matrix < 0.02);
data_matrix(data_matrix_zero) = 0;

%Set up window/peak finding parameters

Tmin = 15; 
window_width = 0.01/3;
window = (0.005)/2;
binnumber = linspace(0,2,round(2/window))+window_width;
loc_diff = zeros(round(2/window), 1);

%Bin the peaks for the histogram and make the plot
    for i = 1:round(2/window)
        [~,Loc_1] = findpeaks(data_matrix,'MINPEAKHEIGHT',((window_width-window)+window*i),'MINPEAKDIST',Tmin);
        [~,Loc_2] = findpeaks(data_matrix,'MINPEAKHEIGHT',((window_width-window)+window*(i+1)),'MINPEAKDIST',Tmin);
        loc_diff(i,1) = length(Loc_1)-length(Loc_2);
    end
end

