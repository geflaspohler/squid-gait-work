function [ running_average ] = OBDA_window_finder( A_x, MAX, MIN )
% Determine the OBDA window size
MIN = 1;
MAX = 250;
LEN = length(A_x);
value = zeros(ceil(LEN/MIN), 1);
running_average = zeros(MAX-MIN,1);

for window_size = MIN:MAX
    for segment = 1:floor(LEN/window_size)
       lower = ((segment-1) * window_size)+1;
       upper = segment * window_size;
       MEAN = mean(abs(A_x(lower:upper)));
       value(segment) = MEAN;
       %figure; plot(A_x(lower:upper));
       
    end
    %figure; plot(value);
    running_average(window_size,1) = mean(value);
    %figure;plot(running_average);
    value = zeros(ceil(LEN/window_size), 1);
end
%Determined to be 50 samples or ~2 seconds
end

