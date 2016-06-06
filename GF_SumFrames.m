function [ frame_totals, peak_totals] = Sum_frames( frame_size, test_matrix)
%Sum_frames: Iterates through a data_matrix input and finds the
%average value over the frame.  Returns frame_totals, a matrix containing 
%length(data_matrix)/frame_size sum values
    
    sample_step = frame_size;
    frame_totals = zeros(length(test_matrix),1);
    peak_totals = zeros(length(test_matrix),1);
    stop = floor(length(test_matrix)/sample_step);
    

    for iter=1:stop
        sub_sample = test_matrix((iter-1)*sample_step+1:iter*sample_step);

        X = max(findpeaks(abs(sub_sample)));
        if(isempty(X))
           frame_totals(iter*sample_step) = 0; 
        else
            frame_totals(iter*sample_step) = X; 
        end
    end

end
