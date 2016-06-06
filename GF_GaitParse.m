function [average_wait, average_peaks, finning_matrix, swimming_matrix, wait_matrix, count_matrix] = Gait_parse( data_matrix, resting_mag)

bin_len = 25;
%[~,std_matrix] = findpeaks(data_matrix,'MINPEAKDIST',bin_len);
std_matrix = zeros(ceil(length(data_matrix)/bin_len)+5, 1);

time_a = ((1:length(data_matrix))/25)';
h = figure;
set(h,'Color',[1,1,1], 'Position', [1 1 840 330]);
plot(time_a,data_matrix,'k')

hold on
%plot(time_a, 0.1, 'green');
fin_count = 0;
wait_count = 0;
count_pure_fin = 0;

state = 0;
next_state = 0;

NOT_STARTED = 0;
MAYBE_FINNING = 1;
FINNING = 2;
MAYBE_WAITING = 3;
WAITING = 4;
MAYBE_DONE_WAITING = 5;
WAIT_COUNTER = 0;
HOLD_COUNTER = 0;

samples = [1:length(data_matrix)]';
data_matrix = [data_matrix samples];
hold_matrix = [];
hold_wait_matrix = [];
finning_matrix = [];
swimming_matrix = [];
%movement_matrix = [];
%movement_index = 1;
wait_matrix = zeros(length(data_matrix),1);
wait_index = 1;
count_matrix = zeros(length(data_matrix),1);
count_index = 1;

for i=1:bin_len:length(data_matrix)
    if i+bin_len < length(data_matrix)
        submatrix = data_matrix(i:i+bin_len-1)';
    else
        submatrix= data_matrix(i:length(data_matrix))';
    end
    
    std_matrix(i, 1) = std(submatrix);
    %plot(time_a(i),std_matrix(i,1),'.r');
    
    if state == NOT_STARTED
        %text(time_a(i),data_matrix(i,1),'F','BackgroundColor', [.7 .9 .7]);
        hold_matrix = [hold_matrix; submatrix];
        %plot(time_a(i),data_matrix(i,1), '.r', 'markersize', 20);
        if std_matrix(i,1) > resting_mag
            next_state = MAYBE_FINNING;
        else
            next_state = NOT_STARTED;
        end
    elseif state == MAYBE_FINNING
        hold_matrix = [hold_matrix; submatrix];
        %text(time_a(i),data_matrix(i,1),'F','BackgroundColor', [.7 .9 .7]);
        %fin_count = fin_count + 1;
        %plot(time_a(i),data_matrix(i,1), '.r' ,'markersize', 20);
        if std_matrix(i,1) > resting_mag
            next_state = FINNING;
        else
            next_state = MAYBE_DONE_WAITING;
        end
    elseif state == FINNING
        hold_matrix = [hold_matrix; submatrix];
        fin_count = fin_count + 1;
        %plot(time_a(i),data_matrix(i,1), '.r','markersize', 20);
        %hold_matrix(hold_index:hold_index+bin_len, 1) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 1);
        %hold_matrix(hold_index:hold_index+bin_len, 2) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 2);
        if std_matrix(i,1) > resting_mag
            next_state = FINNING;
        else
            next_state = MAYBE_WAITING;
        end
    elseif state == MAYBE_WAITING
        hold_wait_matrix = [hold_wait_matrix; submatrix];
        fin_count = fin_count + 1;
        %plot(time_a(i),data_matrix(i,1), '.g');
        %text(time_a(i),data_matrix(i,1),'MW','BackgroundColor', [.7 .9 .7]);
        %hold_matrix(hold_index:hold_index+bin_len, 1) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 1);
        %hold_matrix(hold_index:hold_index+bin_len, 2) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 2);
        if std_matrix(i,1) > resting_mag
            next_state = FINNING;
        elseif WAIT_COUNTER < 1
            WAIT_COUNTER = WAIT_COUNTER + 1;
            next_state = MAYBE_WAITING;
            wait_count = 0;
        else
            next_state = WAITING;
            WAIT_COUNTER = 0;
            count_matrix(count_index) = fin_count;
            count_index = count_index + 1;
            fin_count = 0;
        end   
    elseif state == WAITING
        hold_wait_matrix = [hold_wait_matrix; submatrix];
        %plot(time_a(i),data_matrix(i,1), '.g');
        %text(time_a(i),data_matrix(i,1),'W','BackgroundColor', [.7 .9 .7]);
        wait_count = wait_count + 1;
        %hold_matrix(hold_index:hold_index+bin_len, 1) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 1);
        %hold_matrix(hold_index:hold_index+bin_len, 2) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 2);
        
        if std_matrix(i,1) > resting_mag
            next_state = MAYBE_DONE_WAITING;
%         elseif HOLD_COUNTER < 1
%             HOLD_COUNTER = HOLD_COUNTER + 1;
%             next_state = WAITING;
        else
            next_state = WAITING;
            HOLD_COUNTER = 0;
        end       
    elseif state == MAYBE_DONE_WAITING
        hold_wait_matrix = [hold_wait_matrix; submatrix];
        %plot(time_a(i),data_matrix(i,1), '.g');
        %hold_matrix(hold_index:hold_index+bin_len, 1) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 1);
        %hold_matrix(hold_index:hold_index+bin_len, 2) = data_matrix(std_matrix(i):std_matrix(i)+bin_len, 2);
        if std_matrix(i,1) > resting_mag
            if HOLD_COUNTER > 1
                next_state = MAYBE_FINNING;
                wait_index = wait_index + 1;     
                wait_matrix(wait_index) = wait_count;
                if wait_count > 3
                    %text(time_a(i),data_matrix(i,1),strcat('W',int2str(wait_count)),'BackgroundColor', [.7 .2 .7]);
                    plot(time_a(i),data_matrix(i,1), '.r', 'markersize', 20);
                    swimming_matrix = [swimming_matrix; hold_wait_matrix];
                    finning_matrix = [finning_matrix; hold_matrix];
                    hold_matrix = [];
                    hold_wait_matrix = [];
                else
                    count_pure_fin = count_pure_fin + 1;
                    finning_matrix = [finning_matrix; hold_matrix; hold_wait_matrix];
                    hold_matrix = [];
                    hold_wait_matrix = [];
                    %text(time_a(i),data_matrix(i,1),strcat('W',int2str(wait_count)),'BackgroundColor', [.7 .9 .7]);
                end
                wait_count = 0;
            else
               next_state = MAYBE_DONE_WAITING;
               HOLD_COUNTER = HOLD_COUNTER + 1;
            end
        else
            next_state = WAITING;
        end         
    end
    state = next_state;
end

if wait_count > 3
    %text(time_a(i),data_matrix(i,1),strcat('W',int2str(wait_count)),'BackgroundColor', [.7 .2 .7]);
    plot(time_a(i),data_matrix(i,1), '.r', 'markersize', 20);
    swimming_matrix = [swimming_matrix; hold_wait_matrix];
    finning_matrix = [finning_matrix; hold_matrix];
else

    finning_matrix = [finning_matrix; hold_matrix; hold_wait_matrix];
    %text(time_a(i),data_matrix(i,1),strcat('W',int2str(wait_count)),'BackgroundColor', [.7 .9 .7]);
end  

average_wait = sum(wait_matrix)/(wait_index-1);
average_peaks = sum(count_matrix)/(count_index-1);

