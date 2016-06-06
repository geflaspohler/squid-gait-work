function [A_all, P_all, R_all, H_all, pressure, M_all] = GF_ConcatFiles( file_base, file_term, sample_number, starting_number)
%Read in and concatonate variables
A_all = [];
P_all = [];
R_all = [];
H_all = [];
M_all = [];
pressure = [];

file_base_less10 = strcat(file_base, '00');
file_base_less100 = strcat(file_base, '0');
file_base_less1000 = file_base;

file_names = cell(sample_number, 1);
for i=starting_number:starting_number+sample_number-1
    if i < 10
        file_names{i} = strcat(file_base_less10, num2str(i), file_term)
    elseif i < 100
        file_names{i} = strcat(file_base_less100, num2str(i), file_term)
    else
        file_names{i} = strcat(file_base_less1000, num2str(i), file_term)
    end
end

for i=starting_number:starting_number+sample_number-1
    file_string = strcat(file_names{i});
    load(file_string);
    
    A_all = [A_all; A];
    P_all = [P_all; pitch];
    R_all = [R_all; roll];
    H_all = [H_all; head];
    M_all = [M_all; M];
    pressure = [pressure; p];
end
end

