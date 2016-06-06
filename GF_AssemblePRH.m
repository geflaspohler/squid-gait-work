%% Assemble data matrixes from files
% %User file input
file_prefix = 'data';
file_suffix = 'prh200.mat';
file_number = 1;
starting_number = 1;

% M_S = 10;
% S_S = 30;
% M_E = 12;
% S_E = 30;


%Load files and concatante into large matrixes A_all, P_all, R_all, H_all
%Concat_files(file_prefix, file_suffix, file_number);
[A_all, P_all, R_all, H_all, pressure, M_all] = GF_ConcatFiles(file_prefix, file_suffix, file_number, starting_number);
% 
% %Decimate data to make for faster processing
A_dec = decdc(A_all,1);
P_dec = decdc(P_all,1);
R_dec = decdc(R_all,1);
H_dec = decdc(H_all,1);
M_dec = decdc(M_all, 1);
pressure = decdc(pressure,1);

save flips_decimated A_dec pressure P_dec R_dec M_dec pressure file_prefix