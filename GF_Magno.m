load PIV_decimated
%%
time = 1:length(M_dec);
M_x = M_dec(:, 1);
M_y = M_dec(:, 2);
M_z = M_dec(:, 3);
time = time./200;
M_mag = sqrt(M_x.^2+M_y.^2+M_z.^2);
figure; plot(time, M_dec); 
hold on; plot(time, M_mag, 'y');

%%
