p = 0
r = 90
y = 0

p = (p/360)*2*pi;
r = (r/360)*2*pi;
y = (y/360)*2*pi;

H = [[cos(h) -sin(h) 0];[sin(h) cos(h) 0];[0 0 1]];
P = [[cos(p) 0 -sin(p)];[0 1 0];[sin(p) 0 cos(p)]];
R = [[1 0 0];[0 cos(r) -sin(r)];[0 sin(r) cos(r)]];


rotation_body_intertal = [[cos(y)*cos(p) cos(y)*sin(r)*sin(p)-cos(r)*sin(y) sin(r)*sin(y)+cos(r)*cos(y)*sin(p)];
                         [cos(y)*sin(r)*sin(p)-cos(r)*sin(y)  cos(r)*cos(y)+sin(r)*sin(y)*sin(p) cos(p)*sin(r)];
                         [sin(r)*sin(y)+cos(r)*cos(y)*sin(p) cos(r)*sin(y)*sin(p)-cos(y)*sin(r) cos(r)*cos(p)]]
                     
body_accel = [0 10 0]'                     
                   
interal_accel = (H*P*R)*body_accel
inter = rotation_body_intertal*body_accel

%%
figure;
prh = [pitch roll head];
plot(prh)

%%
A_with = zeros(length(A), 3);

for i=1:length(A)
    p = P_dec(i);
    r = R_dec(i);
    y = head(i);
    
    H = [[cos(h) -sin(h) 0];[sin(h) cos(h) 0];[0 0 1]];
    P = [[cos(p) 0 -sin(p)];[0 1 0];[sin(p) 0 cos(p)]];
    R = [[1 0 0];[0 cos(r) -sin(r)];[0 sin(r) cos(r)]];
    
    Q = P;
    
    rotation_nav_animal = Q;
    
    rotation_nav_animal = [[cos(y)*cos(p) cos(y)*sin(p)*sin(r)-cos(r)*sin(y) sin(r)*sin(y)+cos(r)*cos(y)*sin(p)];
                         [cos(y)*sin(r)*sin(p)-cos(r)*sin(y)  cos(r)*cos(y)+sin(r)*sin(y)*sin(p) cos(p)*sin(r)];
                         [sin(r)*sin(y)+cos(r)*cos(y)*sin(p) cos(r)*sin(y)*sin(p)-cos(y)*sin(r) cos(r)*cos(p)]];
      
    %rotation_animal_nav = rotation_nav_animal';
    rotation_animal_nav = Q';
    body_accel = [A(i,1) A(i,2) A(i,3)]';                        
    intertial_accel = rotation_animal_nav*body_accel;
    
    A_with(i,1) = intertial_accel(1);
    A_with(i,2) = intertial_accel(2);
    A_with(i,3) = intertial_accel(3);

    
end

%%
%figure;
%subplot(211)
%plot(A_dec);
%A_dec_diff = diff(A_dec);
%subplot(212)
%plot(A_dec_diff);
% figure;
% plot(P_dec);
% figure;
% plot(R_dec);
% figure;
% plot(H_dec);


A_with = zeros(length(A_dec), 3);


i_start = 90000;
i_end = 500000;
for i=1:length(A)
    p = P_dec(i);
    r = R_dec(i);
    h = H_dec(i);
    
    H = [[cos(h) -sin(h) 0];[sin(h) cos(h) 0];[0 0 1]];
    P = [[cos(p) 0 -sin(p)];[0 1 0];[sin(p) 0 cos(p)]];
    R = [[1 0 0];[0 cos(r) -sin(r)];[0 sin(r) cos(r)]];
    
    Q = P;
    
    rotation_nav_animal = Q;
    
    rotation_nav_animal = [[cos(y)*cos(p) cos(y)*sin(p)*sin(r)-cos(r)*sin(y) sin(r)*sin(y)+cos(r)*cos(y)*sin(p)];
                         [cos(y)*sin(r)*sin(p)-cos(r)*sin(y)  cos(r)*cos(y)+sin(r)*sin(y)*sin(p) cos(p)*sin(r)];
                         [sin(r)*sin(y)+cos(r)*cos(y)*sin(p) cos(r)*sin(y)*sin(p)-cos(y)*sin(r) cos(r)*cos(p)]];
      
    %rotation_animal_nav = rotation_nav_animal';
    rotation_animal_nav = Q';
    body_accel = [A(i,1) A(i,2) A(i,3)]';                        
    intertial_accel = rotation_animal_nav*body_accel;
    
    A_with(i,1) = intertial_accel(1);
    A_with(i,2) = intertial_accel(2);
    A_with(i,3) = intertial_accel(3);

    
end
%% 
figure;
%subplot(211)
% figure;plot(A_dec);
% plot(A_dec(i_start:i_end,:))
A_x = A_dec(i_start:i_end,1);
A_y = A_dec(i_start:i_end,2);
A_z = A_dec(i_start:i_end,3);
A_sub =[A_x A_y A_z];


% subplot(212)
% figure;
% plot(A_with(i_start:i_end, :))
A_with_x = A_with(i_start:i_end,1);
A_with_y = A_with(i_start:i_end,2);
A_with_z = A_with(i_start:i_end,3)-1;
A_with_small = [A_with_x A_with_y A_with_z];

% A_mag = (A_x.^2+A_y.^2+A_z.^2).^0.5-1;
% A_diff = diff(A_mag);
% figure; plot(A_mag)
% figure; plot(A_diff)

figure;plot(A_sub)
figure;plot(A_with_small)
%figure; plot(A_x)
%figure; plot(A_z, 'red')
mean(A_z)
mean(A_with_z)
figure; plot(A_with_z, 'red'), axis([0 10e5 -1 1.5])
figure; plot(A_with_x)



%figure;
% subplot(211)
% plot(A_dec(i_start:i_end,:))
% subplot(212)
%plot(A_with_small(i_start:i_end, :))


% figure;
% subplot(221)
% plot(A_with(i_start:i_end, :))
% subplot(222)
% plot(A_with_x(i_start:i_end, :) ,'blue')
% subplot(223)
% plot(A_with_y(i_start:i_end, :), 'green')
% subplot(224)
% plot(A_with_z(i_start:i_end, :), 'red')
%%

save 
%A_with = zeros(length(A_dec)); 
% for i=1:2
%     p = P_dec(i);
%     r = R_dec(i);
%     h = H_dec(i);
% 
%     rotation_nav_animal = [[cos(y)*cos(p) cos(y)*sin(r)*sin(p)-cos(r)*sin(y) sin(r)*sin(y)+cos(r)*cos(y)*sin(p)];
%                          [cos(y)*sin(r)*sin(p)-cos(r)*sin(y)  cos(r)*cos(y)+sin(r)*sin(y)*sin(p) cos(p)*sin(r)];
%                          [sin(r)*sin(y)+cos(r)*cos(y)*sin(p) cos(r)*sin(y)*sin(p)-cos(y)*sin(r) cos(r)*cos(p)]];
%                      
%     rotation_animal_nav = rotation_nav_animal';
%     body_accel = [A_dec(i,1) A_dec(i,2) A_dec(i,3)]';                    
%                    
%     interal_accel = rotation_animal_nav*body_accel;
%     A_with(i,1) = interal_accel(1);
%     A_with(i,2) = interal_accel(2);
%     A_with(i,3) = interal_accel(3);
% end
% 
% plot(A_with);
                     