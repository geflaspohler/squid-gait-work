function [ interal_accel ] = Rotation_correct( A, pitch, roll)

pitch_roll = [pitch roll];

figure;
plot(pitch_roll);

'Click a single point to calibrate tag orientation on animal'
[x,y] = ginput(1);
close all

pitch_user = pitch_roll(floor(x),1);
roll_user = pitch_roll(floor(x),2);

p = -pitch_user;
r = -roll_user;

P = [[cos(p) 0 -sin(p)];[0 1 0];[sin(p) 0 cos(p)]];
R = [[1 0 0];[0 cos(r) -sin(r)];[0 sin(r) cos(r)]];
    
Q = P*R;
    
rotation_nav_animal = Q;                    
                   
interal_accel = A*(P*R);

figure;
subplot(211)
plot(A), title('Before rotation correction');
subplot(212)
plot(interal_accel), title('After rotation correction');

end

