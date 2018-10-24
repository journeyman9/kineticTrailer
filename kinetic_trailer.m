%% Kinetic Trailer Control 
% Journey McDowell (c) 2018

clear; close all; clc;

%% Nominal Parameters
m1 = 8994.4; %[kg]
J1 = 53314.6; %[kg-m2]
a1 = 2.5359; %[m]
b1 = 3.1977; %[m]
L1 = a1 + b1; %[m]
h1 = 2.9691; %[m]
C1 = 399010; %[N/rad] 2 tires
C2 = 1031136; %[N/rad] 8 tires

m2 = 27256.8; %[kg]
J2 = 729225; %[kg-m2]
a2 = 6.2941; %[m]
L2 = 12.192; %[m]
b2 = L2 - a2; %[m]
C3 = 1057244; %[N/rad] 4 tires
v = 10; %[m/s] 2.2352

L1_star = a1 + h1; %[m]
e1 = L1_star - L1; %[m]

orientation = 'up'; % right for horizontal, up for vertical, left for pi, and down for 3pi/2

%% x = [y_d1, psi_d1, theta_d, theta, d1, psi_1]
M = [(m1 + m2)      -m2 * (a2 + h1)          -m2 * a2           0       0       0;
     -m2 * h1    J1 + m2*h1*(h1+a2)          m2 * h1 * a2       0       0       0;  
     -m2 * a2    J2 + m2*a2*(a2+h1)         J2 + m2*a2.^2       0       0       0;
         0               0                       0              1       0       0;
         0               0                       0              0       1       0;
         0               0                       0              0       0       1];
     
K = (-1 / v) * [(C1+C2+C3)              (C1*a1 - C2*b1 - C3*(h1+L2) + (m1+m2)*v.^2)          (-C3*L2)     (-C3*v)   0       0;
                (C1*a1 - C2*b1 - C3*h1) (C1*a1.^2 + C2*b1.^2 + C3*h1*(h1+L2) - m2*h1*v.^2)   (C3*h1*L2)  (C3*h1*v)  0       0;
                (-C3*L2)                     (C3*L2*(h1+L2) - m2*a2*v.^2)                    (C3*L2.^2)  (C3*L2*v)  0       0;
                   0                                       0                                    -v           0      0       0;
                  -v                                       0                                     0           0      0      -v.^2;
                   0                                      -v                                     0           0      0       0];

N = [C1; C1*a1; 0; 0; 0; 0];

A = M \ K;
B = M \ N;
% C = [0 0 0 1 0 0;
%      0 0 0 0 1 0;
%      0 0 0 0 0 1];
C = eye(6);
D = zeros(6, 1);

sys = ss(A, B, C, D);

%% Controllability
controllability = rank(ctrb(A, B));
if controllability == size(A, 1)
    disp('Controllable')
else
    disp('NOT Controllable')
end

%% Observability
observability = rank(obsv(A, C));
if observability == size(A, 1)
    disp('Observable')
    else
        disp('NOT Observable')
end

%% Asymptotically stable if eigenvalues all negative
A_ev = eig(A);
if A_ev < 0
    disp('Asymptotically stable from eigenvalues')
    elseif A_ev <= 0
        disp('Marginally stable from eigenvalues')
    else
        disp('NOT asymptotically stable from eigenvalues')
end

%% LQR Gains
% G = [0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      0 0 0 1 0 0;
%      0 0 0 0 1 0;
%      0 0 0 0 0 1];
% H = zeros(6, 1);
% rho = 1;
% R = 1;
% Q = [0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      0 0 0 1 0 0;
%      0 0 0 0 1 0;
%      0 0 0 0 0 1];
% 
% QQ = G'*Q*G;
% RR = H'*Q*H + rho*R;
% NN = G'*Q*H;
% QN = eye(2); %Match dimension of state
% RN = eye(2); %Match dimension of output
% Bbar = B;
% 
% [K S e] = lqry(sys, QQ, RR, NN);

%% Feedforward
track_vector = csvread('t_oval.txt');
s = track_vector(:, 5);
t = abs(s / v);
curv = [t track_vector(:, 3)];
if v < 0
    yaw_tractor = [t track_vector(:, 4)-pi];
else
    yaw_tractor = [t track_vector(:, 4)];
end
theta_r = [t 0*track_vector(:, 4)];
y_r = [t track_vector(:, 2)];
x_r = [t track_vector(:, 1)];

sim_time = t(end, 1);

%% Simulink
y_IC = 0;
% x = [y_d1, psi_d1, theta_d, theta, d1, psi_1]
switch orientation
    case 'right'
        trailerIC = [track_vector(1,1)-y_IC*sin(0), track_vector(1, 2)+y_IC*cos(0)]; %x_t y_t
        tractorIC = [trailerIC(1) + (a2+h1), trailerIC(2)]; 
        ICs = [0; deg2rad(0); deg2rad(0); deg2rad(0); y_IC; deg2rad(0)]; 
    case 'up'
        trailerIC = [track_vector(1,1)-y_IC*sin(pi/2), track_vector(1, 2)+y_IC*cos(pi/2)]; %x_t y_t
        tractorIC = [trailerIC(1), trailerIC(2) + (a2+h1)];
        ICs = [0; deg2rad(0); deg2rad(0); deg2rad(0); y_IC; deg2rad(90)];
    case 'left'
        trailerIC = [track_vector(1,1)-y_IC*sin(pi), track_vector(1, 2)+y_IC*cos(pi)]; %x_t y_t
        tractorIC = [trailerIC(1) - (a2+h1), trailerIC(2)]; 
        ICs = [0; deg2rad(0); deg2rad(0); deg2rad(0); y_IC; deg2rad(180)];
    case 'down'
        trailerIC = [track_vector(1,1)-y_IC*sin(3*pi/2), track_vector(1, 2)+y_IC*cos(3*pi/2)]; %x_t y_t
        tractorIC = [trailerIC(1), trailerIC(2) - (a2+h1)];
        ICs = [0; deg2rad(0); deg2rad(0); deg2rad(0); y_IC; deg2rad(270)]; 
end

sim('trailer_kinetic.slx')

% %e = [theta_e, d1_e, psi1_e]
% theta_e = error(:, 1);
% d1_e = error(:, 2);
% psi1_e = error(:, 3);

%% Jack-knife check 
hitch_angle = odometry(:, 8);
for terminal_index = 1:length(hitch_angle)
    if hitch_angle(terminal_index) > deg2rad(60)
        fprintf('Jackknifed!\n')
        break
    elseif hitch_angle(terminal_index) < deg2rad(-60)
        fprintf('Jackknifed!\n')
        break
    else
        continue
    end
end

tractor_x = odometry(1:terminal_index, 7);
tractor_y = odometry(1:terminal_index, 6);
trailer_x = odometry(1:terminal_index, 5);
trailer_y = odometry(1:terminal_index, 4);
psi_tractor = odometry(1:terminal_index, 1);
psi_trailer = odometry(1:terminal_index, 3);

%% Plots
% figure
% ax1 = subplot(3, 1, 1);
% plot(tout, rad2deg(theta_e))
% hold on
% plot(tout, 0*linspace(0, length(tout), length(tout))', '--r')
% hold off
% ylabel('\theta [{\circ}]')
% ax2 = subplot(3, 1, 2);
% plot(tout, rad2deg(d1_e))
% hold on
% plot(tout, 0*linspace(0, length(tout), length(tout))', '--r')
% hold off
% ylabel('d_{1e} [m]')
% ax3 = subplot(3, 1, 3);
% plot(tout, psi1_e)
% hold on
% plot(tout, 0*linspace(0, length(tout), length(tout))', '--r')
% hold off
% ylabel('psi_{1e} [{\circ}]')
% 
% xlabel('time [s]')
% legend('response', 'desired')
% movegui('west')
% linkaxes([ax1 ax2, ax3], 'x')

figure
hold on
plot(track_vector(:, 1), track_vector(:, 2), '--r')
plot(trailer_x, trailer_y, 'b') % trailer
plot(tractor_x, tractor_y, 'g') % tractor

plot(trailer_x(1), trailer_y(1), 'ob')
plot(tractor_x(1), tractor_y(1), 'og')
plot(trailer_x(end), trailer_y(end), 'xb')
plot(tractor_x(end), tractor_y(end), 'xg')
axis square
axis equal
xlabel('Position in x [m]')
ylabel('Position in y [m]')
legend('desired path', 'trailer path', 'tractor path')
movegui('east')
hold off

%% Animation
H_c = L2 / 3;
H_t = L2 / 3;

time = 0:.01:tout(terminal_index);
tractor_x = interp1(tout(1:terminal_index), tractor_x, time);
tractor_y = interp1(tout(1:terminal_index), tractor_y, time);
trailer_x = interp1(tout(1:terminal_index), trailer_x, time);
trailer_y = interp1(tout(1:terminal_index), trailer_y, time);
psi_tractor = interp1(tout(1:terminal_index), psi_tractor, time);
psi_trailer = interp1(tout(1:terminal_index), psi_trailer, time);

DCM = @(ang) [cos(ang) -sin(ang) 0;
              sin(ang)  cos(ang) 0;
                0         0      1];

% homogenous transformation
center = @(x, y) [1 0 x;
                  0 1 y;
                  0 0 1];
figure
for i = 1:length(time)
    plot(track_vector(:, 1), track_vector(:, 2), '--r')
    hold on
    
    ang0 = psi_trailer(i);
    ang1 = psi_tractor(i);
    
    % trailer ccw pts starting with top right -- cg
    x_trail = [trailer_x(i)+a2 trailer_x(i)-b2 trailer_x(i)-b2 trailer_x(i)+a2 trailer_x(i)+a2]; 
    y_trail = [trailer_y(i)+H_t/2 trailer_y(i)+H_t/2 trailer_y(i)-H_t/2 trailer_y(i)-H_t/2 trailer_y(i)+H_t/2];
    corners_trail = zeros(5, 3);
    for j = 1:length(x_trail)
        corners_trail(j, 1:3) = center(trailer_x(i), trailer_y(i)) * DCM(ang0) * center(-trailer_x(i), -trailer_y(i)) * [x_trail(j); y_trail(j); 1];
    end
    plot(corners_trail(:, 1), corners_trail(:, 2), 'b-', 'LineWidth', 2)
    
    % tractor ccw pts starting with top right -- cg
    x_trac = [tractor_x(i)+a1 tractor_x(i)-b1 tractor_x(i)-b1 tractor_x(i)+a1 tractor_x(i)+a1]; 
    y_trac = [tractor_y(i)+H_c/2 tractor_y(i)+H_c/2 tractor_y(i)-H_c/2 tractor_y(i)-H_c/2 tractor_y(i)+H_c/2];
    corners_trac = zeros(5, 3);
    for j = 1:length(x_trac)
        corners_trac(j, 1:3) = center(tractor_x(i), tractor_y(i)) * DCM(ang1) * center(-tractor_x(i), -tractor_y(i)) * [x_trac(j); y_trac(j); 1];
    end
    plot(corners_trac(:, 1), corners_trac(:, 2), 'g-', 'LineWidth', 2)
    
    % cg
    plot(trailer_x(i), trailer_y(i), 'b+')
    plot(tractor_x(i), tractor_y(i), 'g+')
    
    % hitch point (should be the same for both)
    hitch_trail = center(trailer_x(i), trailer_y(i)) * DCM(ang0) * center(-trailer_x(i), -trailer_y(i)) * [trailer_x(i)+a2; trailer_y(i); 1];
    plot(hitch_trail(1), hitch_trail(2), 'b*')
    
    hitch_trac = center(tractor_x(i), tractor_y(i)) * DCM(ang1) * center(-tractor_x(i), -tractor_y(i)) * [tractor_x(i)-h1; tractor_y(i); 1];
    plot(hitch_trac(1), hitch_trac(2), 'g*')
    
    xlim([trailer_x(i)-25 trailer_x(i)+25])
    ylim([ trailer_y(i)-25 trailer_y(i)+25])
    xlabel('Position in x [m]')
    ylabel('Position in y [m]')
    drawnow
    hold off
end

