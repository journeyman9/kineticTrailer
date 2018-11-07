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
v = -4.5; %[m/s] 2.2352

L1_star = a1 + h1; %[m]
e1 = L1_star - L1; %[m]

%% x = [y_d1, psi_d1, psi_d2, psi_1, psi_2, y2] where d_d2 = y_d1 - (h1+a2)\psi_d_1 - a2\theta_d + v\psi_1 (with -v \theta)
M = [(m1 + m2)     -m2 * h1          -m2 * a2          0       0    0;
     -m2 * h1    J1 + m2*h1.^2       m2 * h1 * a2      0       0    0;  
     -m2 * a2      m2*h1*a2         J2 + m2*a2.^2      0       0    0;
         0               0              0              1       0    0;
         0               0              0              0       1    0;
         0               0              0              0       0    1];
     
S = (-1 / v) * [      (C1+C2+C3)                  (C1*a1 - C2*b1 - C3*h1 + (m1+m2)*v.^2)      -C3*L2     C3*v      -C3*v        0;
                (C1*a1 - C2*b1 - C3*h1)       (C1*a1.^2 + C2*b1.^2 + C3*h1.^2 - m2*h1*v.^2)  C3*h1*L2   -C3*h1*v    C3*h1*v     0;
                       -C3*L2                     (C3*L2*h1 - m2*a2*v.^2)                    C3*L2.^2  -C3*L2*v  C3*L2*v        0;
                   0                                      -v                                     0         0           0        0;
                   0                                       0                                    -v         0           0        0;
                  -v                                      -h1*(-v)                              -a2*(-v)  -v*(-v)     2*v*(-v)  0];
                  
E = [C1; C1*a1; 0; 0; 0; 0];

A = M \ S;
B = M \ E;
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
steer_max = 45; %[degrees]
 
G = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
H = zeros(3, 1);
rho = 1;
R = 1;
Q = eye(3);
% R = 1 / (deg2rad(steer_max).^2);
% Q = [1/(deg2rad(5).^2)       0                       0;
%      0                   1/(deg2rad(5).^2)           0;
%      0                        0                1/(1.^2)];

QQ = G'*Q*G;
RR = H'*Q*H + rho*R;
NN = G'*Q*H;
QN = eye(2); %Match dimension of state
RN = eye(2); %Match dimension of output
Bbar = B;

[K S e] = lqr(sys, QQ, RR, NN);
% [est, L, P] = kalman(ss(A, [B Bbar], C, zeros(2,2)), QN, RN);

%% Set Point Control
Q_sp = [A, B; G, H];
[n, n] = size(A);
[l, p] = size(G); %Number of controlled outputs
m = 1; %Number of process inputs, or just inputs
MM = pinv(Q_sp); %psuedo inverse because matrix is not square
F = MM(1:n, end-l+1:end);
N = MM(end-m+1:end, end-l+1:end);

%% Trajectory Generation
track_vector = csvread('t_horizontal.txt');
if v < 0
    fprintf('Going Backwards!')
    track_vector(:, 4) = track_vector(:, 4) + pi;
end

hitch_max = 90; %[degrees]

%% Simulink
y_IC = 0;
psi_2_IC = deg2rad(0) + track_vector(1, 4);
hitch_IC = deg2rad(0);

look_ahead = 0; %indices

% x = [y_d1, psi_d1, theta_d, psi_1, psi_2, y2]
psi_1_IC = -hitch_IC + psi_2_IC;

trailerIC = [track_vector(1, 1)-y_IC*sin(psi_2_IC)+b2*cos(psi_2_IC), track_vector(1, 2)+y_IC*cos(psi_2_IC)+b2*sin(psi_2_IC)]; %x2, y2
tractorIC = [trailerIC(1)+a2*cos(psi_2_IC)+h1*cos(psi_1_IC), trailerIC(2)+a2*sin(psi_2_IC)+h1*sin(psi_1_IC)]; %x1, y1
ICs = [0; deg2rad(0); deg2rad(0); psi_1_IC; psi_2_IC; y_IC];

sim('trailer_kinetic.slx')

% x = [y_d1, psi_d1, psi_d2, psi_1, psi_2, y2]
psi_1_e = error(:, 4);
psi_2_e = error(:, 5);
d1_e = error(:, 6);

% state
psi_1 = state(:, 4);
psi_2 = state(:, 5);
d1 = state(:, 6);

%% Jack-knife check 
hitch_angle = odometry(:, 8);

for terminal_index = 1:length(hitch_angle)
    if hitch_angle(terminal_index) > deg2rad(hitch_max)
        fprintf('Jackknifed! theta = %4.2f \n', rad2deg(hitch_angle(terminal_index)))
        break
    elseif hitch_angle(terminal_index) < deg2rad(-hitch_max)
        fprintf('Jackknifed! theta = %4.2f \n', rad2deg(hitch_angle(terminal_index)))
        break
    else
        continue
    end
end

%% Goal check
if goal(end) == 1
    fprintf('GOAL with d = %4.2f m and psi = %4.2f degrees\n', d_goal(end), rad2deg(psi_goal(end)))
else
    [minimum, best_index] = min(d_goal(1:terminal_index));
    fprintf('TIMES UP. Closest: d = %4.2f m and psi = %4.2f degrees\n', minimum, rad2deg(psi_goal(best_index)))
end

tractor_x = odometry(1:terminal_index, 7);
tractor_y = odometry(1:terminal_index, 6);
trailer_x = odometry(1:terminal_index, 5);
trailer_y = odometry(1:terminal_index, 4);
psi_tractor = odometry(1:terminal_index, 1);
psi_trailer = odometry(1:terminal_index, 3);

%% Plots
figure
ax1 = subplot(3, 1, 1);
plot(tout, rad2deg(psi_1_e))
hold on
plot(tout, 0*linspace(0, length(tout), length(tout))', '--r')
line([tout(terminal_index) tout(terminal_index)], [max(rad2deg(psi_1_e)) min(rad2deg(psi_1_e))],'Color','red')
hold off
ylabel('\psi_{1e} [{\circ}]')
ax2 = subplot(3, 1, 2);
plot(tout, rad2deg(psi_2_e))
hold on
plot(tout, 0*linspace(0, length(tout), length(tout))', '--r')
line([tout(terminal_index) tout(terminal_index)], [max(rad2deg(psi_2_e)) min(rad2deg(psi_2_e))],'Color','red')
hold off
ylabel('\psi_{2e} [{\circ}]')
ax3 = subplot(3, 1, 3);
plot(tout, d1_e)
hold on
plot(tout, 0*linspace(0, length(tout), length(tout))', '--r')
line([tout(terminal_index) tout(terminal_index)], [max(d1_e) min(d1_e)],'Color','red')
hold off
ylabel('d_{1e} [m]')


xlabel('time [s]')
legend('response', 'desired')
movegui('west')
linkaxes([ax1 ax2, ax3], 'x')

figure
ax1 = subplot(3, 1, 1);
plot(tout, rad2deg(psi_1))
hold on
plot(tout, rad2deg(r(:, 1)), '--r')
line([tout(terminal_index) tout(terminal_index)], [max(rad2deg(psi_1)) min(rad2deg(psi_1))],'Color','red')
hold off
ylabel('\psi_{1} [{\circ}]')
ax2 = subplot(3, 1, 2);
plot(tout, rad2deg(psi_2))
hold on
plot(tout, rad2deg(r(:, 2)), '--r')
line([tout(terminal_index) tout(terminal_index)], [max(rad2deg(psi_2)) min(rad2deg(psi_2))],'Color','red')
hold off
ylabel('\psi_{2} [{\circ}]')
ax3 = subplot(3, 1, 3);
plot(tout, d1)
hold on
plot(tout, r(:, 3), '--r')
line([tout(terminal_index) tout(terminal_index)], [max(d1) min(d1)],'Color','red')
hold off
ylabel('d_{1} [m]')
xlabel('time [s]')
legend('response', 'desired')
movegui('south')
linkaxes([ax1 ax2, ax3], 'x')

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
%     frames(i) = getframe(gcf);
end

% video = VideoWriter('x_same.avi', 'Motion JPEG AVI');
% video.Quality = 50;
% open(video)
% writeVideo(video, frames);
% close(video)
