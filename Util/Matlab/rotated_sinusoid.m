clc; clear all; close all;

% % Script to check the sinusoid path

% Time settings
end_pt = 1/4;

dt = 100;
t = linspace(0,1,dt);

% Sine wave
amp = 1;
freq = dt/end_pt;
ph = 0;
y = amp*sin(2*pi*freq*t + ph);

th = deg2rad(45);

R = [cos(th) -sin(th);
    sin(th) cos(th)];

x = [t; y];
traj = R*x;

plot(traj(1,:), traj(2,:));
% plot(t,y);
grid on;
axis([-5 5 -5 5]);
axis square;