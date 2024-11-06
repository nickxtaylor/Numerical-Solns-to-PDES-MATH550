clc; clear all; close all;
% Define system parameters
% Measured parameters
m = 0.39; % Mass in kg (replace with your measured value)
k = 24; % Spring constant in N/m (replace with your measured value)
zeta = 0.05; % Damping ratio (replace with your calculated value)
A = 0.2; % Initial amplitude of oscillation (in meters)
% Define the phase shift (in radians) - helps shift sine wave phase horizontally
phi = pi/3; % Example: 45-degree phase shift (adjust as needed)
% Define the experimental data file name - put in same folder as this code
filename = 'Waveform_lvm_compatibility_test.lvm'; % Replace with your file name in same folder as code (or full path)
% Calculate damping coefficient
c = 2 * zeta * sqrt(k * m); % Damping coefficient in Ns/m
% Define the spring-mass-damper system
spring_mass_damper = @(t, x) [x(2); (-c/m)*x(2) - (k/m)*x(1)];
% Calculate natural and damped natural frequencies
omega_n = sqrt(k/m); % Natural frequency in rad/s
omega_d = omega_n * sqrt(1 - zeta^2); % Damped natural frequency
% Define the initial conditions with phase shift - enables horizontal shift of data
x0 = [A * cos(phi); -A * omega_d * sin(phi)]; % Initial displacement and velocity
% Time span for simulation
tspan = [0 5]; % Based on 5 second experimental data duration
% Solve the system of equations using ode45
[t, x] = ode45(spring_mass_damper, tspan, x0);
% Plot the theoretical response
figure(1);
plot(t, x(:,1), 'b');
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Theoretical Spring-Mass-Damper Response');
legend('Theoretical');
% Open the file
fileID = fopen(filename, 'r');
% Initialize variables to store data
displacement_exp = [];
% Read through the lvm file line by line
while ~feof(fileID)
line = fgetl(fileID);
% Skip header lines (if any) or lines that don't contain data
if startsWith(line, '#') || isempty(line)
continue;
end
% Read numerical data from each line
data = textscan(line, '%f %f', 'Delimiter', '\t');
% If data is valid, append to arrays
if ~isempty(data{1})
displacement_exp = [displacement_exp; data{2}]; % Second column: displacement
end
end
% Close the file
fclose(fileID);
% Define time vector for 5 second experiment
time_exp = linspace(0,5,length(displacement_exp));
% Plot the experimental data
figure(2);
plot(time_exp, displacement_exp, 'r');
hold on;
plot(t, x(:,1), 'b');
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Experimental Data from .lvm File');
legend('Experimental Displacement');
hold off;
%% MATLAB Code End