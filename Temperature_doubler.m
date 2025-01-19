%% Temperature Dynamics with Opposite Thermal Diodes
%  "Temperature Doubler"

% Simulation parameters
T0 = 300; % Initial temperature of the system (K)
C1 = 5; % Thermal capacitance of node 1 (J/K)
C2 = 5; % Thermal capacitance of node 2 (J/K)
R_forward = 0.2; % Thermal resistance in forward bias (K/W)
R_reverse = 10; % Reduced high thermal resistance in reverse bias (K/W)
R_breakdown = 0.4; % Thermal resistance in reverse breakdown (K/W)
T_breakdown = 50; % Reverse temperature difference for breakdown (K)
dt = 0.01; % Time step (s)
total_time = 100; % Total simulation time (s)

% Parameters for TR(t)
TR_amplitude = 50; % Amplitude of TR oscillation (K)
TR_frequency = 0.05; % Frequency of TR oscillation (Hz)

% Initializing variables
time = 0:dt:total_time;
T1 = zeros(size(time)); % Temperature at node 1 (K)
T2 = zeros(size(time)); % Temperature at node 2 (K)
TR = T0 + TR_amplitude * sin(2 * pi * TR_frequency * time); % External oscillation
T1(1) = T0; % Initial temperature at node 1
T2(1) = T0; % Initial temperature at node 2

% Simulation loop
for i = 2:length(time)
    % Calculate temperature differences for each path
    delta_T1 = TR(i-1) - T1(i-1); % Path TR -> T1
    delta_T2 = T2(i-1) - TR(i-1); % Path T2 -> TR

    % Determine thermal resistance for TR -> T1 path (positive diode)
    if delta_T1 > 0
        R1 = R_forward; % Forward bias
    elseif abs(delta_T1) < T_breakdown
        R1 = R_reverse; % Reverse bias (no breakdown)
    else
        R1 = R_breakdown; % Reverse breakdown
    end

    % Determine thermal resistance for TR -> T2 path (negative diode)
    if delta_T2 > 0
        R2 = R_forward; % Forward bias
    elseif abs(delta_T2) < T_breakdown
        R2 = R_reverse; % Reverse bias (no breakdown)
    else
        R2 = R_breakdown; % Reverse breakdown
    end

    % Update temperatures using thermal circuit equations
    dT1 = (TR(i-1) - T1(i-1)) / (R1 * C1); % TR to T1
    dT2 = (TR(i-1) - T2(i-1)) / (R2 * C2); % TR to T2

    T1(i) = T1(i-1) + dT1 * dt; % Update T1
    T2(i) = T2(i-1) + dT2 * dt; % Update T2

end

% Calculate T1 - T2 difference
delta_T12 = T1 - T2;

% Plot results
figure;
plot(time, T1, 'r-', 'LineWidth', 2); % Node 1 temperature
hold on;
plot(time, T2, 'b-', 'LineWidth', 2); % Node 2 temperature
plot(time, TR, 'k--', 'LineWidth', 1.5); % External temperature TR(t)
plot(time, delta_T12, 'g-', 'LineWidth', 2); % T1 - T2 difference
hold off;
grid on;
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('T1 (Node 1)', 'T2 (Node 2)', 'TR(t)', 'T1 - T2');
title('Temperature Dynamics with Opposite Thermal Diodes');
