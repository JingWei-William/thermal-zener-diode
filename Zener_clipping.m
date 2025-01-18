%% Temperature Clipping Circuit with a Thermal Zener Diode
%  Assuming Reverse Blocking at Each Step

% Simulation parameters
T0 = 300; % Base temperature (K)
T_breakdown = 40; % Reverse temperature difference for breakdown (K)
R_source = 1; % Source thermal resistance (K/W)
R_forward = 0.2; % Thermal resistance in forward bias (K/W)
R_reverse_high = 10; % Thermal resistance in reverse bias before breakdown (K/W)
R_reverse_breakdown = 0.8; % Thermal resistance in reverse bias after breakdown (K/W)
dt = 0.01; % Time step (s)
total_time = 10; % Total simulation time (s)

% Parameters for input sinusoidal temperature wave
TR_amplitude = 50; % Amplitude of sinusoidal temperature wave (K)
TR_frequency = 1; % Frequency of sinusoidal temperature wave (Hz)

% Time and input temperature wave
time = 0:dt:total_time;
TR = T0 + TR_amplitude * sin(2 * pi * TR_frequency * time); % Input sinusoidal temperature wave

% Initialize output temperature
T_out = zeros(size(time)); % Output temperature

% Simulation loop
for i = 1:length(time)
    % Step 1: Assume initial state (reverse high resistance)
    R_diode = R_reverse_high;

    % Step 2: Calculate initial temperature difference across the diode
    delta_T_source = TR(i) - T0;
    delta_T_diode = delta_T_source * R_diode / (R_source + R_diode);

    % Step 3: Re-evaluate state based on calculated delta_T_diode
    if delta_T_diode > 0
        % Reverse heat flow
        if delta_T_diode >= T_breakdown
            % Reverse breakdown
            R_diode = R_reverse_breakdown;
        else
            % Reverse blocking
            R_diode = R_reverse_high;
        end
    else
        % Forward heat flow
        R_diode = R_forward;
    end

    % Step 4: Recalculate temperature difference with updated resistance
    delta_T_diode = delta_T_source * R_diode / (R_source + R_diode);

    % Step 5: Calculate output temperature
    T_out(i) = T0 + delta_T_diode;
end

% Plot results
figure;
plot(time, TR, 'b-', 'LineWidth', 2); % Input temperature wave
hold on;
plot(time, T_out, 'r-', 'LineWidth', 2); % Output temperature
yline(T0 + TR_amplitude, 'k--', 'LineWidth', 1.5); % Upper boundary of input wave
yline(T0 - TR_amplitude, 'k--', 'LineWidth', 1.5); % Lower boundary of input wave
hold off;
grid on;
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('Input Temperature Wave', 'Output Temperature', 'Upper Boundary', 'Lower Boundary');
title('Heat Clipping Circuit with a Thermal Zener Diodes');
