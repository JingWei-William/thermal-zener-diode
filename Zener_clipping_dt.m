%% Temperature Clipping Circuit with a Thermal Zener Diode
%  Using Previous Step's Temperature Difference

% Simulation parameters
T0 = 300; % Base temperature (K)
T_breakdown = 40; % Reverse temperature difference for breakdown (K)
R_source = 1; % Source thermal resistance (K/W)
R_forward = 0.2; % Thermal resistance in forward bias (K/W)
R_reverse = 10; % Thermal resistance in reverse bias before breakdown (K/W)
R_breakdown = 0.4; % Thermal resistance in reverse bias after breakdown (K/W)
dt = 0.01; % Time step (s)
total_time = 5; % Total simulation time (s)

% Parameters for input sinusoidal temperature wave
TR_amplitude = 50; % Amplitude of sinusoidal temperature wave (K)
TR_frequency = 1; % Frequency of sinusoidal temperature wave (Hz)

% Time and input temperature wave
time = 0:dt:total_time;
TR = T0 + TR_amplitude * sin(2 * pi * TR_frequency * time); % Input sinusoidal temperature wave

% Initialize output temperature and previous temperature difference
T_out = zeros(size(time)); % Output temperature
delta_T_diode_prev = 0; % Initial temperature difference across the diode (at time step 0)

% Simulation loop
for i = 1:length(time)
    % Step 1: Calculate temperature difference across the diode from previous step
    delta_T_source = TR(i) - T0;
    
    % Step 2: Determine diode's state based on previous temperature difference
    if delta_T_diode_prev >= T_breakdown
        % Reverse breakdown
        R_diode = R_breakdown;
    elseif delta_T_source < 0
        % Forward heat flow
        R_diode = R_forward;
    else
        % Reverse blocking
        R_diode = R_reverse;
    end

    % Step 3: Calculate temperature difference across the diode with current resistance
    delta_T_diode = delta_T_source * R_diode / (R_source + R_diode);

    % Step 4: Calculate output temperature
    T_out(i) = T0 + delta_T_diode;

    % Step 5: Store the current temperature difference for the next time step
    delta_T_diode_prev = delta_T_diode;
end

% Plot results
figure;
plot(time, TR, 'b-', 'LineWidth', 2); % Input temperature wave
hold on;
plot(time, T_out, 'r-', 'LineWidth', 2); % Output temperature
yline(T0 + T_breakdown, 'k--', 'LineWidth', 1.5); % Temperature difference for breakdown
hold off;
grid on;
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('Input Temperature Wave', 'Output Temperature', 'T breakdown');
title('Temperature Clipping Circuit with a Thermal Zener Diode');
