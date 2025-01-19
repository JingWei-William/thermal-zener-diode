%% Temperature Clipping Circuit with Double Thermal Zener Diodes
%  "Poor Man's Square Wave Generator"

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

% Initialize output temperature and previous temperature differences
T_out = zeros(size(time)); % Output temperature
delta_T_diode1_prev = 0; % Initial temperature difference across diode 1 (at time step 0)
delta_T_diode2_prev = 0; % Initial temperature difference across diode 2 (at time step 0)

% Simulation loop
for i = 1:length(time)
    % Step 1: Calculate the current temperature difference across the diodes
    delta_T_source = TR(i) - T0;

    % Step 2: Determine the state of diode 1(positive Zener diode) based on previous temperature difference
    if delta_T_diode1_prev <= -T_breakdown
        % Negative Zener Diode: Reverse breakdown
        R_diode1 = R_breakdown;
    elseif delta_T_source >= 0
        % Negative Zener Diode: Forward heat flow
        R_diode1 = R_forward;
    else
        % Negative Zener Diode: Reverse blocking
        R_diode1 = R_reverse;
    end
    
    % Step 3: Determine the state of diode 2(negative Zener Diode) based on previous temperature difference
    if delta_T_diode2_prev >= T_breakdown
        % Positive Zener Diode: Reverse breakdown
        R_diode2 = R_breakdown;
    elseif delta_T_source <= 0
        % Positive Zener Diode: Forward heat flow
        R_diode2 = R_forward;
    else
        % Positive Zener Diode: Reverse blocking
        R_diode2 = R_reverse;
    end

    % Step 4: Calculate the temperature across both diodes using voltage division
    delta_T_diode1 = delta_T_source * R_diode1 / (R_source + R_diode1 + R_diode2);
    delta_T_diode2 = delta_T_source * R_diode2 / (R_source + R_diode1 + R_diode2);
    
    % Step 5: Output temperature is the sum of both diodes' temperature drops
    T_out(i) = T0 + delta_T_diode1 + delta_T_diode2;

    % Step 6: Store the current temperature difference for the next time step
    delta_T_diode1_prev = delta_T_diode1; % Update previous step temperature difference for diode 1
    delta_T_diode2_prev = delta_T_diode2; % Update previous step temperature difference for diode 2
end

% Plot results
figure;
plot(time, TR, 'b-', 'LineWidth', 2); % Input temperature wave
hold on;
plot(time, T_out, 'r-', 'LineWidth', 2); % Output temperature (Square wave)
yline(T0 + T_breakdown, 'k--', 'LineWidth', 1.5); % Upper boundary of Zener diode 1
yline(T0 - T_breakdown, 'k--', 'LineWidth', 1.5); % Lower boundary of Zener diode 2
hold off;
grid on;
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('Input Temperature Wave', 'Output Temperature (Square Wave)', 'Upper Boundary', 'Lower Boundary');
title('Temperature Clipping Circuit with Double Thermal Zener Diodes');
