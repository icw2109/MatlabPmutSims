% Focused Ultrasound through Skull with Heating and Cavitation Simulation
clear;

% safety sim for pmut, still needs checking through, work, etc.
% Parameters
c0 = 1500;                 % sound speed in tissue (m/s)
rho0 = 1000;               % density in tissue (kg/m^3)
c_skull = 3000;            % sound speed in skull (m/s)
rho_skull = 1900;          % density in skull (kg/m^3)
skull_thickness = 5e-3;    % skull thickness (m)
alpha_coeff = 0.002;       % absorption coefficient (dB/MHz^2)
alpha_power = 2;           % absorption power law exponent
frequency = 650e3;         % ultrasound frequency (Hz)
focal_distance = 50e-3;    % focal distance (m)
num_elements_x = 20;       % number of elements in x direction
num_elements_y = 20;       % number of elements in y direction
element_width = 280e-6;    % transducer element width (m)
element_spacing = 400e-6;  % transducer element spacing (m)
curvature_effect = 0.1;    % effective focal distance correction
source_amp = 1e6;          % source amplitude (Pa)

% Grid parameters
dx = 180e-6;               % grid spacing (m)
Nx = 500;                  % grid points in x direction
Ny = 500;                  % grid points in y direction
dt = 40e-9;                % time step (s)
kgrid = kWaveGrid(Nx, dx, Ny, dx);


medium.sound_speed = c0 * ones(Nx, Ny);    % sound speed in tissue
medium.density = rho0 * ones(Nx, Ny);      % density in tissue
medium.alpha_coeff = alpha_coeff;          % absorption coefficient
medium.alpha_power = alpha_power;          % absorption power law
medium.BonA = 6;                           % nonlinearity parameter

% Add skull layer properties
skull_end = Ny - round(20e-3 / dx);        % end of skull
skull_start = skull_end - round(skull_thickness / dx); % start of skull
medium.sound_speed(:, skull_start:skull_end) = c_skull; % skull sound speed
medium.density(:, skull_start:skull_end) = rho_skull;   % skull density

Nt = 3000;    % number of time steps
kgrid.setTime(Nt, dt);

% Create transducer array
array_width_x = (num_elements_x - 1) * element_spacing + element_width;
array_width_y = (num_elements_y - 1) * element_spacing + element_width;
element_size = max(1, round(element_width / dx));
gap_size = max(0, round((element_spacing - element_width) / dx));

% Initialize source mask for the transducer
source.p_mask = zeros(Nx, Ny);
start_index_x = round((Nx - array_width_x / dx) / 2);
start_index_y = Ny - 20 - num_elements_y * (element_size + gap_size);

for i = 1:num_elements_x
    for j = 1:num_elements_y
        idx_x = start_index_x + (i-1) * (element_size + gap_size);
        idx_y = start_index_y + (j-1) * (element_size + gap_size);
        source.p_mask(idx_x:idx_x+element_size-1, idx_y:idx_y+element_size-1) = 1;
    end
end

% Calculate delays and apodization
[element_pos_x, element_pos_y] = meshgrid(((1:num_elements_x) - (num_elements_x+1)/2) * element_spacing, ...
                                          ((1:num_elements_y) - (num_elements_y+1)/2) * element_spacing);

effective_focal_distance = focal_distance - curvature_effect;
delays = (sqrt((Ny*dx - effective_focal_distance).^2 + element_pos_x.^2 + element_pos_y.^2) ...
         - (Ny*dx - effective_focal_distance)) / c0;
apodization = tukeywin(num_elements_x, 0.3) * tukeywin(num_elements_y, 0.3)';

% Create source signals for each element
source.p = zeros(sum(source.p_mask(:)), length(kgrid.t_array));
element_index = 1;
for i = 1:num_elements_x
    for j = 1:num_elements_y
        idx_x = start_index_x + (i-1) * (element_size + gap_size);
        idx_y = start_index_y + (j-1) * (element_size + gap_size);
        num_element_points = sum(sum(source.p_mask(idx_x:idx_x+element_size-1, idx_y:idx_y+element_size-1)));
        element_signals = source_amp * apodization(i,j) * sin(2 * pi * frequency * (kgrid.t_array - delays(i, j)));
        source.p(element_index:element_index+num_element_points-1, :) = repmat(element_signals, num_element_points, 1);
        element_index = element_index + num_element_points;
    end
end

sensor.mask = ones(Nx, Ny);
sensor.record = {'p', 'p_final'};

% Simulation parameters
num_pulses = 5;
pulse_interval = 1e-3;  % 1 ms between pulses
temperature_rise = zeros(Nx, Ny);
cavitation_index = zeros(Nx, Ny);
input_args = {'PMLSize', 10, 'PlotPML', false, 'PlotSim', false, 'DataCast', 'single'};

% Run simulation for multiple pulses
for pulse = 1:num_pulses
    % Run k-Wave simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    % Calculate temperature rise (simplified model)
    specific_heat = 3600;  % J/(kg·K) for tissue
    absorption = medium.alpha_coeff * (frequency*1e-6)^medium.alpha_power;  
    intensity = mean(sensor_data.p.^2, 2) / (medium.density(1) * medium.sound_speed(1));
    delta_T = absorption * intensity * kgrid.dt * Nt / (medium.density(1) * specific_heat);
    temperature_rise = temperature_rise + reshape(delta_T, Nx, Ny);
    
    % Calculate cavitation index 
    p_neg = min(sensor_data.p, [], 2);
    cavitation_index = cavitation_index + reshape(abs(p_neg) / 1e6, Nx, Ny);
    
    % Pause between pulses (simulated delay)
    pause(0.01);  % Reduced for faster simulation
end

% Visualization of results
figure;

% Pressure field visualization
subplot(2,2,1);
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, flipud(sensor_data.p_final));
colormap(hot);
colorbar;
ylabel('Lateral position (mm)');
xlabel('Axial position (mm)');
title('Pressure Field at Final Time Step');
hold on;
yline((Ny*dx - focal_distance)*1e3, 'r--', 'Target Focal Point');

% Central axis pressure visualization
subplot(2,2,2);
central_axis = squeeze(sensor_data.p_final(round(Nx/2), :));
central_axis = flipud(central_axis);
y_vec = kgrid.y_vec;
y_vec = linspace(min(y_vec), max(y_vec), length(central_axis));
plot(y_vec*1e3, central_axis);
hold on;
xline((Ny*dx - focal_distance)*1e3, 'r--', 'Target Focal Point');
xlabel('Axial position (mm)');
ylabel('Pressure (Pa)');
title('Pressure Along Central Axis');

% Temperature rise visualization
subplot(2,2,3);
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, flipud(temperature_rise));
colormap(hot);
colorbar;
ylabel('Lateral position (mm)');
xlabel('Axial position (mm)');
title('Cumulative Temperature Rise (K)');

% Cavitation index visualization
subplot(2,2,4);
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, flipud(cavitation_index));
colormap(hot);
colorbar;
ylabel('Lateral position (mm)');
xlabel('Axial position (mm)');
title('Cumulative Cavitation Index');

% Display calculated focal point information
[max_pressure, max_index] = max(central_axis);
calculated_focal_point = (Ny*dx - y_vec(max_index))*1e3;
disp(['Calculated focal point location: ', num2str(calculated_focal_point), ' mm']);
disp(['Target focal point: ', num2str(focal_distance*1e3), ' mm']);
focal_point_difference = abs(calculated_focal_point - focal_distance*1e3);
disp(['Difference between calculated and target focal points: ', num2str(focal_point_difference), ' mm']);
