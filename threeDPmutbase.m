clear;

% No skull or brain; however, this is the 3d model. 

% This is one of the models I used to estimate performance of
% the PMUT for effects on brain. This is a baseline 3d model without
% MRI data. This model although only 150 lines of code was extremely
% difficult to construct as it is simulating many PMUT elements with
% modified performance constraints as opposed to one singular PMUT element
% Happy to talk though this code as well as challanges in writing in,
% and setting it up for performance. I wrote this entirely from scratch as
% well. 

% Grid parameters
Nx = 160; Ny = 128; Nz = 160;  % Increase grid size to accommodate the array
dx = 1e-3; dy = 1e-3; dz = 1e-3;  % Grid spacing in meters
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

c_medium = 1540;  % [m/s] Sound speed (e.g., soft tissue or water)
rho_medium = 1000;  % [kg/m^3] Density (e.g., soft tissue or water)

medium.sound_speed = c_medium * ones(Nx, Ny, Nz);  % Homogeneous medium
medium.density = rho_medium * ones(Nx, Ny, Nz);  % Homogeneous medium

% Time array
kgrid.t_array = makeTime(kgrid, max(medium.sound_speed(:)));

% Source parameters
source_f0 = 1e6;  % Source frequency [Hz]
source_amp = 1e6; % Source amplitude [Pa]
source_cycles = 3;

% Create the transducer array
array_size = 20;  
element_spacing = 2;  % Element spacing in grid points

% Define the diameter of each element in grid points (0.18 mm -> 0.18/1e-3 = 0.18)
element_diameter = 180e-6;  % Diameter of 180 µm (0.18 mm)
element_radius = element_diameter / 2 / dx;  % Radius in grid points

% Calculate the starting positions for the array
array_pos_x = Nx/2 - (array_size-1)*element_spacing/2;  % Ensure array is centered in x
array_pos_z = Nz/2 - (array_size-1)*element_spacing/2;  % Ensure array is centered in z
array_pos = [array_pos_x, 1, array_pos_z];

% Create array mask with circular elements
source.p_mask = zeros(Nx, Ny, Nz);
for ix = 1:array_size
    for iz = 1:array_size
        % Calculate element center positions
        x_center = array_pos(1) + (ix-1)*element_spacing;
        z_center = array_pos(3) + (iz-1)*element_spacing;
        
        % Define circular element in x-z plane (y = 1)
        for x = floor(x_center-element_radius):ceil(x_center+element_radius)
            for z = floor(z_center-element_radius):ceil(z_center+element_radius)
                % Check if within the radius of the element
                if (x - x_center)^2 + (z - z_center)^2 <= element_radius^2
                    if x > 0 && x <= Nx && z > 0 && z <= Nz
                        source.p_mask(x, 1, z) = 1;
                    end
                end
            end
        end
    end
end

% Calculate element positions
[elem_x, elem_y, elem_z] = ind2sub(size(source.p_mask), find(source.p_mask));
elem_pos = [elem_x, elem_y, elem_z] .* [dx, dy, dz];

% Define the focus point 50 mm away from the array (along the z-axis)
focus_point_x = Nx/2 * dx;  % Centered in x
focus_point_y = Ny/2 * dy;  % Centered in y
focus_point_z = 50e-3;      % 50 mm away in the z-axis

focus_point = [focus_point_x, focus_point_y, focus_point_z];

% Calculate delays to focus the beam at the focus_point (50 mm away)
distances_to_focus = sqrt(sum((elem_pos - focus_point).^2, 2));
max_distance_focus = max(distances_to_focus);

% Get sound speed for each element's position (assume 2D slice at y = 1)
element_sound_speed = medium.sound_speed(sub2ind(size(medium.sound_speed), elem_x, elem_y, elem_z));

% Perform element-wise division to calculate delays
focusing_delays = (max_distance_focus - distances_to_focus) ./ element_sound_speed;

% Create the source signal
source_signal = source_amp * toneBurst(1/kgrid.dt, source_f0, source_cycles);

% Apply the delays for focusing at the center
num_elements = sum(source.p_mask(:));
source.p = zeros(num_elements, length(kgrid.t_array));

for i = 1:num_elements
    % Calculate the shift in samples for focusing
    shift_samples_focus = round(focusing_delays(i) / kgrid.dt);
    
    % Create a time array for this element, incorporating the focusing delay
    t_element_focus = (0:length(kgrid.t_array)-1) * kgrid.dt - focusing_delays(i);
    
    % Use interpolation to create the delayed signal
    source.p(i, :) = interp1((-shift_samples_focus:length(source_signal)-1)*kgrid.dt, ...
                             [zeros(1, shift_samples_focus), source_signal], ...
                             t_element_focus, ...
                             'linear', 0);
end

% Define the sensor mask
sensor.mask = zeros(Nx, Ny, Nz);

% Define the sensor mask at the center slice of the grid
sensor_x_start = round(Nx/4);
sensor_x_end = round(3*Nx/4);
sensor_z_start = round(Nz/4);
sensor_z_end = round(3*Nz/4);
sensor_y = round(Ny/2);  

sensor.mask(sensor_x_start:sensor_x_end, sensor_y, sensor_z_start:sensor_z_end) = 1;

sensor.record = {'p'};  % Record the pressure at the sensor locations

input_args = {'PMLInside', false, 'PlotPML', false};

sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

p_max_full = zeros(Nx, Ny, Nz);  % Full grid of zeros

% Get the indices of the sensor mask where recording took place
sensor_indices = find(sensor.mask);

% Extract the recorded pressure data from sensor_data.p
if isfield(sensor_data, 'p')
    % Calculate the maximum pressure over time at each sensor point
    p_max_full(sensor_indices) = max(sensor_data.p, [], 2);  % Take the max over time for each sensor point
    
    % Visualize results for the focused beam
    figure;
    imagesc(squeeze(p_max_full(:, Ny/2, :)));  % Plot the x-z plane at y = Ny/2
    title('Maximum Pressure Focused at 50mm (x-z plane at y = Ny/2)');
    xlabel('z'); ylabel('x');
    colorbar;
    axis image;

    % Display maximum pressure
    max_pressure = max(p_max_full(:));
    fprintf('Maximum pressure: %.2f MPa\n', max_pressure/1e6);
else
    error('The sensor_data structure does not contain the expected field "p".');
end

