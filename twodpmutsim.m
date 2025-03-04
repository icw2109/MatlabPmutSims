clear;


c0 = 1500; rho0 = 1000; c_skull = 3000; rho_skull = 1900;
skull_thickness = 5e-3; alpha_coeff = 0.002; alpha_power = 2;
num_elements_x = 20; num_elements_y = 20;
element_width = 280e-6;
element_spacing = 400e-6; 
frequency = 650e3; focal_distance = 50e-3;

dx = 180e-6;
Nx = 500;    
Ny = 500;   
dt = 40e-9;  
kgrid = kWaveGrid(Nx, dx, Ny, dx);

% Correct medium fields
medium.sound_speed = c0 * ones(Nx, Ny);
medium.density = rho0 * ones(Nx, Ny);
medium.alpha_coeff = alpha_coeff;
medium.alpha_power = alpha_power;


skull_end = Ny - round(20e-3 / dx);
skull_start = skull_end - round(skull_thickness / dx);
medium.sound_speed(:, skull_start:skull_end) = c_skull;
medium.density(:, skull_start:skull_end) = rho_skull;

Nt = 3000;
kgrid.setTime(Nt, dt);

array_width_x = (num_elements_x - 1) * element_spacing + element_width;
array_width_y = (num_elements_y - 1) * element_spacing + element_width;
element_size = max(1, round(element_width / dx));
gap_size = max(0, round((element_spacing - element_width) / dx));

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

[element_pos_x, element_pos_y] = meshgrid(((1:num_elements_x) - (num_elements_x+1)/2) * element_spacing, ...
                                          ((1:num_elements_y) - (num_elements_y+1)/2) * element_spacing);

curvature_effect = 0.1; 
effective_focal_distance = focal_distance - curvature_effect;

delays = (sqrt((Ny*dx - effective_focal_distance).^2 + element_pos_x.^2 + element_pos_y.^2) ...
         - (Ny*dx - effective_focal_distance)) / c0;

apodization = tukeywin(num_elements_x, 0.3) * tukeywin(num_elements_y, 0.3)';

source_amp = 1e6;
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
sensor.record = {'p_final'};

input_args = {'PMLSize', 10, 'PlotPML', false, 'PlotSim', false, 'DataCast', 'single'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});


figure;
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, flipud(sensor_data.p_final));
colormap(hot);
colorbar;
ylabel('Lateral position (mm)');
xlabel('Axial position (mm)');
title('Pressure field at final time step');
hold on;
yline((Ny*dx - focal_distance)*1e3, 'r--', 'Target focal point');
legend('Pressure field', 'Target focal point');

% Pressure along central axis
central_axis = squeeze(sensor_data.p_final(round(Nx/2), :));
central_axis = flipud(central_axis);
y_vec = kgrid.y_vec;
y_vec = linspace(min(y_vec), max(y_vec), length(central_axis));

figure;
plot(y_vec*1e3, central_axis);
hold on;
xline((Ny*dx - focal_distance)*1e3, 'r--', 'Target focal point');
xlabel('Axial position (mm)');
ylabel('Pressure (Pa)');
title('Pressure along central axis');
legend('Pressure', 'Target focal point');


[max_pressure, max_index] = max(central_axis);
calculated_focal_point = (Ny*dx - y_vec(max_index))*1e3;
disp(['Calculated focal point location: ', num2str(calculated_focal_point), ' mm']);
disp(['Target focal point: ', num2str(focal_distance*1e3), ' mm']);
focal_point_difference = abs(calculated_focal_point - focal_distance*1e3);
disp(['Difference between calculated and target focal points: ', num2str(focal_point_difference), ' mm']);
