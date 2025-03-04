import gdspy
import time
import os

# This file generates the GDS file for our custome PMUT chip, each of the 15 layers must comply with design rules given by the foundry, as well as 
# simulated chip performance based on performance specifications and material properties done using finite element analysis. Each layer cannot conflict with any other 
# layer or the chip will fail. 
# To OPEN GDS file you need K layout or related software, K layout is free if interested.
# I have also included screenshots of what this program produces when ran. Happy to talk in detail
# about design choices, constraints and how program is structured.
# I created this program myself from scratch.  


# Define the layers for metal connections and additional rectangles
layers = {
    "ZR": (1, 0),  # Zero Alignment
    "BS": (2, 0),  # Backside etch
    "BE": (3, 0),  # Bottom Electrode (platinum)
    "TV": (13, 0),  # Top Electrode via
    "FE": (5, 0),  # Ferroelectric (PZT etch)
    "ILD": (8, 0),  # ILD Via
    "M1": (9, 0),  # Metal One (Cr/Au wire bond pad)
    "TR": (11, 0),  # Trench for dicing lanes
    "TE": (7, 0),  # Top Electrode (platinum)
    "BEVIA": (19, 0),  # Vias to BE through PZT
    "M1_CU": (17, 0),  # Metal One Copper (Cr/Cu/Au solder pad)
    "IB": (30, 0),  # Ibeam1 mechanical support
    "IB2": (31, 0),  # Ibeam2 mechanical support
    "FLD": (52, 0),  # ILD Field
    "DICE": (54, 0)  # Dicing Lane
}

# Define the PMUT element configuration
class PMUTArrayConfig:
    def __init__(
        self, spacing_x, spacing_y, num_pmuts, pmut_diameter, connection_width,
        pad_size, pad_spacing, pad_x_offset, pad_y_offset, extra_pad_y_offset_bottom, extra_pad_y_offset_top,
        connection_x_offset, connection_y_offset, bottom_connection_position, top_connection_position,
        bottom_direction, top_direction
    ):
        self.spacing_x = spacing_x
        self.spacing_y = spacing_y
        self.num_pmuts = num_pmuts
        self.pmut_diameter = pmut_diameter  # Diameter of PMUT elements
        self.connection_width = connection_width  # Connection width for solder pads
        self.pad_size = pad_size  # Size of the solder pads
        self.pad_spacing = pad_spacing  # Spacing between solder pads
        self.pad_x_offset = pad_x_offset  # X-offset for solder pad placement
        self.pad_y_offset = pad_y_offset  # Y-offset for solder pad placement
        self.extra_pad_y_offset_bottom = extra_pad_y_offset_bottom  # Y-offset for extra bottom rectangles
        self.extra_pad_y_offset_top = extra_pad_y_offset_top  # Y-offset for extra top rectangles
        self.connection_x_offset = connection_x_offset  # X-offset for connection placement
        self.connection_y_offset = connection_y_offset  # Y-offset for connection placement
        self.bottom_connection_position = bottom_connection_position  # (x, y) starting position for bottom connections
        self.top_connection_position = top_connection_position  # (x, y) starting position for top connections
        self.bottom_direction = bottom_direction  # Direction for bottom connections ('+x', '-x', '+y', '-y')
        self.top_direction = top_direction  # Direction for top connections ('+x', '-x', '+y', '-y')

# Base configuration for PMUT array
base_spacing_x = 390  # Pitch in X direction (center-to-center distance)
base_spacing_y = 390  # Pitch in Y direction (center-to-center distance)
base_num_pmuts = 20  # Number of PMUTs per group (should be 20x20)
base_pmut_diameter = 280  # Diameter of each PMUT element
base_connection_width = 8  # Width of the connections to solder pads
base_pad_size = (760, 700)  # Size of the solder pads
base_pad_spacing = 2400  # Spacing between solder pads
base_pad_x_offset = 100  # X-offset for solder pad placement
base_pad_y_offset = 100  # Y-offset for solder pad placement
base_extra_pad_y_offset_bottom = 900  # Y-offset for extra bottom rectangles
base_extra_pad_y_offset_top = 900  # Y-offset for extra top rectangles
base_connection_x_offset = 10000  # X-offset for connection placement
base_connection_y_offset = 10000  # Y-offset for connection placement
base_bottom_connection_position = (1000, 1000)  # Starting (x, y) for bottom connections
base_top_connection_position = (2000, 1000)  # Starting (x, y) for top connections
bottom_direction = '-y'  # Direction for bottom connections
top_direction = '+y'  # Direction for top connections

# Initialize GDS library and top cell
gdsii_lib = gdspy.GdsLibrary()
cell_name = f'METAL_CONNECTIONS_{int(time.time())}'
top_cell = gdsii_lib.new_cell(cell_name)

def add_rectangles_with_xy_control(group_positions, intra_group_spacing, config, scaling_factor, top_cell):
    """Adds sets of rectangles to the PMUT array with XY control for each group."""
    pad_width, pad_height = config.pad_size

    for group, (x_start, y_start) in enumerate(group_positions):
        group_x_start = x_start * scaling_factor
        group_y_start = y_start * scaling_factor

        # Use the M1 layer for metallization
        layer = layers["M1"]

        for col in range(3):  # Modified to create only 3 columns
            x_pos = group_x_start + col * intra_group_spacing * scaling_factor + config.pad_x_offset * scaling_factor
            y_pos = group_y_start + config.pad_y_offset * scaling_factor

            # Create rectangles on the M1 layer
            rectangle = gdspy.Rectangle(
                (x_pos - pad_width / 2 * scaling_factor, y_pos - pad_height / 2 * scaling_factor),
                (x_pos + pad_width / 2 * scaling_factor, y_pos + pad_height / 2 * scaling_factor),
                layer=layer[0], datatype=layer[1]
            )
            top_cell.add(rectangle)

def add_lines_with_horizontal_top(x_start, y_start, line_count, line_length, line_spacing, line_width, line_layer, top_cell, scaling_factor, horizontal_length=80, vertical_offset=10):
    """Add a set of lines with specified parameters, including a horizontal line slightly below the top."""
    for i in range(line_count):
        x_pos = x_start * scaling_factor + i * line_spacing * scaling_factor
        
        # Create the vertical line, reducing its length by the vertical_offset
        line = gdspy.Path(line_width * scaling_factor, (x_pos, y_start * scaling_factor))
        line.segment((line_length - vertical_offset) * scaling_factor, direction='+y', layer=line_layer[0], datatype=line_layer[1])
        
        # Add a horizontal line starting 10 µm below the top of the original line length
        line.segment(horizontal_length * scaling_factor, direction='+x', layer=line_layer[0], datatype=line_layer[1])
        
        top_cell.add(line)

# Replace existing calls with the updated function
# Add two sets of lines, with a horizontal line 10 µm below the top of the vertical line
add_lines_with_horizontal_top(37830, 37820, 20, 7900, 400, 30, layers["BE"], top_cell, 1)
add_lines_with_horizontal_top(38160, 37980, 20, 7900, 400, 20, layers["TE"], top_cell, 1)

def create_backside_circles(x_start, y_start, num_pmuts, spacing_x, spacing_y, pmut_diameter, top_cell, scaling_factor):
    """Create backside etch circles for each PMUT element."""
    # Create a 20x20 array of circles
    for row in range(num_pmuts):
        for col in range(num_pmuts):
            x_pos = x_start + col * spacing_x * scaling_factor
            y_pos = y_start + row * spacing_y * scaling_factor
            backside_circle = gdspy.Round(
                (x_pos, y_pos), 
                pmut_diameter / 2 * scaling_factor, 
                number_of_points=64, 
                layer=layers["BS"][0], 
                datatype=layers["BS"][1]
            )
            top_cell.add(backside_circle)


def add_electrode_rectangles(via_positions, top_electrode_size, bottom_electrode_size, scaling_factor, top_cell, 
                             horizontal_offset=0, gap=10, x_start_offset=0, y_start_offset=0, 
                             exclude_bottom_two_rows=False, scale_factor=1.14, is_top_electrode=False):
    """
    Add top or bottom electrode rectangles centered around the vias with adjustable horizontal offset and gap.
    The function determines whether to place top or bottom electrodes based on the `is_top_electrode` flag.
    """
    
    electrode_size = top_electrode_size if is_top_electrode else bottom_electrode_size
    width, height = electrode_size
    
    for x_pos, y_pos in via_positions:
        adjusted_x_pos = x_pos + horizontal_offset * scaling_factor + x_start_offset
        adjusted_y_pos = y_pos + y_start_offset
        
        if is_top_electrode:
            # Top electrode rectangle
            x_start_top = adjusted_x_pos - (width / 2) * scaling_factor
            y_start_top = adjusted_y_pos - (height / 2) * scaling_factor
            electrode_rect = gdspy.Rectangle(
                (x_start_top, y_start_top),
                (x_start_top + width * scaling_factor - gap * scaling_factor, 
                 y_start_top + height * scaling_factor),
                layer=layers["TE"][0], datatype=layers["TE"][1]
            )
        else:
            # Bottom electrode rectangle
            x_start_bottom = adjusted_x_pos - width * scaling_factor * scale_factor / 2
            y_start_bottom = adjusted_y_pos - (height / 2) * scaling_factor
            electrode_rect = gdspy.Rectangle(
                (x_start_bottom, y_start_bottom),
                (x_start_bottom + width * scaling_factor * scale_factor, 
                 y_start_bottom + height * scaling_factor),
                layer=layers["BE"][0], datatype=layers["BE"][1]
            )
        
        top_cell.add(electrode_rect)


def add_vias(x_start, y_start, num_vias, via_size, via_spacing, row_spacing, top_cell, scaling_factor, layer):
    """Add vias with specified parameters, centered for each electrode."""
    via_positions = []
    for row in range(2):  # Two rows of vias
        for i in range(num_vias):
            x_pos = x_start * scaling_factor + (i + 0.5) * via_spacing * scaling_factor - (via_size / 2) * scaling_factor
            y_pos = y_start * scaling_factor + (row + 0.5) * row_spacing * scaling_factor - (via_size / 2) * scaling_factor
            via = gdspy.Rectangle(
                (x_pos, y_pos),
                (x_pos + via_size * scaling_factor, y_pos + via_size * scaling_factor),
                layer=layer[0], datatype=layer[1]
            )
            top_cell.add(via)
            via_positions.append((x_pos + via_size / 2 * scaling_factor, y_pos + via_size / 2 * scaling_factor))
    return via_positions

def connect_vias_to_pads(via_positions, connection_params, line_layer, top_cell, scaling_factor):
    """Connect each via to solder pads with vertical lines using the provided connection parameters."""
    y_offset = 80  # Raise Y start by 80 μm
    for (x_pos, y_pos), (line_width, line_length, direction, offset) in zip(via_positions, connection_params):
        if direction == '-y':
            start_y = y_pos + offset * scaling_factor + 20 * scaling_factor + y_offset  # Start 20 um above the bottom via
        else:
            start_y = y_pos - offset * scaling_factor - y_offset  # Start 10 um below the top via

        connection_line = gdspy.Path(line_width * scaling_factor, (x_pos, start_y))
        connection_line.segment(line_length * scaling_factor, direction=direction, layer=line_layer[0], datatype=line_layer[1])
        top_cell.add(connection_line)

def add_connection_squares(via_positions, square_dimensions, layer, top_cell, scaling_factor):
    """Add customizable connection squares connecting the vias."""
    square_width, square_height = square_dimensions
    for x_pos, y_pos in via_positions:
        square = gdspy.Rectangle(
            (x_pos - square_width / 2 * scaling_factor, y_pos - square_height / 2 * scaling_factor),
            (x_pos + square_width / 2 * scaling_factor, y_pos + square_height / 2 * scaling_factor),
            layer=layer[0], datatype=layer[1]
        )
        top_cell.add(square)

def add_electrodes(x_start, y_start, num_pmuts, spacing_x, spacing_y, top_electrode_radius, bottom_electrode_radius, top_cell, scaling_factor):
    """Add top and bottom electrodes to each PMUT element with sufficient overlap."""
    for row in range(num_pmuts):
        for col in range(num_pmuts):
            x_pos = x_start + col * spacing_x * scaling_factor
            y_pos = y_start + row * spacing_y * scaling_factor
            
            # Top electrode
            top_electrode = gdspy.Round(
                (x_pos, y_pos), 
                top_electrode_radius * scaling_factor, 
                number_of_points=64, 
                layer=layers["TE"][0], 
                datatype=layers["TE"][1]
            )
            top_cell.add(top_electrode)
            
            # Bottom electrode with specified larger radius
            bottom_electrode = gdspy.Round(
                (x_pos, y_pos), 
                bottom_electrode_radius * scaling_factor, 
                number_of_points=64, 
                layer=layers["BE"][0], 
                datatype=layers["BE"][1]
            )
            top_cell.add(bottom_electrode)

def add_ibeam_ring(x_start, y_start, num_pmuts, spacing_x, spacing_y, ibeam_ring_inner_radius, ibeam_ring_outer_radius, top_cell, scaling_factor):
    """Add IBEAM ring for mechanical support around each PMUT element, ensuring separation from the metal layers."""
    for row in range(num_pmuts):
        for col in range(num_pmuts):
            x_pos = x_start + col * spacing_x * scaling_factor
            y_pos = y_start + row * spacing_y * scaling_factor
            ibeam_ring = gdspy.Round(
                (x_pos, y_pos), 
                ibeam_ring_outer_radius * scaling_factor, 
                inner_radius=ibeam_ring_inner_radius * scaling_factor, 
                number_of_points=64, 
                layer=layers["IB"][0], 
                datatype=layers["IB"][1]
            )
            top_cell.add(ibeam_ring)

def add_electrode_connection_lines(
    x_start, y_start, num_pmuts, spacing_x, spacing_y, 
    top_electrode_radius, bottom_electrode_radius, 
    bottom_connection_length, top_connection_length, 
    top_cell, scaling_factor):
    """Add connection lines for electrodes with specified widths and adjustable lengths."""
    for row in range(num_pmuts):
        for col in range(num_pmuts):
            x_pos = x_start + col * spacing_x * scaling_factor
            y_pos = y_start + row * spacing_y * scaling_factor
            
            # Bottom left electrode connection line (30 µm wide) combined with bottom electrode
            bottom_left_connection = gdspy.Path(
                30 * scaling_factor, (x_pos - 0.5 * bottom_electrode_radius * scaling_factor, y_pos)
            )
            bottom_left_connection.segment(
                bottom_connection_length * scaling_factor, direction='-x', layer=layers["BE"][0], datatype=layers["BE"][1]
            )
            top_cell.add(bottom_left_connection)

            # Top right electrode connection line (20 µm wide) combined with top electrode
            top_right_connection = gdspy.Path(
                20 * scaling_factor, (x_pos + 0.5 * top_electrode_radius * scaling_factor, y_pos)
            )
            top_right_connection.segment(
                top_connection_length * scaling_factor, direction='+x', layer=layers["TE"][0], datatype=layers["TE"][1]
            )
            top_cell.add(top_right_connection)

def add_fe_square(x_start, y_start, width, height, top_cell, scaling_factor):
    """Add a ferroelectric (FE) square with specified dimensions."""
    fe_square = gdspy.Rectangle(
        (x_start * scaling_factor, y_start * scaling_factor),
        ((x_start + width) * scaling_factor, (y_start + height) * scaling_factor),
        layer=layers["FE"][0], datatype=layers["FE"][1]
    )
    top_cell.add(fe_square)

def add_dicing_lanes(x_start, y_start, lane_width, lane_height, trench_offset, top_cell, scaling_factor):
    """Add dicing lanes using TRENCH and BACKSIDE layers with hollow centers."""
    # Define the outer rectangle for the dicing lane
    outer_rect = gdspy.Rectangle(
        (x_start * scaling_factor, y_start * scaling_factor),
        ((x_start + lane_width) * scaling_factor, (y_start + lane_height) * scaling_factor),
        layer=layers["TR"][0], datatype=layers["TR"][1]
    )

    # Define the inner rectangle for the hollow center
    inner_rect = gdspy.Rectangle(
        ((x_start + trench_offset) * scaling_factor, (y_start + trench_offset) * scaling_factor),
        ((x_start + lane_width - trench_offset) * scaling_factor, (y_start + lane_height - trench_offset) * scaling_factor)
    )

    # Create the hollow dicing lane by subtracting the inner rectangle from the outer rectangle
    hollow_lane = gdspy.boolean(outer_rect, inner_rect, 'not', layer=layers["TR"][0], datatype=layers["TR"][1])
    top_cell.add(hollow_lane)

    # Create an identical structure for the backside layer
    backside_hollow_lane = gdspy.boolean(outer_rect, inner_rect, 'not', layer=layers["BS"][0], datatype=layers["BS"][1])
    top_cell.add(backside_hollow_lane)

def add_dicing_lanes_with_fld(x_start, y_start, lane_width, lane_height, trench_offset, top_cell, scaling_factor):
    """Add dicing lanes using BACKSIDE and TRENCH layers with hollow centers and overlay with FLD."""
    # Define the outer rectangle for the dicing lane (TR layer)
    outer_rect_tr = gdspy.Rectangle(
        (x_start * scaling_factor, y_start * scaling_factor),
        ((x_start + lane_width) * scaling_factor, (y_start + lane_height) * scaling_factor),
        layer=layers["TR"][0], datatype=layers["TR"][1]
    )

    # Define the inner rectangle for the hollow center (common for both TR and FLD layers)
    inner_rect = gdspy.Rectangle(
        ((x_start + trench_offset) * scaling_factor, (y_start + trench_offset) * scaling_factor),
        ((x_start + lane_width - trench_offset) * scaling_factor, (y_start + lane_height - trench_offset) * scaling_factor)
    )

    # Create the hollow dicing lane for the TR layer
    hollow_lane_tr = gdspy.boolean(outer_rect_tr, inner_rect, 'not', layer=layers["TR"][0], datatype=layers["TR"][1])
    top_cell.add(hollow_lane_tr)

    # Create the hollow dicing lane for the FLD layer
    outer_rect_fld = gdspy.Rectangle(
        (x_start * scaling_factor, y_start * scaling_factor),
        ((x_start + lane_width) * scaling_factor, (y_start + lane_height) * scaling_factor),
        layer=layers["FLD"][0], datatype=layers["FLD"][1]
    )

    hollow_lane_fld = gdspy.boolean(outer_rect_fld, inner_rect, 'not', layer=layers["FLD"][0], datatype=layers["FLD"][1])
    top_cell.add(hollow_lane_fld)

def add_extension_lines(xy_positions, extension_length, line_width, line_layer, top_cell, scaling_factor, direction='+y'):
    """Add extension lines to specific bond pads with customizable length and direction."""
    for x_pos, y_pos in xy_positions:
        if direction in ('+y', '-y'):
            extension_line = gdspy.Path(line_width * scaling_factor, (x_pos * scaling_factor, y_pos * scaling_factor))
            extension_line.segment(extension_length * scaling_factor, direction=direction, layer=line_layer[0], datatype=line_layer[1])
        elif direction in ('+x', '-x'):
            extension_line = gdspy.Path(line_width * scaling_factor, (x_pos * scaling_factor, y_pos * scaling_factor))
            extension_line.segment(extension_length * scaling_factor, direction=direction, layer=line_layer[0], datatype=line_layer[1])
        top_cell.add(extension_line)

def add_bond_pad(x_pos, y_pos, pad_size, pad_layer, top_cell, scaling_factor):
    """Add a single bond pad at a specified position."""
    pad_width, pad_height = pad_size
    bond_pad = gdspy.Rectangle(
        (x_pos - pad_width / 2 * scaling_factor, y_pos - pad_height / 2 * scaling_factor),
        (x_pos + pad_width / 2 * scaling_factor, y_pos + pad_height / 2 * scaling_factor),
        layer=pad_layer[0], datatype=pad_layer[1]
    )
    top_cell.add(bond_pad)

# Define the base position for a single array at the center of the wafer
wafer_center = (38000, 38000)  # Center of a 100mm wafer in micrometers
scaling_factor = 1  # Set scaling factor to 1 for actual size

# Define XY positions for each group of rectangles
group_positions = [
    (38250, 37000),  # Position for the first group   
    (38920, 35000),  # Position for the second group
    (39490, 33000),  # Position for the third group
    (38280, 46500),  # Position for the fourth group
    (38950, 48500),  # Position for the fifth group
    (39620, 50500)   # Position for the sixth group
]

# Spacing control within each group
intra_group_spacing = 2500  # Increased spacing in micrometers between rectangles within a group

# Create the PMUT array configuration
pmut_config = PMUTArrayConfig(
    base_spacing_x, base_spacing_y, base_num_pmuts, base_pmut_diameter,
    base_connection_width, base_pad_size, base_pad_spacing, base_pad_x_offset, base_pad_y_offset,
    base_extra_pad_y_offset_bottom, base_extra_pad_y_offset_top,
    base_connection_x_offset, base_connection_y_offset, base_bottom_connection_position, base_top_connection_position,
    bottom_direction, top_direction
)

# Add rectangles to the design with XY control
add_rectangles_with_xy_control(
    group_positions, intra_group_spacing, pmut_config, scaling_factor, top_cell
)

# Define parameters for the lines
line_count_20 = 20  # Number of lines

# Length of each line for different sets in micrometers
line_length_1 = 7900
line_length_2 = 7900

# Width of each line for different sets in micrometers
line_width_1 = 20
line_width_2 = 30

# Spacing between each line for different sets in micrometers
line_spacing_1 = 400
line_spacing_2 = 400

# Position for each set of lines
x_start_1 = 37830
y_start_1 = 37820

x_start_2 = 38160
y_start_2 = 37980

# Add two sets of lines, M1 as part of bottom electrode and M2 as part of top electrode
# add_lines(x_start_1, y_start_1, line_count_20, line_length_1, line_spacing_1, line_width_1, layers["BE"], top_cell, scaling_factor)
# add_lines(x_start_2, y_start_2, line_count_20, line_length_2, line_spacing_2, line_width_2, layers["TE"], top_cell, scaling_factor)

# Add backside circles
create_backside_circles(
    wafer_center[0], wafer_center[1],
    base_num_pmuts, base_spacing_x, base_spacing_y, base_pmut_diameter, top_cell, scaling_factor
)

# Define positions for the vias
# Adjust the start positions to match the desired locations
x_start_vias_bottom = 37890  # Starting x-position of the bottom electrode vias
y_start_vias_bottom = 37530  # Starting y-position of the bottom electrode vias

x_start_vias_top = 38090  # Starting x-position of the top electrode vias
y_start_vias_top = 45790  # Starting y-position of the top electrode vias

# Add vias for bottom electrode
via_size = 20  # Size of each via square
via_spacing_bottom = 800  # Adjusted spacing for two rows
row_spacing_bottom = 150  # Spacing between rows of vias for bottom
num_vias_bottom = 10  # Number of vias per row for bottom, each covering one electrode

via_positions_bottom = add_vias(
    x_start_vias_bottom, y_start_vias_bottom, num_vias_bottom, via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["BEVIA"]
)

# Add vias for top electrode
via_spacing_top = 800  # Adjusted spacing for two rows
row_spacing_top = 150  # Spacing between rows of vias for top
num_vias_top = 10  # Number of vias per row for top, each covering one electrode

via_positions_top = add_vias(
    x_start_vias_top, y_start_vias_top, num_vias_top, via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["TV"]
)

# Define connection parameters for bottom electrode to bond pads
connection_params_bottom = [
    (20, 500, '-y', 20),
    (20, 2500, '-y', 25),
    (20, 4600, '-y', 30),
    (20, 500, '-y', 15),
    (20, 2500, '-y', 10),
    (20, 4800, '-y', 5),
    (20, 500, '-y', 25),
    (20, 2500, '-y', 20),
    (20, 4600, '-y', 30),
    (20, 800, '-y', 15),
]

# Define connection parameters for top electrode to bond pads
connection_params_top = [
    (20, 500, '+y', 20),
    (20, 2800, '+y', 25),
    (20, 4550, '+y', 30),
    (20, 500, '+y', 15),
    (20, 2800, '+y', 10),
    (20, 4550, '+y', 5),
    (20, 500, '+y', 25),
    (20, 2800, '+y', 20),
    (20, 4550, '+y', 30),
    (20, 800, '+y', 15),
]

# Connect vias to the solder pads for the bottom electrode with customizable parameters
connect_vias_to_pads(via_positions_bottom, connection_params_bottom, line_layer=layers["M1"], top_cell=top_cell, scaling_factor=scaling_factor)

# Connect vias to the solder pads for the top electrode with customizable parameters
connect_vias_to_pads(via_positions_top, connection_params_top, line_layer=layers["M1"], top_cell=top_cell, scaling_factor=scaling_factor)

# Add ILD vias for bottom electrode
ild_via_size = 10  # Size of each ILD via square

add_vias(
    x_start_vias_bottom, y_start_vias_bottom, num_vias_bottom, ild_via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["ILD"]
)

# Add ILD vias for top electrode
add_vias(
    x_start_vias_top, y_start_vias_top, num_vias_top, ild_via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["ILD"]
)

# Add an additional set of vias 30 µm to the right of the current vias

# New positions for the additional vias
additional_x_offset = 30  # Additional offset for the new set of vias
x_start_vias_bottom_new = x_start_vias_bottom + additional_x_offset
x_start_vias_top_new = x_start_vias_top + additional_x_offset

# Add additional vias for bottom electrode
additional_via_positions_bottom = add_vias(
    x_start_vias_bottom_new, y_start_vias_bottom, num_vias_bottom, via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["BEVIA"]
)

# Add additional vias for top electrode
additional_via_positions_top = add_vias(
    x_start_vias_top_new, y_start_vias_top, num_vias_top, via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["TV"]
)

# Add ILD vias for the additional bottom electrode
add_vias(
    x_start_vias_bottom_new, y_start_vias_bottom, num_vias_bottom, ild_via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["ILD"]
)

# Add ILD vias for the additional top electrode
add_vias(
    x_start_vias_top_new, y_start_vias_top, num_vias_top, ild_via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["ILD"]
)

new_additional_x_offset = 30  # Another additional offset for the new set of vias
x_start_vias_bottom_new_2 = x_start_vias_bottom_new + new_additional_x_offset
x_start_vias_top_new_2 = x_start_vias_top_new + new_additional_x_offset

# Add new additional vias for bottom electrode
new_additional_via_positions_bottom = add_vias(
    x_start_vias_bottom_new_2, y_start_vias_bottom, num_vias_bottom, via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["BEVIA"]
)

# Add new additional vias for top electrode
new_additional_via_positions_top = add_vias(
    x_start_vias_top_new_2, y_start_vias_top, num_vias_top, via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["TV"]
)

# Add ILD vias for the new additional bottom electrode
add_vias(
    x_start_vias_bottom_new_2, y_start_vias_bottom, num_vias_bottom, ild_via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["ILD"]
)

# Add ILD vias for the new additional top electrode
add_vias(
    x_start_vias_top_new_2, y_start_vias_top, num_vias_top, ild_via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["ILD"]
)


#########################################################################

new_additional_x_offset = -60  # Additional offset for the new set of vias
new_additional_y_offset = 30  # Y offset for the new set of vias

# Calculate the new positions for the additional vias
x_start_vias_bottom_new_3 = x_start_vias_bottom_new_2 + new_additional_x_offset
y_start_vias_bottom_new_3 = y_start_vias_bottom + new_additional_y_offset

x_start_vias_top_new_3 = x_start_vias_top_new_2 + new_additional_x_offset
y_start_vias_top_new_3 = y_start_vias_top + new_additional_y_offset

# Add new additional vias for bottom electrode
new_additional_via_positions_bottom_3 = add_vias(
    x_start_vias_bottom_new_3, y_start_vias_bottom_new_3, num_vias_bottom, via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["BEVIA"]
)

# Add new additional vias for top electrode
new_additional_via_positions_top_3 = add_vias(
    x_start_vias_top_new_3, y_start_vias_top_new_3, num_vias_top, via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["TV"]
)

# Add ILD vias for the new additional bottom electrode
add_vias(
    x_start_vias_bottom_new_3, y_start_vias_bottom_new_3, num_vias_bottom, ild_via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["ILD"]
)

# Add ILD vias for the new additional top electrode
add_vias(
    x_start_vias_top_new_3, y_start_vias_top_new_3, num_vias_top, ild_via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["ILD"]
)

####################################################

new_additional_x_offset = -30  # Additional offset for the new set of vias
new_additional_y_offset = 30  # Y offset for the new set of vias

# Calculate the new positions for the additional vias
x_start_vias_bottom_new_3 = x_start_vias_bottom_new_2 + new_additional_x_offset
y_start_vias_bottom_new_3 = y_start_vias_bottom + new_additional_y_offset

x_start_vias_top_new_3 = x_start_vias_top_new_2 + new_additional_x_offset
y_start_vias_top_new_3 = y_start_vias_top + new_additional_y_offset

# Add new additional vias for bottom electrode
new_additional_via_positions_bottom_3 = add_vias(
    x_start_vias_bottom_new_3, y_start_vias_bottom_new_3, num_vias_bottom, via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["BEVIA"]
)

# Add new additional vias for top electrode
new_additional_via_positions_top_3 = add_vias(
    x_start_vias_top_new_3, y_start_vias_top_new_3, num_vias_top, via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["TV"]
)

# Add ILD vias for the new additional bottom electrode
add_vias(
    x_start_vias_bottom_new_3, y_start_vias_bottom_new_3, num_vias_bottom, ild_via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["ILD"]
)

# Add ILD vias for the new additional top electrode
add_vias(
    x_start_vias_top_new_3, y_start_vias_top_new_3, num_vias_top, ild_via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["ILD"]
)

######################################################

new_additional_x_offset = 0  # Additional offset for the new set of vias
new_additional_y_offset = 30  # Y offset for the new set of vias

# Calculate the new positions for the additional vias
x_start_vias_bottom_new_3 = x_start_vias_bottom_new_2 + new_additional_x_offset
y_start_vias_bottom_new_3 = y_start_vias_bottom + new_additional_y_offset

x_start_vias_top_new_3 = x_start_vias_top_new_2 + new_additional_x_offset
y_start_vias_top_new_3 = y_start_vias_top + new_additional_y_offset

# Add new additional vias for bottom electrode
new_additional_via_positions_bottom_3 = add_vias(
    x_start_vias_bottom_new_3, y_start_vias_bottom_new_3, num_vias_bottom, via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["BEVIA"]
)

# Add new additional vias for top electrode
new_additional_via_positions_top_3 = add_vias(
    x_start_vias_top_new_3, y_start_vias_top_new_3, num_vias_top, via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["TV"]
)

# Add ILD vias for the new additional bottom electrode
add_vias(
    x_start_vias_bottom_new_3, y_start_vias_bottom_new_3, num_vias_bottom, ild_via_size, via_spacing_bottom, row_spacing_bottom, top_cell, scaling_factor, layers["ILD"]
)

# Add ILD vias for the new additional top electrode
add_vias(
    x_start_vias_top_new_3, y_start_vias_top_new_3, num_vias_top, ild_via_size, via_spacing_top, row_spacing_top, top_cell, scaling_factor, layers["ILD"]
)


##########################################################
# Add customizable connection squares for additional vias
connection_square_dimensions = (100, 100)  # Width and height of the connection squares
add_connection_squares(additional_via_positions_bottom, connection_square_dimensions, layers["M1"], top_cell, scaling_factor)
add_connection_squares(additional_via_positions_top, connection_square_dimensions, layers["M1"], top_cell, scaling_factor)

# Define adjustable sizes for the electrode rectangles
top_electrode_size = (650, 160)  # Width and height for the top electrode rectangles, with increased gap
bottom_electrode_size = (660, 170)  # Width and height for the bottom electrode rectangles, with increased gap

# Add electrode rectangles centered around vias with horizontal offset
horizontal_offset = -100  # Example horizontal offset to move electrodes left or right
# Add bottom electrode rectangles centered around bottom electrode vias
add_electrode_rectangles(via_positions_bottom, top_electrode_size, bottom_electrode_size, scaling_factor, top_cell, 
                         horizontal_offset, gap=10, is_top_electrode=False)

# Add top electrode rectangles centered around top electrode vias
add_electrode_rectangles(via_positions_top, top_electrode_size, bottom_electrode_size, scaling_factor, top_cell, 
                         horizontal_offset, gap=10, is_top_electrode=True)

# Add additional electrode rectangles for new vias
#add_electrode_rectangles(additional_via_positions_bottom, top_electrode_size, bottom_electrode_size, scaling_factor, top_cell, horizontal_offset, gap=10, exclude_bottom_two_rows=True)
#add_electrode_rectangles(additional_via_positions_top, top_electrode_size, bottom_electrode_size, scaling_factor, top_cell, horizontal_offset, gap=10)

# Define adjustable radii for the electrode circles
top_electrode_radius = 125  # Radius of the top electrodes
bottom_electrode_radius = 130  # Radius of the bottom electrodes

# Add electrodes
add_electrodes(
    wafer_center[0], wafer_center[1],
    base_num_pmuts, base_spacing_x, base_spacing_y,
    top_electrode_radius=top_electrode_radius, bottom_electrode_radius=bottom_electrode_radius,  # Adjustable radii
    top_cell=top_cell, scaling_factor=scaling_factor
)

# Add IBEAM rings
add_ibeam_ring(
    wafer_center[0], wafer_center[1],
    base_num_pmuts, base_spacing_x, base_spacing_y,
    ibeam_ring_inner_radius=135, ibeam_ring_outer_radius=0,  # Example radii for IBEAM rings
    top_cell=top_cell, scaling_factor=scaling_factor
)

# Add electrode connection lines with adjustable lengths
bottom_connection_length = 90  # Adjustable length for the bottom connection line
top_connection_length = 90     # Adjustable length for the top connection line

add_electrode_connection_lines(
    wafer_center[0], wafer_center[1],
    base_num_pmuts, base_spacing_x, base_spacing_y,
    top_electrode_radius=top_electrode_radius, bottom_electrode_radius=bottom_electrode_radius,  # Use radii for connection line positioning
    bottom_connection_length=bottom_connection_length,  # Adjustable length
    top_connection_length=top_connection_length,        # Adjustable length
    top_cell=top_cell, scaling_factor=scaling_factor
)

# Add FE square
fe_square_x = 37790  # X position of the bottom-left corner of the FE square
fe_square_y = 37470  # Y position of the bottom-left corner of the FE square
fe_square_width = 8230  # Width of the FE square
fe_square_height = 8750  # Height of the FE square

add_fe_square(fe_square_x, fe_square_y, fe_square_width, fe_square_height, top_cell, scaling_factor)

# Add dicing lanes
dicing_lane_x = 37000  # Example X position for the bottom-left corner of the dicing lane
dicing_lane_y = 30700  # Example Y position for the bottom-left corner of the dicing lane
dicing_lane_width = 9200  # Width of the dicing lane
dicing_lane_height = 21500  # Height of the dicing lane
trench_offset = 50  # Offset for the trench to create a hollow center with 50 µm thickness

add_dicing_lanes(dicing_lane_x, dicing_lane_y, dicing_lane_width, dicing_lane_height, trench_offset, top_cell, scaling_factor)

# Use the updated function to add dicing lanes with FLD
add_dicing_lanes_with_fld(dicing_lane_x, dicing_lane_y, dicing_lane_width, dicing_lane_height, trench_offset, top_cell, scaling_factor)

# Define bond pad extension positions with XY control
bond_pad_offset_bottom = (12600, -2000)  # Offset for the extra bond pad at the bottom
bond_pad_offset_top = (12800, 1800)     # Offset for the extra bond pad at the top

# Add two additional bond pads
extra_bond_pad_size = (700, 700)  # Size of the bond pads

# Calculate the positions for the two additional bond pads
# For the bottom electrode bond pad at (20, 800, '-y', 15)
bottom_bond_pad_x = via_positions_bottom[-1][0]  # X position
bottom_bond_pad_y = via_positions_bottom[-1][1] - connection_params_bottom[-1][1]  # Start at the end of the last bottom connection

# For the top electrode bond pad at (20, 800, '+y', 15)
top_bond_pad_x = via_positions_top[-1][0]  # X position
top_bond_pad_y = via_positions_top[-1][1] + connection_params_top[-1][1]  # Start at the end of the last top connection

# Add bond pads at these positions
add_bond_pad(bottom_bond_pad_x, bottom_bond_pad_y, extra_bond_pad_size, layers["M1"], top_cell, scaling_factor)
add_bond_pad(top_bond_pad_x, top_bond_pad_y, extra_bond_pad_size, layers["M1"], top_cell, scaling_factor)

# Define a new cell to place multiple instances of the original design
final_cell_name = f'FINAL_LAYOUT_{int(time.time())}'
final_cell = gdsii_lib.new_cell(final_cell_name)

# Add references to the top cell at different positions
def add_references_to_final_cell(base_position, offsets, top_cell, final_cell):
    """Add references of the top cell to the final cell at specified offsets."""
    for offset_x, offset_y in offsets:
        ref_position = (base_position[0] + offset_x, base_position[1] + offset_y)
        cell_ref = gdspy.CellReference(top_cell, ref_position)
        final_cell.add(cell_ref)

########################################################################################


# Define the arrays to exclude in the final layout (0-based index)
arrays_to_exclude = [0, 7, 8, 9, 14, 15, 16]  # Adjust these indices based on your requirements

# Define positions for the references
base_position = (0, 0)
layout_width = 9200  # Width of a single layout instance
layout_height = 21500  # Height of a single layout instance

# Calculate the offsets for all possible arrays
offsets = []

# Determine the number of layouts per row (adjust as necessary)
layouts_per_row = 8  # This will give us 24, and we need 23, so we'll remove one later

# Add offsets for each row and column
for row in range(3):
    for col in range(layouts_per_row):
        offsets.append((col * (layout_width  -50), row * (layout_height - 50)))

# Remove the last offset to ensure we only have 23 instances
offsets = offsets[:23]

# Exclude specific arrays based on the indices in `arrays_to_exclude`
offsets = [offset for i, offset in enumerate(offsets) if i not in arrays_to_exclude]

# Add references to the final cell using the modified offsets list
add_references_to_final_cell(base_position, offsets, top_cell, final_cell)


###########################################################################################

output_directory = os.path.expanduser("~/Desktop")
if not os.path.isdir(output_directory):
    os.makedirs(output_directory)

output_path = os.path.join(output_directory, f'{final_cell_name}.gds')
gdsii_lib.write_gds(output_path)

def add_ibeam_connections(via_positions, connection_params, line_layer, top_cell, scaling_factor, line_width=50):
    """Add IBEAM connections that mimic the metal bond pad connections but are 50 µm wide."""
    for (x_pos, y_pos), (line_width_original, line_length, direction, offset) in zip(via_positions, connection_params):
        if direction == '-y':
            start_y = y_pos + offset * scaling_factor + 20 * scaling_factor  # Adjust start for bottom connections
        else:
            start_y = y_pos - offset * scaling_factor  # Adjust start for top connections

        ibeam_line = gdspy.Path(line_width * scaling_factor, (x_pos, start_y))
        ibeam_line.segment(line_length * scaling_factor, direction=direction, layer=line_layer[0], datatype=line_layer[1])
        top_cell.add(ibeam_line)

# Define the IBEAM layer
ibeam_layer = layers["IB"]

# Add IBEAM connections mimicking the bottom metal bond pad connections
# add_ibeam_connections(via_positions_bottom, connection_params_bottom, ibeam_layer, top_cell, scaling_factor, line_width=50)

# Add IBEAM connections mimicking the top metal bond pad connections
# add_ibeam_connections(via_positions_top, connection_params_top, ibeam_layer, top_cell, scaling_factor, line_width=50)

# Save the updated GDS file with IBEAM connections included
gdsii_lib.write_gds(output_path)


print(f"GDS file '{output_path}' updated with IBEAM connections successfully.")

