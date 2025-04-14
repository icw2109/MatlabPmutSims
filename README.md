# Focused Ultrasound Simulation with PMUT Arrays

This project contains MATLAB code to simulate focused ultrasound propagation using 2D and 3D PMUT arrays through soft tissue and skull. It includes pressure field visualization, beamforming, skull modeling, and safety evaluations such as heating and cavitation effects.

---

## File Overview

### 1. `simulation_3D_focus.m`
- **Type:** 3D simulation (no skull)
- **Description:** 
  - Simulates beam focusing in a homogeneous medium using a full 3D PMUT array.
  - Uses time-delay focusing on a 20×20 circular-element grid.
  - Outputs a pressure field and verifies the focal spot at 50 mm.

### 2. `simulation_with_skull.m`
- **Type:** 2D simulation (with skull)
- **Description:**
  - Models ultrasound propagation through a skull layer with heating and cavitation analysis over multiple pulses.
  - Outputs temperature rise maps, cavitation index, and the final pressure field.

### 3. `simulation_with_skull_simplefocus.m`
- **Type:** Simplified 2D skull model
- **Description:**
  - Focuses on the single-pulse pressure field and focal point accuracy.
  - Does not calculate heating or cavitation effects.

---

## Key Simulation Parameters

- **Frequency:** 650 kHz  
- **Focal Distance:** 50 mm  
- **Skull Thickness:** 5 mm  
- **Element Spacing:** 400 μm  
- **Element Width:** 280 μm  
- **Grid Resolution:** 180 μm  
- **Absorption Coefficient:** 0.002 dB/MHz²  
- **Nonlinearity (B/A):** 6  
- **Time Step:** 40 ns  
- **Medium Model:** Homogeneous soft tissue with a layered skull

---

## How to Run

1. Download and install the [k-Wave MATLAB toolbox](http://www.k-wave.org/).
2. Add it to your MATLAB path.
3. Run any of the scripts:
   ```matlab
   run('simulation_3D_focus.m')
   run('simulation_with_skull.m')
   run('simulation_with_skull_simplefocus.m')

