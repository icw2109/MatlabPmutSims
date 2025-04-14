%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Focused Ultrasound Simulation with PMUT Arrays
%
% This project contains MATLAB code to simulate focused ultrasound
% propagation using 2D and 3D PMUT arrays through soft tissue and skull.
% It includes pressure field visualization, beamforming, skull modeling,
% and safety evaluations such as heating and cavitation effects.
%
% -------------------------------------------------------------------------
%  FILE OVERVIEW
%
% 1. simulation_3D_focus.m
%    - 3D simulation (no skull)
%    - Simulates beam focusing in homogeneous medium using a full 3D PMUT array
%    - Uses time-delay focusing on a 20x20 circular-element grid
%    - Outputs pressure field and verifies focal spot at 50 mm
%
% 2. simulation_with_skull.m
%    - 2D simulation (with skull)
%    - Models ultrasound passing through a skull layer with heating and
%      cavitation analysis over multiple pulses
%    - Outputs temperature rise map, cavitation index, and final pressure field
%
% 3. simulation_with_skull_simplefocus.m
%    - Simplified 2D skull model
%    - Focuses on single-pulse pressure field and focal point accuracy
%    - Does not calculate heating or cavitation
%
% -------------------------------------------------------------------------
% 🧪 KEY SIMULATION PARAMETERS
%
% Frequency             : 650 kHz
% Focal Distance        : 50 mm
% Skull Thickness       : 5 mm
% Element Spacing       : 400 µm
% Element Width         : 280 µm
% Grid Resolution       : 180 µm
% Absorption Coefficient: 0.002 dB/MHz^2
% Nonlinearity (B/A)    : 6
% Time Step             : 40 ns
% Medium Model          : Homogeneous soft tissue + layered skull
%
% -------------------------------------------------------------------------
% ▶ HOW TO RUN
%
% 1. Download and install the k-Wave MATLAB toolbox: http://www.k-wave.org/
% 2. Add it to your MATLAB path.
% 3. Run any of the scripts:
%
%    >> run('simulation_3D_focus.m')
%    >> run('simulation_with_skull.m')
%    >> run('simulation_with_skull_simplefocus.m')
%
% -------------------------------------------------------------------------
% OUTPUTS
%
% - 3D or 2D pressure field plots
% - Pressure along central axis
% - Focal spot validation (target vs. actual)
% - Optional:
%     • Temperature rise map
%     • Cavitation index heatmap
%
% -------------------------------------------------------------------------
%  APPLICATIONS
%
% - Transcranial Focused Ultrasound (tFUS)
% - Neuromodulation studies
% - Therapeutic ultrasound safety analysis
% - PMUT array performance benchmarking
%
% -------------------------------------------------------------------------
%  NOTES FROM DEVELOPER
%
% This codebase was written entirely from scratch and was designed to handle
% dynamic focusing using spatial delays across a 2D/3D PMUT matrix. The model
% simulates realistic skull interaction and cumulative energy effects over time.
%
% Challenges included:
% - Multi-element focusing without built-in transducer objects
% - Delay interpolation and array geometry alignment
% - Large-scale 3D memory optimization
%
% Feel free to reach out to collaborate or discuss extensions (e.g., MRI-based
% head models, frequency modulation, transducer optimization, etc).
%
% -------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
