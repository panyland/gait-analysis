# MATLAB scripts for analyzing human gait using motion capture (Vicon) and IMU data.

## Toe_trajectory_and_minTC_fromVicon.m
Analyzes 3D motion capture data to extract toe clearance during the swing phase. Computes and visualizes minimum toe clearance (minTC) for left and right feet using .c3d files.

-Uses BTK toolkit 
-Detects toe-off and heel-strike events 
-Normalizes and plots toe clearance trajectories

## Gait_Analysis_with_IMU.m
Processes IMU data to extract gait features such as heel strikes, stride lengths, vertical foot displacement, and symmetry measures.

-Works with multi-sensor IMU data
-Filters and transforms acceleration data
-Computes stride metrics and autocorrelation
-Includes visualizations of gait dynamics
