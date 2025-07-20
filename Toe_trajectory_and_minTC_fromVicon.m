clc;
clear;
close all;

% Add BTK path and read the C3D file
addpath('../btk');
c3dfile = 'F:\Vicon_data_full sets\236_24082022\236\236_24082022 level2new.c3d';
btkData = btkReadAcquisition(c3dfile);

% Extract the marker data
Markers = btkGetMarkers(btkData);

% Extract the frequency of data collection
f = btkGetPointFrequency(btkData);

%% Filtering and initial steps
[B, A] = butter(4, 6 / (f / 2), 'low');

% Left foot markers
filtLTOE = filtfilt(B, A, Markers.LTOE);
filtLHEE = filtfilt(B, A, Markers.LHEE);
filtLASI = filtfilt(B, A, Markers.LASI);

% Right foot markers
filtRTOE = filtfilt(B, A, Markers.RTOE);
filtRHEE = filtfilt(B, A, Markers.RHEE);
filtRASI = filtfilt(B, A, Markers.RASI);

% Initial Position for Normalization
initial_toe_position_L = filtLTOE(1, 3);  % 3 for the z-coordinate
initial_toe_position_R = filtRTOE(1, 3);

%% Calculate toe clearances for both feet
time_per_frame = 1/f;

[all_minTCs_L, all_minTCs_normalized_L, all_step_cycles_L] = computeToeClearances(filtLTOE, filtLHEE, filtLASI, initial_toe_position_L, time_per_frame);
[all_minTCs_R, all_minTCs_normalized_R, all_step_cycles_R] = computeToeClearances(filtRTOE, filtRHEE, filtRASI, initial_toe_position_R, time_per_frame);

% Compute mean and standard deviation
mean_minTCs_L = mean(all_minTCs_L(:, 2));
std_minTCs_L = std(all_minTCs_L(:, 2));

mean_minTCs_normalized_L = mean(all_minTCs_normalized_L(:, 2));
std_minTCs_normalized_L = std(all_minTCs_normalized_L(:, 2));

mean_minTCs_R = mean(all_minTCs_R(:, 2));
std_minTCs_R = std(all_minTCs_R(:, 2));

mean_minTCs_normalized_R = mean(all_minTCs_normalized_R(:, 2));
std_minTCs_normalized_R = std(all_minTCs_normalized_R(:, 2));

%% Plot the toe clearance trajectories for each step cycle with percentage of step cycle

% Left foot
figure;
hold on;
for i = 1:length(all_step_cycles_L)
    percentage_of_step_cycle = linspace(0, 100, length(all_step_cycles_L{i}));
    plot(percentage_of_step_cycle, all_step_cycles_L{i}, 'LineWidth', 1.5);
end
xlabel('Percentage of Swing Phase (Toe-Off to Heel-Strike)');
ylabel('Normalized Toe Clearance (mm)');
title('Left Toe Clearance Trajectories for Each Swing Phase');
grid on;
hold off;

% Right foot
figure;
hold on;
for i = 1:length(all_step_cycles_R)
    percentage_of_step_cycle = linspace(0, 100, length(all_step_cycles_R{i}));
    plot(percentage_of_step_cycle, all_step_cycles_R{i}, 'LineWidth', 1.5);
end
xlabel('Percentage of Swing Phase (Toe-Off to Heel-Strike)');
ylabel('Normalized Toe Clearance (mm)');
title('Right Toe Clearance Trajectories for Each Swing Phase');
grid on;
hold off;


%% Function to compute toe clearances
function [all_minTCs, all_minTCs_normalized, all_step_cycles] = computeToeClearances(filtTOE, filtHEE, filtASI, initial_toe_position, time_per_frame)

    % Distances to hip markers to identify heel strike and toe-off
    distance_TOE = filtTOE(:, 2) - filtASI(:, 2); % 2 for the y-coordinate
    distance_HEE = filtHEE(:, 2) - filtASI(:, 2);
    maxima_points = islocalmax(distance_HEE);
    minima_points = islocalmin(distance_TOE);

    % Switch if necessary (y-axis direction in the coordinate system)
    heel_strike_frames = find(minima_points);
    toe_off_frames = find(maxima_points);

    all_toe_off_frames = toe_off_frames;
    all_heel_strike_frames = heel_strike_frames;

    all_minTCs = [];
    all_minTCs_normalized = [];
    all_step_cycles = {};

    for i = 1:length(toe_off_frames)
        next_heel_strike_frame = find(heel_strike_frames > toe_off_frames(i), 1);
        
        if isempty(next_heel_strike_frame)
            break;
        end
        
        segment_frames = toe_off_frames(i):heel_strike_frames(next_heel_strike_frame);
        segment_toe = filtTOE(segment_frames, :);
        
        z_segment_toe = segment_toe(:, 3);
        z_segment_toe_normalized = z_segment_toe - initial_toe_position;

        maxTC_points = islocalmax(z_segment_toe);
        maxTC_frames = find(maxTC_points);
        minTC_points = islocalmin(z_segment_toe);
        minTC_frames = find(minTC_points);

        for j = 1:length(maxTC_frames) - 1
            interval_min_frames = minTC_frames(minTC_frames > maxTC_frames(j) & minTC_frames < maxTC_frames(j+1));
            if isempty(interval_min_frames)
                continue;
            end

            minTC = min(z_segment_toe(interval_min_frames));
            minTC_normalized = min(z_segment_toe_normalized(interval_min_frames));

            minTC_frame = interval_min_frames(z_segment_toe(interval_min_frames) == minTC);
            minTC_time = (segment_frames(1) + minTC_frame - 1) * time_per_frame;
            
            % Time stamp (1. column) and minTC (2. column) for every step
            all_minTCs = [all_minTCs; minTC_time, minTC];
            all_minTCs_normalized = [all_minTCs_normalized; minTC_time, minTC_normalized];
        end

        all_step_cycles{end+1} = z_segment_toe_normalized;
    end

    return
end
