clc
clear
close all;

%% Import and setup
sessionData = importSession('222_6mwt');

samplePeriod = 1 / 200; 
[sessionData, time] = resampleSession(sessionData, samplePeriod);

numberOfSamples = length(time);

% IMU numbers (from list in sessionData)
imu_numbers = struct('left_shank', 9, 'left_foot', 4, 'right_shank', 8, 'right_foot', 1, 'hip', 6);

sides = fieldnames(imu_numbers);
nSides = numel(sides);

% Initialize variables for acceleration, gyroscope, and rotations
acc = struct();
gyro = struct();
rotMat = struct();
euler = struct();
acc_filtered = struct();

% Design the Butterworth filter
n = 4; 
wn = 0.1;
[b, a] = butter(n, wn, 'low');

for i = 1:nSides
    side = sides{i};
    imu_nro = imu_numbers.(side);

    % Quaternions to rotation matrices
    quaternion = sessionData.(sessionData.deviceNames{imu_nro}).quaternion.vector;
    rotMat.(side) = quatern2rotMat(quaternion);

    % Rotation matrices to euler angles
    euler.(side) = manualRotm2eul(rotMat.(side));

    % Extract and process accelerometer data
    acc.(side).X = sessionData.(sessionData.deviceNames{imu_nro}).sensors.accelerometerX * 9.81;
    acc.(side).Y = sessionData.(sessionData.deviceNames{imu_nro}).sensors.accelerometerY * 9.81;
    acc.(side).Z = sessionData.(sessionData.deviceNames{imu_nro}).sensors.accelerometerZ * 9.81;
    acc_vector = [acc.(side).X, acc.(side).Y, acc.(side).Z];
    
    % Filter acceleration
    for j = 1:3
        acc_filtered.(side)(:,j) = filtfilt(b, a, acc_vector(:,j));
    end
    
    % Extract gyroscope data
    gyro.(side).X = sessionData.(sessionData.deviceNames{imu_nro}).sensors.gyroscopeX;
    gyro.(side).Y = sessionData.(sessionData.deviceNames{imu_nro}).sensors.gyroscopeY;
    gyro.(side).Z = sessionData.(sessionData.deviceNames{imu_nro}).sensors.gyroscopeZ;
    gyro_vector = [gyro.(side).X, gyro.(side).Y, gyro.(side).Z];

    % Rotate back to global frame
    acc.(side) = applyRotation(acc_filtered.(side), rotMat.(side));
end

%% Detect heelstrikes

% Parameters for peak detection
minPeakHeight = 12;  
minPeakDistance = 100;  
maxPeakWidth = 20; 

% Find peaks for left shank vertical acceleration
[left_peaks, left_locs] = findpeaks(acc_filtered.left_shank(:,1), ...
                                    'MinPeakHeight', minPeakHeight, ...
                                    'MinPeakDistance', minPeakDistance, ...
                                    'MaxPeakWidth', maxPeakWidth);

% Find peaks for right shank vertical acceleration
[right_peaks, right_locs] = findpeaks(acc_filtered.right_shank(:,1), ...
                                     'MinPeakHeight', minPeakHeight, ...
                                     'MinPeakDistance', minPeakDistance,...
                                     'MaxPeakWidth', maxPeakWidth);

% Plot shank accelerations with heelstrikes
figure;
subplot(2,1,1);
plot(acc_filtered.left_shank(:,1)); hold on;
plot(left_locs, left_peaks, 'ro');
title('Left Shank Vertical Acceleration with Detected Heel Strikes');
xlabel('Sample Number');
ylabel('Acceleration (m/s^2)');

subplot(2,1,2);
plot(acc_filtered.right_shank(:,1)); hold on;
plot(right_locs, right_peaks, 'ro');
title('Right Shank Vertical Acceleration with Detected Heel Strikes');
xlabel('Sample Number');
ylabel('Acceleration (m/s^2)');
hold off;

% Plotting acceleration for the Left Foot
figure;
subplot(3, 1, 1);  
plot(time, acc_filtered.left_foot(:, 1), 'r'); 
title('Left Foot Acceleration - X Axis');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

subplot(3, 1, 2); 
plot(time, acc_filtered.left_foot(:, 2), 'g');  
title('Left Foot Acceleration - Y Axis');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

subplot(3, 1, 3);  
plot(time, acc_filtered.left_foot(:, 3), 'b'); 
title('Left Foot Acceleration - Z Axis');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

% Plotting acceleration for the Right Foot
figure;
subplot(3, 1, 1); 
plot(time, acc_filtered.right_foot(:, 1), 'r');  
title('Right Foot Acceleration - X Axis');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

subplot(3, 1, 2);  
plot(time, acc_filtered.right_foot(:, 2), 'g');  
title('Right Foot Acceleration - Y Axis');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

subplot(3, 1, 3); 
plot(time, acc_filtered.right_foot(:, 3), 'b');  
title('Right Foot Acceleration - Z Axis');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

%%  Stride length

% Calculate the number of complete strides for each foot
nStridesLeft = length(left_locs) - 1;
nStridesRight = length(right_locs) - 1;

% Initialize arrays to store stride lengths
stride_lengths_left = zeros(nStridesLeft, 1);
stride_lengths_right = zeros(nStridesRight, 1);

% Calculate stride lengths for the left foot
for i = 1:nStridesLeft
    startIndex = left_locs(i);
    endIndex = left_locs(i + 1) - 1;
    strideAccel = acc_filtered.left_foot(startIndex:endIndex, 1);  % Assuming X-axis is forward
    strideAccel = strideAccel - mean(strideAccel);  
    strideVel = cumtrapz(strideAccel) * samplePeriod;
    strideDisp = cumtrapz(strideVel) * samplePeriod;
    stride_lengths_left(i) = abs(strideDisp(end));
end

% Calculate stride lengths for the right foot
for i = 1:nStridesRight
    startIndex = right_locs(i);
    endIndex = right_locs(i + 1) - 1;
    strideAccel = acc_filtered.right_foot(startIndex:endIndex, 1);  % Assuming X-axis is forward
    strideAccel = strideAccel - mean(strideAccel);  
    strideVel = cumtrapz(strideAccel) * samplePeriod;
    strideDisp = cumtrapz(strideVel) * samplePeriod;
    stride_lengths_right(i) = abs(strideDisp(end));
end

% Display the calculated stride lengths
disp('Left Foot Stride Lengths (m):');
disp(stride_lengths_left);
disp('Right Foot Stride Lengths (m):');
disp(stride_lengths_right);

%% Vertical displacement of the foot during a swing

offset_left = 7.475;
offset_right = 7.247;
acc_filtered.left_foot(:, 3) = acc_filtered.left_foot(:, 3) - offset_left;
acc_filtered.right_foot(:, 3) = acc_filtered.right_foot(:, 3) - offset_right;

% Calculate the number of complete strides for each foot
nStridesLeft = length(left_locs) - 1;
nStridesRight = length(right_locs) - 1;

% Initialize arrays to store maximum vertical displacements
max_vertical_displacements_left = zeros(nStridesLeft, 1);
max_vertical_displacements_right = zeros(nStridesRight, 1);

% Swing phase starts slightly after the heel strike; define the offset
swing_phase_offset = 20; 

% Initialize figure for plotting all trajectories
figure;
% Subplot for the left foot
subplot(1, 2, 1);
hold on;

% Process each stride for the left foot
for i = 1:nStridesLeft
    startIndex = left_locs(i) + swing_phase_offset;
    endIndex = left_locs(i + 1) - 1;
    verticalAccel = acc_filtered.left_foot(startIndex:endIndex, 3);  % Assuming Z-axis is vertical
    verticalAccel = verticalAccel - mean(verticalAccel);  

    % Double integration to get vertical displacement
    verticalVel = cumtrapz(verticalAccel) * samplePeriod;
    verticalDisp = cumtrapz(verticalVel) * samplePeriod;

    % Store the maximum vertical displacement for this stride
    max_vertical_displacements_left(i) = max(abs(verticalDisp)); 

    % Plotting the vertical trajectory for this stride
    plot(linspace(0, 100, length(verticalDisp)), abs(verticalDisp), 'DisplayName', sprintf('Left Stride %d', i));
end
title('Left Foot Vertical Trajectories');
xlabel('Percentage of Swing Phase');
ylabel('Vertical Displacement (m)');
legend show;
hold off;

% Subplot for the right foot
subplot(1, 2, 2);
hold on;

% Process each stride for the right foot
for i = 1:nStridesRight
    startIndex = right_locs(i) + swing_phase_offset;
    endIndex = right_locs(i + 1) - 1;
    verticalAccel = acc_filtered.right_foot(startIndex:endIndex, 3);  % Assuming Z-axis is vertical
    verticalAccel = verticalAccel - mean(verticalAccel);  

    % Double integration to get vertical displacement
    verticalVel = cumtrapz(verticalAccel) * samplePeriod;
    verticalDisp = cumtrapz(verticalVel) * samplePeriod;

    % Store the maximum vertical displacement for this stride
    max_vertical_displacements_right(i) = max(abs(verticalDisp));  

    % Plotting the vertical trajectory for this stride
    plot(linspace(0, 100, length(verticalDisp)), abs(verticalDisp), 'DisplayName', sprintf('Right Stride %d', i));
end
title('Right Foot Vertical Trajectories');
xlabel('Percentage of Swing Phase');
ylabel('Vertical Displacement (m)');
legend show;
hold off;

% Display the maximum vertical displacements
disp('Maximum Vertical Displacements for Left Foot (m):');
disp(max_vertical_displacements_left);
disp('Maximum Vertical Displacements for Right Foot (m):');
disp(max_vertical_displacements_right);

%% Compute autocorrelation and find the first two peaks to the right of zero lag

verticalAccelLeft = (acc_filtered.left_foot(:, 3) - mean(acc_filtered.left_foot(:, 3))) / std(acc_filtered.left_foot(:, 3));
verticalAccelRight = (acc_filtered.right_foot(:, 3) - mean(acc_filtered.right_foot(:, 3))) / std(acc_filtered.right_foot(:, 3));
verticalAccelHip = (acc_filtered.hip(:, 1) - mean(acc_filtered.hip(:, 1))) / std(acc_filtered.hip(:, 1));

% Compute autocorrelation
[autocorrLeft, lagsLeft] = xcorr(verticalAccelLeft, 'coeff');
[autocorrRight, lagsRight] = xcorr(verticalAccelRight, 'coeff');
[autocorrHip, lagsHip] = xcorr(verticalAccelHip, 'coeff');

%% Find first two peaks to the right of zero lag

[firstPeakHip, secondPeakHip] = findFirstTwoPeaks(autocorrHip, lagsHip);

disp('First two peaks (hip):');
disp(firstPeakHip);
disp(secondPeakHip);

% Plot autocorrelation functions with marked peaks
figure;
subplot(3, 1, 1);
plot(lagsLeft, autocorrLeft); hold on;
title('Autocorrelation of Left Foot Vertical Acceleration');
xlabel('Lag');
ylabel('Autocorrelation');
legend('show');
grid on;

subplot(3, 1, 2);
plot(lagsRight, autocorrRight); hold on;
title('Autocorrelation of Right Foot Vertical Acceleration');
xlabel('Lag');
ylabel('Autocorrelation');
legend('show');
grid on;

subplot(3, 1, 3);
plot(lagsHip, autocorrHip); hold on;
plot(firstPeakHip(2), firstPeakHip(1), 'ro', 'MarkerSize', 8, 'DisplayName', 'First Peak');
plot(secondPeakHip(2), secondPeakHip(1), 'go', 'MarkerSize', 8, 'DisplayName', 'Second Peak');
title('Autocorrelation of Hip Vertical Acceleration');
xlabel('Lag');
ylabel('Autocorrelation');
legend('show');
grid on;

%%

% Angular rate data vectors
left_foot_angular_rate = gyro.left_foot.Y; 
right_foot_angular_rate = gyro.right_foot.Y;

% Calculate cross-correlation
[ccor, lags] = xcorr(left_foot_angular_rate, right_foot_angular_rate, 'coeff');

% Find the maximum correlation and its index
[maxCorr, maxIndex] = max(ccor);

% Display the maximum correlation coefficient
disp(['Maximum Cross-Correlation Coefficient: ', num2str(maxCorr)]);

% Plot the cross-correlation function
figure;
plot(lags, ccor);
title('Cross-Correlation of Left and Right Foot Angular Rates');
xlabel('Lag');
ylabel('Cross-Correlation');
grid on;

%% Functions

function eulerAnglesDeg = manualRotm2eul(rotMatrices)
    % Initialize the output matrix
    N = size(rotMatrices, 3); 
    eulerAnglesDeg = zeros(N, 3); % Each row will hold [yaw, pitch, roll] for a time point
    
    for i = 1:N
        R = rotMatrices(:,:,i); 

        % Assuming the rotation matrix is in ZYX format (yaw-pitch-roll)
        if R(3,1) < 1
            if R(3,1) > -1
                % Non-singular case
                pitch = asin(-R(3,1));
                yaw = atan2(R(3,2)/cos(pitch), R(3,3)/cos(pitch));
                roll = atan2(R(2,1)/cos(pitch), R(1,1)/cos(pitch));
            else
                % Gimbal lock at the lower limit
                pitch = pi / 2;
                yaw = atan2(R(1,2), R(1,3));
                roll = 0;
            end
        else
            % Gimbal lock at the upper limit
            pitch = -pi / 2;
            yaw = atan2(-R(1,2), -R(1,3));
            roll = 0;
        end

        % Convert radians to degrees
        eulerAnglesDeg(i,:) = rad2deg([yaw, pitch, roll]);
    end
end

function globalData = applyRotation(localData, rotationMatrices)
    % Angle of incline
    angle_deg = -5;  
    angle_rad = deg2rad(angle_deg);

    % Define the rotation matrix for the Y-axis incline
    Ry = [cos(angle_rad) 0 sin(angle_rad); 0 1 0; -sin(angle_rad) 0 cos(angle_rad)];

    % Initialize the matrix to store global data
    globalData = zeros(size(localData));

    % Apply the rotation to each sample
    for i = 1:size(localData, 1)
        % Multiply the per-sample rotation matrix by the incline adjustment matrix
        adjustedRotationMatrix = Ry * rotationMatrices(:, :, i);
        
        % Now apply the adjusted rotation matrix to the local data
        globalData(i, :) = (adjustedRotationMatrix * localData(i, :)');
    end
end

function [firstPeak, secondPeak] = findFirstTwoPeaks(autocorr, lags)
    % Find the index of zero lag
    zeroLagIndex = find(lags == 0);
    
    % Extract the positive lag part of the autocorrelation
    positiveLags = lags(zeroLagIndex+1:end);
    positiveAutocorr = autocorr(zeroLagIndex+1:end);
    
    % Find peaks in the positive lag part
    [peaks, locs] = findpeaks(positiveAutocorr);
    
    % Ensure there are at least two peaks
    if length(peaks) < 2
        error('Less than two peaks found in the autocorrelation.');
    end
    
    % Get the first two peaks
    firstPeak = [peaks(1), positiveLags(locs(1))];
    secondPeak = [peaks(2), positiveLags(locs(2))];
end
