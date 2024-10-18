function [params] = radonVbarDir(inpDat, params, plotFlag)
%% Apply Radon to define the velocity search bounds (vBounds) to + or -
theta = -85:85;  % the range of angles we're considering

% Check if plotFlag is provided, if not, set it to 0 (no plotting)
if nargin < 2
    plotFlag = 0;
end

% Convert time to numeric 'datenum' format for plotting
numeric_time = inpDat.tGrid(1,:) / (3600.0 * 24) + datenum(1970,1,1);
plot_time = datetime(numeric_time, 'ConvertFrom', 'datenum');  % Use datetime for better interpretation

% Perform Radon transform on all angles between -85 and +85 degrees
[R, xp] = radon(inpDat.rawGrid, theta);

% Exclude energy from stationary objects (~0 degrees = shore normal)
excluded_angles = (theta >= -5 & theta <= 5);  % Exclude shore-normal angles
R(:, excluded_angles) = 0;

% Find peaks in the Radon transform to identify the dominant angle

maxR = max(R, [], 1);
[peaks, locations] = islocalmax(maxR);

% Plot the results if plotFlag is set
if plotFlag
    figure();
    % Plot the intensity data (raw pixel data)
    pcolor(numeric_time, inpDat.yGrid(:,1), inpDat.rawGrid);
    hold on;
    shading flat;
    colormap gray;
    title('Detected Oblique Foam Traces');
    xlabel('Time');
    ylabel('Y-Position (m)');
    datetick('x', 'HH:MM:SS', 'keeplimits');  % Convert numeric time back to readable format
    hold off;

    figure();
    plot(theta, maxR, 'LineWidth', 2);
    hold on;
    scatter(theta(peaks), maxR(peaks), 'r*');  % Plot detected peaks
    xlabel('Angle (degrees)');
    ylabel('Max Radon Energy');
    title('Radon Transform - Detected Peaks');
    grid on;
end

% Detect the angle of the highest peak
detected_angle = theta(locations(peaks == max(peaks)));  % Angle of the highest peak

neg_angle_energy = sum(peak1d(find(theta >= -60 & theta >= -5)));  % Sum energy over all negative angles
pos_angle_energy = sum(peak1d(find(theta >= 5 & theta >= 60)));  % Sum energy over all positive angles

% Interpretation of the detected angle for setting velocity bounds (vBounds)
if neg_angle_energy > pos_angle_energy
    disp('The foam is drifting south');
    params.vBounds = [-2.5 -0.01];  % Adjust bounds for southward motion
elseif neg_angle_energy < pos_angle_energy
    disp('The foam is drifting north');
    params.vBounds = [0.01 2.5];  % Adjust bounds for northward motion
else
    disp('The dominant angle is within the excluded range.');
end