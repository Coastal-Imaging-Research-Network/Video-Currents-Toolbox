%% Updated version of "vidCurrentsGenManual" 
%   where the original function is defined:
% videoCurrentOut = videoCurrentGen(stack, time, xy, vB, fkB, Twin, Tstep, plotFlag);
%   & the inputs are as follows:
% stack: data [
% time: time vector (matlab datenum) with size [Mx1s]
% xy: xyzAll [Nx3] & x-search locations
% vB: (vBounds) [minV maxV], units m/s or vector of desired velocity steps,
%       to set this to empty, [], to use defaults [-3 3] ***These defaults
%       will change****
% fkB: (fkBounds) [fmin fmax kmin kmax], vector of frequency and wavenumber
%      bounds energy out of side of these bounds will be set to 0.  Useful
%      to eliminate some of the wave contamination that leaks in.  Set this
%      to empty, [], to use defaults.
% Twin: the time length of the FFT window (in points)
% For 2 Hz data, Twin = 128 will yield a 64s average current
% Tstep:
% plotFlag: optional, if true (~=0) will display a running plot of
%      the data processing

%% Part 1: Set Up

% params definitions moved to a new file 

% fileSearchPath is the folder within which the function vBarRawFile will
% search for the "raw" data file 
fileSearchPath = ("C:\Users\13emo3\OneDrive - Queen's University\PhD\Chapter 4 vBar\Chapter 4 Data\FTP Downloads");

[T, RAW, XYZ, CAM] = loadVbarRawFile(fileSearchPath, params.searchDate, params.searchX);

% time is loaded as epoch, so convert it to datetime for plotting
% interpretation
mtime = (T/(3600*24)+datenum(1970,1,1))';
dTime = datetime(mtime, 'ConvertFrom', 'datenum');

%% Part 2: "sorted"
% Sort data & discard unneeded data
% Get the sorted indices based on Y values
Y = XYZ(:,2);
[~, sorted_indices] = sort(Y);

% Reorder RAW and XYZ based on the sorted indices
sorted.RAW = RAW(:, sorted_indices);
sorted.XYZ = XYZ(sorted_indices, :);
sorted.CAM = CAM(sorted_indices, :); 

if sorted.XYZ(1,2) < params.yLims(1) || sorted.XYZ(end,2) > params.yLims(2)
    iy1 = find(sorted.XYZ(:,2) >= params.yLims(1), 1, 'first');
    iy2 = find(sorted.XYZ(:,2) <= params.yLims(2), 1, 'last');
    sorted.XYZ = sorted.XYZ(iy1:iy2,:);
    sorted.RAW = sorted.RAW(:,iy1:iy2);
end

%clearvars XYZ RAW iy1 iy2 sorted_indices

%% Part 3: "unik"
% Remove Redundant x-y Locations
% Find unique XY locations and average corresponding RAW data
[~, uInd, ~] = unique(sorted.XYZ(:, 1:2), 'rows', 'stable');
uXY = sorted.XYZ(uInd, 1:2);
unik.CAM = sorted.CAM(uInd); 
uRAW = nan(length(T), length(uXY));

for iX = 1:length(uXY)
    indices = ismember(sorted.XYZ(:, 1:2), uXY(iX, :), 'rows');
    uRAW(:, iX) = mean(sorted.RAW(:, indices), 2);
end

unik.RAW = uRAW;
unik.XY = uXY;

%clearvars uInd uXY uRAW indices

%% Part 4: "gridded"
% Project Input Variables onto Regular Grids
% Define the search domain
% Create grid space to project RAW onto
[tq, yq] = meshgrid(T, y);
gridded.stack = griddata(unik.XY(:, 2), T, unik.RAW, yq, tq);
gridded.xy(:,2) = y;
gridded.xy(:,1) = params.searchX;

% ** what i'm doing with CAM index here could be simplified
F = scatteredInterpolant(XYZ(:,1), XYZ(:,2), double(CAM));
gridded.CAM = round(F(gridded.xy(:,1), gridded.xy(:,2)));

bad = isnan(gridded.stack);
if sum(bad, 'all') > 0
    gridded.stackOg = gridded.stack; % Save original as a different name
    gridded.stack = fillmissing(gridded.stack, 'nearest'); % Fill NaNs with mean
end

% clearvars unik bad

%% Part 5: Apply Radon
% This section limits the velocity search bounds (vBounds) to + or - 

% This part could be optional -- sometimes the direction is known

% Define the range of angles for the Radon transform
theta = 0:179;

% Perform the Radon transform
[R, xp] = radon(gridded.stack, theta);

% Exclude unwanted angle ranges
excluded_angles = (theta >= 85 & theta <= 95) | (theta >= 175 | theta <= 5);
R(:, excluded_angles) = 0;

% Find the angle with the maximum value in the Radon transform
[max_val, max_idx] = max(R(:));
[max_row, max_col] = ind2sub(size(R), max_idx);
dominant_angle = theta(max_col);

% Interpret the dominant angle
if dominant_angle > 95 && dominant_angle < 175
    disp('The foam is drifting north');
    params.vBounds = [0.01 2];                                 % <--- *see note in "vidCurrentsParams.m"
elseif dominant_angle > 5 && dominant_angle < 85
    disp('The foam is drifting south');
    params.vBounds = [-2 -0.01];
else
    disp('The dominant angle is within the excluded range.');
end

%% Part 6: Run the main function 

% Initialize the video currents table
vc150 = table();
count = 1;

% I'm not certain on how to do this here, but I'm trying to identify the
% locations at which the camera changes so that the defined tiles don't
% include camera seams 

% Calculate the difference between adjacent elements
diff_CAM = diff(gridded.CAM);

% Find the indices where the difference is non-zero
change_indices = find(diff_CAM ~= 0) + 1;

% Maybe there should be one "tileCenters" for each camera? 
tileCenters = y(1):10:y(end); % Define tile centers with a step of 10 meters

for iCenter = 1:length(tileCenters)
    % the below commented out stuff needs to be re-written (with the cam
    % lims) 
    % 
    % % Calculate the start and end indices for the current tile
    % centerY = tileCenters(iCenter);
    % startY = centerY - currentTileSize / 2;
    % endY = centerY + currentTileSize / 2;
    % 
    % % Find the corresponding indices in the y array
    % i1 = find(y >= startY, 1);
    % i2 = find(y <= endY, 1, 'last');

    % Ensure indices are within y-search bounds
    if isempty(i1) || isempty(i2)
        continue;
    end

    % Run video-current-toolbox
    vC = videoCurrentGen(gridded.stack(i1:i2,:)', params.mtime, gridded.xy(i1:i2,:), params.vBounds, params.fkBounds, params.tWindow, params.tStep, params.plotFlag);
    vC.ci1 = vC.ci(:, 1);
    vC.ci2 = vC.ci(:, 2);
    vC = rmfield(vC, "ci");

    % Save vC to the table 
    vc150.vC(count) = vC;
    vc150.y(count) = centerY;
    count = count + 1;
end


% clearvars -except vc175 gridded
%% Part 7a: Filter & Combine Results 
% New method of calculating single meanV from timeseries 

for i = 1:height(vc150)
    vc150.wV(i) = wmean(vc150.vC(i).meanV, 1./vc150.vC(i).stdV, 'omitnan');
end

%% Part 7b: Filter & Combine Results 
% Old method, as described in paragraph [47] in Chickadel et al., 2003

for i = 1:height(vc150)
    % Criteria 1:
    % Significance of fit from the model skill must be > 90%
    iPass1 = find(vc150.vC(i).prob >= 0.9);
    pass1{i} = iPass1;
    totPass1(i) = numel(iPass1);

    % Criteria 2:
    % 95ci range must be <= 0.2 m/s
    iPass2 = find(vc150.vC(i).cispan <= 0.2);
    pass2{i} = iPass2;
    totPass2(i) = numel(iPass2);

    % Criteria 3:
    % I_range must be > 40
    iPass3 = find(vc150.vC(i).QCspan > 40);
    pass3{i} = iPass3;
    totPass3(i) = numel(iPass3);
    allPass{i} = intersect(intersect(iPass1, iPass2), iPass3);

    vc150.mV_c03(i) = mean(vc150.vC(i).meanV(allPass{i}));
end
