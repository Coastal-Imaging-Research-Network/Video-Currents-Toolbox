%% My own version of the vBar algorithm "videoCurrentsDemo.m"
%   where the function is defined:
% videoCurrentOut = videoCurrentGen(stack, time, xy, vB, fkB, Twin, Tstep, plotFlag);
%   & the inputs are as follows:
% stack: data [
% time: time vector (matlab datenum) with size [Mx1s]
% xy: xyzAll [Nx3] & x-search locations
% vB: (vBounds) [minV maxV], units m/s or vector of desired velocity steps,
%       to set this to empty, [], to use defaults [-3 3]
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

% Define search folder paths
c = genpath("C:\Users\13emo3\OneDrive - Queen's University\PhD\Chapter 4 vBar"); addpath(c);
q = genpath("C:\Users\13emo3\OneDrive - Queen's University\PhD\Matlab Functions"); addpath(q);

%%
params.yLims = [0 1100]; % Maximum & minimum y-values you're looking at
params.searchX = 200; % (m)
params.delY = 0.2;
y = params.yLims(1):params.delY:params.yLims(2);
params.yCam = 565;

params.searchDate = datetime(2017, 10, 01, 16, 00, 00);
%fileSearchPath = ("/Users/eloraoades/Library/CloudStorage/OneDrive-Queen'sUniversity/PhD/Chapter 3 vBar/Chapter 3 Data/FTP Downloads/");
fileSearchPath = ("C:\Users\13emo3\OneDrive - Queen's University\PhD\Chapter 4 vBar\Chapter 4 Data\FTP Downloads");

[T, RAW, XYZ, CAM] = loadVbarRawFile(fileSearchPath, params.searchDate, params.searchX);

params.mtime = (T/(3600*24)+datenum(1970,1,1))';
dTime = datetime(params.mtime, 'ConvertFrom', 'datenum');

%% Part 2: "sorted"
% Sort data & discard unneeded data
% Get the sorted indices based on Y values
Y = XYZ(:,2);
[~, sorted_indices] = sort(Y);

% Reorder RAW and XYZ based on the sorted indices
sorted.RAW = RAW(:, sorted_indices);
sorted.XYZ = XYZ(sorted_indices, :);

if sorted.XYZ(1,2) < params.yLims(1) || sorted.XYZ(end,2) > params.yLims(2)
    iy1 = find(sorted.XYZ(:,2) >= params.yLims(1), 1, 'first');
    iy2 = find(sorted.XYZ(:,2) <= params.yLims(2), 1, 'last');
    sorted.XYZ = sorted.XYZ(iy1:iy2,:);
    sorted.RAW = sorted.RAW(:,iy1:iy2);
end

%clearvars XYZ RAW iy1 iy2 sorted_indices
% Part 3: "unik"
% Remove Redundant x-y Locations
% Find unique XY locations and average corresponding RAW data
[~, uInd, ~] = unique(sorted.XYZ(:, 1:2), 'rows', 'stable');
uXY = sorted.XYZ(uInd, 1:2);
uRAW = nan(length(T), length(uXY));

for iX = 1:length(uXY)
    indices = ismember(sorted.XYZ(:, 1:2), uXY(iX, :), 'rows');
    uRAW(:, iX) = mean(sorted.RAW(:, indices), 2);
end

unik.RAW = uRAW;
unik.XY = uXY;

%clearvars uInd uXY uRAW indices
% Part 4: "gridded"
% Project Input Variables onto Regular Grids
% Define the search domain
% Create grid space to project RAW onto
[tq, yq] = meshgrid(T, y);
gridded.stack = griddata(unik.XY(:, 2), T, unik.RAW, yq, tq);
gridded.xy(:,2) = y;
gridded.xy(:,1) = params.searchX;

bad = isnan(gridded.stack);
if sum(bad, 'all') > 0
    gridded.stackOg = gridded.stack; % Save original as a different name
    gridded.stack = fillmissing(gridded.stack, 'nearest'); % Fill NaNs with mean
end

%clearvars unik bad
%% Part 3: Apply Radon
% Define the range of angles for the Radon transform
theta = 0:179;
Iz = gridded.stack;

% Perform the Radon transform
[R, xp] = radon(Iz, theta);

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
    params.vB = [0.01 2];
elseif dominant_angle > 5 && dominant_angle < 85
    disp('The foam is drifting south');
    params.vB = [-2 -0.01];
else
    disp('The dominant angle is within the excluded range.');
end

%% Define videoCurrentGen input parameters

% Set Twin: the time length of the FFT window (in points)
% For 2 Hz data, Twin = 128 will yield a 64s average current
params.Twin = 128;

% Set Tstep: time length to step the window (in points)
% For 2 Hz data, Tstep = 64 will yield a current estimate every 32 s
params.Tstep = 64;

% Set fkBounds = [fmin fmax kmin kmax], vector of frequency and wavenumber
%      bounds energy out of side of these bounds will be set to 0.  Useful
%      to eliminate some of the wave contamination that leaks in.  Set this
%      to empty, [], to use defaults.

params.vB = [-2 0]; 
params.fkB = [];

%params.vB = [-1.5 -0.075]; 
%params.fkB = [0.015 0.5 0.02 0.5];
params.plotFlag = 0;

% start/end indices for tile sections with varying sizes
% Define start/end indices for tile sections with varying sizes
yLength = y(end) - y(1);
tileCenters = y(1):10:y(end); % Define tile centers with a step of 10 meters

% Initialize the video currents table
vc150 = table();
count = 1;

for iCenter = 1:length(tileCenters)
    % Calculate the distance from the camera location
    %distanceFromCamera = ceil(abs(tileCenters(iCenter) - params.yCam)/10)*10;
    currentTileSize = 20; 
    % Determine the appropriate tile size
    %currentTileSize = round(distanceFromCamera/10, -1);
    %if currentTileSize < 10
    %    currentTileSize = 10;
    %end
    % Calculate the start and end indices for the current tile
    centerY = tileCenters(iCenter);
    startY = centerY - currentTileSize / 2;
    endY = centerY + currentTileSize / 2;

    % Find the corresponding indices in the y array
    i1 = find(y >= startY, 1);
    i2 = find(y <= endY, 1, 'last');

    % Ensure indices are within bounds
    if isempty(i1) || isempty(i2)
        continue;
    end
    stack = gridded.stack(i1:i2,:)';

    % Run video-current-toolbox
    vC = videoCurrentGen(stack, params.mtime, gridded.xy(i1:i2,:), params.vB, params.fkB, params.Twin, params.Tstep, params.plotFlag);
    vC.ci1 = vC.ci(:, 1);
    vC.ci2 = vC.ci(:, 2);
    vC = rmfield(vC, "ci");
    %vC.dTime = datetime(vC.t, 'ConvertFrom', 'datenum');

    % Calculate mean and standard deviation of meanV
    vC.mean_meanV = mean(vC.meanV);
    vC.std_meanV = std(vC.meanV);

    % Identify indices where meanV values are more than 3 standard deviations from the mean
    vC.outlierIndices = abs(vC.meanV - vC.mean_meanV) > 2 * vC.std_meanV;

    % List of all fields in the structure
    %fields = fieldnames(vC);
    % Remove data at these indices for all fields
    %for i = 1:numel(fields)
    %    vC.(fields{i})(outlierIndices) = [];
    %end

    % Save vC to the table
    vc150.vC(count) = vC;
    vc150.y(count) = centerY;
    vc150.tileSize(count) = currentTileSize;

    count = count + 1;
end

%vc150 = sortrows(vc150, 2); % Sort by y-coordinate    %%TAKE OUT

% clearvars -except vc175 gridded
% Weighted meanV
for i = 1:height(vc150)
    try
        vc150.wV(i) = wmean(vc150.vC(i).meanV, 1./vc150.vC(i).stdV, 'omitnan');
        vc150.snr(i) = wmean(vc150.vC(i).meanV, abs(vc150.vC(i).SNR), 'omitnan');
        vc150.mV(i) = mean(vc150.vC(i).meanV, "omitnan"); 
    catch
    end
end

% Post-Alg Filtering
% From paragraph [47] in Chickadel et al., 2003

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
