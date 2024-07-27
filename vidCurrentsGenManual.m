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
paramsFile = ("vidCurrentsParams");
eval(paramsFile);

% fileSearchPath is the folder within which the function vBarRawFile will
% search for the "raw" data file
fileSearchPath = ("C:\Users\13emo3\OneDrive - Queen's University\PhD\Chapter 4 vBar\Chapter 4 Data\FTP Downloads");

[T, RAW, XYZ] = loadVbarRawFile(fileSearchPath, params.searchDate, params.searchX);

% time is loaded as epoch, so convert it to datetime for plotting
% interpretation
mtime = (T/(3600*24)+datenum(1970,1,1))';
dTime = datetime(mtime, 'ConvertFrom', 'datenum');

%% Part 2: Sort & Grid

vidCurrTable = table(); 

%create "list" of cameras.
for i = 1:max(double(CAM))
    %For each camera, find extent of XY and RAW. Store variable
    fieldNameI = sprintf('cam%1.0d', i);
    camList.(fieldNameI).XY = XYZ(CAM == i, 1:2);
    camList.(fieldNameI).RAW = RAW(:,CAM == i);
    % temp vars for iteration of loop
    tempXY = camList.(fieldNameI).XY;
    tempRAW = camList.(fieldNameI).RAW;

    minY = round(min(tempXY(:,2)));
    maxY = round(max(tempXY(:,2)));

    % % Use Y limits of camera (and T on x-axis) to create grid
    % yBounds = minY : params.delY : maxY;
    % [tGrid, yGrid] = meshgrid(T, yBounds);
    % 
    % % Define a new variable called "inputData" which stores the gridded
    % % data that will be used as input for "videoCurrentGen"
    % inputDat.(fieldNameI).tGrid = tGrid;
    % inputDat.(fieldNameI).yGrid = yGrid;
    % 
    % % Project 'RAW' data for this camera onto the grid
    % inputDat.(fieldNameI).rawGrid = griddata(tempXY(:,2), T, double(tempRAW), yGrid, tGrid, 'natural');
    
    % Define the center points of the analysis windows
    inputDat.yCentres = minY+(params.tileSize/2):params.tileSize:(maxY-params.tileSize/2); 

    
    % % Apply Radon to define the velocity search bounds (vBounds) to + or -
    % theta = 0:179;  % the range of angles we're considering
    % % Note that given the data configuration, theta represents the angle
    % % *from which* energy is coming
    % 
    % % if vBounds are undefined (in params file),
    % if isempty(params.vBounds)
    %     % Perform Radon on all angles between 0 to 180 (theta),
    %     [R, xp] = radon(interpRAW, theta);
    % 
    %     % Exclude energy from unwanted angle ranges (where ~90 deg = shore
    %     % normal (stationary objects) & ~180 deg = shore parallel (inf speed)
    %     excluded_angles = (theta >= 85 & theta <= 95) | (theta >= 175 | theta <= 5);
    %     R(:, excluded_angles) = 0;
    % 
    %     % The angle with the maximum value from the Radon transform is
    %     % associated with the direction of the foam drift
    %     [max_val, max_idx] = max(R(:));
    %     [max_row, max_col] = ind2sub(size(R), max_idx);
    %     dominant_angle = theta(max_col);
    % 
    %     % Interpret the dominant angle:
    %     if dominant_angle > 95 && dominant_angle < 175
    %         disp('The foam is drifting north');
    %         params.vBounds = [0.01 2.5];                                 % <--- *see note in "vidCurrentsParams.m"
    %     elseif adominant_angle > 5 && dominant_angle < 85
    %         disp('The foam is drifting south');
    %         params.vBounds = [-2.5 -0.01];
    %     else
    %         disp('The dominant angle is within the excluded range.');
    %     end
end


%% Part 2: Run the main function 

% Step 1: Initialize the table with y-values and cam-values 
vc150 = table();

for iCenter = 1:length(tileCenters)

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
