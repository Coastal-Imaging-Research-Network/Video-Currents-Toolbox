function [inpDat] = prepDataForInput(sampleStack,params)

% Grab stack: grayscale image of timestack, size [MxN]
data = sampleStack.RAW;

% Grab time: time line (starting from zero) of stack, size [1xM]
time = sampleStack.T;

% Grab xy: x,y position [Nx2] of each pixel in a dimension (equally spaced)
xyz = sampleStack.XYZ;

cam = sampleStack.CAM; 

% create "list" of cameras.
for j = 1:max(double(sampleStack.CAM))
    %For each camera, find extent of XY and RAW. Store variable
    fieldNameI = sprintf('cam%1.0d', j);
    camList.(fieldNameI).XY = xyz(cam == j, 1:2);
    camList.(fieldNameI).RAW = data(:,cam == j);
    % temp vars for iteration of loop
    tempXY = camList.(fieldNameI).XY;
    tempRAW = camList.(fieldNameI).RAW;

    minY = max(round(min(tempXY(:,2))), min(params.yLims));
    maxY = min(round(max(tempXY(:,2))), max(params.yLims));

    % Use Y limits of camera (and T on x-axis) to create grid
    yBounds = minY : params.delY : maxY;
    [tGrid, yGrid] = meshgrid(time, yBounds);

    % Define a new variable called "inputData" which stores the gridded
    % data that will be used as input for "videoCurrentGen"
    inpDat.(fieldNameI).tGrid = tGrid;
    inpDat.(fieldNameI).yGrid = yGrid;

    % Project 'RAW' data for this camera onto the grid
    inpDat.(fieldNameI).rawGrid = griddata(tempXY(:,2), time, double(tempRAW), yGrid, tGrid, 'natural');

    % Define the center points of the analysis windows
    inpDat.(fieldNameI).yCentres = minY+(params.tileSize/2):params.tileSize:(maxY-params.tileSize/2);
end

end