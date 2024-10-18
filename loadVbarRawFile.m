function [sampleStack, timex] = loadVbarRawFile(fileSearchPath, searchDate, transects)
% Function to load vBar and timex files for specified transects on a specific date
%
% fileSearchPath: master directory to search
% searchDate: the date to match (in datetime format)
% transects: array of x-locations to match (e.g., [125, 150, 175, 200, 225])

% Convert the target date to the required format (e.g., 'Tue.Oct.17')
targetDateStr = datestr(searchDate, 'ddd.mmm.dd');
targetTimeStr1 = datestr(searchDate, 'HH_MM');  % e.g., '11_59'
targetTimeStr2 = datestr(searchDate+minutes(1), 'HH_MM');

% Generate the path for all subfolders
searchPath = genpath(fileSearchPath);

% Split the search path into individual folders
folders = strsplit(searchPath, pathsep);

% Initialize an empty array to hold the .mat file information
matFileNames = [];

% Loop through each folder and search for .mat files
for i = 1:length(folders)
    folder = folders{i};
    if ~isempty(folder)
        % Get the .mat files in the current folder
        files = dir(fullfile(folder, '*.mat'));
        % Append the files to the matFileNames array
        matFileNames = [matFileNames; files];
    end
end

% Initialize variables to store the files for vbar and timex
vbarFiles = cell(length(transects), 1);  % Initialize cell array for each transect
timexFiles = [];

for i = 1:length(matFileNames)
    fileName = matFileNames(i).name;
    if contains(fileName, targetDateStr) && contains(fileName, targetTimeStr2) && ...
            contains(fileName, '.timex.merge')
        timexFiles = [timexFiles; matFileNames(i)];
    end

    % Loop through the transects
    for j = 1:length(transects)
        transect = transects(j);
        transectStr = sprintf('.vbar%d', transect);  % Create the string for the current transect (e.g., .vbar125)
        % Check if the file name matches the desired pattern for vbar files
        if contains(fileName, targetDateStr) && contains(fileName, targetTimeStr1) && ...
                contains(fileName, 'vbar') && ...
                contains(fileName, transectStr)  % Check the current transect
            vbarFiles{j} = [vbarFiles{j}; matFileNames(i)];  % Store vbar files for each transect
        end
    end
end

% Load the specified variables from the filtered files

% There should only be one timex
if ~isempty(timexFiles)
    timex = load(fullfile(timexFiles.folder, timexFiles.name), 'x', 'y', 'Ip');
else
    error('No timex file found.');
end

% Initialize sampleStack for holding vBar data
sampleStack = cell(length(transects), 1);

% Loop through vBar files for each transect and load data into struct array
for j = 1:length(transects)
    for i = 1:length(vbarFiles{j})
        filePath = fullfile(vbarFiles{j}(i).folder, vbarFiles{j}(i).name);
        disp(['Loading file: ', filePath]);
        % Load the variables and store them in the sampleStack struct
        loadedData = load(filePath, 'XYZ', 'T', 'RAW', 'CAM');

        % Assign the loaded data to the sampleStack struct at the j-th position
        sampleStack{j}.XYZ = loadedData.XYZ;
        sampleStack{j}.T = loadedData.T;
        sampleStack{j}.RAW = loadedData.RAW;
        sampleStack{j}.CAM = loadedData.CAM;
    end
end

end