function [T, RAW, XYZ, CAM] = loadVbarRawFile(fileSearchPath, searchDate, searchX)
% masterFolder: the master directory to search
    % targetDate: the date to match (datetime format)
    % xLocation: the x-location to match (double)

    % Convert the target date to the required format (e.g., 'Tue.Oct.17')
    targetDateStr = datestr(searchDate, 'ddd.mmm.dd');
    targetTimeStr1 = datestr(searchDate, 'HH_MM');  % e.g., '11_59'
    targetTimeStr2 = datestr(searchDate-minutes(1), 'HH_MM'); 

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

    % Filter files based on the provided date, time, x-location, and type 'vbar'
    filteredFiles = [];
    for i = 1:length(matFileNames)
        fileName = matFileNames(i).name;
        % Check if the file name matches the desired pattern
        if contains(fileName, targetDateStr) && contains(fileName, targetTimeStr1) && ...
           contains(fileName, ['vbar', num2str(searchX)]) && ...
           contains(fileName, '.argus02b.cx.vbar')
            filteredFiles = [filteredFiles; matFileNames(i)];
        elseif contains(fileName, targetDateStr) && contains(fileName, targetTimeStr2) && ... 
                contains(fileName, ['vbar', num2str(searchX)]) && ...
                contains(fileName, '.argus02b.cx.vbar')
            filteredFiles = [filteredFiles; matFileNames(i)];
        end
    end

    % Load the specified variables from the filtered files
    for i = 1:length(filteredFiles)
        filePath = fullfile(filteredFiles(i).folder, filteredFiles(i).name);
        disp(['Loading file: ', filePath]);
        data = load(filePath, 'XYZ', 'T', 'RAW', 'CAM');

        % Check if the variables exist and then process them as needed
        if isfield(data, 'XYZ') && isfield(data, 'T') && isfield(data, 'RAW') && isfield(data, 'CAM')
            XYZ = data.XYZ;
            T = data.T;
            RAW = data.RAW;
            CAM = data.CAM;
            % Process the variables as needed
            disp('Variables loaded: XYZ, T, RAW, CAM');
        else
            disp('One or more variables not found in the file.');
        end
    end
end
