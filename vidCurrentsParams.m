%% Define videoCurrentGen input parameters

% Define input parameters <-- Make this a separate file 
%       (This is currently all in FRF local coordinates)
params.yLims = [0 1100];                                                     % Maximum & minimum y-values you're looking at
params.transects = [125, 150, 175, 200, 225];                                % x-transects
params.delY = 0.2;                                                           % ground pixel spacing in y <-- ** RELEVANT FOR VB
y = params.yLims(1):params.delY:params.yLims(2);                             % define the y-grid
params.yCam = 565;                                                           % y-location of the camera (m)
params.tileSize = 20; % m 
params.searchDate = datetime(2017, 10, 01, 15, 59, 00);                       % This date corresponds to the demo data

                    % *** NEW VARIABLE -- image sampling frequency ***
                    %       So that params.tWindow * params.fSample = number of points in the window
                    %       but tWindow & tStep can be defined in seconds 
params.fSample = 2; % Hz 

% NOW IN SECONDS 
params.tWindow = 64; % s                                                        % The time length of the FFT window 
params.tStep = 32;   % s                                                        % Temporal resolution output (i.e. 32 = 1 estimate every 32 seconds) 

params.vBounds = [];  % m/s                                                 % if empty, vBounds are defined by the Radon step                                                
params.radonCamNum = 1;                                                     % select camera to use y bounds for radon filter
params.fkBounds = [0.015 0.5 0.02 0.5];                                     % [fmin fmax kmin kmax], vector of search frequency and wavenumber windows

params.plotFlag = 0;


% NOTES: 
% vBounds min should be just greater than 0, to eliminate stationary objects
% from analysis
% 
% vBounds max is determined by the dy & dt spacing: 
%   if the spatial resolution is 0.2 m, 
%   f_N(spatial) = (1/2)*(1/0.2m) = 2.5 cycles/m
%   if the temporal resolution is 2 Hz, 
%   f_N(temporal) = (1/2)*2 Hz = 1 Hz 
% therefore, v_max = 2.5 cycles/m * 1 Hz = 2.5 m/s 