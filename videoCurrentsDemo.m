% Description:
% This script applies the videoCurrentGen code to compute optical currents
% for an example video-currents stack structure.
%
% Here, the expected fields of the stack structure are:
% inst: 
%        type: 'line'
%        xyz: [x0 y0 z; x1 y1 z], where (x0,y0,z) and (x1,y1,z) are the
%              endpoints of pixel line, e.g. a 20-m alongshore line:
%              xyz: [115 640 0; 115 660 0], size [3x2]
%        name: 'vBar115'
% dn: time vector (matlab datenum) with size [Mx1s]
% xyzAll: vector with x,y,z coordinates of points along pixel line
%        instrument, size [Nx3]
% data: grayscale image of timestack for one linear pixel instrument,
%        with size [MxN]
% Note on dimensions:
%        M is the length of the whole record, e.g., for a 15 min video,
%           subsampled at 2 Hz, M~1801
%           The timestack will be analyzed in windows Twin (e.g., 64 s)
%           for the video-current estimation in videoCurrentGen
%        N is the length of pixel array given by points (x,y,z), e.g.,
%           for a 20 m pixel array with .2 m spacing, N~101
% Note on instruments:
%        If there are multiple pixel instruments, adapt this script to loop
%        through stackstruct(1), stackstruct(2), etc.
%
% (This should be the same form as the Aerielle video demo output from the
%  2017 bootcamp.)

%% Inputs

% Load example stack structure
stackStruct = load('CIRN_example_stackstruct');

% Grab stack: grayscale image of timestack, size [MxN]
stack = stackStruct.data;

% Grab time: time line (starting from zero) of stack, size [1xM]
time = stackStruct.dn;

% Grab xy: x,y position [Nx2] of each pixel in a dimension (equally spaced)
xyz = stackStruct.xyzAll;
xy = xyz(:,1:2);

% Set Twin: the time length of the FFT window (in points)
% For 2 Hz data, Twin = 128 will yield a 64s average current
Twin = 128; 

% Set Tstep: time length to step the window (in points)
% For 2 Hz data, Tstep = 64 will yield a current estimate every 32 s
Tstep = 64; 

% Set vBounds: [minV maxV], units m/s or vector of desired velocity steps,
% Set this to empty, [], to use defaults [-3 3]
vB = [];

% Set fkBounds = [fmin fmax kmin kmax], vector of frequency and wavenumber
%      bounds energy out of side of these bounds will be set to 0.  Useful
%      to eliminate some of the wave contamination that leaks in.  Set this
%      to empty, [], to use defaults.
fkB = [];

% Set {plotFlag}: optional, if true (~=0) will display a running plot of
%      the data processing
plotFlag = 1;

%% Run the videoCurrentGen code

% Run code
videoCurrentOut = videoCurrentGen(stack, time, xy, vB, fkB, ...
    Twin, Tstep, plotFlag);

%  OUTPUT fields in videoCurrentOut returned:
%    meanV - video current estimate of mean current for each time window
%    t - time index for meanV
%    ci - the 95% conf. interval around meanV
%    cispan -  the width of ci
%    prob - the probability of the model fit
%    QCspan - the 95th percentile minus the 50th percentile of the timestack
%             histogram, used to measure the amount of video "texture"
%    stdV - the width (std. dev.) of the energy in velocity spectrum
%    vAngle - orientation of the pixel array (radians)
