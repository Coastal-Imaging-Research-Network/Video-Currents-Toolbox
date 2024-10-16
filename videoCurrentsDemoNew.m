% Description: 



%% Inputs

% First, let's customize the input parameters file
edit vidCurrentsParams

% And load them 
run vidCurrentsParams.m

% If you'd like, you can use the search function to find the demo data 
% This comes in handy when you have folders with oodles of files 
fileSearchPath = ("D:\Elora PhD\GitHub\Video-Currents-Toolbox\DemoData");
sampleStack = loadVbarRawFile(fileSearchPath, params.searchDate, params.searchX);

% Or, uncomment the following line to load example stack structure directly 
% sampleStack = load('1506873540.Sun.Oct.01_15_59_00.GMT.2017.argus02b.cx.vbar125.mat');

% Time is initially defined as epoch, so let's convert it to datetime for 
% easier figure interpretation
params.mtime = (sampleStack.T/(3600*24)+datenum(1970,1,1))';
params.dTime = datetime(params.mtime, 'ConvertFrom', 'datenum');

% Extract the number of cameras in your data 
params.numCams = max(sampleStack.CAM, [], 'all');

% Sort the camera data & prep it for input
[inpDat] = prepDataForInput(sampleStack,params);  % <-- This line might take awhile! be patient :) 

% If you aren't sure which direction the current is heading (north/ south),
% leave the params.vBounds empty & use radonVbarDir: 
%       when plotFlag == 1, the following function will also output a figure that
%       you can verify & visualize the direction of foam propagation
plotFlag = 1; 
if isempty(params.vBounds)
    [params] = radonVbarDir(inpDat, params, plotFlag); 
end

% Run the videoCurrentGen code
[vcTable] = vcTableGen(inpDat, params);


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
