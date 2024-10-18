% Description:



%% Inputs

% First, let's customize the input parameters file
edit vidCurrentsParams

% And load them
run vidCurrentsParams.m

% If you'd like, you can use the search function to find the demo data
% This comes in handy when you have folders with oodles of files
fileSearchPath = ("D:\Elora PhD\GitHub\Video-Currents-Toolbox\DemoData");
fileSearchPath = 'D:\Argus Downloads';
[sampleStack, timex] = loadVbarRawFile(fileSearchPath, params.searchDate, params.transects);

% Or, uncomment the following line to load example stack structure directly
% sampleStack = load('1506873540.Sun.Oct.01_15_59_00.GMT.2017.argus02b.cx.vbar125.mat');

% Time is initially defined as epoch, so let's convert it to datetime for
% easier figure interpretation
params.mtime = (sampleStack.T/(3600*24)+datenum(1970,1,1))';
params.dTime = datetime(params.mtime, 'ConvertFrom', 'datenum');

% Extract the number of cameras in your data
params.numCams = max(sampleStack.CAM, [], 'all');

% Sort the camera data & prep it for input
for i = 1:length(params.transect)
    [inpDat.(params.transects(i))] = prepDataForInput(sampleStack{i}, params);  % <-- This line might take awhile! be patient :)
end 

% If you aren't sure which direction the current is heading (north/ south),
% leave the params.vBounds empty & use radonVbarDir:
%       when plotFlag == 1, the following function will also output a figure that
%       you can verify & visualize the direction of foam propagation
plotFlag = 1;
if isempty(params.vBounds)
    % select y bounds for radon to test
    radonCam = sprintf('cam%1.0d', params.radonCamNum);
    [params] = radonVbarDir(inpDat.(radonCam), params, plotFlag);
end

% Run the videoCurrentGen code
[vcTable] = vcTableGen(inpDat, params);


%% Plot the Data

figure();
tcolor(timex.x, timex.y, timex.Ip, 'corners'); shading flat; hold on;
axis tight equal; set(gca, 'Layer', 'top', 'FontName', 'Cambria', 'FontSize', 14, 'box', 'on');
scatter(vcTable.x, vcTable.y, 20, vcTable.wV, 'o', 'filled', 'MarkerEdgeColor', 'k');
xlabel('x (m)'); ylabel('y (m)');
ylim([params.yLims(1) params.yLims(2)])
c = colorbar();
colormap(redblue);

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
