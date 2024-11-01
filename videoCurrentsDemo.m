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
params.mtime = (sampleStack{1}.T/(3600*24)+datenum(1970,1,1))';
params.dTime = datetime(params.mtime, 'ConvertFrom', 'datenum');

% Extract the number of cameras in your data
params.numCams = max(sampleStack{1}.CAM, [], 'all');

% Sort the camera data & prep it for input
for i = 1:length(params.transects)
    fieldNameI = sprintf('x%1.0d', params.transects(i));
    [inpDat.(fieldNameI)] = prepDataForInput(sampleStack{i}, params);  % <-- This line might take awhile! be patient :)
end 

% If you aren't sure which direction the current is heading (north/ south),
% leave the params.vBounds empty & use radonVbarDir:
%       when plotFlag == 1, the following function will also output a figure that
%       you can verify & visualize the direction of foam propagation
plotFlag = 1;
if isempty(params.vBounds)
    % select y bounds for radon to test
    radonCam = sprintf('cam%1.0d', params.radonCamNum);
    [params] = radonVbarDir(inpDat.(fieldNameI).(radonCam), params, plotFlag);
end

% Run the videoCurrentGen code
[vcTable125] = vcTableGen(inpDat.x125, params, params.transects(1));
[vcTable150] = vcTableGen(inpDat.x150, params, params.transects(2)); 
[vcTable175] = vcTableGen(inpDat.x175, params, params.transects(3)); 
[vcTable200] = vcTableGen(inpDat.x200, params, params.transects(4)); 
[vcTable225] = vcTableGen(inpDat.x225, params, params.transects(5)); 

%% Plot the Data

figure();
tcolor(timex.x, timex.y, timex.Ip, 'corners'); shading flat; hold on;
axis tight equal; set(gca, 'Layer', 'top', 'FontName', 'Cambria', 'FontSize', 14, 'box', 'on');
scatter(vcTable125.x, vcTable125.y, 20, vcTable125.wV, 'o', 'filled', 'MarkerEdgeColor', 'k');
scatter(vcTable150.x, vcTable150.y, 20, vcTable150.wV, 'o', 'filled', 'MarkerEdgeColor', 'k');
scatter(vcTable175.x, vcTable175.y, 20, vcTable175.wV, 'o', 'filled', 'MarkerEdgeColor', 'k');
scatter(vcTable200.x, vcTable200.y, 20, vcTable200.wV, 'o', 'filled', 'MarkerEdgeColor', 'k');
scatter(vcTable225.x, vcTable225.y, 20, vcTable225.wV, 'o', 'filled', 'MarkerEdgeColor', 'k');

xlabel('x (m)'); ylabel('y (m)');
ylim([params.yLims(1) params.yLims(2)])
c = colorbar();
colormap(parula);

%% 
figure(); 
T = tiledlayout(1, 5); 

nexttile(); 
scatter(vcTable125.wV, vcTable125.y, 'o', 'filled'); 
xlabel('v_y (m/s)'); ylabel('y-position (m)');
title(sprintf('x = %d', vcTable125.x(1))); 
xlim([floor(min(vcTable125.wV)) ceil(max(vcTable125.wV))]);
ylim([min(params.yLims) max(params.yLims)])
set(gca, 'FontName', 'Cambria', 'FontSize', 12, 'box', 'on');

nexttile(); 
scatter(vcTable150.wV, vcTable150.y, 'o', 'filled'); 
xlabel('v_y (m/s)');
title(sprintf('x = %d', vcTable150.x(1))); 
xlim([floor(min(vcTable150.wV)) ceil(max(vcTable150.wV))]);
ylim([min(params.yLims) max(params.yLims)])
set(gca, 'FontName', 'Cambria', 'FontSize', 12, 'box', 'on');

nexttile(); 
scatter(vcTable175.wV, vcTable175.y, 'o', 'filled'); 
xlabel('v_y (m/s)');
title(sprintf('x = %d', vcTable175.x(1))); 
xlim([floor(min(vcTable175.wV)) ceil(max(vcTable175.wV))]);
ylim([min(params.yLims) max(params.yLims)])
set(gca, 'FontName', 'Cambria', 'FontSize', 12, 'box', 'on');

nexttile(); 
scatter(vcTable200.wV, vcTable200.y, 'o', 'filled'); 
xlabel('v_y (m/s)');
title(sprintf('x = %d', vcTable200.x(1))); 
xlim([floor(min(vcTable200.wV)) ceil(max(vcTable200.wV))]);
ylim([min(params.yLims) max(params.yLims)])
set(gca, 'FontName', 'Cambria', 'FontSize', 12, 'box', 'on');

nexttile(); 
scatter(vcTable225.wV, vcTable225.y, 'o', 'filled'); 
xlabel('v_y (m/s)');
title(sprintf('x = %d', vcTable225.x(1))); 
xlim([floor(min(vcTable225.wV)) ceil(max(vcTable225.wV))]);
ylim([min(params.yLims) max(params.yLims)])
set(gca, 'FontName', 'Cambria', 'FontSize', 12, 'box', 'on');
