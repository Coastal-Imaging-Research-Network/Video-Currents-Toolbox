function [params] = radonVbarDir(inpDat, params, plotFlag)
%% Apply Radon to define the velocity search bounds (vBounds) to + or -
theta = 0:179;  % the range of angles we're considering
% Note that given the data configuration, theta represents the angle
% *from which* energy is coming

% Check if plotflag is provided, if not, set it to 0 (no plotting)
if nargin < 2
    plotFlag = 0;
end

%select y bounds for radon to test
radonCam = sprintf('cam%1.0d', params.radonCamNum);

plotTimeGrid = datetime(inpDat.(radonCam).tGrid(1,:)/(3600.0*24) + datenum(1970,1,1), 'ConvertFrom', 'datenum');

% Perform Radon on all angles between 0 to 180 (theta),
[R, xp] = radon(inpDat.(radonCam).rawGrid, theta);

% Exclude energy from unwanted angle ranges (where ~90 deg = shore
% normal (stationary objects) & ~180 deg = shore parallel (inf speed)
excluded_angles = (theta >= 85 & theta <= 95) | (theta >= 175 | theta <= 5);
R(:, excluded_angles) = 0;

% The angle with the maximum value from the Radon transform is
% associated with the direction of the foam drift
[~, max_idx] = max(R(:));
[~, max_col] = ind2sub(size(R), max_idx);
dominant_angle = theta(max_col);

imshow(R,[],'Xdata',theta,'Ydata',xp,'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(gca,hot), colorbar

% Interpret the dominant angle:
if dominant_angle > 95 && dominant_angle < 175
    disp('The foam is drifting south');
    params.vBounds = [-2.5 -0.01];                                 % <--- *see note in "vidCurrentsParams.m"
elseif dominant_angle > 5 && dominant_angle < 85
    disp('The foam is drifting north');
    params.vBounds = [0.01 2.5];
else
    disp('The dominant angle is within the excluded range.');
end

% Plot check if plotflag is set to 1
if plotFlag == 1    
    figure();
    pcolor(plotTimeGrid, inpDat.(radonCam).yGrid, inpDat.(radonCam).rawGrid);
    title(sprintf('%s', plotTimeGrid(1))); 
    shading flat;
    colormap(gray);
    xlabel('Time');
    ylabel('y-position (m)');
    title('Figure Check');
end

end