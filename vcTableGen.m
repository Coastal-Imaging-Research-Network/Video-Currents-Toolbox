function [vcTable] = vcTableGen(inputDataX, params, xLocation)

% Step 1: Initialize the table with y-values and cam-values
%vcTable = table([], struct([]), [], 'VariableNames', {'y', 'vC', 'wV'});
vcTable = table(); 
count = 1; 

for j = 1:params.numCams
    fieldNameJ = sprintf('cam%1.0d', j);
    
    minY = min(inputDataX.(fieldNameJ).yGrid, [],  'all'); 
    maxY = max(inputDataX.(fieldNameJ).yGrid, [], 'all'); 

    % Define the center points of the analysis windows
    inputDataX.(fieldNameJ).yCenters = (minY+params.tileSize/2):params.tileSize:(maxY-params.tileSize/2);

    for k = 1:length(inputDataX.(fieldNameJ).yCenters)
        y = inputDataX.(fieldNameJ).yGrid(:,1); 
        
        y1 = inputDataX.(fieldNameJ).yCenters(k) - params.tileSize/2;    % start point 
        y2 = inputDataX.(fieldNameJ).yCenters(k) + params.tileSize/2;    % end point
        
        i1 = find(y == y1, 1, 'first'); 
        i2 = find(y == y2, 1, 'first'); 

        % Run video-current-toolbox
        % stack, time, xy, vBounds, fkBounds, Twin, Tstep {plotFlag})
        stack = inputDataX.(fieldNameJ).rawGrid(i1:i2,:)'; 
        xy = inputDataX.(fieldNameJ).yGrid(i1:i2,1); 
        vC = videoCurrentGen(stack, params.mtime, xy, ...
                params.vBounds, params.fkBounds, params.tWindow, params.tStep, params.plotFlag);

        % Save vC to the table
        vcTable.x(count) = xLocation; 
        vcTable.y(count) = inputDataX.(fieldNameJ).yCenters(k);       % midpoint
        vcTable.vC{count} = vC; 
        indexKeep = find((vcTable.vC{count}.SNR > 5)&(vcTable.vC{count}.prob>0.8)&(vcTable.vC{count}.stdI>5));
        vcTable.wV(count) = wmean(vC.meanV(indexKeep), 1./vC.stdV(indexKeep), 'omitnan');
        vcTable.camNum(count) = j; 
        count = count+1; 
    end
end

vcTable = sortrows(vcTable, 'y');   
end