function [vcTable] = vcTableGen(inputDat, params)

% Step 1: Initialize the table with y-values and cam-values
%vcTable = table([], struct([]), [], 'VariableNames', {'y', 'vC', 'wV'});
vcTable = table(); 
count = 1; 

for j = 1:params.numCams
    fieldNameJ = sprintf('cam%1.0d', j);
    for k = 1:length(inputDat.(fieldNameJ).yCentres)
        y = inputDat.(fieldNameJ).yGrid(:,1); 
        
        y1 = inputDat.(fieldNameJ).yCentres(k) - params.tileSize/2;    % start point 
        y2 = inputDat.(fieldNameJ).yCentres(k) + params.tileSize/2;    % end point
        
        i1 = find(y == y1, 1, 'first'); 
        i2 = find(y == y2, 1, 'first'); 

        % Ensure indices are within y-search bounds
        if isempty(i1) || isempty(i2)
            continue;
        end

        % Run video-current-toolbox
        % stack, time, xy, vBounds, fkBounds, Twin, Tstep {plotFlag})
        stack = inputDat.(fieldNameJ).rawGrid(i1:i2,:)'; 
        xy = inputDat.(fieldNameJ).yGrid(i1:i2,1); 
        vC = videoCurrentGen(stack, params.mtime, xy, ...
                params.vBounds, params.fkBounds, params.tWindow, params.tStep, params.plotFlag);

        % Save vC to the table
        vcTable.y{count} = inputDat.(fieldNameJ).yCentres(k);       % midpoint
        vcTable.vC{count} = vC; 
        vcTable.wV{count} = wmean(vC.meanV, 1./vC.stdV, 'omitnan');
        %newRow = {y, vC, wV};
        %vcTable(2,:) = [vcTable; newRow];
        count = count+1; 
    end
end

vcTable = sortrows(vcTable, 1);
   
end