# Video-Currents-Toolbox
This repository contains code and documentation to quantify longshore surf-zone currents from video imagery.

The code is written in MATLAB, and requires version *****

Based on the original work of: Chickadel, C.C., Holman, R.A. and Freilich, M.H., 2003. An optical technique for the measurement of longshore currents. Journal of Geophysical Research: Oceans, 108(C11).

Toolbox Input: This code has been written for vbar pixel arrays with naming convention: 1506873540.Sun.Oct.01_15_59_00.GMT.2017.argus02b.cx.vbar125.mat, for example. Input file data should contain (1) pixel data (2) time (3) camera numbers (4) xyz position. 
Toolbox Output: This code will generate a table for each vbar transect listed in params.transects, containing 'x', 'y', 'vC', and 'wV'. 'vC' is a  structure variable that contains the timeseries analysis variables generated for the window, and 'wV' is the representative weighted velocity over the timeseries.  
More detail on data requirements and output structure is provided in the user manual document. 

The "videoCurrentsDemoNew.m" code will run you through the program starting from parameter selection ("vidCurrentsParams.m" file). 

If longshore current direction is unknown, params.vBounds = [], and a direction will be interpreted from the data using the radonVbarDir.m function.

"tcolor.m", "wmean.m", "redblue.m" are support functions. 