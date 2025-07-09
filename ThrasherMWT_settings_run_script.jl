#Settings

#path to Directory containing MWT Output Directories of the experiment
expt_dir = raw"Example\MWT_Output";

#path to Directory where you want to save the analysis output
analysis_dir = raw"Example\analysis";

wormsperwell = 4;

bin_duration = 60;#seconds

negativebaseline = false;

normalize_wl = false;

MWT_output_timestamp_offset = true;

roundpixelsize = false;

sdp_per_worm = false;

#ROIs for wells
    #all x coordinates are meassured at the left edge of the well
    #all y coordinates are meassured at the top edge of the well

    #number of rows
        nrows = 5;
        
    #number of columns
        ncolumns = 7;
    
    #diameter of well
        d = 510;

    #top left well
        xtl = 110;
        ytl = 190;

    #top right well
        xtr = 3242;
        ytr = 228;

    #bottom left well
        xbl = 82;
        ybl = 2276;

#################################################
#DO NOT CHANGE BELOW THIS LINE
#################################################

import ThrasherMWT
if length(ARGS) > 0
    @time ThrasherMWT.old(ARGS[1]);
else
    oldscript = false;
    @time ThrasherMWT.main(oldscript, expt_dir, analysis_dir, wormsperwell, bin_duration, negativebaseline, normalize_wl, MWT_output_timestamp_offset, roundpixelsize, sdp_per_worm, nrows, ncolumns, d, xtl, ytl, xtr, ytr, xbl, ybl);
end
