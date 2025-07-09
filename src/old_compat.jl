function old(old_script_path)
    oldscript = true;
    s = read(old_script_path,String);
    
    starti = findnext("\"",s,findnext("parent_dir",s,1).stop).stop+1;
    stopi = findnext("\";",s,starti).start-1;
    parent_splitdrive = splitdrive(s[starti:stopi]);
    parent_dir = parent_splitdrive[1] * join(split(parent_splitdrive[2],"\\\\"),"\\");
    
    starti = findnext("\"",s,findnext("experiment_name",s,1).stop).stop+1;
    stopi = findnext("\";",s,starti).start-1;
    experiment_name = s[starti:stopi];
    
    expt_dir = joinpath(parent_dir,"1_data",experiment_name);
    analysis_dir = joinpath(parent_dir,"2_analysis",experiment_name);
    
    starti = findnext("=",s,findnext("wormsperwell",s,1).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    wormsperwell = parse(Int64,s[starti:stopi]);
    
    bin_duration = 60;#seconds
    
    negativebaseline = false;
    
    starti = findnext("=",s,findnext("normalize_wl",s,1).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    normalize_wl = parse(Bool,s[starti:stopi]);
    
    MWT_output_timestamp_offset = true;
    
    roundpixelsize = true;
    
    sdp_per_worm = true;
    
    nrows = 4;
    
    ncolumns = 5;
    
    dowi = findnext("#diameter of well",s,1).stop;
    starti = findnext("=",s,findnext("d",s,dowi).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    d = parse(Float64,s[starti:stopi]);
    
    starti = findnext("=",s,findnext("xtl",s,dowi).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    xtl = parse(Float64,s[starti:stopi]);
    
    starti = findnext("=",s,findnext("ytl",s,dowi).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    ytl = parse(Float64,s[starti:stopi]);
    
    starti = findnext("=",s,findnext("xtr",s,dowi).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    xtr = parse(Float64,s[starti:stopi]);
    
    starti = findnext("=",s,findnext("ytr",s,dowi).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    ytr = parse(Float64,s[starti:stopi]);
    
    starti = findnext("=",s,findnext("xbl",s,dowi).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    xbl = parse(Float64,s[starti:stopi]);
    
    starti = findnext("=",s,findnext("ybl",s,dowi).stop).stop+1;
    stopi = findnext(";",s,starti).start-1;
    ybl = parse(Float64,s[starti:stopi]);
    
    ThrasherMWT.main(oldscript, expt_dir, analysis_dir, wormsperwell, bin_duration, negativebaseline, normalize_wl, MWT_output_timestamp_offset, roundpixelsize, sdp_per_worm, nrows, ncolumns, d, xtl, ytl, xtr, ytr, xbl, ybl);
    return nothing
end
