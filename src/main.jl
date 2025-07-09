function main(oldscript, expt_dir, analysis_dir, wormsperwell, bin_duration, negativebaseline, nwl, MWT_output_timestamp_offset, roundpixelsize, sdp_per_worm, nrows, ncolumns, d, xtl, ytl, xtr, ytr, xbl, ybl)
    #setup
    if !isdir(analysis_dir) mkpath(analysis_dir) end
    speed_output = joinpath(analysis_dir,"speed.txt");
    nworms_output = joinpath(analysis_dir,"nworms.txt");
    
    global normalize_wl = nwl;
    
    pixelsize = 9/d;#9mm well
    if roundpixelsize
        pixelsize = round(pixelsize,digits=4);
    end
    
    recordings = readdir(expt_dir, join=true);
    
    if (length(recordings) == 1) & negativebaseline
        negativebaseline == false;
        println("Only 1 recording, so adjusted negativebaseline = false");
    end
    
    println("\nLoading MWT Output Files ...")
    if !negativebaseline
        global clew = load_recording(recordings[1]);#baseline
        global all_ftime = get_all_ftime(recordings[1]);
        if length(recordings) > 1
            current_ftime_offset = 0;
            dt1 = getmwtdatetime(recordings[1]);
            for ri = 2:length(recordings)
                if MWT_output_timestamp_offset
                    current_ftime_offset = Dates.value(getmwtdatetime(recordings[ri]) - dt1)/1000;#seconds
                else
                    current_ftime_offset += get_rec_dur(recordings[ri-1]);
                end
                global clew = vcat(clew, load_recording(recordings[ri],ftime_offset = current_ftime_offset));
                global all_ftime = vcat(all_ftime, get_all_ftime(recordings[ri]) .+ current_ftime_offset);
            end
        end
    else
        if MWT_output_timestamp_offset
            dt1 = getmwtdatetime(recordings[1]);
            dt2 = getmwtdatetime(recordings[2]);
            baseline_ftime_offset = -Dates.value(dt2-dt1)/1000;#seconds
        else
            baseline_ftime_offset = -get_rec_dur(recordings[1]);
        end
        global clew = load_recording(recordings[1],ftime_offset = baseline_ftime_offset);#baseline
        global all_ftime = get_all_ftime(recordings[1]);
        global all_ftime = all_ftime .+ baseline_ftime_offset;
        global clew = vcat(clew, load_recording(recordings[2]));#first post-baseline recording
        global all_ftime = vcat(all_ftime, get_all_ftime(recordings[2]));
        if length(recordings) > 2
            current_ftime_offset = 0;
            for ri = 3:length(recordings)
                if MWT_output_timestamp_offset
                    current_ftime_offset = Dates.value(getmwtdatetime(recordings[ri]) - dt2)/1000;#seconds
                else
                    current_ftime_offset += get_rec_dur(recordings[ri-1]);
                end
                global clew = vcat(clew, load_recording(recordings[ri],ftime_offset = current_ftime_offset));
                global all_ftime = vcat(all_ftime, get_all_ftime(recordings[ri]) .+ current_ftime_offset);
            end
        end
    end
    println("Done!                                ");
    println("                                      ");
    println("Analyzing ...");
    
    ##
    
    #remove worms that have less than 2 frames
    ix = findall(map(w->length(w.fnum),clew) .>= 2);
    clew = clew[ix];
    
    #remove worms that are observed for less than 10 seconds
    ix = findall(map(w->w.ftime[end]-w.ftime[1],clew) .>= 10);
    clew = clew[ix];
    
    ##
    
    if oldscript
        old_setwells(nrows, ncolumns, d, xtl, ytl, xtr, ytr, xbl, ybl);
        
        clew = find_wells_and_sdp(clew);
        
        old_place_all();
    end
    
    if !oldscript
        setwells(nrows, ncolumns, d, xtl, ytl, xtr, ytr, xbl, ybl);
                
        if sdp_per_worm
            sdp_perworm!(clew);
        else
            sdp_perclew!(clew);
        end
        
        place_all();
    end
    
    println("  Calculating Speed ...");
    all_wormsperwell = ones(Int64,length(all_ftime))*wormsperwell;
    for wi in eachindex(wells)
        wells[wi] ./= maximum([nworms[wi] all_wormsperwell],dims=2);
        if wormsperwell == 0#and no worms observed in well at this timepoint, set to NaN
            ix = findall(wells[wi] .== Inf);
            wells[wi][ix] .= NaN;
        end 
    end
    tbins = [reverse(collect(0:-bin_duration:all_ftime[1])[1:end-1]); collect(bin_duration:bin_duration:all_ftime[end])];
    wellbins = zeros(Float64,length(tbins),length(wells));
    nwormbins = zeros(Float64,length(tbins),length(wells));
    for wi in eachindex(eachcol(wellbins))
        for ti in eachindex(eachrow(wellbins))
            wellbins[ti,wi] = sum(wells[wi][(all_ftime .> tbins[ti] - bin_duration) .& (all_ftime .<= tbins[ti])])/bin_duration;
            nwormbins[ti,wi] = mean(nworms[wi][(all_ftime .> tbins[ti] - bin_duration) .& (all_ftime .<= tbins[ti])]);
        end
    end
    
    #set data as NaN if no recorded data at that timepoint
    if MWT_output_timestamp_offset
        for ti in eachindex(tbins)
            if ti == 1
                if !any(all_ftime .<= tbins[ti])
                    wellbins[ti,:] .= NaN;
                end
                continue
            end
            if !any((all_ftime .> tbins[ti-1]) .& (all_ftime .<= tbins[ti]))
                wellbins[ti,:] .= NaN;
            end
        end
    end
    
    if !normalize_wl
        wellbins .*= pixelsize;#adjust for scale
    end
    
    if negativebaseline
        tbins = (tbins .- bin_duration/2) ./ bin_duration;#adjust time by half a bin
    else
        tbins = tbins ./ bin_duration;
    end
    println("    Progress: 100.0%");
    println("  Done!");
    println("Done!");
    println();
    println("Saving ...");
    println("  speed.txt");
    println("  nworms.txt");
    writedlm(speed_output, [tbins wellbins], ' ');
    writedlm(nworms_output, [tbins nwormbins], ' ');
    println("Done!\n");
    return nothing
end
