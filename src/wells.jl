
function old_setwells(nrows, ncolumns, d, xtl, ytl, xtr, ytr, xbl, ybl)
    r = round(Int64,d/2);
    x = xtl+r;#long axis (horizontal)
    y = ytl+r;#short axis (vertical)
    dx = round((xtr-xtl)/(ncolumns-1));
    dy = round((ybl-ytl)/(nrows-1));
    xe = round((xbl-xtl)/(ncolumns-1));#error
    ye = round((ytr-ytl)/(nrows-1));#error
    #should have been:
    #xe = round((xbl-xtl)/(nrows-1));
    #ye = round((ytr-ytl)/(ncolumns-1));
     
    plate = Array{Float64}(undef,nrows,ncolumns,2);
    for i in eachindex(eachcol(plate[:,:,1]))
        plate[:,i,1] = collect(0:nrows-1);
    end
    for i in eachindex(eachrow(plate[:,:,2]))
        plate[i,:,2] = collect(0:ncolumns-1);
    end

    global well_centers = [plate[:,:,1][:] plate[:,:,2][:]];

    for i in eachindex(eachrow(well_centers))
        well_centers[i,:] = [y+dy*well_centers[i,1]+ye*well_centers[i,2], x+dx*well_centers[i,2]+xe*well_centers[i,1]];
    end

    global well_radius = r;
    
    nwells = nrows*ncolumns;
    global wells = Array{Array{Float64}}(undef,nwells);
    global nworms = Array{Array{Int64}}(undef,nwells);
    for wi in eachindex(wells)
        wells[wi] = zeros(Float64,length(all_ftime));
        nworms[wi] = zeros(Int64,length(all_ftime));
    end

    return nothing
end

function setwells(nrows, ncolumns, d, xtl, ytl, xtr, ytr, xbl, ybl)
    r = d/2;
    well1x = xtl+r;#long axis (horizontal)
    well1y = ytl+r;#short axis (vertical)
    dx = vmagnitude([xtl,ytl] .- [xtr,ytl])/(ncolumns-1);
    dy = vmagnitude([xtl,ytl] .- [xbl,ybl])/(nrows-1);
    a = -mean([tan((ytl-ytr)/(xtr-xtl)), tan((xbl-xtl)/(ybl-ytl))]);
    
    plate = Array{Float64}(undef,nrows,ncolumns,2);
    for i in eachindex(eachcol(plate[:,:,1]))
        plate[:,i,1] = collect(0:nrows-1);
    end
    for i in eachindex(eachrow(plate[:,:,2]))
        plate[i,:,2] = collect(0:ncolumns-1);
    end

    global well_centers = [plate[:,:,1][:] plate[:,:,2][:]];

    for i in eachindex(eachrow(well_centers))
        wellx = well1x+dx*well_centers[i,2];
        welly = well1y+dy*well_centers[i,1];
        dwellx = wellx - well1x;
        dwelly = welly - well1y;
        
        rotatedx = dwellx * cos(a) - dwelly * sin(a);
        rotatedy = dwellx * sin(a) + dwelly * cos(a);
        
        well_centers[i,:] = [well1y+rotatedy, well1x+rotatedx];#switch x and y to match data
    end

    global well_radius = r;
    
    nwells = nrows*ncolumns;
    global wells = Array{Array{Float64}}(undef,nwells);
    global nworms = Array{Array{Int64}}(undef,nwells);
    for wi in eachindex(wells)
        wells[wi] = zeros(Float64,length(all_ftime));
        nworms[wi] = zeros(Int64,length(all_ftime));
    end

    return nothing
end

##

#Reference global variables

#finds well that worm's mean position is in
function whichwell(i)
    dc = vmagnitude.(eachrow(well_centers .- repeat(mean(clew[i].fpos,dims=1),size(well_centers,1),1)));
    return minimum(dc) <= well_radius ? findfirst(dc .== minimum(dc)) : nothing
end


#=
#finds well that worm is in the entire time, else nothing
function whichwell(i)
    dc = Array{Float64}(undef,length(clew[i].ftime),length(wells));
    for dci in eachindex(eachrow(dc))
        dc[dci,:] = vmagnitude.(eachrow(well_centers .- repeat(clew[i].fpos[dci,:]',size(well_centers,1),1)));
    end
    selectedwells = map(dci -> findfirst(dc[dci,:] .== minimum(dc[dci,:])),eachindex(eachrow(dc)));
    if all(selectedwells .== selectedwells[1])
        return maximum(dc[:, selectedwells[1]]) <= well_radius ? selectedwells[1] : nothing
    else
        return nothing
    end
end
=#

#finds well that worm's spent any time in
function inwell(wi,i)
    dc = vmagnitude.(eachrow(clew[i].fpos .- repeat(well_centers[wi,:]',size(clew[i].fpos,1),1)));
    return any(dc .<= well_radius)
end

#for old_scripts
function find_wells_and_sdp(clew)
    println("  Calculating Significant Movement ...\n")
    for wormi in eachindex(clew)
        if wormi%10 == 0
            print("\u1b[1F");
            print("    Progress: ");
            printstyled(string(round(wormi/length(clew)*100,digits=1),"%"), color=:cyan, bold=true);
            println();
        end
        wellix = findall(inwell.(collect(1:size(well_centers,1)),wormi));
        if length(wellix) == 0
            continue
        end
        
        if length(wellix) == 1
            dc = vmagnitude.(eachrow(clew[wormi].fpos .- repeat(well_centers[wellix[1],:]',size(clew[wormi].fpos,1),1)));
            validix = findall(dc .<= well_radius);
            if (clew[wormi].ftime[validix[end]] - clew[wormi].ftime[validix[1]]) < 10
                continue
            end
            if length(validix) == length(dc)
                sdp_perworm!(clew[wormi]);
                clew[wormi].wellnum = wellix[1];
            else
                tempworm = copy(clew[wormi]);
                tempworm.fnum = clew[wormi].fnum[validix];
                tempworm.ftime = clew[wormi].ftime[validix];
                tempworm.fpos = clew[wormi].fpos[validix,:];
                tempworm.dp = clew[wormi].dp[validix];
                sdp_perworm!(tempworm);
                clew[wormi].dp[validix] = tempworm.dp;
                clew[wormi].dp_threshold = tempworm.dp_threshold;
                
                valid_ranges = map((x, y) -> validix[x]:validix[y], [[1]; findall((validix[2:end] .- validix[1:end-1]) .> 1) .+ 1], [findall((validix[2:end] .- validix[1:end-1]) .> 1); [length(validix)]]);
            
                
                starti = findfirst(all_ftime .== clew[wormi].ftime[validix[1]]);
                endi = findfirst(all_ftime .== clew[wormi].ftime[validix[end]]);
                fnum = clew[wormi].fnum[validix[1]];
                fpos = Array{Float64}(undef,length(starti:endi),2);
                fill!(fpos,NaN);
                newworm = Worm(collect(fnum:fnum+(endi-starti)), all_ftime[starti:endi], fpos, clew[wormi].spinelength, zeros(length(starti:endi)), clew[wormi].dp_threshold, wellix[1]);
                
                for vr in valid_ranges
                    #if length(vr) > 1
                        #push!(clew, Worm(clew[wormi].fnum[vr], clew[wormi].ftime[vr], clew[wormi].fpos[vr,:], clew[wormi].spinelength, clew[wormi].dp[vr], clew[wormi].dp_threshold, wellix[1]));
                        #clew[end].dp[1] = 0;
                        starti = findfirst(newworm.ftime .== clew[wormi].ftime[vr][1]);
                        endi = starti + (vr.stop - vr.start);
                        newworm.fnum[starti:endi] = clew[wormi].fnum[vr];
                        newworm.ftime[starti:endi] = clew[wormi].ftime[vr];
                        newworm.fpos[starti:endi,:] = clew[wormi].fpos[vr,:];
                        newworm.dp[starti:endi] = clew[wormi].dp[vr];
                    #end
                end
                clew[wormi] = newworm;
            end
        else
            for welli in wellix
                dc = vmagnitude.(eachrow(clew[wormi].fpos .- repeat(well_centers[welli,:]',size(clew[wormi].fpos,1),1)));
                validix = findall(dc .<= well_radius);
                if (clew[wormi].ftime[validix[end]] - clew[wormi].ftime[validix[1]]) < 10
                    continue
                end
                
                tempworm = copy(clew[wormi]);
                tempworm.fnum = clew[wormi].fnum[validix];
                tempworm.ftime = clew[wormi].ftime[validix];
                tempworm.fpos = clew[wormi].fpos[validix,:];
                tempworm.dp = clew[wormi].dp[validix];
                sdp_perworm!(tempworm);
                clew[wormi].dp[validix] = tempworm.dp;
                clew[wormi].dp_threshold = tempworm.dp_threshold;
                
                valid_ranges = map((x, y) -> validix[x]:validix[y], [[1]; findall((validix[2:end] .- validix[1:end-1]) .> 1) .+ 1], [findall((validix[2:end] .- validix[1:end-1]) .> 1); [length(validix)]]);
                
                starti = findfirst(all_ftime .== clew[wormi].ftime[validix[1]]);
                endi = findfirst(all_ftime .== clew[wormi].ftime[validix[end]]);
                fnum = clew[wormi].fnum[validix[1]];
                fpos = Array{Float64}(undef,length(starti:endi),2);
                fill!(fpos,NaN);
                newworm = Worm(collect(fnum:fnum+(endi-starti)), all_ftime[starti:endi], fpos,clew[wormi].spinelength, zeros(length(starti:endi)), clew[wormi].dp_threshold, welli);
                
                for vr in valid_ranges
                    #if length(vr) > 1
                        #push!(clew, Worm(clew[wormi].fnum[vr], clew[wormi].ftime[vr], clew[wormi].fpos[vr,:], clew[wormi].spinelength, clew[wormi].dp[vr], clew[wormi].dp_threshold, welli));
                        #clew[end].dp[1] = 0;
                        starti = findfirst(newworm.ftime .== clew[wormi].ftime[vr][1]);
                        endi = starti + (vr.stop - vr.start);
                        newworm.fnum[starti:endi] = clew[wormi].fnum[vr];
                        newworm.ftime[starti:endi] = clew[wormi].ftime[vr];
                        newworm.fpos[starti:endi,:] = clew[wormi].fpos[vr,:];
                        newworm.dp[starti:endi] = clew[wormi].dp[vr];
                    #end
                end
                push!(clew, newworm);
            end
        end
    end
    #remove worms that are not in a well and worms split into multiple paths
    validix = findall(map(w -> !isnothing(w.wellnum), clew));
    clew = clew[validix];
    print("\u1b[1F");
    println("    Progress: 100.0%");
    println("  Done!");
    return clew
end

function place_data(wi, i)
    startix = findfirst(all_ftime .== clew[i].ftime[1]);
    stopix = findfirst(all_ftime .== clew[i].ftime[end]);
    global wells[wi][startix:stopix] .+= clew[i].dp;
    global nworms[wi][startix:stopix] .+= 1;
    return nothing
end

#uses whichwell to place
function place_all()
    println("  Sorting Worms into Wells ...\n");
    for wormi in eachindex(clew)
        if wormi%10 == 0
            print("\u1b[1F");
            print("    Progress: ");
            printstyled(string(round(wormi/length(clew)*100,digits=1),"%"), color=:cyan, bold=true);
            println();
        end
        welli = whichwell(wormi);
        clew[wormi].wellnum = welli;
        if isnothing(welli)# | (length(ftime[i]) < 30)#Put in filters here?
            continue
        end
        place_data(welli,wormi);
    end
    print("\u1b[1F");
    println("    Progress: 100.0%");
    println("  Done!");
    return nothing
end

#uses Worm.wellnum to place
function old_place_all()
    println("  Sorting Worms into Wells ...\n");
    for wormi in eachindex(clew)
        if wormi%10 == 0
            print("\u1b[1F");
            print("    Progress: ");
            printstyled(string(round(wormi/length(clew)*100,digits=1),"%"), color=:cyan, bold=true);
            println();
        end
        place_data(clew[wormi].wellnum, wormi);
    end
    print("\u1b[1F");
    println("    Progress: 100.0%");
    println("  Done!");
    return nothing
end
