function get_dp_threshold(clew)
    all_nframes = map(w->length(w.ftime)-1,clew);
    all_dp = Array{Float64}(undef,sum(all_nframes));
    currenti = 1;
    for wi in eachindex(clew)
        all_dp[currenti:(currenti+all_nframes[wi]-1)] = vmagnitude.(eachrow(clew[wi].fpos[2:end,:] .- clew[wi].fpos[1:end-1,:])) ./ (normalize_wl ? clew[wi].spinelength : 1);#normalize to wormlength
        currenti += all_nframes[wi];
    end
    return median(all_dp) + 4 * (quantile(all_dp,.75) - median(all_dp))
end

function sdp_perclew!(clew)
    println("  Calculating Significant Movement ...\n")
    dp_threshold = get_dp_threshold(clew);
    for wi in eachindex(clew)
        if wi%100 == 0
            print("\u1b[1F");
            print("    Progress: ");
            printstyled(string(round(wi/length(clew)*100,digits=1),"%"), color=:cyan, bold=true);
            println();
        end
        currentp = clew[wi].fpos[1,:];
        for fi in eachindex(clew[wi].dp)
            if fi == 1
                clew[wi].dp[fi] = 0;
                continue
            end
            clew[wi].dp[fi] = vmagnitude(clew[wi].fpos[fi,:] .- currentp) ./ (normalize_wl ? clew[wi].spinelength : 1);#normalize to wormlength
            if clew[wi].dp[fi] > dp_threshold
                currentp = clew[wi].fpos[fi,:];
            else
                clew[wi].dp[fi] = 0;
            end
        end
    end
    print("\u1b[1F");
    println("    Progress: 100.0%");
    println("  Done!");
    return dp_threshold
end

function sdp_perworm!(clew::Array{Worm})
    println("  Calculating Significant Movement ...\n")
    for wi in eachindex(clew)
        if wi%100 == 0
            print("\u1b[1F");
            print("    Progress: ");
            printstyled(string(round(wi/length(clew)*100,digits=1),"%"), color=:cyan, bold=true);
            println();
        end
        dp = vmagnitude.(eachrow(clew[wi].fpos[2:end,:] .- clew[wi].fpos[1:end-1,:])) ./ (normalize_wl ? clew[wi].spinelength : 1);
        clew[wi].dp_threshold = median(dp) + 4 * (quantile(dp,.75) - median(dp));
        currentp = clew[wi].fpos[1,:];
        for fi in eachindex(clew[wi].dp)
            if fi == 1
                clew[wi].dp[fi] = 0;
                continue
            end
            clew[wi].dp[fi] = vmagnitude(clew[wi].fpos[fi,:] .- currentp) ./ (normalize_wl ? clew[wi].spinelength : 1);#normalize to wormlength
            if clew[wi].dp[fi] > clew[wi].dp_threshold
                currentp = clew[wi].fpos[fi,:];
            else
                clew[wi].dp[fi] = 0;
            end
        end
    end
    print("\u1b[1F");
    println("    Progress: 100.0%");
    println("  Done!");
    return nothing
end

function sdp_perworm!(worm::Worm)
    dp = vmagnitude.(eachrow(worm.fpos[2:end,:] .- worm.fpos[1:end-1,:])) ./ (normalize_wl ? worm.spinelength : 1);
    worm.dp_threshold = median(dp) + 4 * (quantile(dp,.75) - median(dp));
    currentp = worm.fpos[1,:];
    for fi in eachindex(worm.dp)
        if fi == 1
            worm.dp[fi] = 0;
            continue
        end
        worm.dp[fi] = vmagnitude(worm.fpos[fi,:] .- currentp) ./ (normalize_wl ? worm.spinelength : 1);#normalize to wormlength
        if worm.dp[fi] > worm.dp_threshold
            currentp = worm.fpos[fi,:];
        else
            worm.dp[fi] = 0;
        end
    end
    return nothing
end
