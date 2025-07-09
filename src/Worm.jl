mutable struct Worm
    fnum
    ftime
    fpos
    spinelength
    dp
    dp_threshold
    wellnum
end
function Worm()
    return Worm(Array{Int64}(undef,0), Array{Float64}(undef,0), Array{Float64}(undef,0), nothing, Array{Float64}(undef,0), nothing, nothing)
end
function Worm(totalframes)
    fnum = Array{Int64}(undef,totalframes);
    ftime = Array{Float64}(undef,totalframes);
    fpos = Array{Float64}(undef,totalframes,2);
    spinelength = nothing;
    dp = Array{Float64}(undef,totalframes);
    dp_threshold = nothing;
    wellnum = nothing;
    return Worm(fnum, ftime, fpos, spinelength, dp, dp_threshold,wellnum)
end
function Worm(worm,nframes)
    return Worm(worm.fnum[1:nframes], worm.ftime[1:nframes], worm.fpos[1:nframes,:], worm.spinelength,worm.dp[1:nframes],worm.dp_threshold,worm.wellnum) 
end

import Base.copy
function copy(worm::Worm)
    return Worm(worm.fnum, worm.ftime, worm.fpos, worm.spinelength, worm.dp, worm.dp_threshold, worm.wellnum)
end
