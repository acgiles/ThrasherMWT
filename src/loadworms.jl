mutable struct MarkedBA
    ba::Array{UInt8, 1}
    i::Int64
end

import Base.readuntil
function readuntil(mba::MarkedBA, b::UInt8; keep=false)
    #Get the parse value of the string in ba from starti to the next b
    ba = mba.ba;
    starti = mba.i;
    endi = mba.i;
    if ba[endi] == b
        mba.i += 1;
        return keep ? String(ba[starti:endi]) : ""
    end
    try
        while ba[endi+1] != b
            endi += 1;
        end
    catch e
        if typeof(e) == BoundsError
            mba.i = length(mba.ba);
            return nothing
        else
            error(e)
        end
    end
    mba.i = endi+2;
    return keep ? String(ba[starti:endi+1]) : String(ba[starti:endi])
end

function find_spinelength(mba)
    spinesequence = Array{Int64}(undef,22);
    for i in eachindex(spinesequence)
        spinesequence[i] = parse(Int64,readuntil(mba,0x20));#UInt8(' ') = 0x20
    end
    spine = Matrix(reshape(spinesequence,2,Int64(length(spinesequence)/2))');
    return sum(map(r -> vmagnitude(r), eachrow(spine[2:end,:] .- spine[1:end-1,:])))
end

function get_spinestart(sub_ba,startix)
    #find the ix that has a ' ', with a '%' before it, but not a second '%' before it.
    #UInt8(' ') = 0x20
    #UInt8('%') = 0x25
    #if not nothing, add the startix to make it relative to entire ba
    ix = findfirst((sub_ba .== 0x20) .& ([[0x00]; sub_ba[1:end-1]] .== 0x25) .& ([[0x00,0x00]; sub_ba[1:end-2]] .!= 0x25));
    return isnothing(ix) ? nothing : ix + startix;
end

getfnum!(mba) = parse(Int64, readuntil(mba, 0x20));#UInt8(' ') = 0x20
getftime!(mba) = parse(Float64, readuntil(mba, 0x20));#UInt8(' ') = 0x20
getfpos!(mba) = parse(Float64, readuntil(mba, 0x20));#UInt8(' ') = 0x20

function advance!(mba)
    readuntil(mba, 0x20);#UInt8(' ') = 0x20
    return nothing
end

function getworm!(ba, worm; ftime_offset=0)
    fnum = worm.fnum;
    ftime = worm.ftime;
    fpos = worm.fpos;
    linestarts = [[0]; findall(ba .== 0x0a)] .+ 1;#UInt8('\n') = 0x0a
    
    #removes lines without spines
    linestops = [linestarts[2:end] .- 1;[length(ba)]];
    spinestarts = map((startix,stopix)->get_spinestart(ba[startix:stopix],startix),linestarts,linestops);
    linestarts = linestarts[.!isnothing.(spinestarts)];
    spinestarts = spinestarts[.!isnothing.(spinestarts)];
    spinelengths = Array{Float64}(undef,length(spinestarts));
    if length(linestarts) == 0
        return 0
    end
    
    mba = MarkedBA(ba,1);
    for i in eachindex(linestarts)
        try
        mba.i = linestarts[i];
        fnum[i] = getfnum!(mba);
        ftime[i] = getftime!(mba) + ftime_offset;
        advance!(mba);
        fpos[i,1] = getfpos!(mba);
        fpos[i,2] = getfpos!(mba);
        if normalize_wl
            mba.i = spinestarts[i];
            spinelengths[i] = find_spinelength(mba);
        end
        catch e
            println(i)
            println(String(mba.ba))
            foreach(println,linestarts)
            error(e)
        end
    end
    if normalize_wl
        worm.spinelength = median(spinelengths);
    end
    return length(linestarts)
end

function getworms(ba, totalframes; ftime_offset=0)
    #indices for worm id, and start/stop of data rows
    #for io position it is ix-1;
    wids = [[3];findall((ba[1:end-1] .== 0x0a) .& (ba[2:end] .== 0x25)).+3];#UInt8('\n') = 0x0a and UInt8('%') = 0x25
    wstops = [wids[2:end] .- 5; [length(ba)-2]];
    wstarts = map(i->findfirst(ba[wids[i]:wstops[i]] .== 0x0a),eachindex(wids));
    ix = findall(.!isnothing.(wstarts));
    wids = wids[ix];
    wstops = wstops[ix];
    wstarts = wstarts[ix] .+ wids;
    clew = Array{Worm}(undef,length(wids));
    worm = Worm(totalframes);
    progress_string = "";
    println();
    for i in eachindex(wids)
        if i % 100 == 0
            print("\u1b[1F");
            progress_string = string("      ", round(i/length(wids)*100, digits=2),"%    ");
            printstyled(progress_string, color=:cyan, bold=true);
            println(" ");
        end
        nframes = getworm!(ba[wstarts[i]:wstops[i]], worm; ftime_offset = ftime_offset);
        if nframes > 0
            clew[i] = Worm(worm,nframes);
        else
            clew[i] = Worm();
        end
    end
    print("\u1b[1F");
    return clew
end

function load_recording(mwt_output; ftime_offset=0)
    println("  Working on ",splitpath(mwt_output)[end]," ...");
    fs = readdir(mwt_output,join=true);
    summary_file = fs[endswith.(fs,".summary")][1];
    totalframes = sum(read(summary_file) .== 0x0a);
    blobs = fs[endswith.(fs,".blobs")];
    clew = Array{Worm}(undef,0);
    for bi in eachindex(blobs)
        print("    ", splitpath(blobs[bi])[end], " ");
        printstyled(string(bi,"/",length(blobs),"\n"), color=:cyan, bold=true);
        new_clew = getworms(read(blobs[bi]),totalframes, ftime_offset = ftime_offset);
        clew = vcat(clew, new_clew);
        print("\u1b[1F");
    end
    print("\u1b[1F");
    println("  Finished ",splitpath(mwt_output)[end]," ", length(blobs), "/", length(blobs), " blobs files");
    return clew
end

function get_rec_dur(recording)
    fs = readdir(recording);
    summary_file = joinpath(recording,fs[findfirst(endswith.(fs,".summary"))]);
    ba = read(summary_file);
    offset = findall(ba .== 0x0a)[end-1];
    io = IOBuffer(ba);
    seek(io,offset);
    readuntil(io,' ');
    return parse(Float64,readuntil(io,' '));
end

function get_all_ftime(recording)
    fs = readdir(recording);
    summary_file = joinpath(recording,fs[findfirst(endswith.(fs,".summary"))]);
    ba = read(summary_file);
    lineix = [[0];findall(ba .== 0x0a)[1:end-1]];
    io = IOBuffer(ba);
    all_ftime = Array{Float64}(undef,length(lineix));
    for ti in eachindex(all_ftime)
        seek(io,lineix[ti]);
        readuntil(io,' ');
        all_ftime[ti] = parse(Float64,readuntil(io,' '));
    end
    return all_ftime
end

function getmwtdatetime(recording)
    recordingstring = splitdir(recording)[2];
    return DateTime(parse(Int64,recordingstring[1:4]),parse(Int64,recordingstring[5:6]),parse(Int64,recordingstring[7:8]),parse(Int64,recordingstring[10:11]),parse(Int64,recordingstring[12:13]),parse(Int64,recordingstring[14:15]))
end
