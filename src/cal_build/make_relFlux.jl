



f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
df = DataFrame(read(f["$(parg["tele"])/$(mjd)/exposures"]))
close(f)
function get_1d_name_partial(expid)
    return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal=true) * ".jld2"
end
all1Da = get_1d_name_partial.(expid_list)

function get_and_save_relFlux(fname)
    absthrpt, relthrpt, bitmsk_relthrpt = get_relFlux(fname)
    outfname = replace(replace(fname, "apred"=>"$(cal_type)_flats"),"ar1Dcal" => "$(cal_type)Flux")
    if !ispath(outfname)
        mkpath(outfname)
    end
    jldsave(outfname; absthrpt, relthrpt, bitmsk_relthrpt)
end