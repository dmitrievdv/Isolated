function get_tess_sectors(star_ra, star_dec; product = "SPOC")
    resp = HTTP.get("https://mast.stsci.edu/tesscut/api/v0.1/sector?ra=$star_ra&dec=$star_dec&radius=1m&product=$product")
    JSON3.read(resp.body).results
end

function get_tess_cutouts(star_ra, star_dec, width, height; star_name = "star", product = "SPOC")
    nospace_star_name = replace(star_name, " " => "_")
    mkpath("$star_directory/$nospace_star_name")
    url = "https://mast.stsci.edu/tesscut/"
    api = "api/v0.1/astrocut?ra=$star_ra&dec=$star_dec&y=$height&x=$width"
    try 
        HTTP.download(url*api, "$star_directory/$nospace_star_name/$nospace_star_name.zip")
    catch e
        if isa(e, HTTP.StatusError)
            println(e.response["location"])
            HTTP.download(e.response["location"], "$star_directory/$nospace_star_name/$nospace_star_name.zip")
        else
            throw(e)
        end
    end
        #, "$nospace_star_name/$nospace_star_name.zip")
    # archive = ZipReader(mmap(open("$nospace_star_name/$nospace_star_name.zip")))
    # entry_names = zip_names(archive)
    # for entry_name in entry_names
    #     data = zip_readentry(archive, entry_name)
    #     open("$nospace_star_name/$entry_name", "w") do io 
    #         write(io, data)
    #     end
    # end
end