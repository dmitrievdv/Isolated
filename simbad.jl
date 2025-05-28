using VirtualObservatory
using DataFrames
using HTTP
using ZipArchives
using Mmap
using BufferedStreams
using JSON3
using CSV

include("lightcurves.jl")

function get_star_coords(star_name)
    df_simbad = DataFrame(execute(TAPService(:simbad), "select top 5 * FROM ident JOIN basic ON ident.oidref = basic.oid where id = '$star_name'"))

    return df_simbad.ra[1], df_simbad.dec[1]
end

function check_star_for_isolation(star_name, box_width, box_height, mag_threshold; kwargs...)
    df_simbad = DataFrame(execute(TAPService(:simbad), "select top 5 * FROM ident JOIN basic ON ident.oidref = basic.oid where id = '$star_name'"))

    star_ra = df_simbad.ra[1]
    star_dec = df_simbad.dec[1]

    check_star_for_isolation(star_ra, star_dec, box_width, box_height, mag_threshold; kwargs...)
end

function check_star_for_isolation(star_ra, star_dec, box_width, box_height, mag_threshold; gaia = "dr3")
    df_gaia = DataFrame(execute(TAPService(:gaia), "select * from gaia$gaia.gaia_source where CONTAINS(POINT('ICRS', ra, dec),"*
                                              " BOX('ICRS', $star_ra, $star_dec, $box_width, $box_height)) = 1"))

    gaia_dist = @. sqrt((cos(star_dec/180*π)*(df_gaia.ra - star_ra))^2 + (df_gaia.dec - star_dec)^2)
    df_gaia[!, "dist"] = gaia_dist
    df_gaia[:, ["dist", "designation", "ra", "dec", "phot_rp_mean_mag"]]
    sort!(df_gaia, "dist")[:, ["dist", "designation", "ra", "dec", "phot_rp_mean_mag"]]

    star_g_mag = df_gaia.phot_g_mean_mag[1]
    bright_stars_mag = df_gaia.phot_g_mean_mag[df_gaia.phot_g_mean_mag .< star_g_mag + mag_threshold]

    return length(bright_stars_mag) - 1 == 0
end

function get_tess_sectors(star_ra, star_dec; product = "SPOC")
    resp = HTTP.get("https://mast.stsci.edu/tesscut/api/v0.1/sector?ra=$star_ra&dec=$star_dec&radius=1m&product=$product")
    JSON3.read(resp.body).results
end

function get_tess_cutouts(star_ra, star_dec, width, height; star_name = "star", product = "SPOC")
    nospace_star_name = replace(star_name, " " => "_")
    mkpath("$nospace_star_name")
    url = "https://mast.stsci.edu/tesscut/"
    api = "api/v0.1/astrocut?ra=$star_ra&dec=$star_dec&y=$height&x=$width"
    HTTP.download(url*api, "$nospace_star_name/$nospace_star_name.zip")
    # archive = ZipReader(mmap(open("$nospace_star_name/$nospace_star_name.zip")))
    # entry_names = zip_names(archive)
    # for entry_name in entry_names
    #     data = zip_readentry(archive, entry_name)
    #     open("$nospace_star_name/$entry_name", "w") do io 
    #         write(io, data)
    #     end
    # end
end

function extract_tess_cutouts(star_name = "star")
    nospace_star_name = replace(star_name, " " => "_")
    archive = ZipReader(mmap(open("$nospace_star_name/$nospace_star_name.zip")))
    entry_names = zip_names(archive)
    for entry_name in entry_names
        data = zip_readentry(archive, entry_name)
        open("$nospace_star_name/$entry_name", "w") do io 
            write(io, data)
        end
    end
end

function get_simbad_young_stars(mag_min, mag_max)
    df_simbad_sp = DataFrame(execute(TAPService(:simbad), """select oid, otype, main_id, ra, dec, sptype, V, R, "year",  ids
from basic 
join allfluxes on oid = allfluxes.oidref 
join ids on oid = ids.oidref
join (select oidref, sptype, "year"
      from mesSpT join REF on mesSpT.bibcode = REF.bibcode order by oidref,"year" desc) as sptyear on oid = sptyear.oidref
where (otype = 'TT*' or otype = 'Ae*' or otype = 'Or*') and V < $mag_min and V > $mag_max"""))

    n_rows = nrow(df_simbad_sp)
    df_simbad = DataFrame()

    old_main_id = ""
    for i = 1:n_rows
        main_id = df_simbad_sp[i, :main_id]
        if main_id != old_main_id
            old_main_id = main_id
            push!(df_simbad, df_simbad_sp[i, :])
        end
    end

    n_stars = nrow(df_simbad)
    df_simbad[!, "star_name"] = df_simbad[:, "main_id"]

    for i_star = 1:n_stars
        star_names = split(df_simbad[i_star, "ids"], "|")
        main_name = df_simbad[i_star, "star_name"]
        for star_name in star_names
            if star_name[1:2] == "V*"
                main_name = star_name[4:end]
            end
        end
        first_word = split(main_name)[1]
        if first_word[end] == '*'
            main_name = join(split(main_name)[2:end], " ")
        end
        main_name = join(split(main_name), " ")
        df_simbad[i_star, "star_name"] = main_name
    end
    df_simbad
end

function get_lc(star_name)
    nospace_star_name = replace(star_name, " " => "_")
    archive = ZipReader(mmap(open("$nospace_star_name/$nospace_star_name.zip")))
    entry_names = zip_names(archive)

    jd_all = Float64[]
    flux_all = Float64[]

    for entry_name in entry_names
        data = zip_readentry(archive, entry_name)
        fits = FITS(data)
        sector = parse(Int, entry_name[7:10])

        jd, flux = getlc(nospace_star_name, fits, sector = sector)
        jd, flux = uniformlc(jd, flux)

        jd, flux = cleanlc(jd, flux)

        mag = -2.5*log10.(flux) .+ 20.44

        open("$nospace_star_name/$sector-lc.dat", "w") do io
            println(io, "#jd flux")
            for n = 1:length(jd)
                @printf(io, "%15.6f %15.6e %15.6f\n", jd[n], flux[n], mag[n])
            end
        end
        if length(jd) > 0
            int_dates = collect(ceil(Int, jd[1]/5)*5:5:floor(Int, jd[end]/5)*5)
            string_dates = Dates.format.(julian2datetime.(int_dates), "d u Y")
            plt = plot(jd, tessmag.(flux), xticks = (int_dates, string_dates), label = false, rightmargin = 15px, yflip = true, ylabel = "TESS magnitude")
            # savefig(plt, "plots/$star-$sector.pdf")
            # savefig(plt, "plots/$star-$sector.png")
            savefig(plt, "$nospace_star_name/$sector-lc.pdf")
            savefig(plt, "$nospace_star_name/$sector.png")
            jd_all, flux_all = mergelc(jd_all, flux_all, jd, flux)
        end
    end
    jd_all, flux_all
end

df_stars = get_simbad_young_stars(12, 8)
CSV.write("young_simbad.csv", df_stars)
df_stars = CSV.read("young_simbad.csv", DataFrame)



df_isolated = begin
    df_isolated = DataFrame()
    n_stars = nrow(df_stars)
    
    for i_star = 1:n_stars
        star_name = df_stars.star_name[i_star]
        star_ra = df_stars.ra[i_star]
        star_dec = df_stars.dec[i_star]

        is_isolated = check_star_for_isolation(star_ra, star_dec, 2/60, 2/60, 4, gaia = "dr2")
        if is_isolated
            println("$star_name is isolated")
            push!(df_isolated, df_stars[i_star, :])
        else
            println("$star_name is not isolated")
        end
    end
    df_isolated
end


begin
    i_star = 105
    star_name = df_stars.star_name[i_star]
    star_ra = df_stars.ra[i_star]
    star_dec = df_stars.dec[i_star]

    is_isolated = check_star_for_isolation(star_ra, star_dec, 4/60, 4/60, 4, gaia = "dr2")
    
    if is_isolated
        println("$star_name is isolated")
        tess_sectors = get_tess_sectors(star_ra, star_dec)
        if length(tess_sectors) > 0
            println("$star_name: $(length(tess_sectors)) TESS sectors")
            get_tess_cutouts(star_ra, star_dec, 15, 15, star_name = star_name)
            get_lc(star_name)
        else
            println("$star_name: No TESS sectors")
        end
    else
        println("$star_name isn't isolated")
    end
end

# star_ra, star_dec = get_star_coords("TW Hya")
# tess_sectors = get_tess_sectors(star_ra, star_dec)
# get_tess_cutouts(star_ra, star_dec, 15, 15, star_name = "TW Hya")
# check_star_for_isolation.(["GM Aur", "TW Hya"], 2/60, 2/60, 3)

