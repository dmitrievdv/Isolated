function get_star_gaia_data(star_name; gaia = "dr3")
    df_simbad = DataFrame(execute(TAPService(:simbad), "select top 5 oid, main_id, ids 
    FROM ident 
    JOIN basic ON ident.oidref = basic.oid
    JOIN ids on ident.oidref = ids.oidref
    where id = '$star_name'"))

    gaia_regex = Regex("gaia $gaia (?<id>[0-9]+)", "i")
    gaia_id = parse(Int, match(gaia_regex, df_simbad.ids[1])[:id])

    df_gaia = DataFrame(execute(TAPService(:gaia), "select top 1 * from gaia$gaia.gaia_source where source_id = $gaia_id"))
    return df_gaia[1,:]
end

function get_star_coords(star_name)
    df_simbad = DataFrame(execute(TAPService(:simbad), "select top 5 * FROM ident JOIN basic ON ident.oidref = basic.oid where id = '$star_name'"))

    return df_simbad.ra[1], df_simbad.dec[1]
end

function get_star_coords_gaia(star_name; gaia = "dr3")
    df_simbad = DataFrame(execute(TAPService(:simbad), "select top 5 oid, main_id, ids 
    FROM ident 
    JOIN basic ON ident.oidref = basic.oid
    JOIN ids on ident.oidref = ids.oidref
    where id = '$star_name'"))

    gaia_regex = Regex("gaia $gaia (?<id>[0-9]+)", "i")
    gaia_id = parse(Int, match(gaia_regex, df_simbad.ids[1])[:id])

    df_gaia = DataFrame(execute(TAPService(:gaia), "select top 1 source_id, ra, dec from gaia$gaia.gaia_source where source_id = $gaia_id"))
    return df_gaia.ra[1], df_gaia.dec[1]
end

function brightest_near_relative_mag(star_ra, star_dec, box_width, box_height; gaia = "dr3")
    df_gaia = DataFrame(execute(TAPService(:gaia), "select * from gaia$gaia.gaia_source where CONTAINS(POINT('ICRS', ra, dec),"*
                                              " BOX('ICRS', $star_ra, $star_dec, $box_width, $box_height)) = 1"))

    gaia_dist = @. sqrt((cos(star_dec/180*π)*(df_gaia.ra - star_ra))^2 + (df_gaia.dec - star_dec)^2)
    df_gaia[!, "dist"] = gaia_dist
    df_gaia[:, ["dist", "designation", "ra", "dec", "phot_g_mean_mag"]]
    sort!(df_gaia, "dist")

    near_g_mag_sorted = filter(x -> !ismissing(x), df_gaia[2:end, :phot_g_mean_mag])
    
    star_g_mag = try
        if ismissing(df_gaia[1, :phot_g_mean_mag])
            throw(MissingException("missing phot_g_mean_mag for the star"))
        end
        df_gaia[1, :phot_g_mean_mag]
    catch e
        showerror(stdout, e)
        0.0
    end

    brightest_near_g_mag = try 
        minimum(near_g_mag_sorted)
    catch e
        showerror(stdout, e)
        star_g_mag
    end

    return brightest_near_g_mag - star_g_mag
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

    star_rp_mag = df_gaia.phot_rp_mean_mag[1]
    bright_stars_mag = df_gaia.phot_rp_mean_mag[df_gaia.phot_rp_mean_mag .< star_rp_mag + mag_threshold]

    return length(bright_stars_mag) - 1 == 0
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

function check_for_isolation(df_stars, box_size, mag_thres; gaia = "dr2")
    df_isolated = DataFrame()
    n_stars = nrow(df_stars)
    
    for i_star = 1:n_stars
        star_name = df_stars.star_name[i_star]
        star_ra = df_stars.ra[i_star]
        star_dec = df_stars.dec[i_star]

        is_isolated = check_star_for_isolation(star_ra, star_dec, box_size/60, box_size/60, mag_thres, gaia = gaia)
        if is_isolated
            println("$star_name is isolated")
            push!(df_isolated, df_stars[i_star, :])
        else
            println("$star_name is not isolated")
        end
    end
    df_isolated
end