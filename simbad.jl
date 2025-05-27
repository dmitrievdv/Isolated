using VirtualObservatory
using DataFrames
using HTTP
using JSON3



function star_coords(star_name)
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
    df_gaia[:, ["dist", "designation", "ra", "dec", "phot_g_mean_mag"]]
    sort!(df_gaia, "dist")[:, ["dist", "designation", "ra", "dec", "phot_g_mean_mag"]]

    star_g_mag = df_gaia.phot_g_mean_mag[1]
    bright_stars_mag = df_gaia.phot_g_mean_mag[df_gaia.phot_g_mean_mag .< star_g_mag + mag_threshold]

    return length(bright_stars_mag) - 1 == 0
end

function get_tess_sectors(star_ra, star_dec; product = "SPOC")
    resp = HTTP.get("https://mast.stsci.edu/tesscut/api/v0.1/sector?ra=$star_ra&dec=$star_dec&radius=1m&product=$product")
    JSON3.read(resp.body)
end

function get_tess_cutouts(star_ra, star_dec, width, height; star_name = "star", product = "SPOC")
    mkpath("$star_name")
    HTTP.download("https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=$star_ra&dec=$star_dec&y=$height&x=$width", "$star_name/$star_name.zip")
end

star_ra, star_dec = star_coords("T Cha")

# tess_sectors = get_tess_sectors(star_ra, star_dec)
get_tess_cutouts(star_ra, star_dec, 15, 15)
# check_star_for_isolation.(["GM Aur", "TW Hya"], 2/60, 2/60, 3)

