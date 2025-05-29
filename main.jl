using VirtualObservatory
using DataFrames
using HTTP
using ZipArchives
using Mmap
# using BufferedStreams
using JSON3
using CSV
using Photometry

star_directory = "stars_julia"

include("lightcurves.jl")
include("databases.jl")
include("tess-queries.jl")

function extract_tess_cutouts(star_name = "star")
    nospace_star_name = replace(star_name, " " => "_")
    archive = ZipReader(mmap(open("$star_directory/$nospace_star_name/$nospace_star_name.zip")))
    entry_names = zip_names(archive)
    for entry_name in entry_names
        data = zip_readentry(archive, entry_name)
        open("$star_directory/$nospace_star_name/$entry_name", "w") do io 
            write(io, data)
        end
    end
end

function get_sectors(star_name)
    nospace_star_name = replace(star_name, " " => "_")
    if !isfile("$star_directory/$nospace_star_name/$nospace_star_name.zip")
        println("No tesscut file is found. Download it using get_tess_sectors(...).")
        return Int[]
    end

    archive = ZipReader(mmap(open("$star_directory/$nospace_star_name/$nospace_star_name.zip")))
    entry_names = zip_names(archive)
    n_entries = zip_nentries(archive)
    sectors = zeros(Int, n_entries)
    for i_entry = 1:n_entries
        data = zip_readentry(archive, entry_names[i_entry])
        fits = FITS(data)
        sectors[i_entry] = read_key(fits[1], "SECTOR")[1]
    end
    return sectors
end

function get_tesscut_lc(fits)
    println("lol")
end

function get_lcs(star_name)
    nospace_star_name = replace(star_name, " " => "_")
    archive = ZipReader(mmap(open("$star_directory/$nospace_star_name/$nospace_star_name.zip")))
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

        open("$star_directory/$nospace_star_name/$sector-lc.dat", "w") do io
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
            savefig(plt, "$star_directory/$nospace_star_name/$sector-lc.pdf")
            savefig(plt, "$star_directory/$nospace_star_name/$sector-lc.png")
            jd_all, flux_all = mergelc(jd_all, flux_all, jd, flux)
        end
    end
    jd_all, flux_all
end

# df_stars = get_simbad_young_stars(12, 8)
# CSV.write("young_simbad.csv", df_stars)
df_stars = CSV.read("isolated_224.csv", DataFrame)



# df_isolated = begin
#     df_isolated = DataFrame()
#     n_stars = nrow(df_stars)
    
#     for i_star = 1:n_stars
#         star_name = df_stars.star_name[i_star]
#         star_ra = df_stars.ra[i_star]
#         star_dec = df_stars.dec[i_star]

#         is_isolated = check_star_for_isolation(star_ra, star_dec, 2/60, 2/60, 4, gaia = "dr2")
#         if is_isolated
#             println("$star_name is isolated")
#             push!(df_isolated, df_stars[i_star, :])
#         else
#             println("$star_name is not isolated")
#         end
#     end
#     df_isolated
# end


# begin
#     i_star = 60
#     star_name = df_stars.star_name[i_star]
#     star_ra = df_stars.ra[i_star]
#     star_dec = df_stars.dec[i_star]

#     is_isolated = check_star_for_isolation(star_ra, star_dec, 2/60, 2/60, 4, gaia = "dr2")
    
#     if is_isolated
#         println("$star_name is isolated")
#         tess_sectors = get_tess_sectors(star_ra, star_dec)
#         if length(tess_sectors) > 0
#             println("$star_name: $(length(tess_sectors)) TESS sectors")
#             get_tess_cutouts(star_ra, star_dec, 15, 15, star_name = star_name)
#             get_lcs(star_name)
#         else
#             println("$star_name: No TESS sectors")
#         end
#     else
#         println("$star_name isn't isolated")
#     end
# end

# star_ra, star_dec = get_star_coords("TW Hya")
# tess_sectors = get_tess_sectors(star_ra, star_dec)
# get_tess_cutouts(star_ra, star_dec, 15, 15, star_name = "TW Hya")
# check_star_for_isolation.(["GM Aur", "TW Hya"], 2/60, 2/60, 3)

