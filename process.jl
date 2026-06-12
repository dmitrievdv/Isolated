include("main.jl")
include("process-func.jl")

star_names = strip.(replace.(readdir(star_directory), "_" => " "))

df_lc_stats = DataFrame(star_name = String[], sector = Int[], periodicity = Float64[], asymmetry = Float64[])
cut_size = 15

for star_name in star_names
    print(star_name, " ")
    sectors = try
        find_tess_sectors(star_name, tess_max_sectors)
    catch e
        print("\n")
        continue
    end
    print(", sectors: ")
    for sector in sectors
        try 
            df_lc = load_light_curve(star_name, sector, cut_size)
            periodicity = calc_periodicity(df_lc)
            asymmetry = calc_asymmetry(df_lc)
            push!(df_lc_stats, [star_name, sector, periodicity, asymmetry])
        catch e
            continue
        end
        print(sector, " ")
    end
    print("\n")
end

CSV.write("lc_stats.csv", df_lc_stats)
fig = Figure()
ax = Axis(fig[1,1], yreversed = true)

for (i_star, star_name) in enumerate(star_names)
    df_star_lc_stats = df_lc_stats[df_lc_stats.star_name .== star_name, :]
    text!(ax, df_star_lc_stats.periodicity, df_star_lc_stats.asymmetry; text = fill(star_name, nrow(df_star_lc_stats)), 
                color = i_star, colormap = :prism, colorrange = (1,length(star_names)))
end

fig

# df_lc = load_light_curve("CQ Tau", 43, 15)

# jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)
# cleaned_mags = deepcopy(mags)
# clean_flux_sigma!(jds, cleaned_mags, 0.1, 10, 20)
# cleaned_jds, cleaned_mags = delete_nans(jds, cleaned_mags)

# calc_periodicity(df_lc)

# fig = Figure()
# ax_lc = Axis(fig[1,1], yreversed = true)
# lines!(ax_lc, jds, mags)
# lines!(ax_lc, cleaned_jds, cleaned_mags)

# fig