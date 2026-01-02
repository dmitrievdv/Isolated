include("main.jl")

# isolated_df = CSV.read("isolated_2.csv", DataFrame)
# star_names = isolated_df.star_name[isolated_df[:, :max_near_mag] .> 3]
star_name = "CD-30 12558"
cut_size = 30

get_all_data([star_name], cut_size; rewrite_files = false, rewrite_gaia_stars_file = false,
     aperture_radius = 3, Δm_R = 5, day_step = 1, jd_box = 0.1, σ_tol = 5, n_out = 20)

sectors = find_tess_sectors(star_name, tess_max_sectors)

fig = plot_cuts(star_name, sectors[1], cut_size, cut_size; aperture_radius = 3, Δm_R = 5)
