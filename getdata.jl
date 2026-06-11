include("main.jl")

# isolated_df = CSV.read("isolated_2.csv", DataFrame)
# star_names = isolated_df.star_name[isolated_df[:, :max_near_mag] .> 3]
star_name = "T Cha"
cut_size = 45

get_all_data([star_name], cut_size; rewrite_files = true,
     aperture_radius = 2, Δm_R = 5, day_step = 1, jd_box = 0.1, σ_tol = 5, n_out = 20)

sectors = find_tess_sectors(star_name, tess_max_sectors)

sector = sectors[2]
fig, i_cut = plot_cuts(star_name, sector, cut_size, cut_size; Δm_R = 5)
resize!(fig.scene, (1000,1000)); fig

# record(fig, "58.mp4", 400:500; framerate = 10) do t
#      i_cut[] = t
# end