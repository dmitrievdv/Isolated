include("main.jl")

isolated_df = CSV.read("isolated_224.csv", DataFrame)
star_names = isolated_df.star_name

get_all_data(star_names, 15; rewrite_files = false, day_step = 1, jd_box = 0.1, σ_tol = 5, n_out = 20)