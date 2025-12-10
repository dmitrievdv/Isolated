include("main.jl")

# isolated_df = CSV.read("isolated_2.csv", DataFrame)
# star_names = isolated_df.star_name[isolated_df[:, :max_near_mag] .> 3]
star_names = ["CI Tau"]

get_all_data(star_names, 30; rewrite_files = true, day_step = 1, jd_box = 0.1, σ_tol = 5, n_out = 20)