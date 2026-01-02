include("main.jl")

# isolated_df = CSV.read("isolated_2.csv", DataFrame)
# star_names = isolated_df.star_name[isolated_df[:, :max_near_mag] .> 3]
star_name = "19 Cep"
cut_size = 30

get_all_data([star_name], cut_size; rewrite_files = true, aperture_radius = 7, day_step = 1, jd_box = 0.1, σ_tol = 5, n_out = 20)

sectors = find_tess_sectors(star_name, tess_max_sectors)

for sector in sectors
    fig = plot_cuts(star_name, sector, cut_size, cut_size; aperture_radius = 7)
    resize!(fig, 800,800)

    save("$sector.png", fig)
end