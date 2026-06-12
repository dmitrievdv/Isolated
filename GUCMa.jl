include("main.jl")
include("process-func.jl")

star_name = "GU CMa"
cut_size = 45
get_all_data([star_name], cut_size; rewrite_files = false,
     aperture_radius = 2, Δm_R = 5, day_step = 1, jd_box = 0.1, σ_tol = 5, n_out = 20)

sectors = find_tess_sectors(star_name, tess_max_sectors)

periods = zeros(length(sectors))
df_lcs = load_light_curve.(Ref(star_name), sectors, cut_size)
sector_phased_jds = Vector{Float64}[]
sector_cleaned_mags = Vector{Float64}[]

for i_sector in eachindex(sectors)
    df_lc = df_lcs[i_sector]
    period = 2*find_period(df_lc, n_out = 20, σ_tol = 5)
    println(period)
    periods[i_sector] = period
    jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)
    cleaned_mags = deepcopy(mags)
    clean_flux_sigma!(jds, mags, 0.1, 5, 20)
    cleaned_jds, cleaned_mags = delete_nans(jds, mags)

    push!(sector_cleaned_mags, cleaned_mags)
    push!(sector_phased_jds, cleaned_jds .% periods[1])
end

fig = Figure()
ax = Axis(fig[1,1], yreversed = true)
for i_sector in eachindex(sectors)
    min_phase = sector_phased_jds[i_sector][findmax(sector_cleaned_mags[i_sector])[2]]
    println(min_phase, findmax(sector_cleaned_mags)[2])
    scatter!(ax, mod.((sector_phased_jds[i_sector] .- min_phase), periods[1]), sector_cleaned_mags[i_sector], markersize = 5, label = "$(sectors[i_sector])")
end
axislegend(ax)
fig
