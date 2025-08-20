include("main.jl")

star_name = "V501 Aur"
cut_size = 15

jd_box = 0.1; σ_tol = 5; n_out = 20

sectors = find_tess_sectors(star_name, 72)

df_lc = load_light_curve(star_name, sectors[1], cut_size)

for sector in sectors[2:end]
    df_lc_sector = load_light_curve(star_name, sector, cut_size)
    df_lc = vcat(df_lc, df_lc_sector)
end

jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)
cleaned_mags = deepcopy(mags)
clean_flux_sigma!(jds, mags, jd_box, σ_tol, n_out)
cleaned_jds, cleaned_mags = delete_nans(jds, mags)

# plan = LombScargle.plan(cleaned_jds, cleaned_mags .- median(cleaned_mags))
# pgram = lombscargle(plan)

# freq, power = freqpower(pgram)
# periods = findmaxperiod(pgram, [0.01,100])

period = 35.35

phased_jds = jds .% period

fig = Figure()
ax = Axis(fig[1,1], yreversed = true)

scatter!(ax, phased_jds, mags)
fig