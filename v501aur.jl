include("main.jl")
using Polynomials

star_name = "V501 Aur"
cut_size = 30

jd_box = 0.1; σ_tol = 5; n_out = 20

sectors = find_tess_sectors(star_name, 72)

df_lcs = Dict([sector => load_light_curve(star_name, sector, cut_size) for sector in sectors])

df_lc_all = foldl(vcat, [df_lcs[sector] for sector in sort(sectors)])

jds_all, mags_all = delete_nans(get_true_jd.(df_lc_all.MJD), df_lc_all.MAG)
cleaned_mags_all = deepcopy(mags_all)
clean_flux_sigma!(jds_all, mags_all, jd_box, σ_tol, n_out)
cleaned_jds_all, cleaned_mags_all = delete_nans(jds_all, mags_all)

plan = LombScargle.plan(cleaned_jds_all, cleaned_mags_all .- median(cleaned_mags_all))
pgram = lombscargle(plan)

freq, power = freqpower(pgram)
periods = findmaxperiod(pgram, [0.01,100], 0.85)

mean_mag = mean(cleaned_mags_all)

jds_by_sector = Dict([sector => get_true_jd.(df_lc.MJD) for (sector, df_lc) in df_lcs])
mags_by_sector = Dict([sector => df_lc.MAG for (sector, df_lc) in df_lcs])
transits_phases = Vector{Float64}[]
transits_mags = Vector{Float64}[]
transits_sectors = Int[]
left_fits = []
right_fits = []

orbital_period = 68.8333

in_transit_window(ϕ) = ((ϕ ≥ 0.65) & (ϕ ≤ 0.75))
in_full_transit(ϕ) = ((ϕ ≥ 0.67) & (ϕ ≤ 0.73))
in_partial_or_full_transit(ϕ) = (((ϕ ≥ 0.66) & (ϕ ≤ 0.74)))
out_of_transit(ϕ) = (in_transit_window(ϕ) & !in_partial_or_full_transit(ϕ))
in_partial_transit(ϕ) = (!in_full_transit(ϕ) & in_partial_or_full_transit(ϕ))
in_left_partial_transit(ϕ) = in_partial_transit(ϕ) & ϕ ≤ 0.7
in_right_partial_transit(ϕ) = in_partial_transit(ϕ) & ϕ ≥ 0.7
in_left_fit_window(ϕ) = in_transit_window(ϕ) & (ϕ ≤ 0.68) & !in_partial_transit(ϕ)
in_right_fit_window(ϕ) = in_transit_window(ϕ) & (ϕ ≥ 0.72) & !in_partial_transit(ϕ)

for sector in sectors
    jds, mags = delete_nans(get_true_jd.(df_lcs[sector].MJD), df_lcs[sector].MAG)
    cleaned_mags = deepcopy(mags)
    clean_flux_sigma!(jds, mags, jd_box, σ_tol, n_out)
    cleaned_jds, cleaned_mags = delete_nans(jds, mags)

    jds_by_sector[sector] = cleaned_jds; mags_by_sector[sector] = cleaned_mags

    phases = ((cleaned_jds .+ 0.2*orbital_period) .% orbital_period) / orbital_period

    transit_indices = findall(in_transit_window, phases)
    indices_for_fit = findall(out_of_transit, phases)
    full_transit_indices = findall(in_full_transit, phases)
    not_partial_transit_indices = findall(ϕ -> !in_partial_transit(ϕ) & in_transit_window(ϕ), phases)
    left_fit_indices = findall(in_left_fit_window, phases)
    right_fit_indices = findall(in_right_fit_window, phases)

    if isempty(transit_indices)
        continue
    end
    

    transit_mags = cleaned_mags[transit_indices]
    transit_phases = phases[transit_indices]
    not_partial_transit_phases = phases[not_partial_transit_indices]
    not_partial_transit_mags = cleaned_mags[not_partial_transit_indices]
    full_transit_phases = phases[full_transit_indices]
    full_transit_mags = cleaned_mags[full_transit_indices]

    mean_transit_mag = mean(transit_mags)

    poly = fit(phases[indices_for_fit], cleaned_mags[indices_for_fit], 4)
    start_fit_pars = [coeffs(poly); 0.01]

    left_fit_phases = phases[left_fit_indices]; left_fit_mags = cleaned_mags[left_fit_indices]
    left_full_mask = [in_full_transit(ϕ) for ϕ in left_fit_phases]
    left_fitting_function(p) = left_fit_mags - ((p[5]*left_fit_phases .^ 4 + p[4]*left_fit_phases .^ 3 + p[3]*left_fit_phases .^ 2 +
                                p[2]*left_fit_phases .+ p[1]) .+ p[6] .* left_full_mask)                         
    left_fit_pars = optimize(left_fitting_function, start_fit_pars, LevenbergMarquardt()).minimizer 

    right_fit_phases = phases[right_fit_indices]; right_fit_mags = cleaned_mags[right_fit_indices]
    right_full_mask = [in_full_transit(ϕ) for ϕ in right_fit_phases]
    right_fitting_function(p) = right_fit_mags - ((p[5]*right_fit_phases .^ 4 + p[4]*right_fit_phases .^ 3 + p[3]*right_fit_phases .^ 2 +
                                p[2]*right_fit_phases .+ p[1]) .+ p[6] .* right_full_mask)                         
    right_fit_pars = optimize(right_fitting_function, start_fit_pars, LevenbergMarquardt()).minimizer 

    println(left_fit_pars[end], " ", right_fit_pars[end])    
    push!(left_fits, (left_fit_phases, (-left_fitting_function(left_fit_pars) .+ left_fit_mags)))
    push!(right_fits, (right_fit_phases, (-right_fitting_function(right_fit_pars) .+ right_fit_mags)))
    push!(transits_mags, cleaned_mags[not_partial_transit_indices])
    push!(transits_phases, phases[not_partial_transit_indices])
    push!(transits_sectors, sector)
end

n_transits = length(transits_sectors)





# phased_jds = (jds .+ 0.2*period) .% period

fig = Figure()
ax = Axis(fig[1,1], yreversed = true, xticks = WilkinsonTicks(10; k_min = 8, k_max = 20))
# ax_pgram = Axis(fig[2,1])
# xlims!(ax, (-0.0, 1.0))

for i_transit = 1:n_transits
    sector = transits_sectors[i_transit]
    phases = transits_phases[i_transit]; mags = transits_mags[i_transit]
    scatter!(ax, phases, mags, markersize = 2, label = "$sector", color = i_transit, colormap = :tab10, colorrange = (1,10))
    lines!(ax, left_fits[i_transit][1], left_fits[i_transit][2], color = i_transit, colormap = :tab10, colorrange = (1,10))
    lines!(ax, right_fits[i_transit][1], right_fits[i_transit][2], color = i_transit, colormap = :tab10, colorrange = (1,10))
end# lines!(ax_pgram, freq, power)

axislegend(ax)
fig