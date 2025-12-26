include("main.jl")
include("eclipse-light-curve.jl")
import Optim 
using Polynomials

star_name = "V501 Aur"
cut_size = 30

jd_box = 0.1; σ_tol = 5; n_out = 20

sectors = find_tess_sectors(star_name, 72)

df_lcs = Dict([sector => load_light_curve(star_name, sector, cut_size, rewrite_file = false) for sector in sectors])

df_lc_all = foldl(vcat, [df_lcs[sector] for sector in sort(sectors)])

jds_all, mags_all = delete_nans(get_true_jd.(df_lc_all.MJD), df_lc_all.MAG)
cleaned_mags_all = deepcopy(mags_all)
clean_flux_sigma!(jds_all, mags_all, jd_box, σ_tol, n_out)
cleaned_jds_all, cleaned_mags_all = delete_nans(jds_all, mags_all)

plan = LombScargle.plan(cleaned_jds_all, cleaned_mags_all .- median(cleaned_mags_all))
pgram = lombscargle(plan)

freq, power = freqpower(pgram)
periods = findmaxperiod(pgram, [0.01,100], 0.85)

mean_mag = calc_tess_magniude(mean(calc_tess_flux_from_mag.(cleaned_mags_all)))

jds_by_sector = Dict([sector => get_true_jd.(df_lc.MJD) for (sector, df_lc) in df_lcs])
mags_by_sector = Dict([sector => df_lc.MAG for (sector, df_lc) in df_lcs])
transits_phases = Vector{Float64}[]
transits_mags = Vector{Float64}[]
transits_sectors = Int[]
left_fits = []
right_fits = []
left_polys = []
right_polys = []
left_luminosity_relation = Float64[]
right_luminosity_relation = Float64[]

orbital_period = 68.8333

in_transit_window(ϕ) = ((ϕ ≥ 0.65) & (ϕ ≤ 0.75))
in_full_transit(ϕ) = ((ϕ ≥ 0.67) & (ϕ ≤ 0.73))
in_partial_or_full_transit(ϕ) = (((ϕ ≥ 0.66) & (ϕ ≤ 0.74)))
out_of_transit(ϕ) = (in_transit_window(ϕ) & !in_partial_or_full_transit(ϕ))
in_partial_transit(ϕ) = (!in_full_transit(ϕ) & in_partial_or_full_transit(ϕ))
in_left_partial_transit(ϕ) = in_partial_transit(ϕ) & ϕ ≤ 0.7
in_right_partial_transit(ϕ) = in_partial_transit(ϕ) & ϕ ≥ 0.7
in_left_fit_window(ϕ) = in_transit_window(ϕ) & (ϕ ≤ 0.68)# & !in_partial_transit(ϕ)
in_right_fit_window(ϕ) = in_transit_window(ϕ) & (ϕ ≥ 0.72)# & !in_partial_transit(ϕ)

calc_luminosity_relation(Δm, mean_mag, transit_mag) = (10^(0.4*Δm) - 1)*10^(0.4*(mean_mag - transit_mag))

n_poly = 5

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

    mean_transit_mag = calc_tess_magniude(mean(calc_tess_flux_from_mag.(full_transit_mags)))

    poly = fit(phases[indices_for_fit], cleaned_mags[indices_for_fit], n_poly)

    lower = [fill(-Inf, n_poly+1); [0.01, 0.0, 0.68, 10, 50]]
    start_fit_pars = [coeffs(poly); [0.3, 0.001, 0.7, 20, 100]]
    upper = [fill(Inf, n_poly+1); [0.5,0.5,0.72,40,200]]
    inner_optimizer = Optim.LBFGS()

    left_fit_phases = phases[left_fit_indices]; left_fit_mags = cleaned_mags[left_fit_indices]
    left_full_mask = [in_full_transit(ϕ) for ϕ in left_fit_phases]
    #left_fitting_function(p) = left_fit_mags - (Polynomial(p[1:end-1]).(left_fit_phases) .+ p[end] .* left_full_mask)                         
    
    left_fitting_function(p) = left_fit_mags - calc_transit_mag.(Polynomial(p[1:n_poly+1]).(left_fit_phases), p[end], p[end-1], left_fit_phases, p[end-2], p[end-3], p[end-4])
    function lsq_left_fitting_function(p) 
        δ = sum( left_fitting_function(p) .^ 2 )
        # print(p)
        # println(": $δ")
        return δ
    end

    left_fit_pars = Optim.optimize(lsq_left_fitting_function, lower, upper, start_fit_pars, Optim.Fminbox(inner_optimizer)).minimizer 
    println("left $sector: $left_fit_pars")

    right_fit_phases = phases[right_fit_indices]; right_fit_mags = cleaned_mags[right_fit_indices]
    right_full_mask = [in_full_transit(ϕ) for ϕ in right_fit_phases]

    right_fitting_function(p) = right_fit_mags - calc_transit_mag.(Polynomial(p[1:n_poly+1]).(right_fit_phases), p[end], p[end-1], right_fit_phases, p[end-2], p[end-3], p[end-4])
    lsq_right_fitting_function(p) = sum( right_fitting_function(p) .^ 2 )
    # right_fitting_function(p) = right_fit_mags - (Polynomial(p[1:end-1]).(right_fit_phases) .+ p[end] .* right_full_mask)                         
    right_fit_pars = Optim.optimize(lsq_right_fitting_function, lower, upper, start_fit_pars, Optim.Fminbox(inner_optimizer)).minimizer 
    println("right $sector: $right_fit_pars")

    push!(left_luminosity_relation, calc_luminosity_relation(left_fit_pars[end], mean_mag, mean_transit_mag))
    push!(right_luminosity_relation, calc_luminosity_relation(right_fit_pars[end], mean_mag, mean_transit_mag))
    println(sector, " ", left_luminosity_relation[end], " ", right_luminosity_relation[end])    
    push!(left_fits, (left_fit_phases, (-left_fitting_function(left_fit_pars) .+ left_fit_mags)))
    push!(right_fits, (right_fit_phases, (-right_fitting_function(right_fit_pars) .+ right_fit_mags)))

    push!(left_polys, (left_fit_phases, Polynomial(left_fit_pars[1:n_poly+1]).(left_fit_phases)))
    push!(right_polys, (right_fit_phases, Polynomial(right_fit_pars[1:n_poly+1]).(right_fit_phases)))
    push!(transits_mags, [left_fit_mags; right_fit_mags])
    push!(transits_phases, [left_fit_phases; right_fit_phases])
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
    lines!(ax, left_polys[i_transit][1], left_polys[i_transit][2], color = i_transit, colormap = :tab10, colorrange = (1,10), linestyle = :dash)
    lines!(ax, right_fits[i_transit][1], right_fits[i_transit][2], color = i_transit, colormap = :tab10, colorrange = (1,10))
    lines!(ax, right_polys[i_transit][1], right_polys[i_transit][2], color = i_transit, colormap = :tab10, colorrange = (1,10), linestyle = :dash)
end# lines!(ax_pgram, freq, power)

axislegend(ax)
fig