function create_gaia_prf_model(star_name, sector, cut_size, Δm_R = 5)
    create_gaia_prf_model(star_name, sector, cut_size, cut_size, Δm_R)
end

function create_gaia_prf_model(star_name, sector, cut_width, cut_height, Δm_R)
    fits = load_tess_cutouts(star_name, cut_width, cut_height)[sector]

    supersampled_prf = get_tesscut_prf_supersampled(fits)
    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R)
    n_gaia_stars = nrow(gaia_stars_data)
    gaia_data = load_star_gaia_data(star_name)

    star_index = findfirst(s -> s == gaia_data.source_id, gaia_stars_data.source_id)
    star_px = gaia_stars_data.px_x[star_index], gaia_stars_data.px_y[star_index]

    model_cut = zeros(cut_width, cut_height)
    for i_star = 1:n_gaia_stars
        prf_cut = get_prf_cut(supersampled_prf, cut_width, cut_height, gaia_stars_data.px_x[i_star], gaia_stars_data.px_y[i_star] + 0.2)
        model_cut += prf_cut*calc_tess_flux_from_mag(gaia_stars_data.phot_rp_mean_mag[i_star])
    end

    return model_cut
end

function plot_gaia_prf_model(star_name, sector, cut_size, Δm_R = 5)
    fits = load_tess_cutouts(star_name, cut_size, cut_size)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    mjds = read(fits[2], "TIME")

    model = create_gaia_prf_model(star_name, sector, cut_size, Δm_R)

    n_cut = 200
    flux_cut = flux_cuts[:,:,n_cut]
    fig = Figure()
    ax_cut = Axis(fig[1,1], aspect = DataAspect())
    ax_mod = Axis(fig[1,2], aspect = DataAspect())
    ax_rsd = Axis(fig[1,3], aspect = DataAspect())

    heatmap!(ax_cut, log10.(abs.(flux_cut)))
    heatmap!(ax_mod, log10.(abs.(model)))
    heatmap!(ax_rsd, flux_cut - model)

    println(√(mean((flux_cut - model) .^ 2)))

    fig
end