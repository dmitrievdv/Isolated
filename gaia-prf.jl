function create_gaia_prf_model(star_name, sector, cut_size, args...)
    create_gaia_prf_model(star_name, sector, cut_size, cut_size, args...)
end

function create_gaia_prf_model(star_name, sector, cut_width, cut_height, Δm_R)
    create_gaia_prf_model(star_name, sector, cut_width, cut_height, Δm_R, 0.0, 0.0)
end

function create_gaia_prf_model(star_name, sector, cut_width, cut_height, Δm_R, shift_x, shift_y)
    fits = load_tess_cutouts(star_name, cut_width, cut_height)[sector]

    supersampled_prf = get_tesscut_prf_supersampled(fits)
    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R)
    # n_gaia_stars = nrow(gaia_stars_data)
    gaia_data = load_star_gaia_data(star_name)

    # star_index = findfirst(s -> s == gaia_data.source_id, gaia_stars_data.source_id)
    # star_px = gaia_stars_data.px_x[star_index], gaia_stars_data.px_y[star_index]

    stars_x = gaia_stars_data.px_x; stars_y = gaia_stars_data.px_y; stars_m_R = gaia_stars_data.phot_rp_mean_mag

    create_gaia_prf_model(supersampled_prf, cut_width, cut_height, stars_x, stars_y, stars_m_R, shift_x, shift_y)
end

function create_gaia_prf_model(supersampled_prf :: AbstractMatrix, cut_width, cut_height, stars_x, stars_y, stars_m_R, shift_x, shift_y)
    n_gaia_stars = length(stars_m_R)
    
    model_cut = zeros(cut_width, cut_height)
    for i_star = 1:n_gaia_stars
        prf_cut = get_prf_cut(supersampled_prf, cut_width, cut_height, stars_x[i_star] + shift_x, stars_y[i_star] + shift_y)
        model_cut += prf_cut*calc_tess_flux_from_mag(stars_m_R[i_star])
    end

    return model_cut
end

function plot_gaia_prf_model(star_name, sector, cut_size, bkg_mod_q, bkg_cut_q, Δm_R)
    fits = load_tess_cutouts(star_name, cut_size, cut_size)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    mjds = read(fits[2], "TIME")

    # fits = load_tess_cutouts(star_name, cut_width, cut_height)[sector]

    supersampled_prf = get_tesscut_prf_supersampled(fits)
    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R)
    # n_gaia_stars = nrow(gaia_stars_data)
    # gaia_data = load_star_gaia_data(star_name)

    # star_index = findfirst(s -> s == gaia_data.source_id, gaia_stars_data.source_id)
    # star_px = gaia_stars_data.px_x[star_index], gaia_stars_data.px_y[star_index]

    stars_x = gaia_stars_data.px_x; stars_y = gaia_stars_data.px_y; stars_m_R = gaia_stars_data.phot_rp_mean_mag
    model = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, 0.0, 0.0)

    bkg_model_cut = sort(vec(model))[round(Int, n_cuts * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)

    # n_cut = 200
    fig = Figure()
    ax_cut = Axis(fig[0,1], aspect = DataAspect())
    ax_bkgcut = Axis(fig[1,1], aspect = DataAspect())
    ax_mod = Axis(fig[1,2], aspect = DataAspect())
    ax_rsd = Axis(fig[1,3], aspect = DataAspect())

    i_cut = Observable(500)

    cut_slider = Slider(fig[2, 1:3], range = 1:n_cuts, startvalue = 500)
    next_button = Button(fig[3,1], label = "Next", tellwidth = false)
    prev_button = Button(fig[3,3], label = "Prev", tellwidth = false)
    cut_slider_label = Label(fig[3, 2], tellwidth = false, text = @lift @sprintf("MJD = %.4f, i_cut = %d", mjds[$i_cut], $i_cut))

    on(next_button.clicks) do n
        i_cut[] += 1
        i_cut[] = (i_cut[] - 1) % n_cuts + 1
    end

    on(prev_button.clicks) do n
        i_cut[] -= 1
        i_cut[] = (i_cut[] - 1) % n_cuts + 1
    end

    on(cut_slider.value) do val
        i_cut[] = val
    end

    flux_cut = @lift flux_cuts[:,:,$i_cut]

    cut_data = lift(flux_cut) do flux_cut
        log10.(abs.(flux_cut))
    end

    

    bkg_cut_pixels = lift(flux_cut) do flux_cut
        bkg_cut_cut = sort(flux_cut[bkg_mod_pixels])[round(Int, length(bkg_mod_pixels)*bkg_cut_q)]
        bkg_mod_pixels[findall(x -> flux_cut[x] < bkg_cut_cut, bkg_mod_pixels)]
    end

    bkg_cut_scatter = lift(bkg_cut_pixels) do bkg_cut_pixels
        [Point2f(index[1], index[2]) for index in bkg_cut_pixels]
    end

    bkg_cut = lift(bkg_cut_pixels, flux_cut) do bkg_cut_pixels, flux_cut
        # flux_cut = flux_cuts[:,:,i_cut]
        fit_flat_background(flux_cut, bkg_cut_pixels)
    end

    bkg_cut_data = lift(bkg_cut, flux_cut) do bkg_cut, flux_cut
        log10.(abs.(flux_cut - bkg_cut))
    end

    rsd_data = lift(bkg_cut, flux_cut) do bkg_cut, flux_cut
        tmp_rsd = (flux_cut - bkg_cut - model) ./ model .- 1.0
        tmp_rsd[bkg_mod_pixels] .= 0.0
        # tmp_rsd[(flux_cut ./ bkg_cut) .< 1.2] .= 1.0
        tmp_rsd[1:5,:] .= 0.0
        tmp_rsd[cut_size-5+1:cut_size,:] .= 0.0
        tmp_rsd[:,1:5] .= 0.0
        tmp_rsd[:,cut_size-5+1:cut_size] .= 0.0
        tmp_rsd .^ 2
    end

    shift_coords = lift(bkg_cut, flux_cut) do bkg_cut, flux_cut

        function f(shift_coords)
            model = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, shift_coords...)
            tmp_rsd = (flux_cut - bkg_cut - model) ./ model .- 1.0
            tmp_rsd[bkg_mod_pixels] .= 0.0
            # tmp_rsd[(flux_cut ./ bkg_cut) .< 1.2] .= 1.0
            tmp_rsd[1:5,:] .= 0.0
            tmp_rsd[cut_size-5+1:cut_size,:] .= 0.0
            tmp_rsd[:,1:5] .= 0.0
            tmp_rsd[:,cut_size-5+1:cut_size] .= 0.0
            # tmp_rsd .^ 2
            s = sum(tmp_rsd .^ 2)
            println(s)
            return s
        end

        res = Opt.optimize(f, [0.0,0.0])
        println(res.minimizer)
        res.minimizer
    end

    heatmap!(ax_cut, cut_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_cut, bkg_cut_scatter, color = :red, marker = :xcross)
    heatmap!(ax_bkgcut, bkg_cut_data, colorrange = (0, log10(maximum(model))))
    heatmap!(ax_mod, log10.(abs.(model)))
    scatter!(ax_mod, [index[1] for index in bkg_mod_pixels], [index[2] for index in bkg_mod_pixels]; 
            color = :red, marker = :xcross)
    heatmap!(ax_rsd, rsd_data)

    fig
end