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
        add_prf_cut!(model_cut, calc_tess_flux_from_mag(stars_m_R[i_star]), supersampled_prf, cut_width, cut_height, stars_x[i_star] + shift_x, stars_y[i_star] + shift_y)
    end

    return model_cut
end

function create_gaia_prf_model_old(supersampled_prf :: AbstractMatrix, cut_width, cut_height, stars_x, stars_y, stars_m_R, shift_x, shift_y)
    n_gaia_stars = length(stars_m_R)
    
    model_cut = zeros(cut_width, cut_height)
    for i_star = 1:n_gaia_stars
        prf_cut = get_prf_cut(supersampled_prf, cut_width, cut_height, stars_x[i_star] + shift_x, stars_y[i_star] + shift_y)
        model_cut += prf_cut*calc_tess_flux_from_mag(stars_m_R[i_star])
    end

    return model_cut
end

function get_shift_coords(star_name, sector, cut_size, bkg_mod_q, bkg_cut_q, Δm_R)
    fits = load_tess_cutouts(star_name, cut_size, cut_size)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    mjds = read(fits[2], "TIME")

    supersampled_prf = get_tesscut_prf_supersampled(fits)
    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R)

    stars_x = gaia_stars_data.px_x; stars_y = gaia_stars_data.px_y; stars_m_R = gaia_stars_data.phot_rp_mean_mag
    model = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, 0.0, 0.0)

    bkg_model_cut = sort(vec(model))[round(Int, n_cuts * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)
    rsd_window_pixels = filter(findall(x -> true, model)) do x
        (5 < x[1] ≤ cut_size - 5) & (5 < x[2] ≤ cut_size - 5)
    end

    shift_coords = zeros(n_cuts, 2)
    start_shift_coords = [0.0, 0.0]
    for i_cut = 1:n_cuts
        flux_cut = flux_cuts[:,:,i_cut]
        bkg_cut_cut = sort(flux_cut[bkg_mod_pixels])[round(Int, length(bkg_mod_pixels)*bkg_cut_q)]
        bkg_cut_pixels = bkg_mod_pixels[findall(x -> flux_cut[x] < bkg_cut_cut, bkg_mod_pixels)]

        bkg_cut = fit_flat_background(flux_cut, bkg_cut_pixels)

        rsd_fit_pixels = rsd_window_pixels #filter(rsd_window_pixels) do x
            # model[x] > 0.5*bkg_cut[x]
        # end
        
        if length(rsd_fit_pixels) < 30
            continue
        end
        function f(shift_coords)
            model = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, shift_coords...)
            s = 0.0
            for px in rsd_fit_pixels
                s += sqrt(abs((flux_cut[px] - bkg_cut[px]) .* model[px]))
            end
            
            # tmp_rsd[bkg_mod_pixels] .= 0.0
            # # tmp_rsd[(flux_cut ./ bkg_cut) .< 1.2] .= 1.0
            # tmp_rsd[1:5,:] .= 0.0
            # tmp_rsd[cut_size-5+1:cut_size,:] .= 0.0
            # tmp_rsd[:,1:5] .= 0.0
            # tmp_rsd[:,cut_size-5+1:cut_size] .= 0.0
            # # tmp_rsd .^ 2
            # s = sum(tmp_rsd .^ 2)
            # println(s)
            return 1/s
        end

            for px in rsd_fit_pixels
                s += sqrt(abs((flux_cut[px] - bkg_cut[px]) .* model[px]))
            end

            res = Opt.optimize(f, start_shift_coords, Opt.Options(x_abstol = 1e-3))

            # printing using the ANSI escape codes:
            # \e[2K clears the entire current line
            # \e[1G moves the cursor to the first column
            print("\e[2K\e[1G$i_cut from $n_cuts; $(length(rsd_fit_pixels)) $(res.minimizer)",)

            start_shift_coords = res.minimizer
        end
        shift_coords[i_cut, :] = start_shift_coords

        # set_postfix(iter, Coords=@sprintf("%7.4f, %7.4f", shift_coords[i_cut, 1], shift_coords[i_cut, 2]))
    end
    print("\n")
    return shift_coords
end

function calc_shifted_light_curve(star_name, sector, cut_size, Δm_R, aperture_radius, shift_coords)
    fits = load_tess_cutouts(star_name, cut_size, cut_size)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    mjds = read(fits[2], "TIME")

    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R, rewrite_file = false)
    gaia_data = load_star_gaia_data(star_name)

    star_index = findfirst(s -> s == gaia_data.source_id, gaia_stars_data.source_id)
    star_px = gaia_stars_data.px_x[star_index], gaia_stars_data.px_y[star_index]

    prf = get_tesscut_prf_supersampled(fits)
    bkg_pixels = find_background_prf(flux_cuts[:,:,n_cuts÷4], prf, gaia_stars_data.px_x, gaia_stars_data.px_y)

    cuts = [flux_cuts[:,:,i_cut] for i_cut = 1:n_cuts]
    aperture_prf_correction = calc_aperture_prf_correction.(aperture_radius, star_px[1], star_px[2], shift_coords[:,1], shift_coords[:,2], Ref(prf), cut_size)
    phot_flux = aperture_prf_correction .* calc_aperture_photometry_bkg_pixels.(cuts, Ref(bkg_pixels), star_px..., aperture_radius)

    return mjds, phot_flux
end

function test_background_interpolation(star_name, sector, cut_size, bkg_mod_q, Δm_R)
    fits = load_tess_cutouts(star_name, cut_size, cut_size)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    mjds = read(fits[2], "TIME")

    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R)
    stars_x = gaia_stars_data.px_x; stars_y = gaia_stars_data.px_y; stars_m_R = gaia_stars_data.phot_rp_mean_mag

    supersampled_prf = get_tesscut_prf_supersampled(fits)
    model = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, 0.0, 0.0)

    bkg_model_cut = sort(vec(model))[round(Int, (cut_size^2) * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)
    bkg_cut_scatter = [Point2f(index[1], index[2]) for index in bkg_mod_pixels]

    # println(bkg_mod_pixels)

    fig = Figure()
    ax_cut = Axis(fig[1,1], aspect = DataAspect())
    ax_bkg = Axis(fig[1,2], aspect = DataAspect())
    ax_bkgcut = Axis(fig[1,3], aspect = DataAspect())

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

    bkg_cut = lift(flux_cut) do flux_cut
        # flux_cut = flux_cuts[:,:,i_cut]
        fit_flat_background(flux_cut, bkg_mod_pixels)
    end

    bkg_cut_data = lift(bkg_cut) do bkg_cut
        log10.(abs.(bkg_cut))
    end

    flux_bkg_data = lift(bkg_cut) do bkg_cut
        log10.(abs.(flux_cut.val - bkg_cut))
    end

    heatmap!(ax_cut, cut_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_cut, bkg_cut_scatter, color = :red, marker = :xcross)
    heatmap!(ax_bkgcut, flux_bkg_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_bkgcut, bkg_cut_scatter, color = :red, marker = :xcross)
    heatmap!(ax_bkg, bkg_cut_data, colorrange = (0, log10(maximum(model))))

    fig
end

function test_gaia_prf_model_fit(star_name, sector, cut_size, bkg_mod_q, bkg_cut_q, Δm_R)
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
    # model2 = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, 0.0, 0.0)

    # println(model - model2)

    bkg_model_cut = sort(vec(model))[round(Int, n_cuts * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)

    rsd_window_pixels = filter(findall(x -> true, model)) do x
        (5 < x[1] ≤ cut_size - 5) & (5 < x[2] ≤ cut_size - 5)
    end

    fig = Figure()
    ax_mod = Axis(fig[1,1], aspect = DataAspect())
    ax_cut = Axis(fig[1,2], aspect = DataAspect())
    ax_rsd = Axis(fig[1,3], aspect = DataAspect())

    bkg_model_cut = sort(vec(model))[round(Int, n_cuts * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)

    rsd_window_pixels = filter(findall(x -> true, model)) do x
        (5 < x[1] ≤ cut_size - 5) & (5 < x[2] ≤ cut_size - 5)
    end

    i_cut = Observable(500)
    # rsd = Observable(0.0)

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

    control_x_slider = Slider(fig[0, 1], range = -1:0.01:1, startvalue = 0, tellwidth = false)
    control_y_slider = Slider(fig[1, 0], range = -1:0.01:1, startvalue = 0, tellheight = false, horizontal = false)

    shift_coords = lift(control_x_slider.value, control_y_slider.value) do x, y
        x, y
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

    bkg_cut = lift(bkg_cut_pixels) do bkg_cut_pixels
        # flux_cut = flux_cuts[:,:,i_cut]
        fit_flat_background(flux_cut.val, bkg_cut_pixels)
    end

    bkg_cut_data = lift(bkg_cut) do bkg_cut
        log10.(abs.(flux_cut.val - bkg_cut))
    end

   

    shifted_model = lift(shift_coords) do shift_coords
        create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, shift_coords...)
    end

    shifted_model_data = lift(shifted_model) do shifted_model
        log10.(abs.(shifted_model))
    end

     rsd_fit_pixels = lift(bkg_cut, shifted_model) do bkg_cut, shifted_model
        # filter(x -> shifted_model[x] ≥ 0.3*bkg_cut[x], rsd_window_pixels)
        rsd_window_pixels
    end

    rsd_fit_scatter = lift(rsd_fit_pixels) do rsd_fit_pixels
        [Point2f(index[1], index[2]) for index in rsd_fit_pixels]
    end

    shifted_rsd_data = lift(shifted_model, bkg_cut) do model, bkg_cut
        tmp_rsd = ((flux_cut.val - bkg_cut) .* model)
        rsd = zeros(cut_size, cut_size)
        tmp_rsd[bkg_mod_pixels] .= 0.0
        # tmp_rsd[(flux_cut ./ bkg_cut) .< 1.2] .= 1.0
        tmp_rsd[1:5,:] .= 0.0
        tmp_rsd[cut_size-5+1:cut_size,:] .= 0.0
        tmp_rsd[:,1:5] .= 0.0
        tmp_rsd[:,cut_size-5+1:cut_size] .= 0.0
        # tmp_rsd .^ 2
        rsd[rsd_fit_pixels.val] = tmp_rsd[rsd_fit_pixels.val]
        log10.((sqrt.(abs.(rsd))))
    end

    rsd = lift(shifted_rsd_data) do shifted_rsd_data
        1 / (sum(10 .^ shifted_rsd_data) / cut_size^2)
    end

    rsd_label = Label(fig[0,3], tellwidth = false, text = @lift @sprintf("RSD = %.6f", $rsd))

    heatmap!(ax_cut, bkg_cut_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_cut, bkg_cut_scatter, color = :red, marker = :xcross)
    # heatmap!(ax_mod, log10.(abs.(model)))
    heatmap!(ax_mod, shifted_model_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_mod, [index[1] for index in bkg_mod_pixels], [index[2] for index in bkg_mod_pixels]; 
            color = :red, marker = :xcross)
    heatmap!(ax_rsd, shifted_rsd_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_rsd, rsd_fit_scatter; color = :red, marker = :xcross)
    fig
end

function test_background_interpolation(star_name, sector, cut_size, bkg_mod_q, Δm_R)
    fits = load_tess_cutouts(star_name, cut_size, cut_size)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    mjds = read(fits[2], "TIME")

    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R)
    stars_x = gaia_stars_data.px_x; stars_y = gaia_stars_data.px_y; stars_m_R = gaia_stars_data.phot_rp_mean_mag

    supersampled_prf = get_tesscut_prf_supersampled(fits)
    model = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, 0.0, 0.0)

    bkg_model_cut = sort(vec(model))[round(Int, (cut_size^2) * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)
    bkg_cut_scatter = [Point2f(index[1], index[2]) for index in bkg_mod_pixels]

    # println(bkg_mod_pixels)

    fig = Figure()
    ax_cut = Axis(fig[1,1], aspect = DataAspect())
    ax_bkg = Axis(fig[1,2], aspect = DataAspect())
    ax_bkgcut = Axis(fig[1,3], aspect = DataAspect())

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

    bkg_cut = lift(flux_cut) do flux_cut
        # flux_cut = flux_cuts[:,:,i_cut]
        fit_flat_background(flux_cut, bkg_mod_pixels)
    end

    bkg_cut_data = lift(bkg_cut) do bkg_cut
        log10.(abs.(bkg_cut))
    end

    flux_bkg_data = lift(bkg_cut) do bkg_cut
        log10.(abs.(flux_cut.val - bkg_cut))
    end

    heatmap!(ax_cut, cut_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_cut, bkg_cut_scatter, color = :red, marker = :xcross)
    heatmap!(ax_bkgcut, flux_bkg_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_bkgcut, bkg_cut_scatter, color = :red, marker = :xcross)
    heatmap!(ax_bkg, bkg_cut_data, colorrange = (0, log10(maximum(model))))

    fig
end

function test_gaia_prf_model_fit(star_name, sector, cut_size, bkg_mod_q, bkg_cut_q, Δm_R)
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
    # model2 = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, 0.0, 0.0)

    # println(model - model2)

    bkg_model_cut = sort(vec(model))[round(Int, n_cuts * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)

    rsd_window_pixels = filter(findall(x -> true, model)) do x
        (5 < x[1] ≤ cut_size - 5) & (5 < x[2] ≤ cut_size - 5)
    end

    fig = Figure()
    ax_mod = Axis(fig[1,1], aspect = DataAspect())
    ax_cut = Axis(fig[1,2], aspect = DataAspect())
    ax_rsd = Axis(fig[1,3], aspect = DataAspect())

    bkg_model_cut = sort(vec(model))[round(Int, n_cuts * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)

    rsd_window_pixels = filter(findall(x -> true, model)) do x
        (5 < x[1] ≤ cut_size - 5) & (5 < x[2] ≤ cut_size - 5)
    end

    i_cut = Observable(500)
    # rsd = Observable(0.0)

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

    control_x_slider = Slider(fig[0, 1], range = -1:0.01:1, startvalue = 0, tellwidth = false)
    control_y_slider = Slider(fig[1, 0], range = -1:0.01:1, startvalue = 0, tellheight = false, horizontal = false)

    shift_coords = lift(control_x_slider.value, control_y_slider.value) do x, y
        x, y
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

    bkg_cut = lift(bkg_cut_pixels) do bkg_cut_pixels
        # flux_cut = flux_cuts[:,:,i_cut]
        fit_flat_background(flux_cut.val, bkg_cut_pixels)
    end

    bkg_cut_data = lift(bkg_cut) do bkg_cut
        log10.(abs.(flux_cut.val - bkg_cut))
    end

   

    shifted_model = lift(shift_coords) do shift_coords
        create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, shift_coords...)
    end

    shifted_model_data = lift(shifted_model) do shifted_model
        log10.(abs.(shifted_model))
    end

     rsd_fit_pixels = lift(bkg_cut, shifted_model) do bkg_cut, shifted_model
        # filter(x -> shifted_model[x] ≥ 0.3*bkg_cut[x], rsd_window_pixels)
        rsd_window_pixels
    end

    rsd_fit_scatter = lift(rsd_fit_pixels) do rsd_fit_pixels
        [Point2f(index[1], index[2]) for index in rsd_fit_pixels]
    end

    shifted_rsd_data = lift(shifted_model, bkg_cut) do model, bkg_cut
        tmp_rsd = ((flux_cut.val - bkg_cut) .* model)
        rsd = zeros(cut_size, cut_size)
        tmp_rsd[bkg_mod_pixels] .= 0.0
        # tmp_rsd[(flux_cut ./ bkg_cut) .< 1.2] .= 1.0
        tmp_rsd[1:5,:] .= 0.0
        tmp_rsd[cut_size-5+1:cut_size,:] .= 0.0
        tmp_rsd[:,1:5] .= 0.0
        tmp_rsd[:,cut_size-5+1:cut_size] .= 0.0
        # tmp_rsd .^ 2
        rsd[rsd_fit_pixels.val] = tmp_rsd[rsd_fit_pixels.val]
        log10.((sqrt.(abs.(rsd))))
    end

    rsd = lift(shifted_rsd_data) do shifted_rsd_data
        1 / (sum(10 .^ shifted_rsd_data) / cut_size^2)
    end

    rsd_label = Label(fig[0,3], tellwidth = false, text = @lift @sprintf("RSD = %.6f", $rsd))

    heatmap!(ax_cut, bkg_cut_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_cut, bkg_cut_scatter, color = :red, marker = :xcross)
    # heatmap!(ax_mod, log10.(abs.(model)))
    heatmap!(ax_mod, shifted_model_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_mod, [index[1] for index in bkg_mod_pixels], [index[2] for index in bkg_mod_pixels]; 
            color = :red, marker = :xcross)
    heatmap!(ax_rsd, shifted_rsd_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_rsd, rsd_fit_scatter; color = :red, marker = :xcross)
    fig
end

function plot_gaia_prf_model(star_name, sector, cut_size, bkg_mod_q, bkg_cut_q, Δm_R)
    fits = load_tess_cutouts(star_name, cut_size, cut_size)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    bkg_cuts = zeros(cut_size, cut_size, n_cuts)
    shifted_models = zeros(cut_size, cut_size, n_cuts)

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
    # model2 = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, 0.0, 0.0)

    # println(model - model2)

    bkg_model_cut = sort(vec(model))[round(Int, n_cuts * bkg_mod_q)]
    bkg_mod_pixels = findall(x -> x < bkg_model_cut, model)

    rsd_window_pixels = filter(findall(x -> true, model)) do x
        (5 < x[1] ≤ cut_size - 5) & (5 < x[2] ≤ cut_size - 5)
    end

    # n_cut = 200
    fig = Figure()
    ax_cut = Axis(fig[0,1], aspect = DataAspect())
    ax_bkgcut = Axis(fig[1,1], aspect = DataAspect())
    ax_mod = Axis(fig[1,2], aspect = DataAspect())
    ax_rsd = Axis(fig[1,3], aspect = DataAspect())
    ax_shftrsd = Axis(fig[0,3], aspect = DataAspect())
    ax_shftmod = Axis(fig[0,2], aspect = DataAspect())

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

    bkg_cut = lift(bkg_cut_pixels) do bkg_cut_pixels
        # flux_cut = flux_cuts[:,:,i_cut]
        fit_flat_background(flux_cut.val, bkg_cut_pixels)
    end

    bkg_cut_data = lift(bkg_cut) do bkg_cut
        log10.(abs.(flux_cut.val - bkg_cut))
    end

    rsd_fit_pixels = lift(bkg_cut) do bkg_cut
        filter(x -> model[x] ≥ 0.3*bkg_cut[x], rsd_window_pixels)
    end

    rsd_fit_scatter = lift(rsd_fit_pixels) do rsd_fit_pixels
        [Point2f(index[1], index[2]) for index in rsd_fit_pixels]
    end

    rsd_data = lift(bkg_cut) do bkg_cut
        tmp_rsd = (flux_cut.val - bkg_cut - model) ./ model .- 1.0
        tmp_rsd[bkg_mod_pixels] .= 0.0
        # tmp_rsd[(flux_cut ./ bkg_cut) .< 1.2] .= 1.0
        tmp_rsd[1:5,:] .= 0.0
        tmp_rsd[cut_size-5+1:cut_size,:] .= 0.0
        tmp_rsd[:,1:5] .= 0.0
        tmp_rsd[:,cut_size-5+1:cut_size] .= 0.0
        tmp_rsd .^ 2
    end

    shift_coords = lift(bkg_cut) do bkg_cut

        function f(shift_coords)
            model = create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, shift_coords...)
            s = 0.0
            for px in rsd_fit_pixels.val
                s += flux_cut.val[px] .* model[px]
            end
            return 1/s
        end

        res = Opt.optimize(f, [0.0,0.0], Opt.LBFGS(), Opt.Options(x_abstol = 1e-3))
        println(res.minimizer)
        res.minimizer
    end

    shifted_model = lift(shift_coords) do shift_coords
        create_gaia_prf_model(supersampled_prf, cut_size, cut_size, stars_x, stars_y, stars_m_R, shift_coords...)
    end

    shifted_model_data = lift(shifted_model) do shifted_model
        log10.(abs.(shifted_model))
    end

    shifted_rsd_data = lift(shifted_model) do model
        tmp_rsd = (flux_cut.val - bkg_cut.val - model) ./ model .- 1.0
        tmp_rsd[bkg_mod_pixels] .= 0.0
        # tmp_rsd[(flux_cut ./ bkg_cut) .< 1.2] .= 1.0
        tmp_rsd[1:5,:] .= 0.0
        tmp_rsd[cut_size-5+1:cut_size,:] .= 0.0
        tmp_rsd[:,1:5] .= 0.0
        tmp_rsd[:,cut_size-5+1:cut_size] .= 0.0
        tmp_rsd .^ 2
    end

    heatmap!(ax_cut, cut_data, colorrange = (0, log10(maximum(model))))
    scatter!(ax_cut, bkg_cut_scatter, color = :red, marker = :xcross)
    heatmap!(ax_bkgcut, bkg_cut_data, colorrange = (0, log10(maximum(model))))
    heatmap!(ax_mod, log10.(abs.(model)))
    scatter!(ax_mod, [index[1] for index in bkg_mod_pixels], [index[2] for index in bkg_mod_pixels]; 
            color = :red, marker = :xcross)
    heatmap!(ax_shftmod, shifted_model_data)
    heatmap!(ax_rsd, rsd_data)
    scatter!(ax_rsd, rsd_fit_scatter; color = :red, marker = :xcross)
    heatmap!(ax_shftrsd, shifted_rsd_data)
    fig
end