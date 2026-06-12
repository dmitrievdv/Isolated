function calc_noise_rms(df_lc :: DataFrame; jd_box = 0.1, σ_tol = 10, n_out = 20)
    jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)
    cleaned_mags = deepcopy(mags)
    clean_flux_sigma!(jds, mags, jd_box, σ_tol, n_out)
    cleaned_jds, cleaned_mags = delete_nans(jds, mags)

    smooth_fluxs = box_smooth(cleaned_jds, cleaned_mags, jd_box)
    anomalous_fluxs = cleaned_mags .- smooth_fluxs

    noout_anomalous_fluxs = sort(anomalous_fluxs)[1:end-n_out]

    σ = √(varm(noout_anomalous_fluxs, mean(noout_anomalous_fluxs)))
end

function find_period(df_lc :: DataFrame; jd_box = 0.1, σ_tol = 10, n_out = 20)
    jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)
    cleaned_mags = deepcopy(mags)
    clean_flux_sigma!(jds, mags, jd_box, σ_tol, n_out)
    cleaned_jds, cleaned_mags = delete_nans(jds, mags)

    plan = LombScargle.plan(cleaned_jds, cleaned_mags .- median(cleaned_mags))
    pgram = lombscargle(plan)

    freq, power = freqpower(pgram)
    period = findmaxperiod(pgram, [0.08,10])[1]
end

function calc_periodicity(df_lc :: DataFrame; jd_box = 0.1, σ_tol = 10, n_out = 20)
    jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)
    cleaned_mags = deepcopy(mags)
    clean_flux_sigma!(jds, mags, jd_box, σ_tol, n_out)
    cleaned_jds, cleaned_mags = delete_nans(jds, mags)

    plan = LombScargle.plan(cleaned_jds, cleaned_mags .- median(cleaned_mags))
    pgram = lombscargle(plan)

    freq, power = freqpower(pgram)
    period = findmaxperiod(pgram, [0.08,10])[1]
    phased_jds = cleaned_jds .% period

    sampling = find_sampling(cleaned_jds)
    jd_boxes = [0:2*sampling:period;]
    n_boxes = length(jd_boxes)

    n_mags = length(cleaned_mags)
    boxed_mags = zeros(n_mags)

    for i_mag = 1:n_mags
        boxed_mags[i_mag] = median(cleaned_mags[@. abs(phased_jds - phased_jds[i_mag]) < period*0.01])
    end

    mean_mag = mean(cleaned_mags)
    mag_dispersion = √(varm(cleaned_mags, mean_mag))
    mean_phased = mean(cleaned_mags - boxed_mags)
    mag_phased_dispersion = √(varm(cleaned_mags - boxed_mags, mean_phased))

    noise = calc_noise_rms(df_lc)

    periodicity = (mag_phased_dispersion - noise)/(mag_dispersion - noise)
end

function calc_asymmetry(df_lc :: DataFrame; jd_box = 0.1, σ_tol = 10, n_out = 20)
    jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)
    cleaned_mags = deepcopy(mags)
    clean_flux_sigma!(jds, mags, jd_box, σ_tol, n_out)
    cleaned_jds, cleaned_mags = delete_nans(jds, mags)

    mag_sorted_indices = sortperm(cleaned_mags)
    n_mags = length(cleaned_mags)
    n_10 = round(Int, n_mags / 10)

    mean_10 = mean([cleaned_mags[mag_sorted_indices[end-n_10+1:end]]; cleaned_mags[mag_sorted_indices[1:n_10]]])

    median_mag = median(cleaned_mags)
    σ = √(varm(cleaned_mags, median_mag))

    return (mean_10 - median_mag)/σ
end