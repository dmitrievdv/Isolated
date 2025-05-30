# using FITSIO
using Plots
using Plots.PlotMeasures
using Printf
using DelimitedFiles


function findbkg(cut, aperture_size)
    cut_width = size(cut)[1]
    cut_height = size(cut)[2]
    aperture_flux = zeros(cut_width-2aperture_size, cut_height-2aperture_size)
    for i=1+aperture_size:cut_width-aperture_size
        width_aperture = i-aperture_size:i+aperture_size
        for j=1+aperture_size:cut_height-aperture_size
            height_aperture = j-aperture_size:j+aperture_size
            aperture_flux[i-aperture_size,j-aperture_size] = sum(nantobig.(cut[width_aperture, height_aperture]))
        end
    end
    min_index = findmin(aperture_flux)[2] + CartesianIndex(aperture_size,aperture_size)
    return min_index[1]-aperture_size:min_index[1]+aperture_size, min_index[2]-aperture_size:min_index[2]+aperture_size
end

negativetonan(x) = x ≤ 0 ? NaN : x

nantobig(x) = isnan(x) ? 1e100 : x

function findstar(cut, aperture_size, search_box)
    cut_width = size(cut)[1]
    cut_height = size(cut)[2]
    center_px_x = (cut_width+1) ÷ 2
    center_px_y = (cut_height+1) ÷ 2
    aperture_flux = zeros(2search_box + 1, 2search_box + 1)
    no_search_width = center_px_x - 1 - search_box
    if no_search_width < aperture_size
        no_search_width = aperture_size
    end
    for i=1:2search_box + 1
        width_aperture = i+no_search_width-aperture_size:i+no_search_width+aperture_size
        for j=1:2search_box + 1
            height_aperture = j+no_search_width-aperture_size:j+no_search_width+aperture_size
            aperture_flux[i,j] = sum((cut[width_aperture, height_aperture]))
        end
    end
    max_index = findmax(aperture_flux)[2] + CartesianIndex(no_search_width,no_search_width)
    return max_index[1]-aperture_size:max_index[1]+aperture_size, max_index[2]-aperture_size:max_index[2]+aperture_size
end

function getlc(star, fits; sector = 0)
    cuts = negativetonan.(read(fits[2], "FLUX"))
    cuts_shape = size(cuts)
    cut_width = cuts_shape[1]
    cut_height = cuts_shape[2]
    N_cuts = cuts_shape[3]
    N_bkg = N_cuts ÷ 4

    jd = read(fits[2], "TIME") .+ 2457000

    plt = heatmap(log10.(cuts[:,:,N_bkg]'), aspect_ratio = :equal)
    
    

    aperture_size = 2
    center_px_x = (cut_width+1) ÷ 2
    center_px_y = (cut_height+1) ÷ 2

    # star_search_aperture = (center_px_x-2aperture_size:center_px_x+2aperture_size, center_px_y-2aperture_size:center_px_y+2aperture_size)

    raw_flux = zeros(N_cuts)
    bkg_flux = zeros(N_cuts)
    flux = zeros(N_cuts)

    for i_cut = 1:N_cuts
        aperture = findstar(cuts[:,:, i_cut], aperture_size, 2)
        # aperture = aperture[1] .+ (center_px_x - 2aperture_size - 1), aperture[2] .+ (center_px_y - 2aperture_size - 1)
        bkg_aperture = findbkg(cuts[:,:,i_cut], aperture_size)

        raw_flux[i_cut] = sum(cuts[aperture..., i_cut], dims = (1,2))[1]
        bkg_flux[i_cut] = sum(cuts[bkg_aperture..., i_cut], dims = (1,2))[1]

        flux[i_cut] = negativetonan(raw_flux[i_cut] - bkg_flux[i_cut])
    end

    aperture = findstar(cuts[:,:, N_bkg], aperture_size, 2)
    # aperture = aperture[1] .+ (center_px_x - 2aperture_size - 1), aperture[2] .+ (center_px_y - 2aperture_size - 1)
    bkg_aperture = findbkg(cuts[:,:,N_bkg], aperture_size)
    
    # aperture = (center_px_x-2aperture_size:center_px_x+2aperture_size, center_px_y-2aperture_size:center_px_y+2aperture_size)
    # aperture = findstar(cuts[aperture..., N_bkg], aperture_size)
    # aperture = aperture[1] .+ (center_px_x - 2aperture_size - 1), aperture[2] .+ (center_px_y - 2aperture_size - 1)
    # # aperture = (center_px_x-aperture_size:center_px_x+aperture_size, center_px_y-aperture_size:center_px_y+aperture_size)
    # bkg_aperture = findbkg(cuts[:,:,N_bkg], aperture_size)
    # raw_flux = reshape(sum(cuts[aperture..., :], dims = (1,2)), N_cuts)
    # bkg_flux = reshape(sum(cuts[bkg_aperture..., :], dims = (1,2)), N_cuts)
    # flux = negativetonan.(raw_flux .- bkg_flux)

    bkg_aperture_rect_x = [bkg_aperture[1][1] - 0.5, bkg_aperture[1][1] - 0.5, 
                           bkg_aperture[1][end] + 0.5, bkg_aperture[1][end] + 0.5, 
                           bkg_aperture[1][1] - 0.5]

    bkg_aperture_rect_y = [bkg_aperture[2][1] - 0.5, bkg_aperture[2][end] + 0.5, 
                           bkg_aperture[2][end] + 0.5, bkg_aperture[2][1] - 0.5, 
                           bkg_aperture[2][1] - 0.5]

    plot!(plt, bkg_aperture_rect_x, bkg_aperture_rect_y, lw = 2, label = "bkg", rightmargin = 20px)


    aperture_rect_x = [aperture[1][1] - 0.5, aperture[1][1] - 0.5, 
                           aperture[1][end] + 0.5, aperture[1][end] + 0.5, 
                           aperture[1][1] - 0.5]

    aperture_rect_y = [aperture[2][1] - 0.5, aperture[2][end] + 0.5, 
                           aperture[2][end] + 0.5, aperture[2][1] - 0.5, 
                           aperture[2][1] - 0.5]

    plot!(plt, aperture_rect_x, aperture_rect_y, lw = 2, label = "star", rightmargin = 20px)

    savefig(plt, "$star_directory/$star/$sector-map.pdf")
    # savefig(plt, "maps/$star-$sector-map.pdf")

    open("$star_directory/$star/$sector-rawlc.dat", "w") do io
        println(io, "#jd flux raw_flux bkg_flux")
        for n = 1:N_cuts
            @printf(io, "%15.6f %12.3e %12.3e %12.3e\n", jd[n], flux[n], raw_flux[n], bkg_flux[n])
        end
    end

    return Float64.(jd), Float64.(flux)
    # uniform_jd_step = round(minimum(jd[2:end] - jd[1:end-1])*24*60)/24/60
    # jd_start = jd[1]
    # jd_end = jd[end]
end

function cleanlc(jd, flux)
    non_nan_index = 0
    nan_index = 0
    clean_jd = Float64[]
    clean_flux = Float64[]
    N = length(jd)
    for i = 1:N
        # println("$nan_index $non_nan_index")
        if !isnan(flux[i]) & (nan_index == non_nan_index)
            non_nan_index = i
            # if nan_index < non_nan_index
            #     nan_index = i
            # end
        elseif isnan(flux[i])
            nan_index = i
            if (non_nan_index != 0) & ((nan_index - non_nan_index) > 10)
                # println(flux[non_nan_index:nan_index])
                push!(clean_jd, jd[non_nan_index:nan_index]...)
                push!(clean_flux, flux[non_nan_index:nan_index]...)
            end
            non_nan_index = i
        end
    end
    if (non_nan_index != 0) & ((N - non_nan_index) > 10)
        # println(flux[non_nan_index:N])
        push!(clean_jd, jd[non_nan_index:N]...)
        push!(clean_flux, flux[non_nan_index:N]...)
    end
    # println(clean_jd)
    # println(clean_flux)
    if (length(clean_jd) > 0)
        if isnan(clean_jd[end])
            return clean_jd[1:end-1], clean_flux[1:end-1]
        end
    end
        
    return clean_jd, clean_flux
end

function uniformlc(jd, flux, ε_jd = 1.5, ε_flux = 0.1)
    uniform_jd = Float64[]
    uniform_flux = Float64[]
    N_samples = length(jd)
    jd_step = jd[2] - jd[1]
    flux_step = flux[2] - flux[1]
    for i = 2:N_samples
        new_jd_step = jd[i] - jd[i-1]
        new_flux_step = flux[i] - flux[i-1]
        if (new_jd_step/jd_step > ε_jd) | (jd_step/new_jd_step > ε_jd)
            push!(uniform_jd, NaN)
            push!(uniform_flux, NaN)
        elseif (abs(new_flux_step)/flux[i-1] > ε_flux)
            if i < N_samples
                if abs(flux[i+1] - flux[i-1])/flux[i-1] < ε_flux
                    push!(uniform_jd, jd[i])
                    push!(uniform_flux, (flux[i+1] + flux[i-1])/2)
                    continue
                else abs(flux[i+1] - flux[i])/flux[i] < ε_flux
                    push!(uniform_jd, jd[i])
                    push!(uniform_flux, flux[i])
                    continue
                end
            end
            push!(uniform_jd, NaN)
            push!(uniform_flux, NaN)
            # if abs(flux_step + new_flux_step)/flux[i-1] < ε_flux
            #     push!(uniform_jd, jd[i])
            #     push!(uniform_flux, flux[i])
            # end
        else
            push!(uniform_jd, jd[i])
            push!(uniform_flux, flux[i])
            if isnan(uniform_flux[end])
                uniform_jd[end] = NaN
            end
        end
        jd_step = new_jd_step
        flux_step = new_flux_step
    end
    
    return uniform_jd, uniform_flux
end

function tessmag(flux)
    return -2.5*log10(flux) + 20.44
end

function mergelc(jd1, flux1, jd2, flux2)
    N_1 = length(jd1)
    N_2 = length(jd2)
    if N_1 == 0
        return jd2, flux2
    end
    if N_2 == 0
        return jd1, flux1
    end
    N = N_1 + N_2 + 1
    jd = zeros(N)
    flux = zeros(N)
    if jd1[end] < jd2[1]
        jd[1:N_1] = jd1[:]
        flux[1:N_1] = flux1[:]
        jd[N_1+1] = jd2[1]
        flux[N_1+1] = NaN
        jd[N_1+2:N] = jd2[:]
        flux[N_1+2:N] = flux2[:]
    else 
        jd[1:N_2] = jd2[:]
        flux[1:N_2] = flux2[:]
        jd[N_2+1] = jd1[1]
        flux[N_2+1] = NaN
        jd[N_2+2:N] = jd1[:]
        flux[N_2+2:N] = flux1[:]
    end
    return jd, flux
end

function get_lcs(star_name)
    nospace_star_name = get_nospace_star_name(star_name)
    archive = ZipReader(mmap(open("$star_directory/$nospace_star_name/$nospace_star_name.zip")))
    entry_names = zip_names(archive)

    jd_all = Float64[]
    flux_all = Float64[]

    for entry_name in entry_names
        data = zip_readentry(archive, entry_name)
        fits = FITS(data)
        sector = parse(Int, entry_name[7:10])

        jd, flux = getlc(nospace_star_name, fits, sector = sector)
        jd, flux = uniformlc(jd, flux)

        jd, flux = cleanlc(jd, flux)

        mag = -2.5*log10.(flux) .+ 20.44

        open("$star_directory/$nospace_star_name/$sector-lc.dat", "w") do io
            println(io, "#jd flux")
            for n = 1:length(jd)
                @printf(io, "%15.6f %15.6e %15.6f\n", jd[n], flux[n], mag[n])
            end
        end
        if length(jd) > 0
            int_dates = collect(ceil(Int, jd[1]/5)*5:5:floor(Int, jd[end]/5)*5)
            string_dates = Dates.format.(julian2datetime.(int_dates), "d u Y")
            plt = plot(jd, tessmag.(flux), xticks = (int_dates, string_dates), label = false, rightmargin = 15px, yflip = true, ylabel = "TESS magnitude")
            # savefig(plt, "plots/$star-$sector.pdf")
            # savefig(plt, "plots/$star-$sector.png")
            savefig(plt, "$star_directory/$nospace_star_name/$sector-lc.pdf")
            savefig(plt, "$star_directory/$nospace_star_name/$sector-lc.png")
            jd_all, flux_all = mergelc(jd_all, flux_all, jd, flux)
        end
    end
    jd_all, flux_all
end