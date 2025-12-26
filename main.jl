using VirtualObservatory
using DataFrames
using HTTP
using ZipArchives
using Mmap
# using BufferedStreams
# using JSON3
using CSV
using Photometry
using Makie
using GLMakie
using Printf
using FITSIO
using Dates
using LinearAlgebra
using Statistics
using LeastSquaresOptim
using ProgressBars
using Dierckx
using LombScargle
# using Logging

star_directory = "stars_julia"

tess_max_sectors = 96

# include("lightcurves.jl")
include("databases.jl")
include("tess-queries.jl")

get_nospace_star_name(star_name) = replace(star_name, " " => "_")

struct ErrorNoTESScut <: Exception 
    msg :: String
end

function Base.showerror(io :: IO, e :: ErrorNoTESScut)
    # println(io, e)
    print(io, typeof(e), "\n\n")
    println(io, e.msg)
end

struct ErrorTESSWrongSector <: Exception 
    sector_int :: Int
end

function Base.showerror(io :: IO, e :: ErrorTESSWrongSector)
    print(io, typeof(e), "\n\n")
    println(io, "Wrong sector: ", e.sector_int)
end

function extract_tess_cutouts(star_name = "star")
    nospace_star_name = get_nospace_star_name(star_name) 
    if !isfile("$star_directory/$nospace_star_name/$nospace_star_name.zip")
        throw(ErrorNoTESScut("No tesscut file for $star_name is found. Download it using get_tess_cutouts(...)."))    
    else
        archive = ZipReader(mmap(open("$star_directory/$nospace_star_name/$nospace_star_name.zip")))
        entry_names = zip_names(archive)
        for entry_name in entry_names
            data = zip_readentry(archive, entry_name)
            open("$star_directory/$nospace_star_name/$entry_name", "w") do io 
                write(io, data)
            end
        end
    end
end

function get_tess_sectors(star_name; query = false)
    nospace_star_name = get_nospace_star_name(star_name)
    if !isfile("$star_directory/$nospace_star_name/$nospace_star_name.zip") & !query
        throw(ErrorNoTESScut("No tesscut file for $star_name is found. Download it using get_tess_cutouts(...)."))
    elseif query
        star_ra, star_dec = get_star_coords(star_name)
        sector_jsons = get_tess_sectors(star_ra, star_dec)
        n_sectors = length(sector_jsons)
        sectors = zeros(Int, n_sectors)
        for i_sector = 1:n_sectors
            sectors[i_sector] = parse(Int, sector_jsons[i_sector].sector)
        end
        return sectors
    else
        archive = ZipReader(mmap(open("$star_directory/$nospace_star_name/$nospace_star_name.zip")))
        entry_names = zip_names(archive)
        n_entries = zip_nentries(archive)
        sectors = zeros(Int, n_entries)
        for i_entry = 1:n_entries
            data = zip_readentry(archive, entry_names[i_entry])
            fits = FITS(data)
            sectors[i_entry] = read_key(fits[1], "SECTOR")[1]
        end
        return sectors
    end
end

get_px_to_radec_matrix(fits) = [read_key(fits[2], "11PC4")[1] read_key(fits[2], "12PC4")[1]
                                read_key(fits[2], "21PC4")[1] read_key(fits[2], "22PC4")[1]]

get_reference_radec(fits) = [read_key(fits[2], "1CRVL4")[1], read_key(fits[2], "2CRVL4")[1]]

get_reference_px(fits) = [read_key(fits[2], "1CRPX4")[1], read_key(fits[2], "2CRPX4")[1]]

function get_tesscut_corners(fits)
    reference_px = get_reference_px(fits)
    reference_radec = get_reference_radec(fits)
    conversion_matrix_px_to_radec = get_px_to_radec_matrix(fits)

    width, height = size(read(fits[2], "FLUX"))[1:2]
    # width = 1.5*width; height = 1.5*height

    bottom_left = conversion_matrix_px_to_radec * ([0.5,0.5] - reference_px)
    bottom_left = get_true_radec(reference_radec..., bottom_left...)
    bottom_right = conversion_matrix_px_to_radec * ([width+0.5,0.5] - reference_px)
    bottom_right = get_true_radec(reference_radec..., bottom_right...)
    top_left = conversion_matrix_px_to_radec * ([0.5,height+0.5] - reference_px)
    top_left = get_true_radec(reference_radec..., top_left...)
    top_right = conversion_matrix_px_to_radec * ([width+0.5,height+0.5] - reference_px)
    top_right = get_true_radec(reference_radec..., top_right...)

    return [bottom_left, top_left, top_right, bottom_right]
end

function get_gaia_stars_in_poly(corners, mag_threshold; gaia = "dr3")
    fold_corner = cs -> foldl((c1, c2) -> c1*", "*c2, cs)
    corners_string = fold_corner([fold_corner(string.(corner)) for corner in corners])
    df_gaia = DataFrame(execute(TAPService(:gaia), "select * from gaia$gaia.gaia_source where CONTAINS(POINT('ICRS', ra, dec),"*
                                              " POLYGON('ICRS', $corners_string)) = 1 and phot_rp_mean_mag < $mag_threshold"))
end

function get_star_tesscut_fits(star_name, sector)
    nospace_star_name = get_nospace_star_name(star_name)

    sectors = get_tess_sectors(star_name; query = false)

    if !isfile("$star_directory/$nospace_star_name/$nospace_star_name.zip")
        throw(ErrorNoTESScut("No tesscut file for $star_name is found. Download it using get_tess_cutouts(...)."))
    end

    archive = ZipReader(mmap(open("$star_directory/$nospace_star_name/$nospace_star_name.zip")))
    entry_names = zip_names(archive)

    i_sector = findfirst(s -> s == sector, sectors)

    if isa(i_sector, Nothing)
        throw(ErrorTESSWrongSector(sector))
    end

    return FITS(zip_readentry(archive, i_sector))
end

function get_n_min_mean_background(cut_flux, n_min)
    return sum(sort(vec(cut_flux))[1:n_min])/n_min
end

function get_n_min_median_background(cut_flux, n_min)
    return median(sort(vec(cut_flux))[1:n_min])
end

function get_tesscut_lc(fits)
    return fits
end

function calc_aperture_photometry(cut, star_px_x, star_px_y, aperture_radius)
    px_width, px_height = size(cut)
    n_px = px_width*px_height
    bkg = get_n_min_median_background(cut, n_px ÷ 4)
    ap = CircularAperture(star_px_x, star_px_y, aperture_radius)
    return photometry(ap, cut .- bkg)[3]
end

function calc_aperture_photometry_bkg_pixels(cut, bkg_positions, star_px_x, star_px_y, aperture_radius)
    # px_width, px_height = size(cut)
    # n_px = px_width*px_height
    if !isnothing(findfirst(x -> abs(x) < 1e-8, cut))
        return NaN
    end
    bkg_cut = fit_flat_background(cut, bkg_positions)
    ap = CircularAperture(star_px_x, star_px_y, aperture_radius)
    return photometry(ap, cut .- bkg_cut)[3]
end

function get_true_radec(α_0, δ_0, Δα, Δδ)
    α_0_rad = α_0/180*π
    δ_0_rad = δ_0/180*π
    Δα_rad = Δα/180*π
    Δδ_rad = Δδ/180*π

    α̂ = [-sin(α_0_rad), cos(α_0_rad), 0]
    δ̂ = [-cos(α_0_rad)*sin(δ_0_rad), -sin(α_0_rad)*sin(δ_0_rad), cos(δ_0_rad)]

    xyz_0 = [cos(α_0_rad)*cos(δ_0_rad), sin(α_0_rad)*cos(δ_0_rad), sin(δ_0_rad)]
    xyz = xyz_0 + Δα_rad * α̂ + Δδ_rad * δ̂
    xyz = xyz / √(sum(xyz .^ 2))

    δ = asin(xyz[3])/π*180
    α = atan(xyz[2], xyz[1])/π*180

    if α < 0
        return [α + 360, δ]
    else
        return [α, δ]
    end
end

function get_rel_radec(α_0, δ_0, α, δ)
    α_0_rad = α_0/180*π
    δ_0_rad = δ_0/180*π
    α_rad = α/180*π
    δ_rad = δ/180*π

    xyz_0 = [cos(α_0_rad)*cos(δ_0_rad), sin(α_0_rad)*cos(δ_0_rad), sin(δ_0_rad)]
    xyz = [cos(α_rad)*cos(δ_rad), sin(α_rad)*cos(δ_rad), sin(δ_rad)]
    xyz = xyz/(xyz ⋅ xyz_0)
    
    Δxyz = xyz - xyz_0

    α̂ = [-sin(α_0_rad), cos(α_0_rad), 0]
    δ̂ = [-cos(α_0_rad)*sin(δ_0_rad), -sin(α_0_rad)*sin(δ_0_rad), cos(δ_0_rad)]

    Δα = Δxyz ⋅ α̂
    Δδ = Δxyz ⋅ δ̂
    return [Δα, Δδ]/π*180  
end

function get_distance(α_1, δ_1, α_2, δ_2)
    α_1_rad = α_1/180*π
    δ_1_rad = δ_1/180*π
    α_2_rad = α_2/180*π
    δ_2_rad = δ_2/180*π
    xyz_1 = [cos(α_1_rad)*cos(δ_1_rad), sin(α_1_rad)*cos(δ_1_rad), sin(δ_1_rad)]
    xyz_2 = [cos(α_2_rad)*cos(δ_2_rad), sin(α_2_rad)*cos(δ_2_rad), sin(δ_2_rad)]

    return acos(xyz_1 ⋅ xyz_2)/180*π
end

function calc_tess_magnitude(flux)
    return -2.5*log10(abs(flux)) + 20.44
end

function calc_tess_flux_from_mag(mag)
    return 10^(0.4*(20.44 - mag))
end

function elliptic_powexp_PSF(width, height, center_x, center_y, x_σ, y_σ, power, max_int; super_samp = 10)
    model = zeros(width, height)
    d_ss = 1/super_samp
    for x_px = 1:width, y_px = 1:height
        corner_x = x_px - 0.5
        corner_y = y_px - 0.5
        for i_ss_x = 1:super_samp, i_ss_y = 1:super_samp
            x_ss = corner_x + i_ss_x*d_ss - d_ss/2
            y_ss = corner_y + i_ss_y*d_ss - d_ss/2
            model[x_px, y_px] += max_int*exp(-(abs((center_x - x_ss)/x_σ)^(power) + abs((center_y - y_ss)/y_σ)^power))*d_ss^2
        end
    end
    return model
end

function elliptic_gaussian_PSF(width, height, center_x, center_y, x_σ, y_σ, max_int; super_samp = 10)
    model = zeros(width, height)
    d_ss = 1/super_samp
    for x_px = 1:width, y_px = 1:height
        corner_x = x_px - 0.5
        corner_y = y_px - 0.5
        for i_ss_x = 1:super_samp, i_ss_y = 1:super_samp
            x_ss = corner_x + i_ss_x*d_ss - d_ss/2
            y_ss = corner_y + i_ss_y*d_ss - d_ss/2
            model[x_px, y_px] += max_int*exp(-(((center_x - x_ss)/x_σ)^2 + ((center_y - y_ss)/y_σ)^2))*d_ss^2
        end
    end
    return model
end

function gaussian_PSF(width, height, center_x, center_y, σ, max_int; super_samp = 10)
    model = zeros(width, height)
    d_ss = 1/super_samp
    for x_px = 1:width, y_px = 1:height
        corner_x = x_px - 0.5
        corner_y = y_px - 0.5
        for i_ss_x = 1:super_samp, i_ss_y = 1:super_samp
            x_ss = corner_x + i_ss_x*d_ss - d_ss/2
            y_ss = corner_y + i_ss_y*d_ss - d_ss/2
            model[x_px, y_px] += max_int*exp(-((center_x - x_ss)^2 + (center_y - y_ss)^2)/σ^2)*d_ss^2
        end
    end
    return model
end

function fit_stars_gauss(cut_no_bkg, stars_px_x, stars_px_y; super_samp = 10)
    width, height = size(cut_no_bkg)
    n_stars = length(stars_px_x)
    function to_optimize(pars)
        max_ints = pars[1:n_stars]
        y_σ = pars[end]
        x_σ = pars[end-1]
        power = pars[end-2]
        model = deepcopy(cut_no_bkg)
        for i_star in 1:n_stars
            model = model .- elliptic_powexp_PSF(width, height, stars_px_x[i_star], stars_px_y[i_star], x_σ, y_σ, power, abs(max_ints[i_star]); super_samp = super_samp)
        end
        return vec(model)
    end

    start_pars = fill(100.0, n_stars+3)
    start_pars[end-1:end] .= 1
    start_pars[end-2] = 2

    optimize(to_optimize, start_pars, LevenbergMarquardt())
end

function fit_stars_prf(supersampled_prf, cut_no_bkg, stars_px_x, stars_px_y, start_pars)
    width, height = size(cut_no_bkg)
    n_stars = length(stars_px_x)
    prf_cuts = vec.(get_prf_cut.(Ref(supersampled_prf), width, height, stars_px_x, stars_px_y))
    model_start = vec(deepcopy(cut_no_bkg))
    
    function to_optimize!(model, pars)
        model .= model_start
        for i_star = 1:n_stars
            model .= model - prf_cuts[i_star] * abs(pars[i_star])
        end
    end

    optimize!(LeastSquaresProblem(x = start_pars, f! = to_optimize!, output_length = width*height), LevenbergMarquardt())
end

function fit_stars_prf_bkg(supersampled_prf, cut, stars_px_x, stars_px_y)
    width, height = size(cut)
    n_stars = length(stars_px_x)
    prf_cuts = get_prf_cut.(Ref(supersampled_prf), width, height, stars_px_x, stars_px_y)
    function to_optimize(pars)
        model = deepcopy(cut)
        for i_star = 1:n_stars
            model = model .- prf_cuts[i_star] * abs(pars[i_star])
        end

        return vec(model .- pars[end])
    end

    start_pars = fill(100.0, n_stars + 1)

    optimize(to_optimize, start_pars, LevenbergMarquardt())
end

function fit_stars_prf_flat_bkg(supersampled_prf, cut, stars_px_x, stars_px_y, start_pars)
    width, height = size(cut)
    n_stars = length(stars_px_x)
    prf_cuts = vec.(get_prf_cut.(Ref(supersampled_prf), width, height, stars_px_x, stars_px_y))
    model_start = vec(deepcopy(cut))

    function jacobian!(jacob, pars)
        jacob .= 0.0

        for i_star = 1:n_stars
            jacob[:,i_star] -= prf_cuts[i_star]*sign(pars[i_star])
        end

        bkg_plane_normal_length = √(pars[end-3]^2 + pars[end-2]^2 + pars[end-1]^2)

        for x=1:width, y=1:height
            jacob[(x-1)*width + y, end] = -bkg_plane_normal_length/pars[end-1]
            jacob[(x-1)*width + y, end-3] = (x - pars[end-3]/bkg_plane_normal_length*pars[end])/pars[end-1]
            jacob[(x-1)*width + y, end-2] = (y - pars[end-2]/bkg_plane_normal_length*pars[end])/pars[end-1]
            jacob[(x-1)*width + y, end-1] = (pars[end]*bkg_plane_normal_length - pars[end-3]*x - pars[end-2]*y)/pars[end-1]^2 - 1.0/bkg_plane_normal_length*pars[end] 
        end
    end
    
    function to_optimize!(model, pars)
        model .= model_start
        for i_star = 1:n_stars
            model .= model - prf_cuts[i_star] * abs(pars[i_star])
        end

        bkg_plane_normal_length = √(pars[end-3]^2 + pars[end-2]^2 + pars[end-1]^2)

        for x = 1:width, y = 1:height
            model[(x-1)*width + y] = model[(x-1)*width + y] - (pars[end]*bkg_plane_normal_length - pars[end-3]*x - pars[end-2]*y)/pars[end-1]
        end
    end

    optimize!(LeastSquaresProblem(x = start_pars, f! = to_optimize!, g! = jacobian!, output_length = width*height), LevenbergMarquardt(), x_tol = 1e-10, f_tol = 1e-10, g_tol = 1e-10)
end

function get_tess_ffi_coordinates(α, δ, sector)
    ptr = @ccall "./tess_point_c/libtess_stars2px.so.1.0.1".tess_stars2px_sector(α::Float64, δ::Float64, sector :: Int)::Ptr{Cdouble}
    coords = unsafe_wrap(Vector{Float64}, ptr, 2)
    return coords
end

function get_tesscut_ffi_coordinates(cut_fits)
    α_cut, δ_cut = get_reference_radec(cut_fits)
    sector = read_key(cut_fits[1], "SECTOR")[1]

    return get_tess_ffi_coordinates(α_cut, δ_cut, sector)
end

function get_tesscut_prf_supersampled(cut_fits)
    α_cut, δ_cut = get_reference_radec(cut_fits)
    sector = read_key(cut_fits[1], "SECTOR")[1] 
    ccd_x, ccd_y = round.(Int, get_tess_ffi_coordinates(α_cut, δ_cut, sector))

    cam = read_key(cut_fits[1], "CAMERA")[1]
    ccd = read_key(cut_fits[1], "CCD")[1]

    prf_url = "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles"
    prf_dir = "start_s000$(sector > 3 ? 4 : 1)/cam$(cam)_ccd$(ccd)"

    prf_prefix = if sector < 4
        if (cam > 2) | ((cam == 2) & (ccd == 4))
            "tess2018243163601"
        else
            "tess2018243163600"
        end
    else
        if (cam > 3) | ((cam == 3) & (ccd ≥ 2))
            "tess2019107181902"
        elseif (cam == 1) & (ccd == 1)
            "tess2019107181900"
        else
            "tess2019107181901"
        end
    end

    B_prf_row = ((ccd_y - 1) ÷ 512)*512 + 1
    L_prf_col = ((ccd_x - 45) ÷ 512)*512 + 45
    T_prf_row = ((ccd_y - 1) ÷ 512 + 1)*512 + 1
    R_prf_col = ((ccd_x - 45) ÷ 512 + 1)*512 + 45

    B_prf_row -= (B_prf_row > 1500)
    L_prf_col -= (L_prf_col > 1500)
    T_prf_row -= (T_prf_row > 1500)
    R_prf_col -= (R_prf_col > 1500)

    B_prf_row_str = @sprintf "%04i"  B_prf_row
    T_prf_row_str = @sprintf "%04i"  T_prf_row
    L_prf_col_str = @sprintf "%04i"  L_prf_col
    R_prf_col_str = @sprintf "%04i"  R_prf_col

    BL_prf_file_name = "$prf_prefix-prf-$cam-$ccd-row$B_prf_row_str-col$L_prf_col_str.fits"
    BR_prf_file_name = "$prf_prefix-prf-$cam-$ccd-row$B_prf_row_str-col$R_prf_col_str.fits"
    TL_prf_file_name = "$prf_prefix-prf-$cam-$ccd-row$T_prf_row_str-col$L_prf_col_str.fits"
    TR_prf_file_name = "$prf_prefix-prf-$cam-$ccd-row$T_prf_row_str-col$R_prf_col_str.fits"

    mkpath("prf/$prf_dir")

    # println("$prf_url/$prf_dir/$BL_prf_file_name")

    if !isfile("prf/$prf_dir/$BL_prf_file_name")
        HTTP.download("$prf_url/$prf_dir/$BL_prf_file_name", "prf/$prf_dir/$BL_prf_file_name")
    end
    if !isfile("prf/$prf_dir/$BR_prf_file_name")
        HTTP.download("$prf_url/$prf_dir/$BR_prf_file_name", "prf/$prf_dir/$BR_prf_file_name")
    end
    if !isfile("prf/$prf_dir/$TL_prf_file_name")
        HTTP.download("$prf_url/$prf_dir/$TL_prf_file_name", "prf/$prf_dir/$TL_prf_file_name")
    end
    if !isfile("prf/$prf_dir/$TR_prf_file_name")
        HTTP.download("$prf_url/$prf_dir/$TR_prf_file_name", "prf/$prf_dir/$TR_prf_file_name")
    end
    
    BR_prf = read(FITS("prf/$prf_dir/$BR_prf_file_name")[1])
    BL_prf = read(FITS("prf/$prf_dir/$BL_prf_file_name")[1])
    TR_prf = read(FITS("prf/$prf_dir/$TR_prf_file_name")[1])
    TL_prf = read(FITS("prf/$prf_dir/$TL_prf_file_name")[1])

    # Linear Interpolation in a square
    interpolated_prf = ((ccd_y - B_prf_row)*(ccd_x - L_prf_col)*TR_prf 
           -(ccd_y - B_prf_row)*(ccd_x - R_prf_col)*TL_prf
           -(ccd_y - T_prf_row)*(ccd_x - L_prf_col)*BR_prf 
           +(ccd_y - T_prf_row)*(ccd_x - R_prf_col)*BL_prf)/(T_prf_row - B_prf_row)/(R_prf_col - L_prf_col)

    return interpolated_prf
end

function get_prf_cut(prf_supersampled, cut_width, cut_height, x_source, y_source)
    x_source_px = round(Int, x_source)
    y_source_px = round(Int, y_source)

    δx = round(Int, (x_source - x_source_px)*9)
    δy = round(Int, (y_source - y_source_px)*9)

    prf = zeros(cut_width, cut_height)
    supersampled_width, supersampled_height = size(prf_supersampled)

    centered_supersampled_sum = [55,63]

    correct_bound(i) = i ≥ 1 ? (i ≤ 117 ? i : 117) : 1

    supersampled_mask = zeros(supersampled_width, supersampled_height)

    for x_px = 1:cut_width, y_px = 1:cut_height
        x_prf_source = (x_px - x_source)*9 + 59
        y_prf_source = (y_px - y_source)*9 + 59
        # Δx_int = round(Int, Δx)
        # Δy_int = round(Int, Δy)

        # Δx_frac = Δx - Δx_int
        # Δy_frac = Δy - Δy_int

        supersampled_mask .= 1.0

        for x_prf = 1:supersampled_width
            for y_prf = 1:supersampled_height
                if (abs(x_prf - x_prf_source) ≥ 5) | (abs(y_prf - y_prf_source) ≥ 5)
                    supersampled_mask[x_prf, y_prf] = 0.0
                    continue
                end
                if (abs(x_prf - x_prf_source) ≤ 4) & (abs(y_prf - y_prf_source) ≤ 4)
                    supersampled_mask[x_prf, y_prf] = 1.0
                    # println("1.0")
                    continue
                end
                if abs(x_prf - x_prf_source) > 4
                    supersampled_mask[x_prf, y_prf] *= 5 - abs(x_prf - x_prf_source)
                end
                if abs(y_prf - y_prf_source) > 4
                    supersampled_mask[x_prf, y_prf] *= 5 - abs(y_prf - y_prf_source)
                end
            end
        end
        # println(sum(prf_supersampled .* supersampled_mask))
        prf[x_px, y_px] = sum(prf_supersampled .* supersampled_mask)/81
        if prf[x_px, y_px] < 2e-4
            prf[x_px, y_px] = 0.0
        end
    end
    return prf
end

function estimate_tess_ffi_coordinates(cut_fits)
    α_cut, δ_cut = get_reference_radec(cut_fits)
    conversion_matrix_px_to_radec = get_px_to_radec_matrix(cut_fits)
    conversion_matrix_radec_to_px = inv(conversion_matrix_px_to_radec)
    sectors_df = CSV.read("sectors.csv", DataFrame)

    cam = read_key(fits[1], "CAMERA")[1]
    cam_ra_name = "cam$cam"*"_ra"
    cam_dec_name = "cam$cam"*"_dec"

    sector = read_key(fits[1], "SECTOR")[1] 

    prf_url_dir = if sector < 4
        "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0001/"
    else
        "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0004/"
    end

    sector_df_i = findfirst(s -> s == sector, sectors_df.sector[:])
    α_cam = sectors_df[sector_df_i, cam_ra_name]
    δ_cam = sectors_df[sector_df_i, cam_dec_name]

    println("$α_cut $δ_cut, $α_cam $δ_cam")
    println(get_rel_radec(α_cut, δ_cut, α_cam, δ_cam))

    Δx, Δy = conversion_matrix_radec_to_px * get_rel_radec(α_cut, δ_cut, α_cam, δ_cam)
    x = if Δx < 0 
        1 - Δx
    else 
        2048 - Δx
    end

    y = if Δy < 0
        1 - Δy
    else
        2048 - Δy
    end

    prf_row = ((y - 1) ÷ 512)*512 + 1
    prf_col = ((x - 45) ÷ 512)*512 + 45

    prf_row -= (prf_row > 1500)
    prf_col -= (prf_col > 1500)

    prf_row_str = @sprintf "%04i"  prf_row
    prf_col_str = @sprintf "%04i"  prf_col

    ccd = read_key(fits[1], "CCD")[1]

    prf_file_name = "tess2019107181901-prf-$cam-$ccd-row$prf_row_str-col$prf_col_str.fits"

    return prf_file_name
end

function create_gaia_datafiles(star_name; rewrite = false)
    gaia_data_file = "$star_directory/$(get_nospace_star_name(star_name))/gaia_target.csv"
    gaia_data = if !isfile(gaia_data_file) | rewrite
        data = get_star_gaia_data(star_name)
        CSV.write(gaia_data_file, DataFrame(data))
        data
    else
        CSV.read(gaia_data_file, DataFrame)[1, :]
    end


    if isfile("$star_directory/$nospace_star_name/$nospace_star_name.zip")
        sectors = get_tess_sectors(star_name)
        fits = get_star_tesscut_fits(star_name, sectors[1])
        corners = get_tesscut_corners(fits)
        Δm_R = 5
        star_R = gaia_data.phot_rp_mean_mag
        gaia_stars_file = "$star_directory/$(get_nospace_star_name(star_name))/gaia_stars_in_view.csv"

        if !isfile(gaia_stars_file) | rewrite
            data = get_gaia_stars_in_poly(corners, star_R + Δm_R; gaia = "dr2")
            CSV.write(gaia_stars_file, data)
        end
    end
end

function find_background_prf(flux_cut, supersampled_prf, stars_x, stars_y)
    n_stars = length(stars_x)
    start_pars = fill(100.0, n_stars+4)
    cut_width, cut_height = size(flux_cut)
    res = fit_stars_prf_flat_bkg(supersampled_prf, flux_cut, stars_x, stars_y, start_pars)
    PRFs = [abs.(res.minimizer[i_star])*get_prf_cut(supersampled_prf, cut_width, cut_height, stars_x[i_star], stars_y[i_star])  for i_star = 1:n_stars]
    empty_cut = sum([abs.(res.minimizer[i_star])*fill(get_n_min_median_background(PRFs[i_star], cut_width*cut_height ÷ 4), (cut_width, cut_height)) for i_star = 1:n_stars])
    PRF_cut = sum(PRFs)
    median_prf = sort(vec(PRF_cut - empty_cut))[cut_width*cut_height ÷ 2]
    median_prf_px = findall(x -> (median_prf - x) > -1e-8, PRF_cut - empty_cut)

    median_cut = median(flux_cut[median_prf_px])
    final_indeces = findall(x -> (flux_cut[x] - median_cut) < -1e-8, median_prf_px)

    # println(PRF_cut)
    # format = Printf.Format("%8.2f "^15 * "\n")
    # for i = 1:cut_width
    #     Printf.format(stdout, format, PRF_cut[:, cut_height - i +1]...)
    # end
    # println(quartile)
    # for i = 1:cut_width
    #     Printf.format(stdout, format, empty_cut[:, cut_height - i +1]...)
    # end
    return median_prf_px[final_indeces]
end

function fit_flat_background(flux_cut, bkg_positions)
    cut_width, cut_height = size(flux_cut)
    bkg_xs = [index[1] for index in bkg_positions]
    bkg_ys = [index[2] for index in bkg_positions]
    bkg_fluxes = flux_cut[bkg_positions]

    bkg_cut = zeros(cut_width, cut_height)
    if isempty(bkg_positions)
        bkg_cut .= NaN
        return bkg_cut
    end

    function to_optimize(pars)
        normal = √(pars[1]^2 + pars[2]^2 + pars[3]^2)
        return bkg_fluxes - (pars[4]*normal .- pars[1]*bkg_xs - pars[2]*bkg_ys)/pars[3]
    end

    bkg_plane = optimize(to_optimize, [0.0,0.0,1.0,100.0], LevenbergMarquardt()).minimizer

    
    normal = √(bkg_plane[1]^2 + bkg_plane[2]^2 + bkg_plane[3]^2)
    for x = 1:cut_width, y = 1:cut_height
        bkg_cut[x,y] = (bkg_plane[4]*normal - bkg_plane[1]*x - bkg_plane[2]*y)/bkg_plane[3]
    end

    return bkg_cut
end

function load_star_gaia_data(star_name)
    mkpath("$star_directory/$(get_nospace_star_name(star_name))")
    gaia_data_file = "$star_directory/$(get_nospace_star_name(star_name))/gaia_target.csv"
    if !isfile(gaia_data_file) 
        data = get_star_gaia_data(star_name)
        CSV.write(gaia_data_file, DataFrame(data))
        data
    else
        CSV.read(gaia_data_file, DataFrame)[1, :]
    end
end

load_tess_cutouts(star_name, cut_size) = load_tess_cutouts(star_name, cut_size, cut_size)

function load_tess_cutouts(star_name, cut_width, cut_height)
    mkpath("$star_directory/$(get_nospace_star_name(star_name))")
    gaia_data = load_star_gaia_data(star_name)
    ra, dec = gaia_data[[:ra, :dec]]
    nospace_star_name = get_nospace_star_name(star_name)
    cutouts_file = "$star_directory/$nospace_star_name/$(cut_width)x$(cut_height)/$(nospace_star_name).zip"
    if !isfile(cutouts_file)
        mkpath("$star_directory/$nospace_star_name/$(cut_width)x$(cut_height)")
        get_tess_cutouts(ra, dec, cut_width, cut_height; star_name = star_name)
    end

    archive = ZipReader(mmap(open(cutouts_file)))
    entry_names = zip_names(archive)

    all_fits = FITS.(zip_readentry.(Ref(archive), entry_names))
    sectors = [read_key(fits[1], "SECTOR")[1] for fits in all_fits]
    return Dict([key => value for (key, value) in zip(sectors, all_fits)])
end

function load_gaia_stars_in_view_data(star_name, fits, Δm_R)
    gaia_data = load_star_gaia_data(star_name)
    star_R = gaia_data.phot_rp_mean_mag
    corners = get_tesscut_corners(fits)
    cut_width, cut_height = read_key(fits[3], "NAXIS1")[1], read_key(fits[3], "NAXIS2")[1] 
    sector = read_key(fits[1], "SECTOR")[1]
    # println(sector)
    gaia_stars_file = "$star_directory/$(get_nospace_star_name(star_name))/$(cut_width)x$(cut_height)/gaia_stars_in_view_sector_$sector.csv"
    if !isfile(gaia_stars_file) 
        mkpath("$star_directory/$(get_nospace_star_name(star_name))/$(cut_width)x$(cut_height)")
        data = get_gaia_stars_in_poly(corners, star_R + Δm_R)
        reference_px = [read_key(fits[2], "1CRPX4")[1], read_key(fits[2], "2CRPX4")[1]]
        reference_radec = [read_key(fits[2], "1CRVL4")[1], read_key(fits[2], "2CRVL4")[1]]
        conversion_matrix_px_to_radec = [read_key(fits[2], "11PC4")[1] read_key(fits[2], "12PC4")[1]
                                         read_key(fits[2], "21PC4")[1] read_key(fits[2], "22PC4")[1]]

        conversion_matrix_radec_to_px = inv(conversion_matrix_px_to_radec)
        n_stars = nrow(data) 

        stars_x = zeros(n_stars)
        stars_y = zeros(n_stars)
        stars_mag = zeros(n_stars)

        for i_star = 1:n_stars
            star_radec = [data.ra[i_star], data.dec[i_star]]
            delta_radec = get_rel_radec(reference_radec..., star_radec...)
            delta_radec[1] = delta_radec[1]
            stars_mag[i_star] = data.phot_rp_mean_mag[i_star]
            stars_x[i_star], stars_y[i_star] = reference_px + conversion_matrix_radec_to_px * delta_radec
        end

        data[!, :px_x] = stars_x
        data[!, :px_y] = stars_y

        CSV.write(gaia_stars_file, data)
        data
    else
        CSV.read(gaia_stars_file, DataFrame)
    end
end

function aperture_prf_correction(aperture, star_px_x, star_px_y, supersampled_prf, cut_size)
    prf_cut = get_prf_cut(supersampled_prf, cut_size, cut_size, star_px_x, star_px_y)
    prf_flux = sum(prf_cut)
    ap_flux = photometry(aperture, prf_cut)[3]
    println("$ap_flux $prf_flux $(prf_flux/ap_flux)")
    return prf_flux/ap_flux
end

function calc_aperture_prf_correction(aperture_radius :: Real, star_px_x, star_px_y, supersampled_prf, cut_size)
    aperture = CircularAperture(star_px_x, star_px_y, aperture_radius)
    aperture_prf_correction(aperture, star_px_x, star_px_y, supersampled_prf, cut_size)
end

load_light_curve(star_name, sector, cut_size; kwargs...) = load_light_curve(star_name, sector, cut_size, cut_size; kwargs...)

function load_light_curve(star_name, sector, cut_width, cut_height; Δm_R = 5, rewrite_file = false)
    fits = load_tess_cutouts(star_name, cut_width, cut_height)[sector]
    flux_cuts = read(fits[2], "FLUX")
    n_cuts= size(flux_cuts)[3]

    mjds = read(fits[2], "TIME")

    gaia_stars_data = load_gaia_stars_in_view_data(star_name, fits, Δm_R)
    gaia_data = load_star_gaia_data(star_name)

    star_index = findfirst(s -> s == gaia_data.source_id, gaia_stars_data.source_id)
    star_px = gaia_stars_data.px_x[star_index], gaia_stars_data.px_y[star_index]

    light_curve_file = "$star_directory/$(get_nospace_star_name(star_name))/$(cut_width)x$(cut_height)/light_curve_sector_$sector.csv"

    aperture_radius = 5
    
    if !isfile(light_curve_file) | rewrite_file
        mkpath("$star_directory/$(get_nospace_star_name(star_name))/$(cut_width)x$(cut_height)")
        prf = get_tesscut_prf_supersampled(fits)
        bkg_pixels = find_background_prf(flux_cuts[:,:,n_cuts÷4], prf, gaia_stars_data.px_x, gaia_stars_data.px_y)
        # println(gaia_stars_data.px_x)
        # println(gaia_stars_data.px_y)
        # println(bkg_pixels)
        cuts = [flux_cuts[:,:,i_cut] for i_cut = 1:n_cuts]
        aperture_prf_correction = calc_aperture_prf_correction(aperture_radius, star_px..., prf, cut_height)
        phot_flux = aperture_prf_correction * calc_aperture_photometry_bkg_pixels.(cuts, Ref(bkg_pixels), star_px..., aperture_radius)

        lc_df = DataFrame(:MJD => mjds, :FLUX => phot_flux, :MAG => calc_tess_magnitude.(abs.(phot_flux)))
        CSV.write("$star_directory/$(get_nospace_star_name(star_name))/$(cut_width)x$(cut_height)/light_curve_sector_$sector.csv", lc_df)
        lc_df
    else
        CSV.read(light_curve_file, DataFrame)
    end
end

function is_in_sector(α, δ, sector)
    ptr = @ccall "tess_point_c/libtess_stars2px.so.1.0.1".tess_stars2px_sector(α::Float64, δ::Float64, sector::Int)::Ptr{Cdouble}
    px_coord = unsafe_wrap(Vector{Float64}, ptr, 2)
    if (px_coord[1] < 0.0) | (px_coord[2] < 0.0)
        return false
    else
        return true
    end
end

function find_tess_sectors(α, δ, max_sector)
    sectors = Int[]
    for sector in 1:max_sector
        if is_in_sector(α, δ, sector)
            push!(sectors, sector)
        end
    end
    return sectors
end

function find_tess_sectors(star_name, max_sector)
    gaia_df = load_star_gaia_data(star_name)
    α, δ = gaia_df.ra, gaia_df.dec
    sectors = Int[]
    for sector in 1:max_sector
        if is_in_sector(α, δ, sector)
            push!(sectors, sector)
        end
    end
    return sectors
end

function plot_cuts(star_name, sector, cut_width, cut_height)
    cut_size = cut_height
    fits = load_tess_cutouts(star_name, cut_width, cut_height)[sector]
    reference_px = [read_key(fits[2], "1CRPX4")[1], read_key(fits[2], "2CRPX4")[1]]
    reference_radec = [read_key(fits[2], "1CRVL4")[1], read_key(fits[2], "2CRVL4")[1]]
    conversion_matrix_px_to_radec = [read_key(fits[2], "11PC4")[1] read_key(fits[2], "12PC4")[1]
                                 read_key(fits[2], "21PC4")[1] read_key(fits[2], "22PC4")[1]]

    conversion_matrix_radec_to_px = inv(conversion_matrix_px_to_radec)

    flux_cuts = read(fits[2], "FLUX")
    n_cuts = size(flux_cuts)[3]

    fig = Figure()

    ax_cut = Axis(fig[1:2,1:2], title = star_name, aspect = DataAspect())
    xlims!(ax_cut, (0.5, cut_width + 0.5))
    ylims!(ax_cut, (0.5, cut_height + 0.5))
    ax_arrows = Axis(fig[2,0], aspect = DataAspect(), title = @sprintf "1 px = %4.1f\"" norm(conversion_matrix_px_to_radec * [1.0, 0.0])*3600)
    ax_light_curve = Axis(fig[3,0:3], yreversed = true, ylabel = "TESS magnitude", xlabel = "MJD")

    cut_slider = Slider(fig[4, 0:3], range = 1:n_cuts, startvalue = 500)

    next_button = Button(fig[6,0], label = "Next")
    prev_button = Button(fig[6,3], label = "Prev")
    bkg_check= Checkbox(fig[6,1], checked = false, tellwidth = false, halign = :right)
    bkg_check_label = Label(fig[6,2], "Background pixels", halign = :left, tellwidth = false)

    df_lc = load_light_curve(star_name, sector, cut_size)
    phot_flux = df_lc.FLUX; mjds = df_lc.MJD
    df_star = load_star_gaia_data(star_name)
    frame_stars_df = load_gaia_stars_in_view_data(star_name, fits, 5)

    star_index = findfirst(x -> x == df_star.source_id, frame_stars_df.source_id)
    stars_x = frame_stars_df.px_x; stars_y = frame_stars_df.px_y
    star_px = frame_stars_df.px_x[star_index], frame_stars_df.px_y[star_index]

    stars_mag = frame_stars_df.phot_rp_mean_mag

    prf = get_tesscut_prf_supersampled(fits)
    bkg_pixels = find_background_prf(flux_cuts[:,:,n_cuts÷4], prf, stars_x, stars_y)

    max_flux_i_cut = findmax(df_lc.FLUX)[2]
    # xlabel!(ax_cut, @sprintf "1 px = %4.1f\"" norm(conversion_matrix_px_to_radec * [1.0, 0.0])*3600)

    Δδ = conversion_matrix_radec_to_px[:,2]/180
    Δα = conversion_matrix_radec_to_px[:,1]/180

    xlims!(ax_arrows, -1.2, 1.2)
    ylims!(ax_arrows, -1.2, 1.2)

    α_arrow_label = Δα + 0.3*Δδ
    δ_arrow_label = Δδ + 0.3*Δα


    text!(ax_arrows, α_arrow_label..., text = "α", align = (:center, :center))
    text!(ax_arrows, δ_arrow_label..., text = "δ", align = (:center, :center))

    arrows2d!(ax_arrows, [0, 0], [0, 0], [Δα[1], Δδ[1]], [Δα[2], Δδ[2]])
    hidedecorations!(ax_arrows)
    hidespines!(ax_arrows)
    # ax2 = Axis(fig[1,2])

    min_mag = minimum(stars_mag)
    max_mag = maximum(stars_mag)
    min_int_mag = round(Int, min_mag)
    max_int_mag = round(Int, max_mag)
    n_sizes = max_int_mag - min_int_mag + 1
    stars_int_mag = round.(Int, stars_mag)

    sizes_groups_stars_x = [stars_x[stars_int_mag .== int_mag] for int_mag = max_int_mag:-1:min_int_mag]
    sizes_groups_stars_y = [stars_y[stars_int_mag .== int_mag] for int_mag = max_int_mag:-1:min_int_mag]
    sizes = (maximum(stars_mag) .- stars_mag)

    i_cut = Observable(500)

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

    cut_hm_data = lift(i_cut) do i_cut
        bkg_cut = fit_flat_background(flux_cuts[:,:,i_cut], bkg_pixels)
        # bkg = get_n_min_mean_background(flux_cuts[:,:,i_cut], n_px ÷ 4)
        log10.(abs.(flux_cuts[:,:,i_cut] - bkg_cut))
    end



    cut_slider_label = Label(fig[5, 0:3], text = @lift @sprintf("MJD = %.4f, i_cut = %d, mag = %.2f", mjds[$i_cut], $i_cut, df_lc.MAG[$i_cut]))

    # bkg = get_n_min_mean_background(flux_cuts[:,:,500], n_px ÷ 4)
    # bkg = 160
    max_flux_bkg = fit_flat_background(flux_cuts[:,:,max_flux_i_cut], bkg_pixels)
    max_flux_hm_data = flux_cuts[:,:,max_flux_i_cut] - max_flux_bkg
    max_flux_hm_data_no_bkg = flux_cuts[:,:,max_flux_i_cut]

    hm = heatmap!(ax_cut, cut_hm_data, colorrange = (0,max(minimum(max_flux_hm_data), log10(1.5e5))))
    # sc = scatter!(ax_cut, stars_x, stars_y, markersize = 5*sizes, color = :lightgray)
    for i_size = 1:n_sizes
        scatter!(ax_cut, sizes_groups_stars_x[i_size], sizes_groups_stars_y[i_size], 
                                markersize = (25 ÷ n_sizes)*i_size, color = :lightgray, label = string(max_int_mag - i_size + 1))
    end

    scatter!(ax_cut, [index[1] for index in bkg_pixels], [index[2] for index in bkg_pixels]; 
            color = :red, marker = :xcross, alpha = @lift($(bkg_check.checked) ? 1.0 : 0.0))

    scatter!(ax_cut, star_px...; marker = :cross, color = :magenta, label = star_name)
    Colorbar(fig[1:2,3], hm, label = "lg TESS flux")
    Legend(fig[1,0], ax_cut, "GAIA R mag")

    slider_time = @lift mjds[$i_cut]
    slider_flux = @lift phot_flux[$i_cut]

    lines!(ax_light_curve, mjds, calc_tess_magnitude.(phot_flux))
    vlines!(ax_light_curve, slider_time; color = :red)
    scatter!(ax_light_curve, @lift Point2f(mjds[$i_cut], calc_tess_magnitude(phot_flux[$i_cut])); color = :red)


    # Colorbar(fig[1,2], sc)
    fig
end

function delete_nans(jds, fluxs)
    # n_fluxs = length(fluxs)
    cleaned_jds = Float64[]
    cleaned_fluxs = Float64[]

    for (jd, flux) in zip(jds, fluxs)
        if !isnan(flux)
            push!(cleaned_jds, jd)
            push!(cleaned_fluxs, flux)
        end
    end

    cleaned_jds, cleaned_fluxs
end

function box_smooth(jds, fluxs, jd_box)
    n_fluxs = length(fluxs)
    smooth_fluxs = zeros(n_fluxs)
    
    for i_flux = 1:n_fluxs
        jd = jds[i_flux]
        smooth_fluxs[i_flux] = median(fluxs[@. abs.(jds .- jd) < jd_box])
    end
    return smooth_fluxs
end

function clean_flux_sigma!(jds, fluxs, jd_box, σ_tol, n_out)
    n_fluxs = length(fluxs)
    smooth_fluxs = box_smooth(jds, fluxs, jd_box)
    anomalous_fluxs = fluxs .- smooth_fluxs

    noout_anomalous_fluxs = sort(anomalous_fluxs)[1:end-n_out]

    σ = √(varm(noout_anomalous_fluxs, mean(noout_anomalous_fluxs)))

    for i_flux = 1:n_fluxs
        if abs(anomalous_fluxs[i_flux]) > σ * σ_tol
            # println(σ, " ", anomalous_fluxs[i_flux])
            fluxs[i_flux] = NaN 
        end
    end


end

function clean_flux!(mjd, flux; pred_tol = 0.1, jd_tol = 0.1)
    n_flux = length(flux)

    i_flux = 3

    while i_flux <= n_flux
        pred_point = 2flux[i_flux-1] - flux[i_flux-2]

        if (abs(flux[i_flux] - pred_point) > pred_tol)
            # println(mjd[i_flux])
            flux[i_flux] = NaN
            i_flux += 3
        end

        if (mjd[i_flux] - mjd[i_flux-1]) > jd_tol
            flux[i_flux] = NaN
        end
        i_flux += 1
    end

    if abs(flux[1] - 2flux[2] + flux[3]) > pred_tol
        flux[1] = NaN
        flux[2] = NaN
    end

    if abs(flux[end] - 2flux[end-1] + flux[end-2]) > pred_tol
        flux[end] = NaN
        flux[end-1] = NaN
    end
end

function get_true_jd(mjd)
    mjd + 2457000
end

function save_lc_figure(star_name, sector, cut_size; day_step = 2, jd_box = 0.3, σ_tol = 5, n_out = 10)
    df_lc = load_light_curve(star_name, sector, cut_size; rewrite_file = true)
    nospace_star_name = get_nospace_star_name(star_name)

    JDs = get_true_jd.(df_lc.MJD)

    JD_start_tick_pos = (round(Int, JDs[1]) ÷ day_step)*day_step
    JD_end_tick_pos = (round(Int, JDs[end]) ÷ day_step)*day_step
    JD_tick_positions = [JD_start_tick_pos:day_step:JD_end_tick_pos;]
    JD_tick_labels = Dates.format.(julian2datetime.(JD_tick_positions), "dd/mm/yy")
    # date_times = julian2datetime.(JDs)

    fig = Figure(size = (1500,475))
    ax_lc = Axis(fig[1,1], yreversed = true, 
                 xticks = (JD_tick_positions, JD_tick_labels),
                 xticklabelrotation = π/4,
                 limits = (JD_start_tick_pos-day_step/2, JD_end_tick_pos+day_step/2, nothing, nothing),
                 aspect = 4,
                 title = "$star_name, sector $sector")

    # lines!(ax_lc, JDs, df_lc.MAG; linestyle = :dash)
    nonan_jd, nonan_mag = delete_nans(JDs, df_lc.MAG)
    clean_flux_sigma!(nonan_jd, nonan_mag, jd_box, σ_tol, n_out)
    lines!(ax_lc, nonan_jd, nonan_mag)
    
    fig
    mkpath("plots/png/$(cut_size)x$(cut_size)")
    save("plots/png/$(cut_size)x$(cut_size)/$(nospace_star_name)_sector_$(sector).png", fig)
end

function find_sampling(jds)
    return median(jds[2:end] - jds[1:end-1])
end

function find_ACF(jds, fluxs; oversample = 1.5)
    sampling = find_sampling(jds)
    time_length = jds[end] - jds[1]
    lags = collect(sampling:sampling/oversample:time_length)
    n_lags = length(lags)

    times = jds .- jds[1]
    # times_interpolated

    flux_spl = Spline1D(times, fluxs, k = 1)
    interpolated_fluxs = flux_spl.(lags)

    mean_flux = mean(interpolated_fluxs)
    flux_dispersion = sum((interpolated_fluxs .- mean_flux) .^ 2) 
    acf = zeros(n_lags)

    for i_lag = 1:n_lags
        acf[i_lag] = sum((interpolated_fluxs[1:end - i_lag] .- mean_flux) .* (interpolated_fluxs[1+i_lag:end] .- mean_flux))
    end

    return lags, acf / flux_dispersion
end

function get_all_data(star_names, cut_size; rewrite_files = false, save_lc_kwargs...)
    open("get_all_data.log", "w") do log_io
        for star_name in star_names[1:end]
            println("Star $star_name")
            sectors = try 
                find_tess_sectors(star_name, tess_max_sectors)
            catch e
                showerror(log_io, e, catch_backtrace())
                Int[]
            end
            if isempty(sectors)
                continue
            end
            println(log_io, "Resolved in TESS sectors ", join(string.(sectors), ", "))
            for sector in sectors
                print(log_io, "Processing sector $sector: ")
                try
                    load_light_curve(star_name, sector, cut_size; rewrite_file = rewrite_files)
                    print(log_io, "light curve loaded")
                    save_lc_figure(star_name, sector, cut_size; save_lc_kwargs...)
                    print(log_io, ", plot saved\n")
                catch e
                    showerror(log_io, e, catch_backtrace())
                end
            end
            print(log_io, "\n")
            flush(log_io)
        end
    end
end

# df_stars = get_simbad_young_stars(12, 8)
# CSV.write("young_simbad.csv", df_stars)
# df_stars = CSV.read("isolated_224.csv", DataFrame)

# function plot_cutout(star_name, i_sector, i_cut, )


# isolated_df = CSV.read("isolated_334.csv", DataFrame)

# star_names = isolated_df.star_name

# cut_size = 15

# for star_name in star_names[106:end]
#     println(star_name)
# # star_name = "FU Ori"
    
#     sectors = find_tess_sectors(star_name, tess_max_sectors)
#     for sector in sectors
#         load_light_curve(star_name, sector, cut_size; rewrite_file = false)
#         save_lc_figure(star_name, sector, cut_size; day_step = 1, jd_box = 0.1, σ_tol = 10)
#     end
# end

# star_name = "V501 Aur"
# sector = 13

# df_lc = load_light_curve(star_name, sector, 15)

# jds, mags = delete_nans(get_true_jd.(df_lc.MJD), df_lc.MAG)

# fig = Figure()
# ax = Axis(fig[1,1], yreversed = true)
# ax_an = Axis(fig[1,2])

# smoothed_mags = box_smooth(jds, mags, 0.3)
# lines!(ax, jds, mags)
# lines!(ax, jds, smoothed_mags)

# an_mags = mags .- smoothed_mags
# mean_an = median(an_mags)

# n_an = length(an_mags)

# an_mags_percentile = sort(an_mags)[1:end-20]

# σ = √(varm(an_mags_percentile, mean_an))

# scatter!(ax_an, jds, an_mags)
# hlines!(ax_an, [mean_an, mean_an - 10σ, mean_an + 10σ])
# fig
# sampling = find_sampling(df_lc.MJD)

# cleaned_mags = df_lc.MAG
# jds = get_true_jd.(df_lc.MJD)
# time_length = jds[end] - jds[1]
# clean_flux!(jds, cleaned_mags; pred_tol = 0.1)


# cleaned_jds = jds[@. !isnan(cleaned_mags)]
# cleaned_mags = df_lc.MAG[@. !isnan(cleaned_mags)]



# # cleaned_mags = cleaned_mags[@. cleaned_jds > 2459370]
# # cleaned_jds = cleaned_jds[@. cleaned_jds > 2459370]


# n_mags = length(cleaned_mags)


# lags, ACF = find_ACF(cleaned_jds, cleaned_mags)

# plan = LombScargle.plan(cleaned_jds, cleaned_mags .- median(cleaned_mags))
# pgram = lombscargle(plan)

# freq, power = freqpower(pgram)
# period = findmaxperiod(pgram, [0.1,10])[1]
# phased_jds = cleaned_jds .% period

# jd_boxes = [0:2*sampling:period;]
# n_boxes = length(jd_boxes)
# boxed_mags = zeros(n_mags)

# for i_mag = 1:n_mags
#     boxed_mags[i_mag] = median(cleaned_mags[@. abs(phased_jds - phased_jds[i_mag]) < period*0.01])
# end

# mean_mag = mean(cleaned_mags)
# mag_dispersion = varm(cleaned_mags, mean_mag)
# mean_phased = mean(cleaned_mags - boxed_mags)
# mag_phased_dispersion = varm(cleaned_mags - boxed_mags, mean_phased)

# periodicity = mag_phased_dispersion/mag_dispersion
# println(periodicity)

# fig = Figure()
# ax_acf = Axis(fig[1,1], xticks = 0:2:30)
# ax_lc = Axis(fig[2,1], yreversed = true)
# ax_phased = Axis(fig[2,2], yreversed = true)

# pgram_ticks = [0.1:0.1:1..., 1.5, 2, 3, 5, 10]
# get_pgram_tick_label(p) = abs(round(Int,p) - p) > 0.01 ? string(p) : string(round(Int,p))
# pgram_tick_labels = get_pgram_tick_label.(pgram_ticks)
# ax_pgram = Axis(fig[1,2],
#                 xticks = (1 ./ pgram_ticks, pgram_tick_labels))
# xlims!(ax_pgram, 0, 3)

# lines!(ax_acf, lags, ACF)
# lines!(ax_lc, cleaned_jds, cleaned_mags)
# lines!(ax_pgram, freq, power)



# scatter!(ax_phased, phased_jds, cleaned_mags)
# scatter!(ax_phased, phased_jds, boxed_mags, color = :red)
# fig



# plot_cuts(star_name, sectors[1], cut_size, cut_size)
# star_gaia_df = load_star_gaia_data(star_name)
# ra, dec = star_gaia_df.ra, star_gaia_df.dec
# for 

# begin
#     cut = 800
#     prf = get_tesscut_prf_supersampled(fits)
#     bkg = get_n_min_median_background(flux_cuts[:,:,cut], size(flux_cuts)[1]*size(flux_cuts)[2] ÷ 4)
#     start_pars = fill(100.0, n_stars + 4)
#     res = fit_stars_prf_flat_bkg(prf, flux_cuts[:,:,cut], stars_x, stars_y, start_pars)
#     star_flux = res.minimizer[star_index]
#     bkg_plane = res.minimizer[end-3:end]
#     bkg_plane[1:3] /= √(sum(bkg_plane[1:3] .^ 2))

#     cut_width, cut_height = size(flux_cuts)[1:2]

#     bkg_cut = zeros(cut_width, cut_height)
#     for x_px = 1:cut_width, y_px = 1:cut_height
#         bkg_cut[x_px, y_px] = (bkg_plane[end] - bkg_plane[1]*x_px - bkg_plane[2]*y_px)/bkg_plane[3]
#     end
#     # bkg_cut .= bkg

#     PRFs = [get_prf_cut(prf, 15, 15, stars_x[i_star], stars_y[i_star])*abs(res.minimizer[i_star]) for i_star = 1:n_stars]
#     bkg = res.minimizer[end]
#     PRF_cut = sum(PRFs) .+ bkg_cut

#     fig_prf = Figure()
#     ax_prf = Axis(fig_prf[1,1], aspect = DataAspect(), title = "PRF Model")
#     hm_prf = heatmap!(ax_prf, log10.(abs.(PRF_cut .- bkg_cut)), colorrange = (2, 4.5))    
#     Colorbar(fig_prf[1,2], hm_prf, label = "lg TESS Flux")

#     ax_flux = Axis(fig_prf[2,1], aspect = DataAspect(), title = "Observed")
#     hm_flux = heatmap!(ax_flux, log10.(abs.(flux_cuts[:,:,cut] .- bkg_cut)), colorrange = (2, 4.5)) 
#     Colorbar(fig_prf[2,2], hm_flux, label = "lg TESS Flux")
   

#     ax_diff = Axis(fig_prf[1,3], aspect = DataAspect(), title = "(Observed - Model) / star_flux")
#     hm_diff = heatmap!(ax_diff, (flux_cuts[:,:,cut].- PRF_cut) ./ star_flux)
#     Colorbar(fig_prf[1,4], hm_diff)

#     ax_1prf = Axis(fig_prf[2,3], aspect = DataAspect(), title = "PRF")
#     hm_1prf = heatmap!(ax_1prf, log10.(get_prf_cut(prf, cut_width, cut_height, star_px...)))
#     Colorbar(fig_prf[2,4], hm_1prf)

#     for i_size = 1:n_sizes
#         scatter!(ax_flux, sizes_groups_stars_x[i_size], sizes_groups_stars_y[i_size], 
#                                 markersize = (25 ÷ n_sizes)*i_size, color = :lightgray, label = string(max_int_mag - i_size + 1))
#     end
    
#     scatter!(ax_flux, star_px...; marker = :cross, color = :magenta, label = star_name)

#     for i_size = 1:n_sizes
#         scatter!(ax_prf, sizes_groups_stars_x[i_size], sizes_groups_stars_y[i_size], 
#                                 markersize = (25 ÷ n_sizes)*i_size, color = :lightgray, label = string(max_int_mag - i_size + 1))
#     end
    
#     scatter!(ax_prf, star_px...; marker = :cross, color = :magenta, label = star_name)

#     Legend(fig_prf[1,0], ax_flux, "GAIA R mag", tellwidth = false)

#     fig_prf
# end

# begin
#     prf = get_tesscut_prf_supersampled(fits)
#     prf_lc = zeros(n_cuts)
#     n_bright = min(n_stars, 20)
#     bright_indeces = sortperm(stars_mag)[1:n_bright]
#     star_bright_index = findfirst(x -> x == gaia_data.source_id, df_gaia_stars.source_id[bright_indeces])
#     bright_stars_x = stars_x[bright_indeces]
#     bright_stars_y = stars_y[bright_indeces]

#     start_pars = fill(100.0, n_bright + 4)
#     iter = ProgressBar(1:n_cuts)
#     for cut = iter
#         start_pars[1:n_bright] = fill(100.0, n_bright)
#         res = fit_stars_prf_flat_bkg(prf, flux_cuts[:,:,cut], bright_stars_x, bright_stars_y, start_pars)
#         start_pars .= res.minimizer
#         prf_lc[cut] = res.minimizer[star_bright_index]
#         set_postfix(iter, Flux=@sprintf("%.2f", prf_lc[cut]))
#     end
# end

# begin 
#     fig_prflc = Figure()
#     ax_prflc = Axis(fig_prflc[1,1], yreversed = true)
#     lines!(ax_prflc, mjds, calc_tess_magnitude.(phot_flux), label = "Aperture")
#     lines!(ax_prflc, mjds, calc_tess_magnitude.(abs.(prf_lc)), label = "TESS PRF")
#     Legend(fig_prflc[1,1], ax_prflc, tellwidth = false, valign = :top, halign = :right)
#     fig_prflc
# end

# begin
#     fig_app_err = Figure()
#     ax_app_err = Axis(fig_app_err[1,1])
#     xlims!(ax_app_err, (9,12))
#     ax_lc_err = Axis(fig_app_err[1,3], yreversed = true)
#     ylims!(ax_lc_err, (12,9))
#     sc_app_err = scatter!(ax_app_err, calc_tess_magnitude.(abs.(prf_lc)), calc_tess_magnitude.(abs.(prf_lc)) - calc_tess_magnitude.(phot_flux), color = mjds)
#     Colorbar(fig_app_err[1, 2], sc_app_err)
#     lines!(ax_lc_err, mjds, calc_tess_magnitude.(phot_flux), label = "Aperture")
#     lines!(ax_lc_err, mjds, calc_tess_magnitude.(abs.(prf_lc)), label = "TESS PRF")
#     fig_app_err
# end

# df_isolated = begin
#     df_isolated = DataFrame()
#     n_stars = nrow(df_stars)
    
#     for i_star = 1:n_stars
#         star_name = df_stars.star_name[i_star]
#         star_ra = df_stars.ra[i_star]
#         star_dec = df_stars.dec[i_star]

#         is_isolated = check_star_for_isolation(star_ra, star_dec, 2/60, 2/60, 4, gaia = "dr2")
#         if is_isolated
#             println("$star_name is isolated")
#             push!(df_isolated, df_stars[i_star, :])
#         else
#             println("$star_name is not isolated")
#         end
#     end
#     df_isolated
# end


# begin
#     i_star = 60
#     star_name = df_stars.star_name[i_star]
#     star_ra = df_stars.ra[i_star]
#     star_dec = df_stars.dec[i_star]

#     is_isolated = check_star_for_isolation(star_ra, star_dec, 2/60, 2/60, 4, gaia = "dr2")
    
#     if is_isolated
#         println("$star_name is isolated")
#         tess_sectors = get_tess_sectors(star_ra, star_dec)
#         if length(tess_sectors) > 0
#             println("$star_name: $(length(tess_sectors)) TESS sectors")
#             get_tess_cutouts(star_ra, star_dec, 15, 15, star_name = star_name)
#             get_lcs(star_name)
#         else
#             println("$star_name: No TESS sectors")
#         end
#     else
#         println("$star_name isn't isolated")
#     end
# end

# star_ra, star_dec = get_star_coords("TW Hya")
# tess_sectors = get_tess_sectors(star_ra, star_dec)
# get_tess_cutouts(star_ra, star_dec, 15, 15, star_name = "TW Hya")
# check_star_for_isolation.(["GM Aur", "TW Hya"], 2/60, 2/60, 3)

