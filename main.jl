using VirtualObservatory
using DataFrames
using HTTP
using ZipArchives
using Mmap
# using BufferedStreams
using JSON3
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

star_directory = "stars_julia"

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

function get_tesscut_corners(fits)
    reference_px = [read_key(fits[2], "1CRPX4")[1], read_key(fits[2], "2CRPX4")[1]]
    reference_radec = [read_key(fits[2], "1CRVL4")[1], read_key(fits[2], "2CRVL4")[1]]
    cos_dec = cos(reference_radec[2]/180*π)
    conversion_matrix_px_to_radec = [read_key(fits[2], "11PC4")[1] read_key(fits[2], "12PC4")[1]
                                     read_key(fits[2], "21PC4")[1] read_key(fits[2], "22PC4")[1]]

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

function calc_tess_magniude(flux)
    return -2.5*log10(flux) + 20.44
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
        σ = pars[end]
        model = deepcopy(cut_no_bkg)
        for i_star in 1:n_stars
            model = model .- gaussian_PSF(width, height, stars_px_x[i_star], stars_px_y[i_star], σ, abs(max_ints[i_star]); super_samp = super_samp)
        end
        return vec(model)
    end

    start_pars = fill(100.0, n_stars+1)
    start_pars[end] = 3

    optimize(to_optimize, start_pars, LevenbergMarquardt())
end


# df_stars = get_simbad_young_stars(12, 8)
# CSV.write("young_simbad.csv", df_stars)
df_stars = CSV.read("isolated_224.csv", DataFrame)

# function plot_cutout(star_name, i_sector, i_cut, )


begin 
star_name = "RY Lup"

gaia_data = get_star_gaia_data(star_name)
ra, dec = gaia_data[[:ra, :dec]]
nospace_star_name = get_nospace_star_name(star_name)
if !isfile("$star_directory/$nospace_star_name/$nospace_star_name.zip")
    get_tess_cutouts(ra, dec, 15, 15; star_name = star_name)
end
end 

begin
sectors = get_tess_sectors(star_name)
fits = get_star_tesscut_fits(star_name, sectors[1])
flux_cuts = read(fits[2], "FLUX")
n_px = size(flux_cuts)[1]*size(flux_cuts)[2]
n_cuts= size(flux_cuts)[3]

mjds = read(fits[2], "TIME")

reference_px = [read_key(fits[2], "1CRPX4")[1], read_key(fits[2], "2CRPX4")[1]]
reference_radec = [read_key(fits[2], "1CRVL4")[1], read_key(fits[2], "2CRVL4")[1]]
cos_dec = cos(reference_radec[2]/180*π)
conversion_matrix_px_to_radec = [read_key(fits[2], "11PC4")[1] read_key(fits[2], "12PC4")[1]
                                 read_key(fits[2], "21PC4")[1] read_key(fits[2], "22PC4")[1]]

conversion_matrix_radec_to_px = inv(conversion_matrix_px_to_radec)

corners = get_tesscut_corners(fits)
Δm_R = 5
star_R = gaia_data.phot_rp_mean_mag
df_gaia_stars = get_gaia_stars_in_poly(corners, star_R + Δm_R; gaia = "dr2")
n_stars = nrow(df_gaia_stars)

stars_x = zeros(n_stars)
stars_y = zeros(n_stars)
stars_mag = zeros(n_stars)

star_px = reference_px + conversion_matrix_radec_to_px * get_rel_radec(reference_radec..., ra, dec)

for i_star = 1:n_stars
    star_radec = [df_gaia_stars.ra[i_star], df_gaia_stars.dec[i_star]]
    delta_radec = get_rel_radec(reference_radec..., star_radec...)
    delta_radec[1] = delta_radec[1]
    stars_mag[i_star] = df_gaia_stars.phot_rp_mean_mag[i_star]
    stars_x[i_star], stars_y[i_star] = reference_px + conversion_matrix_radec_to_px * delta_radec
end

cuts = [flux_cuts[:,:,i_cut] for i_cut = 1:n_cuts]

phot_flux = calc_aperture_photometry.(cuts, star_px..., 2.5)
end

begin 
fig = Figure()
ax_cut = Axis(fig[1:2,1:2], title = star_name, aspect = DataAspect())
ax_arrows = Axis(fig[2,0], aspect = DataAspect(), title = @sprintf "1 px = %4.1f\"" norm(conversion_matrix_px_to_radec * [1.0, 0.0])*3600)
ax_light_curve = Axis(fig[3,0:3], yreversed = true, ylabel = "TESS magnitude", xlabel = "MJD")

cut_slider = Slider(fig[4, 0:3], range = 1:n_cuts, startvalue = 500)

next_button = Button(fig[6,0:1], label = "Next")
prev_button = Button(fig[6,2:3], label = "Prev")

max_flux_i_cut = findmax(phot_flux)[2]
# xlabel!(ax_cut, @sprintf "1 px = %4.1f\"" norm(conversion_matrix_px_to_radec * [1.0, 0.0])*3600)

Δδ = conversion_matrix_radec_to_px[:,2]/180
Δα = conversion_matrix_radec_to_px[:,1]/180

xlims!(ax_arrows, -1.2, 1.2)
ylims!(ax_arrows, -1.2, 1.2)

α_arrow_label = Δα + 0.3*Δδ
δ_arrow_label = Δδ + 0.3*Δα


text!(ax_arrows, α_arrow_label..., text = "α", align = (:center, :center))
text!(ax_arrows, δ_arrow_label..., text = "δ", align = (:center, :center))

arrows!(ax_arrows, [0, 0], [0, 0], [Δα[1], Δδ[1]], [Δα[2], Δδ[2]])
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
    bkg = get_n_min_mean_background(flux_cuts[:,:,i_cut], n_px ÷ 4)
    log10.(abs.(flux_cuts[:,:,i_cut] .- bkg))
end



cut_slider_label = Label(fig[5, 0:3], text = @lift "MJD = "*string(mjds[$i_cut])*", i_cut = "*string($i_cut))

# bkg = get_n_min_mean_background(flux_cuts[:,:,500], n_px ÷ 4)
# bkg = 160
max_flux_bkg = get_n_min_mean_background(flux_cuts[:,:,max_flux_i_cut], n_px ÷ 4)
max_flux_hm_data = flux_cuts[:,:,max_flux_i_cut] .- max_flux_bkg
max_flux_hm_data_no_bkg = flux_cuts[:,:,max_flux_i_cut]

hm = heatmap!(ax_cut, cut_hm_data, colorrange = (0,max(minimum(max_flux_hm_data), log10(1.5e5))))
# sc = scatter!(ax_cut, stars_x, stars_y, markersize = 5*sizes, color = :lightgray)
for i_size = 1:n_sizes
    scatter!(ax_cut, sizes_groups_stars_x[i_size], sizes_groups_stars_y[i_size], 
                            markersize = (25 ÷ n_sizes)*i_size, color = :lightgray, label = string(max_int_mag - i_size + 1))
end

scatter!(ax_cut, star_px...; marker = :cross, color = :magenta, label = star_name)
Colorbar(fig[1:2,3], hm, label = "lg TESS flux")
Legend(fig[1,0], ax_cut, "GAIA R mag")

slider_time = @lift mjds[$i_cut]
slider_flux = @lift phot_flux[$i_cut]

lines!(ax_light_curve, mjds, calc_tess_magniude.(phot_flux))
vlines!(ax_light_curve, slider_time; color = :red)
scatter!(ax_light_curve, @lift Point2f(mjds[$i_cut], calc_tess_magniude(phot_flux[$i_cut])); color = :red)


# Colorbar(fig[1,2], sc)
fig
end

begin
    fig3d = Figure()
    ax = Axis3(fig3d[1,1], aspect = :equal)
    surface!(ax, max_flux_hm_data_no_bkg)
    # scale!(ax.scene, 1, 1, 5)
    fig3d
end

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

