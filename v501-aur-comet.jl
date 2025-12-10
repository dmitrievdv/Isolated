include("main.jl")

star_name = "V501 Aur"
sector = 43
cut_size = 30

gaia_data = load_star_gaia_data(star_name)
# delta_radec = get_rel_radec(reference_radec..., star_radec...)

fits = load_tess_cutouts(star_name, cut_size)[sector]
flux_cuts = read(fits[2], "FLUX")

conversion_matrix_px_to_radec = [read_key(fits[2], "11PC4")[1] read_key(fits[2], "12PC4")[1]
                                 read_key(fits[2], "21PC4")[1] read_key(fits[2], "22PC4")[1]]

conversion_matrix_radec_to_px = inv(conversion_matrix_px_to_radec)
reference_px = [read_key(fits[2], "1CRPX4")[1], read_key(fits[2], "2CRPX4")[1]]
reference_radec = [read_key(fits[2], "1CRVL4")[1], read_key(fits[2], "2CRVL4")[1]]

comet_df = CSV.read("COBS - Comet OBServation database.csv", DataFrame)
comet_df.RA .= Dates.format.(comet_df.RA, Ref("HH:MM:SS"))

function convert_cobs_angle(cobs_angle :: AbstractString)
    dms_strings = split(cobs_angle, ':')
    dms_nums = parse.(Float64, dms_strings)
    return dms_nums[1] + dms_nums[2]/60 + dms_nums[3]/3600
end

comet_df.RA .= 15 .* convert_cobs_angle.(comet_df.RA)
comet_df.Dec .= convert_cobs_angle.(comet_df.Dec)
delta_radec = [get_rel_radec(gaia_data.ra, gaia_data.dec, ra, dec) for (ra, dec) in zip(comet_df.RA, comet_df.Dec)]
coords = [conversion_matrix_radec_to_px * radec + reference_px for radec in delta_radec]

comet_x = [coord[1] for coord in coords]
comet_y = [coord[2] for coord in coords]

fig = plot_cuts(star_name, sector, cut_size, cut_size)
scatter!(content(fig[1:2,1:2]), comet_x, comet_y, colormap = :tab10, colorrange = (1,10), color = 2)
fig