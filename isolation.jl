include("main.jl")

df_young_stars = CSV.read("young_simbad.csv", DataFrame)
n_young_stars = nrow(df_young_stars)

isolation_box_size = 2
isolation_mag_lim = 3

df_isolated = DataFrame()

brightest_mags = zeros(n_young_stars)

for i_star = 1:n_young_stars
    star_name = df_young_stars.star_name[i_star]
    star_ra = df_young_stars.ra[i_star]
    star_dec = df_young_stars.dec[i_star]
    max_near_mag = try
        brightest_near_relative_mag(star_ra, star_dec, isolation_box_size/60, isolation_box_size/60; gaia = "dr3")
    catch
        0.0
    end 
    brightest_mags[i_star] = max_near_mag
    println("$star_name, $max_near_mag")
end

df_isolated = deepcopy(df_young_stars)
 
df_isolated[!, "max_near_mag"] = brightest_mags
CSV.write("isolated_$isolation_box_size.csv", df_isolated)

