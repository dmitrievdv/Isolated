from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

Simbad.add_votable_fields("G")
res = Simbad.query_object("V 4046 Sgr")

print(res["G"][0])

Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"

star_ra = res["ra"][0]
star_dec = res["dec"][0]

coord = SkyCoord(ra=res["ra"][0], dec=res["dec"][0], unit=(u.degree, u.degree), frame='icrs')

width = u.Quantity(2, u.arcmin)
height = u.Quantity(2, u.arcmin)

r1 = Gaia.query_object_async(coordinate=coord, width=width, height=height)
job = Gaia.launch_job("select * from gaiadr2.gaia_source where CONTAINS(POINT('ICRS', ra, dec), BOX('ICRS', "+
                    str(star_ra)+", "+str(star_dec)+", "+ str(2/60)+", "+ str(2/60)+")) = 1 AND ra IS NOT NULL AND dec IS NOT NULL")
r2 = job.get_results()

r1[["dist", "designation", "ra", "dec", "phot_g_mean_mag"]].pprint(max_lines=12, max_width=130)
r2[["designation", "ra", "dec", "phot_g_mean_mag"]].pprint(max_lines=12, max_width=130)
# central_star_mag = r["phot_g_mean_mag"][0]
# mag_threshold = 3
# bright_indeces = []

# for star_index, star_mag in enumerate(r["phot_g_mean_mag"][1:]):
#     if star_mag < central_star_mag + mag_threshold:
#         bright_indeces.append(star_index + 1)

# r[["dist", "designation", "ra", "dec", "phot_g_mean_mag"]][bright_indeces].pprint(max_lines=12, max_width=130)
# print(r.keys())