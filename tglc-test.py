import os
import numpy as np
import matplotlib.pyplot as plt
from tglc.quick_lc import tglc_lc
from astroquery.mast import Catalogs

target = 'RY Lup'     # TIC ID (preferred, 'TIC 12345678'), Target ID ('TOI 519') or coordinates ('ra dec')
local_directory = f'{target}/'    # directory to save all files
# os.makedirs(local_directory, exist_ok=True)
# tglc_lc(target=target, 
#         local_directory=local_directory, # directory to output everything
#         size=90, # FFI cutsize. Recommand at least 50 or larger for better performance. Cannot exceed 99. 
#                  # Downloading FFI might take longer (or even cause timeouterror) for larger sizes. 
#         save_aper=True, # whether to save 5*5 pixels timeseries of the decontaminated images in fits file primary HDU
#         limit_mag=15, # the TESS magnitude lower limit of stars to output
#         get_all_lc=False, # whether to return all lcs in the region. If False, return the nearest star to the target coordinate
#         first_sector_only=False, # whether to return only lcs from the sector this target was first observed. 
#                                 # If False, return all sectors of the target, but too many sectors could be slow to download.
#         last_sector_only=True, # whether to return only lcs from the sector this target was last observed. 
#         sector=None, # If first_sector_only = True or last_sector_only = True and type(sector) != int, return first or last sector.
#                      # If first(last)_sector_only=False and sector = None, return all observed sectors
#                      # If first(last)_sector_only=False and type(sector) == int, return only selected sector. 
#                      # (Make sure only put observed sectors. All available sectors are printed in the sector table.)
#         prior=None,  # If None, does not allow all field stars to float. SUGGESTED for first use. 
#                      # If float (usually <1), allow field stars to float with a Gaussian prior with the mean 
#                      # at the Gaia predicted value the width of the prior value multiplied on the Gaia predicted value.
#         transient=None, # If generating light curve at an arbitrary RA and Dec, set this to [`name', `RA', `DEC'].
#                         # Make sure target = 'RA Dec'
#         )

import pickle
with open(f'{local_directory}source/source_RY Lup_sector_12.pkl', 'rb') as input_: # replace with each sector's .pkl file if you have more than the first sector
    source = pickle.load(input_)
    print(f'sector = {source.sector}')
epsf = np.load(f'{local_directory}epsf/epsf_RY Lup_sector_12.npy')

f, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow(np.log10(source.flux[0]))
ax1.set_title('FFI cut (log scale)')
ax2.imshow(epsf[0,:23**2].reshape(23,23))
ax2.set_title('ePSF shape')
plt.show()

# from astropy.io import fits
# hdul_s0014 = fits.open(f'{local_directory}lc/hlsp_tglc_tess_ffi_gaiaid-5996151172781298304-s0012-cam1-ccd1_tess_v2_llc.fits')
# hdul_s0014.info()


# # aperture_sequence = hdul_s0014[0].data
# # plt.imshow(aperture_sequence[0])
# # plt.show()

# q_14 = [a and b for a, b in zip(list(hdul_s0014[1].data['TESS_flags'] == 0),
#                                 list(hdul_s0014[1].data['TGLC_flags'] == 0))]
# # filter out bad datapoints from both TESS FFI flags and TGLC flags

# time_14 = hdul_s0014[1].data['time'][q_14]
# psf_flux_14 = hdul_s0014[1].data['psf_flux'][q_14] # raw psf flux
# psf_flux_err_14 = hdul_s0014[1].header['PSF_ERR'] # raw psf flux error
# aper_flux_14 = hdul_s0014[1].data['aperture_flux'][q_14] # raw aper flux
# aper_flux_err_14 = hdul_s0014[1].header['APER_ERR'] # raw aper flux error
# plt.errorbar(time_14, psf_flux_14, psf_flux_err_14, marker = '', label = 'psf')
# plt.errorbar(time_14, aper_flux_14, aper_flux_err_14, marker = '', label = 'aperture')
# plt.title('WW Vul - TESS Sector 14')
# plt.xlabel('TBJD')
# plt.ylabel('Flux e-/s')
# plt.legend()
# plt.show()

