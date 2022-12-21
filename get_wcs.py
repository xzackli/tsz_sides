# %%

import h5py
import numpy as np
from astropy import wcs
import sys
import astropy.units as u

# TAKEN FROM PYSIDES https://gitlab.lam.fr/mbethermin/sides-public-release/-/blob/main/PYSIDES/pysides/make_maps.py
def set_wcs_map(cat, pixel_size):

    #-----------------SET WCS-----------------
    # Coordinate increment at reference point

    ra =  np.asarray(cat["ra"]) * u.deg
    dec = np.asarray(cat["dec"]) * u.deg
    
    ra_mean, dec_mean = np.mean(ra.value) , np.mean(dec.value) 

    # set a first time the wcs
    pix_resol = pixel_size / 3600.

    ra_cen = 0.5 * (ra.max() + ra.min())
    dec_cen = 0.5 * (dec.max() + dec.min())
    delta_ra = ra.max() - ra.min()
    delta_dec = dec.max() - dec.min()

    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra_cen.value, dec_cen.value]
    w.wcs.crpix = [0.5*delta_ra.value / pix_resol, 0.5*delta_dec.value / pix_resol]
    w.wcs.cdelt = [pix_resol, pix_resol]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = [u.deg, u.deg]

    # compute the position in pixel units
    x, y = w.celestial.wcs_world2pix(ra , dec, 0)

    #Offset the central pixel to have all x>=0 and y>=0
    w.wcs.crpix = [0.5*delta_ra.value / pix_resol - np.min(x), 0.5*delta_dec.value / pix_resol  - np.min(y)]

    #recompute x and y in the new WCS
    x, y = w.celestial.wcs_world2pix(ra , dec, 0)

    # list of the position in pixel units: pos[0] = y, pos[1] = x 
    pos = [y , x]
    
    # compute the dimensions of the three axes
    shape = [np.ceil(pos[0].max()), np.ceil(pos[1].max())]
    shape = [i//2*2+1 for i in shape] #force an odd number of pixels to generate better psf
    x_edges = list(np.arange(-0.5, shape[1] + 0.5, 1))
    y_edges = list(np.arange(-0.5, shape[0] + 0.5, 1))

    wcs_dict = {}
    wcs_dict['w'] = w
    wcs_dict['shape'] = shape
    wcs_dict['pos'] = pos
    wcs_dict['x_edges'] = x_edges
    wcs_dict['y_edges'] = y_edges

    return wcs_dict


# %%  CHANGE ME

abhi_file = 'fullcatalog_data_pidminus1_gzipped.h5'
websky_file = 'little_box.h5'
f = h5py.File(abhi_file, 'r')

cat = {
    "halo_mass": np.array(f["halo_mass"]),
    "ra": np.array(f["ra"]),
    "dec": np.array(f["dec"])
}

c = set_wcs_map(cat, 5)  # 30 arcseconds

wcs_info = {
    'cdelt': (c['w'].wcs.cdelt[0], c['w'].wcs.cdelt[1]),
    'crpix': (c['w'].wcs.crpix[0], c['w'].wcs.crpix[1]),
    'crval': (c['w'].wcs.crval[0], c['w'].wcs.crval[1]),
    'shape': [int(s) for s in c["shape"]]
}

import json
with open('wcs.json', 'w') as fp:
    json.dump(wcs_info, fp)


# %%
