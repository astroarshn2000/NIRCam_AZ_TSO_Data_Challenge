""" NIRCam AZ TSO Data Challenge (GJ 436b) Simulations
         ~ Arsh R. Nadkarni (UArizona), 2020 """  

# ----------   Set Environment Variables   ----------

import os
os.environ["MIRAGE_DATA"] = "/home/anadkarni/JWST_Programs/mirage_reference_files/mirage_data/"
os.environ["CRDS_PATH"] = os.path.expandvars("$HOME/crds_cache")
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["PYSYN_CDBS"] = "/home/anadkarni/JWST_Programs/stsynphot_reference_files/grp/hst/cdbs/"
os.environ["WEBBPSF_PATH"] = "/home/anadkarni/JWST_Programs/webbpsf_reference_files/webbpsf-data/"


#  ----------   Import Libraries  ----------

from astropy.io import fits, ascii
from astropy.table import Table
from astropy.visualization import simple_norm, imshow_norm
from astropy import units as u
import batman
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams; rcParams["figure.dpi"] = 300
from matplotlib.ticker import (AutoMinorLocator)
from matplotlib import colors
import matplotlib.cm as cmx
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
import stsynphot as stsyn
from synphot import SourceSpectrum, SpectralElement
from synphot import units
import yaml
from mirage.catalogs.hdf5_catalog import save_tso
from mirage.catalogs.catalog_generator import GrismTSOCatalog, ImagingTSOCatalog, PointSourceCatalog
from mirage.catalogs.catalog_generator import TSO_GRISM_INDEX
from mirage.grism_tso_simulator import GrismTSO
from mirage.imaging_simulator import ImgSim
from mirage.seed_image.catalog_seed_image import Catalog_seed
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator


# ----------   Define Paths to organize Inputs and Outputs  ----------

input_data_path = os.path.abspath('/home/anadkarni/NIRCam_AZ_TSO_Data_Challenge/GJ436_apt_data/')
output_dir = '/home/anadkarni/NIRCam_AZ_TSO_Data_Challenge/GJ436/GJ436_2102/'
output_yaml_dir = os.path.abspath('/home/anadkarni/NIRCam_AZ_TSO_Data_Challenge/GJ436_2102/GJ436_yaml_data/')
ensure_dir_exists(output_yaml_dir)
output_data_dir = os.path.abspath('/home/anadkarni/NIRCam_AZ_TSO_Data_Challenge/GJ436_2102/GJ436_sim_data/')
ensure_dir_exists(output_data_dir)
tsdir = '/home/anadkarni/NIRCam_AZ_TSO_Data_Challenge/GJ436_LC_TS_Params/'


# ----------  Prepare Inputs  ----------

xml_file = os.path.join(input_data_path, 'GJ436.xml')
pointing_file = xml_file.replace('.xml', '.pointing')


# ----------  Stellar Spectrum  ----------

t_eff = 3500  # surface temperature
metallicity = 0.02 # Fe/H
log_g = 5.0  # surface gravity = 182 m/s^2
sp = stsyn.grid_to_spec('ck04models', t_eff, metallicity, log_g)
bp = SpectralElement.from_filter('johnson_k')
vega = SourceSpectrum.from_vega()
sp_norm = sp.normalize(6.073 * units.VEGAMAG, bp, vegaspec=vega)
wavelengths = sp_norm.waveset.to(u.micron)
fluxes = sp_norm(wavelengths, flux_unit='flam')
wavelength_units = 'microns'
flux_units = 'flam'
sed_file = os.path.join(output_dir, 'GJ436_stellar_spectrum.hdf5')
fluxes = [fluxes]
wavelengths = [wavelengths]
with h5py.File(sed_file, "w") as file_obj:
    for i in range(len(fluxes)):
        dset = file_obj.create_dataset(str(i+TSO_GRISM_INDEX), data=[wavelengths[i].value, fluxes[i].value],
                                       dtype='f', compression="gzip", compression_opts=9)
        dset.attrs[u'wavelength_units'] = wavelength_units
        dset.attrs[u'flux_units'] = flux_units


# ----------   Batman Parameters and Light-Curve File  ----------

params = batman.TransitParams()       # object to store transit parameters
params.t0 = 0.04165880981289195       # time of inferior conjunction
params.per = 2.643898012982137        # orbital period
params.rp = 0.0828098015829365        # planet radius (in units of stellar radii)
params.a = 14.656751572765119         # semi-major axis (in units of stellar radii)
params.inc = 86.87843069704299        # orbital inclination (in degrees)
params.ecc = 0.                       # eccentricity
params.w = 90.                        # longitude of periastron (in degrees)
params.fp = 0.0
params.t_secondary = 0.0 
params.limb_dark = "quadratic"        # limb darkening model
params.u = [0.019121089782213305, 0.3124784495689501]  # limb darkening coefficients
times = np.arange(0,0.055,0.00005787)
m = batman.TransitModel(params, times)
flux = m.light_curve(params)
lightcurve_file = os.path.join(output_dir, 'GJ436_lightcurve.hdf5')
contents = {}
contents['1'] = {'times': times,
                 'fluxes': flux}
save_tso(contents, lightcurve_file, time_unit='days')


# ----------  Transmission Spectrum  ----------

in_file = os.path.join(tsdir,'GJ436b_trans_MIRAGE_GRISMR2.txt')
tran_spec_file = in_file.replace('GJ436b_trans_MIRAGE_GRISMR2.txt', 'GJ436b_trans_MIRAGE_GRISMR2_rp_over_rs.txt')
tab = ascii.read(in_file)
new_trans = np.sqrt(tab['Transmission'])
tab['Transmission'] = new_trans
ascii.write(tab, tran_spec_file, overwrite=True)


# ---------- Create Grism TSO Catalog   ----------

grism_tso_catalog = os.path.join(output_dir,'tso_grism_source.cat')
object_ra = +175.55054
object_dec = +26.70307
object_f322w2_mag = 5.980825287207619
grism_cat = GrismTSOCatalog(ra=[object_ra], dec=[object_dec], semimajor_axis=[params.a],
                            orbital_inclination=[params.inc], eccentricity=[params.ecc],
                            orbital_period=[params.per], longitude_of_periastron=[params.w],
                            limb_dark_model=[params.limb_dark], limb_dark_coeffs=[params.u],
                            time_units=['days'], start_time=[np.min(times)],
                            end_time=[np.max(times)], inferior_conj=[params.t0],
                            transmission_spectrum=[tran_spec_file])
grism_cat.add_magnitude_column([object_f322w2_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f322w2')
grism_cat.save(grism_tso_catalog)


# ---------- Create Imaging TSO Catalog   ----------

object_f212n_mag = 6.062028376949228
imaging_tso_catalog = os.path.join(output_dir, 'tso_imaging_source.cat')
tsimg_cat = ImagingTSOCatalog(ra=[object_ra], dec=[object_dec], lightcurve_file=[lightcurve_file])
tsimg_cat.add_magnitude_column([object_f212n_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='wlp4')
tsimg_cat.save(imaging_tso_catalog)


# ---------- Create Input YAML Files for Mirage   ----------

catalogs = {'ROSS-905': {'point_source': None,
                         'tso_imaging_catalog': imaging_tso_catalog,
                         'tso_grism_catalog': grism_tso_catalog,
                        }
           }
background = 'medium'
pav3 = 0.
yam = yaml_generator.SimInput(xml_file, pointing_file, catalogs=catalogs, verbose=True,
                              output_dir=output_yaml_dir, simdata_output_dir=output_data_dir,
                              background=background, roll_angle=pav3,
                              dates=None, datatype='linear, raw', dateobs_for_background=False,
                              reffile_defaults='crds')

yam.use_linearized_darks = True
yam.create_inputs()


# ---------- SIMULATE F322W2 GRISMR TSO ----------

gr_tso_yaml_file = os.path.join(output_yaml_dir, 'jw00042001001_01101_00002_nrca5.yaml')
gr_f322w2 = GrismTSO(gr_tso_yaml_file, SED_file=sed_file, SED_normalizing_catalog_column=None,
                    final_SED_file=None, save_dispersed_seed=True, source_stamps_file=None,
                    extrapolate_SED=True, override_dark=None, disp_seed_filename=None,
                    orders=["+1", "+2"])
gr_f322w2.create()


# ---------- SIMULATE WLP4+CLEAR IMAGING TSO ----------

img_tso_sw_yaml = os.path.join(output_yaml_dir, 'jw00042001001_01101_00001_nrca1.yaml')
img_tso = ImgSim()
img_tso.paramfile = img_tso_sw_yaml
img_tso.create()