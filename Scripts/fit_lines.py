"""Module to manage fits to the lines"""
import copy
import numpy as np
from scipy.constants import speed_of_light
import iminuit
from astropy.io import fits
import h5py

from trident.absorption_spectrum.absorption_line import voigt
# from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit?

SPEED_LIGHT = speed_of_light / 1000.0 # km/s

SPECIES_DICTS = {
    'Lya': {
        'f': [.4160],
        'Gamma':[4.690E8],
        'wavelength':[1215.67]
    },
    'Lyb': {
        'f': [.0791],
        'Gamma':[5.570E7],
        'wavelength':[1025.7222]
    },
    'Lyg': {
        'f': [.0290],
        'Gamma':[1.280E6],
        'wavelength':[972.5367]
    },
    'Lyother': { # Lyman series (excluding Lya, Lyb and Lyg)
        'f': [
            2.440e-05, 2.640e-05, 2.850e-05, 3.090e-05, 3.350e-05, 3.650e-05, 3.980e-05,
            4.360e-05, 4.780e-05, 5.260e-05, 5.800e-05, 6.430e-05, 7.140e-05, 7.970e-05,
            8.920e-05, 1.000e-04, 1.140e-04, 1.290e-04, 1.480e-04, 1.700e-04, 1.970e-04,
            2.290e-04, 2.700e-04, 3.210e-04, 3.850e-04, 4.680e-04, 5.770e-04, 7.220e-04,
            9.210e-04, 1.200e-03, 1.600e-03, 2.210e-03, 3.180e-03, 4.810e-03, 7.790e-03,
            1.390e-02
        ],
        'Gamma': [
            1.220e+02, 1.390e+02, 1.580e+02, 1.810e+02, 2.070e+02, 2.390e+02, 2.760e+02,
            3.200e+02, 3.740e+02, 4.380e+02, 5.160e+02, 6.120e+02, 7.290e+02, 8.740e+02,
            1.060e+03, 1.280e+03, 1.580e+03, 1.950e+03, 2.440e+03, 3.070e+03, 3.920e+03,
            5.070e+03, 6.650e+03, 8.850e+03, 1.200e+04, 1.660e+04, 2.340e+04, 3.390e+04,
            5.060e+04, 7.830e+04, 1.260e+05, 2.140e+05, 3.870e+05, 7.560e+05, 1.640e+06,
            4.120e+06
        ],
        'wavelength': [
            912.3236, 912.3532, 912.3853, 912.4199, 912.4575, 912.4983, 912.5428,
            912.5914, 912.6447, 912.7032, 912.7676, 912.8388, 912.9178, 913.0058,
            913.1041, 913.2145, 913.3391, 913.4802, 913.6411, 913.8256, 914.0385,
            914.2860, 914.5762, 914.9192, 915.3289, 915.8237, 916.4290, 917.1805,
            918.1293, 919.3513, 920.9630, 923.1503, 926.2256, 930.7482, 937.8034,
            949.7429
        ]
    },
    'CI': {
        'f': [
            2.120E+01, 1.240E+01, 4.090E+01, 1.250E+01, 3.940E+01, 9.230E+01, 2.610E+01,
            5.800E+01, 7.160E+01, 1.430E+02
        ],
        'Gamma': [ 
            3.520E+10, 1.950E+10, 6.390E+10, 1.940E+10, 5.510E+10, 1.260E+11, 3.550E+10,
            7.310E+10, 6.540E+10, 1.160E+11
        ],
        'wavelength':[
            1157.9100, 1188.8330, 1193.0310, 1193.9950, 1260.7350, 1277.2450, 1280.1350,
            1328.8340, 1560.3090, 1656.9290
        ] 
    },
    'CII': {
        'f': [3.310E-01, 1.190E-01, 1.180E-01, 1.290E-01, 1.270E-02],
        'Gamma':[2.700E+09, 7.380E+08, 1.460E+09, 2.410E+08, 4.760E+07],
        'wavelength':[903.9620, 1036.3370, 1037.0180, 1334.5320, 1335.6630]
    },
    'CIII': {
        'f': [7.586E-01],
        'Gamma':[1.767E+09],
        'wavelength':[977.0200]
    },
    'CIV': {
        'f': [1.900E-01, 9.520E-02],
        'Gamma':[2.650E+08, 2.640E+08],
        'wavelength':[1548.1870, 1550.7720]
    },
    'NI': {
        'f': [
            2.290E-03, 1.970E-03, 1.290E-02, 2.470E-02, 3.310E-02, 4.000E-03, 1.240E-02,
            7.900E-03, 3.860E-03, 9.690E-04, 1.460E-02, 2.870E-02, 3.870E-03, 1.320E-01,
            8.690E-02
        ],
        'Gamma':[
            1.120E+07, 1.450E+07, 1.900E+08, 1.810E+08, 1.620E+08, 1.950E+07, 5.940E+07,
            5.660E+07, 5.520E+07, 3.520E+06, 1.510E+08, 1.490E+08, 1.820E+07, 4.070E+08,
            4.030E+08
        ],
        'wavelength':[
            952.3030, 952.4150, 953.4150, 953.6550, 953.9700, 954.1040, 963.9900,
            964.6260, 965.0410, 1134.1650, 1134.4150, 1134.9800, 1199.5500, 1200.2230,
            1200.7100
        ]
    },
    'NII': {
        'f': [1.590E-01, 1.110E-01],
        'Gamma':[4.230E+08, 2.100E+08],
        'wavelength':[915.6120, 1083.9900]
    },
    'NIII': {
        'f': [1.230E-01],
        'Gamma':[4.180E+08],
        'wavelength':[989.7990]
    },
    'NIV': {
        'f': [6.109E-01],
        'Gamma':[2.320E+09],
        'wavelength':[765.1470]
    },
    'NV': {
        'f': [1.560E-01, 7.800E-02],
        'Gamma':[3.400E+08, 3.370E+08],
        'wavelength':[1238.8210, 1242.8040]
    },
    'OI': {
        'f': [
            2.190E-02, 1.300E-03, 1.920E-03, 3.060E-03, 1.580E-03, 3.310E-03, 8.460E-03,
            4.640E-02, 9.160E-03, 5.200E-02, 5.180E-02, 5.190E-02
        ],
        'Gamma':[
            1.230E+08, 7.220E+06, 1.060E+07, 1.660E+07, 1.940E+07, 3.860E+07, 5.770E+07,
            2.260E+08, 9.430E+07, 3.410E+08, 2.030E+08, 6.760E+07
        ],
        'wavelength':[
            922.0080, 924.9500, 929.5170, 936.6290, 950.8850, 976.4480, 988.6550, 
            988.7730, 1039.2300, 1302.1680, 1304.8580, 1306.0290
        ]
    },
    'OII': {
        'f': [4.510E-02, 9.010E-02, 1.350E-01],
        'Gamma':[8.670E+08, 8.650E+08, 8.610E+08],
        'wavelength':[832.7583, 833.3303, 834.4654]
    },
    'OIII': {
        'f': [1.340E-01, 4.520E-02, 1.060E-01, 8.770E-02, 8.770E-02],
        'Gamma':[6.060E+08, 1.830E+09, 3.410E+08, 5.990E+08, 5.990E+08],
        'wavelength':[702.3370, 702.8380, 832.9290, 835.2890, 835.2890]
    },
    'OIV': {
        'f': [1.120E-01, 2.240E-01, 2.790E-01, 6.710E-02, 6.690E-02, 1.110E-01, 9.940E-02],
        'Gamma':[1.220E+09, 4.860E+09, 6.060E+09, 1.210E+09, 2.400E+09, 5.950E+08, 7.080E+08],
        'wavelength':[553.3290, 554.0760, 554.5130, 608.3970, 609.8290, 787.7100, 790.1990]
    },
    'OV': {
        'f': [5.122E-01],
        'Gamma':[2.872E+09],
        'wavelength':[629.7320]
    },
    'OVI': {
        'f': [1.330E-01, 6.600E-02],
        'Gamma':[4.160E+08, 4.090E+08],
        'wavelength':[1031.9120, 1037.6130]
    },
    'NeV': {
        'f': [1.100E-01, 8.100E-02, 9.600E-02],
        'Gamma':[7.700E+08, 1.000E+09, 1.400E+09],
        'wavelength':[568.4240, 569.8280, 572.3350]
    },
    'NeVI': {
        'f': [1.400E-01, 1.200E-01],
        'Gamma':[1.500E+09, 1.700E+09],
        'wavelength':[558.6030, 562.8050]
    },
    'NeVII': {
        'f': [5.600E-02],
        'Gamma':[5.800E+08],
        'wavelength':[465.2210]
    },
    'NeVIII': {
        'f': [1.020E-01, 5.020E-02],
        'Gamma':[5.720E+08, 5.500E+08],
        'wavelength':[770.4090, 780.3240]
    },
     'NaIX': {
        'f': [9.190E-02, 4.500E-02],
        'Gamma':[6.600E+08, 6.230E+08],
        'wavelength':[681.7200, 694.1500]
    },
    'MgII': {
        'f': [6.210E-04, 3.510E-04],
        'Gamma':[1.350E+06, 1.520E+06],
        'wavelength':[1239.9253, 1240.3947]
    },
    'MgX': {
        'f': [8.380E-02, 4.070E-02],
        'Gamma':[7.510E+08, 6.950E+08],
        'wavelength':[609.7930, 624.9410]
    },
    'AlII': {
        'f': [1.830E+00],
        'Gamma':[1.460E+09],
        'wavelength':[1670.7874]
    },
    'AlIII': {
        'f': [5.570E-01, 2.770E-01],
        'Gamma':[5.400E+08, 5.330E+08],
        'wavelength':[1854.7164, 1862.7895]
    },
    'SiII': {
        'f': [
            2.000E-01, 1.570E-01, 1.390E-02, 1.390E-02, 2.770E-01, 5.750E-01, 7.370E-01,
            1.500E-01, 1.180E+00, 1.940E-01, 1.090E+00, 1.130E-01, 2.490E-03, 1.970E-03
        ],
        'Gamma':[
            6.810E+08, 7.110E+08, 8.910E+07, 1.770E+08, 6.530E+08, 2.690E+09, 3.450E+09,
            1.400E+09, 2.950E+09, 8.300E+08, 3.040E+09, 4.730E+08, 2.540E+06, 2.650E+06
        ],
        'wavelength':[
            989.8730, 992.6830, 1020.6990, 1023.7000, 1190.4160, 1193.2900, 1194.5000,
            1197.3940, 1260.4220, 1264.7380, 1304.3700, 1309.2760, 1808.0130, 1816.9280]
    },
    'SiIII': {
        'f': [1.630E+00],
        'Gamma':[2.480E+09],
        'wavelength':[1206.5000]
    },
    'SiIV': {
        'f': [2.550E-01],
        'Gamma':[8.630E+08],
        'wavelength':[1402.7700]
    },
    'SiXII': {
        'f': [7.180E-02, 3.420E-02],
        'Gamma':[9.600E+08, 8.410E+08],
        'wavelength':[499.4060, 520.6650]
    },
    'PIV': {
        'f': [1.600E+00],
        'Gamma':[3.940E+09],
        'wavelength':[950.6570]
    },
    'PV': {
        'f': [4.500E-01, 2.210E-01],
        'Gamma':[1.200E+09, 1.160E+09],
        'wavelength':[1117.9770, 1128.0080]
    },
    'SII': {
        'f': [
            4.190E-01, 8.360E-01, 1.250E+00, 2.010E-01, 1.320E-01, 6.530E-02, 6.020E-03,
            1.210E-02, 1.820E-02
        ],
        'Gamma':[
            9.600E+09, 9.550E+09, 9.470E+09, 1.090E+09, 1.060E+09, 1.050E+09, 5.130E+07,
            5.120E+07, 5.100E+07
        ],
        'wavelength':[
            763.6560, 764.4160, 765.6840, 906.8760, 910.4850, 912.7360, 1250.5840,
            1253.8110, 1259.5190
        ]
    },
    'SIII': {
        'f': [
            1.430E+00, 9.170E-01, 1.380E+00, 2.780E-01, 8.500E-01, 4.290E-01, 2.680E-01,
            5.810E-01, 3.780E-01, 3.520E-01, 3.090E-01, 4.410E-02, 1.490E-02, 1.210E-02,
            1.680E-02, 1.020E-02, 3.270E-02, 2.580E-02, 1.920E-02, 2.000E-02
        ],
        'Gamma':[
            6.930E+09, 7.970E+09, 1.420E+10, 3.990E+09, 3.870E+09, 3.500E+09, 1.090E+10,
            7.850E+09, 1.600E+09, 4.450E+09, 6.480E+09, 9.560E+07, 2.880E+08, 7.800E+07,
            6.530E+07, 1.090E+08, 2.090E+08, 4.050E+07, 5.380E+07, 6.620E+07
        ],
        'wavelength':[
            677.7290, 678.4560, 680.6770, 680.9250, 698.7270, 700.1500, 700.2880,
            702.7790, 724.2880, 725.8580, 728.6850, 1012.4950, 1015.5020, 1015.5670,
            1015.7790, 1021.1080, 1021.3230, 1190.2030, 1194.0580, 1200.9660
        ]
    },
    'SIV': {
        'f': [1.130E+00, 2.490E-01, 4.590E-01, 5.970E-01, 1.310E-01, 1.180E-01, 8.500E-02, 4.200E-02],
        'Gamma':[8.690E+09, 1.500E+09, 5.470E+09, 7.080E+09, 3.080E+09, 1.200E+09, 1.700E+09, 1.620E+08],
        'wavelength':[657.3190, 744.9040, 748.3930, 750.2210, 753.7600, 809.6560, 815.9410, 1072.9730]
    },
    'SV': {
        'f': [1.360E+00],
        'Gamma':[4.870E+09],
        'wavelength':[786.4680]
    },
    'SVI': {
        'f': [4.360E-01, 2.150E-01],
        'Gamma':[1.670E+09, 1.610E+09],
        'wavelength':[933.3780, 944.5230]
    },
    'SXIV': {
        'f': [6.350E-02, 2.960E-02],
        'Gamma':[1.210E+09, 9.950E+08],
        'wavelength':[417.6600, 445.7000]
    },
    'ArI': {
        'f': [2.500E-01, 6.090E-02],
        'Gamma':[5.100E+08, 1.190E+08],
        'wavelength':[1048.2200, 1066.6600]
    },
    'ArII': {
        'f': [3.700E-02, 1.800E-01, 1.500E-01],
        'Gamma':[9.500E+08, 2.300E+09, 1.900E+09],
        'wavelength':[718.0899, 723.3606, 725.5485]
    },
    'ArVII': {
        'f': [1.210E+00],
        'Gamma':[7.830E+09],
        'wavelength':[585.7480]
    },
    'CaX': {
        'f': [3.260E-01, 1.600E-01],
        'Gamma':[3.500E+09, 3.200E+09],
        'wavelength':[557.7650, 574.0100]
    },
    'FeII': {
        'f': [
            2.900E-03, 3.300E-02, 3.300E-03, 7.000E-03, 6.100E-03, 2.900E-03, 5.500E-02,
            4.700E-03, 1.300E-02, 2.800E-03, 3.300E-02, 4.500E-03, 2.900E-02, 1.600E-02,
            1.100E-03, 1.900E-02, 5.910E-02
        ],
        'Gamma':[
            2.800E+07, 2.600E+08, 3.200E+07, 4.400E+07, 4.600E+07, 1.700E+07, 3.200E+08,
            3.500E+07, 6.000E+07, 1.600E+07, 2.300E+08, 2.000E+07, 1.900E+08, 1.000E+08,
            6.000E+06, 1.000E+08, 1.910E+08
        ],
        'wavelength':[
            923.8782, 926.2120, 926.8969, 937.6513, 1055.2618, 1062.1531, 1063.1769,
            1063.9718, 1081.8753, 1083.4205, 1096.8769, 1112.0483, 1121.9747, 1125.4476,
            1127.0984, 1143.2257, 1608.4508
        ]
    },
}

MIN_NHI = 1e5
MAX_NHI = 1e24
START_NHI = 1e12
# check values start_NHI, max_NHI and min_NHI
# maxminN or maxmin_N?

for species_dict in SPECIES_DICTS.values():
    species_dict["max_N"] = MAX_NHI
    species_dict["min_N"] = MIN_NHI
    species_dict["init_N"] = START_NHI

ORDER_FITS = [
    "Lya", "Lyb", "CI", "CII", "CIII", "CIV", "NI", "NII", "NIII", "NIV", "NV", "OI", "OII", "OIII", "OIV",
    "OV", "OVI", "NeV", "NeVI", "NeVII", "NeVIII", "NaIX", "MgII", "MgX", "AlII", "AlIII", "SiII", "SiIII",
    "SiIV", "SiXII", "PIV", "PV", "SII", "SIII", "SIV", "SV", "SVI", "SXIV", "ArI", "ArII", "ArIII", "CaX", 
    "FeII", "Lyg", "Lyother"
]

Z_MIN_FACTOR = 0.9
Z_MAX_FACTOR = 1.1
B_MIN_FACTOR = 0.9
B_MAX_FACTOR = 1.1

class LineProfileLeastSquares:
    """This class deals with the line profile fitting.

    It is passed to iminuit and when called it will return the chi2 for a given
    set of parameters

    Methods
    -------
    __init__
    __call__
    compute_line_profile

    Attributes
    ----------
    flux: array
    The flux array

    flux_fit: array
    Current flux fit including the profile of previously fitted lines

    species_dict: dict
    Dictionary containing all relevant parameters needed to create an
    absorption line of a given species (f, Gamma, lambda0)

    wavelength: array
    Wavelength array
    """
    def __init__(self, wavelength, flux, flux_fit, species_dict):
        """Initialize class instance

        Arguments
        ---------
        wavelength: array
        Wavelength array

        flux: array
        The flux array

        flux_fit: array
        Current flux fit including information from previous lines.
        New lines will be added to this

        species_dict: dict
        Dictionary containing all relevant parameters needed to create an
        absorption line of a given species (f, Gamma, lambda0)
        """
        self.flux = flux
        self.flux_fit = flux_fit
        self.species_dict = species_dict
        self.wavelength = wavelength

    def __call__(self, column_density, impact_parameter, redshift):
        """Compute the chi2 of the fit

        Arguments
        ---------
        column_density: float
        The gas column density

        impact_parameter: float
        The absoprtion impact parameter

        redsfhit: float
        The absorption redshift

        Return
        ------
        chi2: float
        The chi2 of this run
        """
        # TODO: add ivar dependence
        chi2_contribution = self.compute_line_profile(
            column_density, impact_parameter, redshift) - self.flux
        return np.sum(chi2_contribution**2)

    def compute_line_profile(self, column_density, impact_parameter, redshift):
        """Calculate the normalized flux for a region of wavelength
        space generated by a set of absorption lines.

        Takes into account the profiles of the lines previously fitted

        Arguments
        ---------
        column_density: float
        The gas column density

        impact_parameter: float
        The absoprtion impact parameter

        redsfhit: float
        The absorption redshift

        Return
        ------
        flux: array
        The predicted flux array
        """
        tau = np.zeros_like(self.wavelength)

        for oscillator_strength, gamma, line_wavelength in zip(
            self.species_dict['f'],
            self.species_dict['Gamma'],
            self.species_dict['wavelength']):
                tau += gen_tau(
                    self.wavelength,
                    column_density,
                    impact_parameter,
                    redshift,
                    oscillator_strength,
                    gamma,
                    line_wavelength)

        flux = self.flux_fit * np.exp(-tau)
        #flux = np.exp(-tau)
        return flux

def gen_tau(wavelength, column_density, impact_parameter, redshift,
            oscillator_strength, gamma, line_wavelength):
    """This calculates a flux distribution for given parameters using the yt
    voigt profile generator"""
    tau_o = (1.4973614E-15 * column_density * oscillator_strength *
             line_wavelength / impact_parameter)
    a = 7.95774715459E-15 * gamma * line_wavelength / impact_parameter
    x = SPEED_LIGHT / impact_parameter * (
        line_wavelength * (1 + redshift) / wavelength - 1)

    H = voigt(a, x)

    tau = tau_o * H

    return tau


def fit_lines(base_name, input_extension, noise, start_z, start_b,
              species_dicts=SPECIES_DICTS, order_fits=ORDER_FITS):
    """
    Fit an absorption-line spectrum into line profiles.

    Fits the spectrum into absorption complexes and iteratively adds and
    optimizes voigt profiles for each complex.

    Arguments
    ---------
    base_name: str
    Base name used to name the outputs

    input_extension: str
    Extension of the file containing the spectra. Valid extensions are
    ".h5" or ".txt"

    noise: float
    Noise level in the spectra (-1 for no noise)

    start_z: float
    Initial guess for the redshift. Redshift range in the fit set from
    Z_MIN_FACTOR * start_z to Z_MAX_FACTOR * start_z

    start_b: float
    Initial guess for the impact parameter. Impact parameter range in the
    fit set from B_MIN_FACTOR * start_b to B_MAX_FACTOR * start_b

    species_dicts: dict of dicts - default: SPECIES_DICTS
    Top level keys should be the names of all the species given
    in order_fits. The entries should be dictionaries containing all
    relevant parameters needed to create an absorption line of a given
    species (f, Gamma, lambda0) as well as max and min values for
    the impact parameter b, redshift z and column density N.

    order_fits: list of str - default: ORDER_FITS
    List of the names of the species in the order that they
    should be fit. Names should correspond to the names of the species
    given in species_dict. (ex: ["lya", "OVI"])

    Return
    ------
    all_species_lines: array
    Array with the parameters of the fitted lines.
    NaNs are placed when encoutering fitting problems

    fit_flux: array
    Fitted flux
    """
    # initialize_output
    dtype = []
    for species in order_fits:
        dtype += [
            (f"{species} N [cm^-2]", float),
            (f"{species} b [km/s]", float),
            (f"{species} zfit", float)
        ]
    all_species_lines = np.array([np.nan], dtype=dtype)

    # load arrays
    if input_extension == ".h5":
        if noise < 0.:
            file = h5py.File(f"{base_name}_spec_nonoise.h5")
        else:
            file = h5py.File(f"{base_name}_spec.h5")
        flux = file['flux'][:]
        wavelength = file["wavelength"][:]
        file.close()
    elif input_extension in [".fits", ".fits.gz"]:
        if noise < 0.:
            hdu = fits.open(f"{base_name}_spec_nonoise{input_extension}")
        else:
            hdu = fits.open(f"{base_name}_spec{input_extension}")
        flux = hdu[1].data["flux"]
        wavelength = hdu[1].data["wavelength"]
        # TODO: read ivar
        hdu.close()
    elif input_extension == ".txt":
        # TODO: read data from txt file
        print(".txt extension not implemented")
        return all_species_lines, np.array([])
    else:
        print("Unknown extension, expected one of .fits, .fits.gz, .txt or .h5")
        return all_species_lines, np.array([])

    # add impact paramenter and redshift guesses
    species_dicts_fit = copy.copy(species_dicts)
    min_z = start_z * Z_MIN_FACTOR
    max_z = start_z * Z_MAX_FACTOR
    min_b = start_b * B_MIN_FACTOR
    max_b = start_b * B_MAX_FACTOR
    for key, species in species_dicts_fit.items():
        species["max_z"] = max_z
        species["min_z"] = min_z
        species["init_z"] = start_z
        species["max_b"] = max_b
        species["min_b"] = min_b
        species["init_b"] = start_b
    # maxminb, maxminN, maxminz or maxmin_b, maxmin_N, maxminz?
    # why use init_z?

    # initialze fit results
    flux_fit = np.ones_like(flux)
    current_chi2 = np.sum((flux_fit - flux)**2)

    # fit all species one at a time in the specified order
    for species in order_fits:

        species_dict = species_dicts_fit.get(species)

        # initialize the fitter class
        leasts_squares = LineProfileLeastSquares(
            wavelength, flux, flux_fit, species_dict)

        # fit line
        minimizer = iminuit.Minuit(
            leasts_squares,
            name=("column_density", "impact_parameter", "redshift"),
            column_density=species_dict["init_N"],
            impact_parameter=species_dict["init_b"],
            redshift=species_dict["init_z"]
        )
        minimizer.limits["column_density"] = (
            species_dict["min_N"], species_dict["max_N"])
        minimizer.limits["impact_parameter"] = (
            species_dict["min_b"], species_dict["max_b"])
        minimizer.limits["redshift"] = (
            species_dict["min_z"], species_dict["max_z"])
        minimizer.errordef = 1.
        minimizer.print_level = 0
        minimizer.migrad()

        # update total fit
        if minimizer.valid and current_chi2 > minimizer.fval:
            minimizer.hesse()

            all_species_lines[f"{species} N [cm^-2]"] = minimizer.values["column_density"]
            all_species_lines[f"{species} N [cm^-2]"] = minimizer.values["impact_parameter"]
            all_species_lines[f"{species} N [cm^-2]"] = minimizer.values["redshift"]

            flux_fit = leasts_squares.compute_line_profile(
                minimizer.values["column_density"],
                minimizer.values["impact_parameter"],
                minimizer.values["redshift"],
            )

            current_chi2 = minimizer.fval

        # fit is now worse, report and skip line
        elif minimizer.valid:
            print(
                f"For spectrum {base_name}, fit on species {species} worsens "
                f"the chi2. Previous chi2 = {current_chi2}."
                f"New chi2 = {minimizer.fval}. Ignoring line")
        # fit did not converge, report and skip line
        else:
            print(
                f"For spectrum {base_name}, fit on species {species} failed. "
                "Ignoring")

    return all_species_lines, flux_fit
