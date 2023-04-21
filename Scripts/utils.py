import yt
import trident
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
# Ignasi: not sure what we are importing here. Not sure if it's used either
# I'm commenting it for the time being
#from yt.mods import *
import math
import h5py
from astropy.io import fits
import numpy as np


# some configuration variables used by run_simple_ray
Z_Solar = 0.02041

speciesDicts = {
    'HI': {
        'name':'HI',
        'f': [.4160],
        'Gamma':[4.690E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':30,
        'init_N':1E14
    },
    'HILy': {
        'name':'HILy',
        'f': [.0791,.0290,.0139,.00779,.00481,.00318,.00221,.0016,.0012,.000921,.000722,.000577,.000468,.000385,.000321,.000270,.000229,.000197,.000170,.000148,.000129,.000114,.0001,.0000892,.0000797,7.140E-5,6.430E-5,5.80E-5,5.26E-5,4.78E-5,4.36E-5,3.98E-5,3.65E-5,3.350000E-05,3.090000E-05,2.850000E-05,2.640000E-05,2.440000E-05],
        'Gamma':[5.570E7,1.280E6,4.120E6,1.640E6,7.560E5,3.870E5,2.140E5,1.260E5,7.830E4,5.060E4,3.390E4,2.340E4,1.660E4,1.200E4,8.850E3,6.650E3,5.070E3,3.920E3,3.070E3,2.440E3,1.950E3,1.580E3,1.280E3,1.060E3,8.740E2,7.290E2,6.120E2,5.160E2,4.380E2,3.740E2,3.20E2,2.760E2,2.390E2,2.070000E+02,1.810000E+02,1.580000E+02,1.390000E+02,1.220000E+02],
        'wavelength':[1025.7222,972.5367,949.7429,937.8034,930.7482,926.2256,923.1503,920.9630,919.3513,918.1293,917.1805,916.4290,915.8237,915.3289,914.9192,914.5762,914.2860,914.0385,913.8256,913.6411,913.4802,913.3391,913.2145,913.1041,913.0058,912.9178,912.8388,912.7676,912.7032,912.6447,912.5914,912.5428,912.4983,912.4575,912.4199,912.3853,912.3532,912.3236],
        'numLines':38,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'CI': {
        'name':'CI',
        'f': [1.430000E-01,7.160000E-02,5.800000E-02,2.610000E-02,9.230000E-02,3.940000E-02,1.250000E-02,4.090000E-02,1.240000E-02,2.120000E-02],
        'Gamma':[1.160000E+08,6.540000E+07,7.310000E+07,3.550000E+07,1.260000E+08,5.510000E+07,1.940000E+07,6.390000E+07,1.950000E+07,3.520000E+07],
        'wavelength':[1656.9290,1560.3090,1328.8340,1280.1350,1277.2450,1260.7350,1193.9950,1193.0310,1188.8330,1157.9100],
        'numLines':10,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'CII': {
        'name':'CII',
        'f': [1.270000E-02,1.290000E-01,1.180000E-01,1.190000E-01,3.310000E-01],
        'Gamma':[4.760000E+07,2.410000E+08,1.460000E+09,7.380000E+08,2.700000E+09],
        'wavelength':[1335.6630,1334.5320,1037.0180,1036.3370,903.9620],
        'numLines':5,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'CIII': {
        'name':'CIII',
        'f': [7.586000E-01],
        'Gamma':[1.767000E+09],
        'wavelength':[977.0200],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'CIV': {
        'name':'CIV',
        'f': [9.520000E-02,1.900000E-01],
        'Gamma':[2.640000E+08,2.650000E+08],
        'wavelength':[1550.7720,1548.1870],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NI': {
        'name':'NI',
        'f': [8.690000E-02,1.320000E-01,3.870000E-03,2.870000E-02,1.460000E-02,9.690000E-04,3.860000E-03,7.900000E-03,1.240000E-02,4.000000E-03,3.310000E-02,2.470000E-02,1.290000E-02,1.970000E-03,2.290000E-03],
        'Gamma':[4.030000E+08,4.070000E+08,1.820000E+07,1.490000E+08,1.510000E+08,3.520000E+06,5.520000E+07,5.660000E+07,5.940000E+07,1.950000E+07,1.620000E+08,1.810000E+08,1.900000E+08,1.450000E+07,1.120000E+07],
        'wavelength':[1200.7100,1200.2230,1199.5500,1134.9800,1134.4150,1134.1650,965.0410,964.6260,963.9900,954.1040,953.970000,953.655000,953.415000,952.415000,952.303000],
        'numLines':15,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NII': {
        'name':'NII',
        'f': [1.110000E-01,1.590000E-01],
        'Gamma':[2.100000E+08,4.230000E+08],
        'wavelength':[1083.990000,915.612000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NIII': {
        'name':'NIII',
        'f': [1.230000E-01],
        'Gamma':[4.180000E+08],
        'wavelength':[989.799000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NIV': {
        'name':'NIV',
        'f': [6.109000E-01],
        'Gamma':[2.320000E+09],
        'wavelength':[765.147000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NV': {
        'name':'NV',
        'f': [7.800000E-02,1.560000E-01],
        'Gamma':[3.370000E+08,3.400000E+08],
        'wavelength':[1242.804000,1238.821000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'OI': {
        'name':'OI',
        'f': [5.190000E-02,5.180000E-02,5.200000E-02,9.160000E-03,4.640000E-02,8.460000E-03,3.310000E-03,1.580000E-03,3.060000E-03,1.920000E-03,1.300000E-03,2.190000E-02],
        'Gamma':[6.760000E+07,2.030000E+08,3.410000E+08,9.430000E+07,2.260000E+08,5.770000E+07,3.860000E+07,1.940000E+07,1.660000E+07,1.060000E+07,7.220000E+06,1.230000E+08],
        'wavelength':[1306.029000,1304.858000,1302.168000,1039.230000,988.773000,988.655000,976.448000,950.885000,936.629000,929.517000,924.950000,922.008000],
        'numLines':12,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'OII': {
        'name':'OII',
        'f': [1.350000E-01,9.010000E-02,4.510000E-02],
        'Gamma':[8.610000E+08,8.650000E+08,8.670000E+08],
        'wavelength':[834.465400,833.330300,832.758300],
        'numLines':3,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'OIII': {
        'name':'OIII',
        'f': [8.770000E-02,8.770000E-02,1.060000E-01,4.520000E-02,1.340000E-01],
        'Gamma':[5.990000E+08,5.990000E+08,3.410000E+08,1.830000E+09,6.060000E+08],
        'wavelength':[835.289000,835.289000,832.929000,702.838000,702.337000],
        'numLines':5,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'OIV': {
        'name':'OIV',
        'f': [9.940000E-02,1.110000E-01,6.690000E-02,6.710000E-02,2.790000E-01,2.240000E-01,1.120000E-01],
        'Gamma':[7.080000E+08,5.950000E+08,2.400000E+09,1.210000E+09,6.060000E+09,4.860000E+09,1.220000E+09],
        'wavelength':[790.199000,787.710000,609.829000,608.397000,554.513000,554.076000,553.329000],
        'numLines':7,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'OV': {
        'name':'OV',
        'f': [5.122000E-01],
        'Gamma':[2.872000E+09],
        'wavelength':[629.732000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'OVI': {
        'name':'OVI',
        'f': [6.600000E-02,1.330000E-01],
        'Gamma':[4.090000E+08,4.160000E+08],
        'wavelength':[1037.613000,1031.912000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NeV': {
        'name':'NeV',
        'f': [9.600000E-02,8.100000E-02,1.100000E-01],
        'Gamma':[1.400000E+09,1.000000E+09,7.700000E+08],
        'wavelength':[572.335000,569.828000,568.424000],
        'numLines':3,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NeVI': {
        'name':'NeVI',
        'f': [1.200000E-01,1.400000E-01],
        'Gamma':[1.700000E+09,1.500000E+09],
        'wavelength':[562.805000,558.603000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NeVII': {
        'name':'NeVII',
        'f': [5.600000E-02],
        'Gamma':[5.800000E+08],
        'wavelength':[465.221000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NeVIII': {
        'name':'NeVIII',
        'f': [5.020000E-02,1.020000E-01],
        'Gamma':[5.500000E+08,5.720000E+08],
        'wavelength':[780.324000,770.409000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'NaIX': {
        'name':'NaIX',
        'f': [4.500000E-02,9.190000E-02],
        'Gamma':[6.230000E+08,6.600000E+08],
        'wavelength':[694.150000,681.720000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'MgII': {
        'name':'MgII',
        'f': [3.510000E-04,6.210000E-04],
        'Gamma':[1.520000E+06,1.350000E+06],
        'wavelength':[1240.394700,1239.925300],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'MgX': {
        'name':'MgX',
        'f': [4.070000E-02,8.380000E-02],
        'Gamma':[6.950000E+08,7.510000E+08],
        'wavelength':[624.941000,609.793000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'AlII': {
        'name':'AlII',
        'f': [1.830000E+00],
        'Gamma':[1.460000E+09],
        'wavelength':[1670.787400],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'AlIII': {
        'name':'AlIII',
        'f': [2.770000E-01,5.570000E-01],
        'Gamma':[5.330000E+08,5.400000E+08],
        'wavelength':[1862.789500,1854.716400],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SiII': {
        'name':'SiII',
        'f': [1.970000E-03,2.490000E-03,1.130000E-01,1.090000E+00,1.940000E-01,1.180000E+00,1.500000E-01,7.370000E-01,5.750000E-01,2.770000E-01,1.390000E-02,1.390000E-02,1.570000E-01,2.000000E-01],
        'Gamma':[2.650000E+06,2.540000E+06,4.730000E+08,3.040000E+09,8.300000E+08,2.950000E+09,1.400000E+09,3.450000E+09,2.690000E+09,6.530000E+08,1.770000E+08,8.910000E+07,7.110000E+08,6.810000E+08],
        'wavelength':[1816.928000,1808.013000,1309.276000,1304.370000,1264.738000,1260.422000,1197.394000,1194.500000,1193.290000,1190.416000,1023.700000,1020.699000,992.683000,989.873000],
        'numLines':14,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SiIII': {
        'name':'SiIII',
        'f': [1.630000E+00],
        'Gamma':[2.480000E+09],
        'wavelength':[1206.500000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SiIV': {
        'name':'SiIV',
        'f': [2.550000E-01],
        'Gamma':[8.630000E+08],
        'wavelength':[1402.770000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SiXII': {
        'name':'SiXII',
        'f': [3.420000E-02,7.180000E-02],
        'Gamma':[8.410000E+08,9.600000E+08],
        'wavelength':[520.665000,499.406000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'PIV': {
        'name':'PIV',
        'f': [1.600000E+00],
        'Gamma':[3.940000E+09],
        'wavelength':[950.657000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'PV': {
        'name':'PV',
        'f': [2.210000E-01,4.500000E-01],
        'Gamma':[1.160000E+09,1.200000E+09],
        'wavelength':[1128.008000,1117.977000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SII': {
        'name':'SII',
        'f': [1.820000E-02,1.210000E-02,6.020000E-03,6.530000E-02,1.320000E-01,2.010000E-01,1.250000E+00,8.360000E-01,4.190000E-01],
        'Gamma':[5.100000E+07,5.120000E+07,5.130000E+07,1.050000E+09,1.060000E+09,1.090000E+09,9.470000E+09,9.550000E+09,9.600000E+09],
        'wavelength':[1259.519000,1253.811000,1250.584000,912.736000,910.485000,906.876000,765.684000,764.416000,763.656000],
        'numLines':9,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SIII': {
        'name':'SIII',
        'f': [2.000000E-02,1.920000E-02,2.580000E-02,3.270000E-02,1.020000E-02,1.680000E-02,1.210000E-02,1.490000E-02,4.410000E-02,3.090000E-01,3.520000E-01,3.780000E-01,5.810000E-01,2.680000E-01,4.290000E-01,8.500000E-01,2.780000E-01,1.380000E+00,9.170000E-01,1.430000E+00],
        'Gamma':[6.620000E+07,5.380000E+07,4.050000E+07,2.090000E+08,1.090000E+08,6.530000E+07,7.800000E+07,2.880000E+08,9.560000E+07,6.480000E+09,4.450000E+09,1.600000E+09,7.850000E+09,1.090000E+10,3.500000E+09,3.870000E+09,3.990000E+09,1.420000E+10,7.970000E+09,6.930000E+09],
        'wavelength':[1200.966000,1194.058000,1190.203000,1021.323000,1021.108000,1015.779000,1015.567000,1015.502000,1012.495000,728.685000,725.858000,724.288000,702.779000,700.288000,700.150000,698.727000,680.925000,680.677000,678.456000,677.729000],
        'numLines':20,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SIV': {
        'name':'SIV',
        'f': [4.200000E-02,8.500000E-02,1.180000E-01,1.310000E-01,5.970000E-01,4.590000E-01,2.490000E-01,1.130000E+00],
        'Gamma':[1.620000E+08,1.700000E+09,1.200000E+09,3.080000E+09,7.080000E+09,5.470000E+09,1.500000E+09,8.690000E+09],
        'wavelength':[1072.973000,815.941000,809.656000,753.760000,750.221000,748.393000,744.904000,657.319000],
        'numLines':8,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SV': {
        'name':'SV',
        'f': [1.360000E+00],
        'Gamma':[4.870000E+09],
        'wavelength':[786.468000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SVI': {
        'name':'SVI',
        'f': [2.150000E-01,4.360000E-01],
        'Gamma':[1.610000E+09,1.670000E+09],
        'wavelength':[944.523000,933.378000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'SXIV': {
        'name':'SXIV',
        'f': [2.960000E-02,6.350000E-02],
        'Gamma':[9.950000E+08,1.210000E+09],
        'wavelength':[445.700000,417.660000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'ArI': {
        'name':'ArI',
        'f': [6.090000E-02,2.500000E-01],
        'Gamma':[1.190000E+08,5.100000E+08],
        'wavelength':[1066.660000,1048.220000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'ArII': {
        'name':'ArII',
        'f': [1.500000E-01,1.800000E-01,3.700000E-02],
        'Gamma':[1.900000E+09,2.300000E+09,9.500000E+08],
        'wavelength':[725.548530,723.360560,718.089852],
        'numLines':3,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'ArVII': {
        'name':'ArVII',
        'f': [1.210000E+00],
        'Gamma':[7.830000E+09],
        'wavelength':[585.748000],
        'numLines':1,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'CaX': {
        'name':'CaX',
        'f': [1.600000E-01,3.260000E-01],
        'Gamma':[3.200000E+09,3.500000E+09],
        'wavelength':[574.010000,557.765000],
        'numLines':2,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },
    'FeII': {
        'name':'FeII',
        'f': [5.910000E-02,1.900000E-02,1.100000E-03,1.600000E-02,2.900000E-02,4.500000E-03,3.300000E-02,2.800000E-03,1.300000E-02,4.700000E-03,5.500000E-02,2.900000E-03,6.100000E-03,7.000000E-03,3.300000E-03,3.300000E-02,2.900000E-03],
        'Gamma':[1.910000E+08,1.000000E+08,6.000000E+06,1.000000E+08,1.900000E+08,2.000000E+07,2.300000E+08,1.600000E+07,6.000000E+07,3.500000E+07,3.200000E+08,1.700000E+07,4.600000E+07,4.400000E+07,3.200000E+07,2.600000E+08,2.800000E+07],
        'wavelength':[1608.450830,1143.225730,1127.098400,1125.447630,1121.974710,1112.048320,1096.876870,1083.420460,1081.875260,1063.971790,1063.176870,1062.153100,1055.261760,937.651330,926.896870,926.212010,923.878210],
        'numLines':17,
        'maxN': 1E17,
        'minN':1E11,
        'maxb': 300,
        'minb':1,
        'maxz': 6,
        'minz':0,
        'init_b':20,
        'init_N':1E12
    },

}
orderFits = ['FeII','CaX','ArVII','ArII','ArI','SXIV','SVI','SV','SIV','SIII','SII','PV','PIV','SiXII','SiIV','SiIII','SiII','AlIII','AlII','MgX','MgII','NaIX','NeVIII','NeVII','NeVI','NeV','OVI','OV','OIV','OIII','OII','OI','NV','NIV','NIII','NII','NI','CIV','CIII','CII','CI','HILy','HI']



def fit_lines(base_name, input_extension, save_fit, noise):
    """Fit line profiles in the spectrum

    Arguments
    ---------
    base_name: str
    Base name used to name the outputs

    input_extension: str
    Extension of the file containing the spectra. Valid extensions are
    ".h5" or ".txt"

    save_fit: bool
    Save the fit results in an h5 file called "{base_name}HIspecFIT.h5"

    noise: float
    Noise level in the spectra (-1 for no noise)

    Return
    ------
    nhi: float
    The column density of the line

    b_kms: float
    The impact parameter of the line, in km/s

    zfit: float
    The fitted redshift
    """
    # load arrays
    if input_extension == ".h5":
        if noise < 0.:
            file = h5py.File(f"{base_name}spec_nonoise.h5")
        else:
            file = h5py.File(f"{base_name}spec.h5")
        flux = file['flux'][:]
        wavelength = file["wavelength"][:]
        file.close()
    elif input_extension in [".fits", ".fits.gz"]:
        if noise < 0.:
            hdu = fits.open(f"{base_name}spec_nonoise{input_extension}")
        else:
            hdu = fits.open(f"{base_name}spec{input_extension}")
        flux = hdu[1].data["flux"]
        wavelength = hdu[1].data["wavelength"]
        hdu.close()
    elif input_extension == ".txt":
        # TODO: read data from txt file
        return 0, 0, 0
    else:
        print("Loading from other than .h5 not implemented")
        return 0, 0, 0

    # do actual fit
    if save_fit:
        output_file = f"{base_name}specFIT.h5"
    else:
        output_file = None

    fitted_lines, fitted_flux = generate_total_fit(
        wavelength,
        flux,
        orderFits,
        speciesDicts,
        maxLength=9000,
        output_file=output_file)

    # parse results
    if fitted_lines['HI']['N'].size > 0:
        nhi = fitted_lines['HI']['N'][0]
        b_kms = fitted_lines['HI']['b'][0]
        zfit = fitted_lines['HI']['z'][0]
    else:
        nhi = 0.
        b_kms = 0.
        zfit = 0.

    return nhi, b_kms, zfit


def initialize_catalogue(catalogue_name, quasar=False, galaxy=False):
    """Initialize the output catalogue.
    User is responsible for closing the file

    Arguments
    ---------
    catalogue_name: str
    The name of the catalogue

    Return
    ------
    catalogue_file: _io.TextIOWrapper
    The open catalogue. User is responsible for closing the file.

    Raise
    -----
    Exception if
    """
    catalogue_file = open(catalogue_name, "w")
    if quasar and galaxy:
        catalogue_file.write(write("; ".join([
            "Name QSO", "Name Galaxy", "z QSO", "z Galaxy", "d",
            "start_shift_x", "start_shift_y", "start_shift_z",
            "end_shift_x", "end_shift_y", "end_shift_z",
            "N [cm^-2]", "Absortion System"])))
    elif quasar:
        catalogue_file.write("; ".join(["Name", "Dataset", "z", "seed"]))
    elif galaxy:
        catalogue_file.write("; ".join([
            "Name", "Dataset", "z", "d",
            "start_shift_x", "start_shift_y", "start_shift_z",
            "end_shift_x", "end_shift_y", "end_shift_z",
            "N [cm^-2]", "b [km/s]", "zfit"]))
    else:
        catalogue_file.close()
        raise IOError("Unkown catalogue to initialize")

    catalogue_file.write("\n")

    return catalogue_file


def load_snapshot(fn, dir=""):
    """Load the simulation snapshot

    Arguments
    ---------
    fn: str
    Name of the snapshot

    dir: str
    Directory where snapshots are kept

    Return
    ------
    ds: yt.data_objects.static_output.Dataset
    The loaded snapshot
    """
    ds = yt.load(f"{dir}{fn}/{fn}")
    ds.add_field(
        ("gas", "metallicity"),
        function=metallicity_e,
        sampling_type="local",
        force_override=True,
        display_name="Z/Z$_{\odot}$",
        take_log=True,
        units="")
    return ds

def metallicity_e(field, data):
    """Compute metallicity"""
    return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

def run_galaxy_snapshot(fn,
                        galaxy_pos,
                        z_min,
                        z_max,
                        z_step,
                        dist_min,
                        dist_max,
                        dist_step,
                        base_name,
                        catalogue_file,
                        starting_n=1):
    """Run a set of simple rays for the specified distances to the center of the galaxy

    Arguments
    ---------
    fn: str
    Name of the snapshot

    galaxy_pos: array
    3D position of the galaxy in the snapshot

    z_min: float
    Minimu redshift

    z_max: float
    Maximum redshift

    z_step: float
    Redshift step

    dist_min: float
    Minimum distance to the center of the galaxy

    dist_max: float
    Maximum distance to the center of the galaxy

    dist_step: float
    Step in the distance coverage (from dist_min to dist_max).

    base_name: str
    Base name used to name the outputs

    catalogue_name: str
    The name of the catalogue
    """
    snapshot_name = fn.split("/")[0]
    ds = load_snapshot(fn)

    n = starting_n
    for z in np.arange(z_min, z_max, z_step):
        #rayo centro
        base_name_alt = base_name.replace(
            "/G_", "/G" + snapshot_name.replace("RD0", "")) + f"_{z:.1f}"
        if snapshot_name == "RD0180":
            # for S II 766 line (765.684000 Angs) as otherwise we get
            # MemoryError: Unable to allocate 6.10 GiB for an array
            start_shift=[-10, -10, -10]
            end_shift=[11, 10, 10]
            catalogue_file.write(run_simple_ray(
                ds,
                z,
                0.0,
                start_shift,
                end_shift,
                snapshot_name,
                galaxy_pos,
                base_name_alt))
        else:
            start_shift=[-10, -10, -10]
            end_shift=[10, 10, 10]
            catalogue_file.write(run_simple_ray(
                ds,
                z,
                0.0,
                start_shift,
                end_shift,
                snapshot_name,
                galaxy_pos,
                base_name_alt))

        # loop over distances
        for dist in np.arange(dist_min, dist_max, dist_step):
            l = math.radians(0.0)
            j = dist * (math.sin(l) / math.sin(math.pi/2.-l))
            start_shift = [-dist, -dist, -j]
            end_shift = [dist, -dist, -j]
            catalogue_file.write(run_simple_ray(
                ds,
                z,
                start_shift,
                end_shift,
                dist,
                snapshot_name,
                galaxy_pos,
                f"{base_name}{n}"))
            n += 1
            # case: angles up to 90
            for i in range(18, 90, 18):
                l = math.radians(i)
                j = dist * (math.sin(l) / math.sin(math.pi/2.-l))
                start_shift = [-dist, -dist, j]
                end_shift = [dist, -dist, -j]
                catalogue_file.write(run_simple_ray(
                    ds,
                    z,
                    dist,
                    start_shift,
                    end_shift,
                    snapshot_name,
                    galaxy_pos,
                    f"{base_name}{n}"))
                n += 1
            # special case: 90
            start_shift = [0, -dist, dist]
            end_shift = [0, -dist, -dist]
            catalogue_file.write(run_simple_ray(
                ds,
                z,
                dist,
                start_shift,
                end_shift,
                snapshot_name,
                galaxy_pos,
                f"{base_name}{n}"))
            n += 1
            # case: angles from 90 to 180
            for i in range(72, 0, -18):
                l = math.radians(i)
                j = dist * (math.sin(l) / math.sin(math.pi/2.-l))
                start_shift = [dist, -dist, j]
                end_shift = [-dist, -dist, -j]
                catalogue_file.write(run_simple_ray(
                    ds,
                    z,
                    dist,
                    start_shift,
                    end_shift,
                    snapshot_name,
                    galaxy_pos,
                    f"{base_name}{n}"))
                n += 1

def run_quasar_snapshot(fn,
                        z_min,
                        z_max,
                        z_step,
                        base_name,
                        catalogue_file,
                        starting_n=1):
    """Run a set of simple rays for the specified distances to the center of the galaxy

    Arguments
    ---------
    fn: str
    Name of the snapshot

    galaxy_pos: array
    3D position of the galaxy in the snapshot

    n: int
    Current simulation number

    z: float
    Redshift

    dist_min: float
    Minimum distance to the center of the galaxy

    dist_max: float
    Maximum distance to the center of the galaxy

    dist_step: float
    Step in the distance coverage (from dist_min to dist_max).

    snapshot_name: str
    Name of the snapshot used (e.g. "RD0196")

    base_name: str
    Base name used to name the outputs

    catalogue_name: str
    The name of the catalogue

    starting_n: int
    First simulation number
    """
    snapshot_name = fn.split("/")[0]
    ds = load_snapshot(fn)

    n = starting_n
    for z in np.arange(z_min, z_max, z_step):
        for j in range (12061943, 12061953, 1):
            ray = trident.make_compound_ray(
                fn,
                simulation_type='Enzo',
                near_redshift=0,
                far_redshift=z,
                max_box_fraction= 1.4,
                lines='all', ftype='gas',
                fields=['density', 'temperature', 'metallicity'],
                solution_filename=f"{base_name}{n}ray.txt",
                data_filename=f"{base_name}{n}ray.h5",
                seed= j)
            spec_gen = trident.SpectrumGenerator(
                lambda_min=3000,
                lambda_max=9000,
                dlambda=0.8)
            spec_gen.make_spectrum(
                ray,
                lines='all',
                output_absorbers_file=f"{base_name}{n}absorbers.txt",
                store_observables=True)
            spec_gen.add_qso_spectrum(emitting_redshift=(i))
            spec_gen.save_spectrum(f"{base_name}{n}spec.h5")
            spec_gen.save_spectrum(f"{base_name}{n}spec.txt")
            name = "QSO_" + str(n)
            catalogue_file.write(
                f"{name}; {snapshot_name}; {i}; {j}; \n")

def run_simple_ray(ds,
                   z,
                   start_shift,
                   end_shift,
                   snapshot_name,
                   galaxy_pos,
                   base_name,
                   noise=None):
    """Run a simple ray from a specified start and end shifts from the centre
    of a galaxy

    Arguments
    ---------
    ds: yt.data_objects.static_output.Dataset
    The loaded snapshot

    z: float
    The redshift of the ray

    start_shift: float, array
    Shift of the starting point of the ray with respect to the galaxy centre
    If a float, all 3 dimensions are equally shifted.
    If an array, it must have size=3 and each dimension will be shifted independently

    end_shift: float, array
    Shift of the ending point of the ray with respect to the galaxy centre.
    If a float, all 3 dimensions are equally shifted.
    If an array, it must have size=3 and each dimension will be shifted independently.

    snapshot_name: str
    Name of the snapshot used (e.g. "RD0196")

    galaxy: array
    3D position of the galaxy in the snapshot

    base_name: str
    Base name used to name the outputs

    noise: float
    The noise to be applied to the spectrum. Ignored.

    Return
    ------
    entry: str
    A string to write in the catalogue
    """
    start = galaxy_pos[:] + start_shift
    end = galaxy_pos[:] + end_shift

    datastart = start*ds.units.kpc
    ray_start=datastart.to('code_length')
    dataend = end*ds.units.kpc
    ray_end=dataend.to('code_length')

    ray = trident.make_simple_ray(
        ds,
        start_position=ray_start,
        end_position=ray_end,
        redshift=z,
        lines='all',
        fields=['density', 'temperature', 'metallicity'],
        data_filename=f"{base_name}ray.h5")
    spec_gen = trident.SpectrumGenerator(
        lambda_min= 3000,
        lambda_max= 9000,
        dlambda=0.8)
    spec_gen.make_spectrum(
        ray,
        lines='all',
        output_absorbers_file=f"{base_name}absorbers.txt",
        store_observables=True)
    spec_gen.save_spectrum(f"{base_name}spec.h5")
    spec_gen.save_spectrum(f"{base_name}spec.txt")
    spec_gen.make_spectrum(
        ray, lines=['H I 1216'],
        output_absorbers_file=f"{base_name}HIabsorbers.txt",
        store_observables=True)
    spec_gen.save_spectrum(f"{base_name}HIspec.h5")
    spec_gen.save_spectrum(f"{base_name}HIspec.txt")
    nhi, b_kms, zfit = fit_lines(f"{base_name}HI")

    name = base_name.split("/")[-1]
    return (
        f"{name}; {snapshot_name}; {z}; "
        f"{start_shift[0]}; {start_shift[1]}; {start_shift[2]}; "
        f"{end_shift[0]}; {end_shift[1]}; {end_shift[2]}; "
        f"{nhi}; {b_kms}; {zfit}\n"
    )

def run_simple_ray_fast(ds,
                        z,
                        start_shift,
                        end_shift,
                        snapshot_name,
                        galaxy_pos,
                        base_name,
                        output_dir,
                        noise):
    """Run a simple ray from a specified start and end shifts from the centre
    of a galaxy

    Arguments
    ---------
    ds: yt.data_objects.static_output.Dataset
    The loaded snapshot

    z: float
    The redshift of the ray

    start_shift: float, array
    Shift of the starting point of the ray with respect to the galaxy centre
    If a float, all 3 dimensions are equally shifted.
    If an array, it must have size=3 and each dimension will be shifted independently

    end_shift: float, array
    Shift of the ending point of the ray with respect to the galaxy centre.
    If a float, all 3 dimensions are equally shifted.
    If an array, it must have size=3 and each dimension will be shifted independently.

    snapshot_name: str
    Name of the snapshot used (e.g. "RD0196")

    galaxy: array
    3D position of the galaxy in the snapshot

    base_name: str
    Base name used to name the outputs

    output_dir: str
    Directory where outputs are saved

    noise: float
    The noise to be applied to the spectrum
    """
    start = galaxy_pos[:] + start_shift
    end = galaxy_pos[:] + end_shift

    datastart = start*ds.units.kpc
    ray_start=datastart.to('code_length')
    dataend = end*ds.units.kpc
    ray_end=dataend.to('code_length')

    ray = trident.make_simple_ray(
        ds,
        start_position=ray_start,
        end_position=ray_end,
        redshift=z,
        lines='all',
        fields=['density', 'temperature', 'metallicity'],
        data_filename=f"{output_dir}{base_name}_ray.h5"
        )
    spec_gen = trident.SpectrumGenerator(
        lambda_min= 3000,
        lambda_max= 9000,
        dlambda=0.8)
    spec_gen.make_spectrum(
        ray,
        lines='all',
        store_observables=True)
    spec_gen.save_spectrum(
        f"{output_dir}{base_name}_spec_nonoise.fits.gz",
        format="FITS")
    if noise > 0.0:
        spec_gen.add_gaussian_noise(noise)
        spec_gen.save_spectrum(
            f"{output_dir}{base_name}_spec.fits.gz",
            format="FITS")
