import yt   
from yt.units import kpc
from yt.units import Mpc
from yt import YTArray
import numpy as np
import trident
from yt.mods import *
fn = 'RD0196/RD0196'
ds = yt.load(fn)
qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_catalog.csv","a")
qc.write(";".join(["Name", "Dataset", "z", "seed"]))
qc.close()

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")
ds.derived_field_list
n = 1
for k in range(250, 255, 5):
    i = k*0.01
    for j in range (12061943, 12061953, 1):
        ray = trident.make_compound_ray(fn, simulation_type='Enzo', near_redshift=0, far_redshift=(i), max_box_fraction= 1.4, lines='all', ftype='gas', fields=['density', 'temperature', 'metallicity'], solution_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO" + str(n) + "ray.txt", data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_" + str(n) + "ray.h5", seed= (j))
        sg = trident.SpectrumGenerator(lambda_min=3000, lambda_max=9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', 	output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.add_qso_spectrum(emitting_redshift=(i))
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.txt')
        name = "QSO_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0196")
        qc.write(";")
        qc.write(str(i))
        qc.write(";")
        qc.write(str(j))
        qc.close()
        n += 1



fn = 'RD0180/RD0180'
ds = yt.load(fn)

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")
ds.derived_field_list
n = 11
for k in range(250, 305, 5):
    i = k*0.01
    for j in range (12061943, 12061953, 1):
        ray = trident.make_compound_ray(fn, simulation_type='Enzo', near_redshift=0, far_redshift=(i), max_box_fraction= 1.4, lines='all', ftype='gas', fields=['density', 'temperature', 'metallicity'], solution_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO" + str(n) + "ray.txt", data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_" + str(n) + "ray.h5", seed= (j))
        sg = trident.SpectrumGenerator(lambda_min=3000, lambda_max=9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', 	output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.add_qso_spectrum(emitting_redshift=(i))
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.txt')
        name = "QSO_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0180")
        qc.write(";")
        qc.write(str(i))
        qc.write(";")
        qc.write(str(j))
        qc.close()
        n += 1



fn = 'RD0165/RD0165'
ds = yt.load(fn)

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")
ds.derived_field_list
n = 121
for k in range(250, 355, 5):
    i = k*0.01
    for j in range (12061943, 12061953, 1):
        ray = trident.make_compound_ray(fn, simulation_type='Enzo', near_redshift=0, far_redshift=(i), max_box_fraction= 1.4, lines='all', ftype='gas', fields=['density', 'temperature', 'metallicity'], solution_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO" + str(n) + "ray.txt", data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_" + str(n) + "ray.h5", seed= (j))
        sg = trident.SpectrumGenerator(lambda_min=3000, lambda_max=9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', 	output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.add_qso_spectrum(emitting_redshift=(i))
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.txt')
        name = "QSO_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0165")
        qc.write(";")
        qc.write(str(i))
        qc.write(";")
        qc.write(str(j))
        qc.close()
        n += 1


fn = 'RD0152/RD0152'
ds = yt.load(fn)

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")
ds.derived_field_list
n = 331
for k in range(250, 405, 5):
    i = k*0.01
    for j in range (12061943, 12061953, 1):
        ray = trident.make_compound_ray(fn, simulation_type='Enzo', near_redshift=0, far_redshift=(i), max_box_fraction= 1.4, lines='all', ftype='gas', fields=['density', 'temperature', 'metallicity'], solution_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO" + str(n) + "ray.txt", data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_" + str(n) + "ray.h5", seed= (j))
        sg = trident.SpectrumGenerator(lambda_min=3000, lambda_max=9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', 	output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.add_qso_spectrum(emitting_redshift=(i))
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_' + str(n) + 'spec.txt')
        name = "QSO_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/QSO_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0152")
        qc.write(";")
        qc.write(str(i))
        qc.write(";")
        qc.write(str(j))
        qc.close()
        n += 1










