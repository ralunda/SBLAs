import yt   
from yt.units import kpc
from yt.units import Mpc
from yt import YTArray
import numpy as np
import trident
from yt.mods import *
import math
import h5py
fn = 'RD0196/RD0196'
ds = yt.load(fn)
qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
qc.write(";".join(["Name", "Dataset", "z", "d", "Ang", "N [cm^-2]", "b [km/s]", "zfit"]))
qc.close()


Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")

galaxy=np.array([11928.448765067616,12822.112568926215,12327.245902774346])

n=1
for k in range (22, 26, 1):
    z = round(k*0.1,1)
    #rayo centro
    start=galaxy[:] - 10
    end=galaxy[:] + [10,10,10]
    u = ds.units
    datastart = start*u.kpc
    ray_start=datastart.to('code_length')
    v = ds.units
    dataend = end*v.kpc
    ray_end=dataend.to('code_length')
    ds.derived_field_list
    line_list = 'all'
    ray = trident.make_simple_ray(ds, 
                              start_position=ray_start,
                              end_position=ray_end,  redshift= (z), 
          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_" + str(z)+ "_ray.h5")
    sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
    sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'absorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'spec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'spec.txt')
    sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'HIabsorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'HIspec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'HIspec.txt')
    f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'HIspec.h5')
    flux = f['flux'][:]
    wavelength = f["wavelength"][:]
    f.close()
    HI_parameters = {'name':'HI',
    'f': [.4164],
    'Gamma':[6.265E8],
    'wavelength':[1215.67],
    'numLines':1,
    'maxN': 1E22, 'minN':1E11,
    'maxb': 300, 'minb':1,
    'maxz': 6, 'minz':0,
    'init_b':30,
    'init_N':1E14}
    speciesDicts = {'HI':HI_parameters}
    from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
    orderFits = ['HI']
    fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G196_' + str(z)+ 'HIspecFIT.h5')
    fitted_lines
    if fitted_lines['HI']['N'].size > 0:
        a=fitted_lines['HI']['N'][0]
        b=fitted_lines['HI']['b'][0]
        c=fitted_lines['HI']['z'][0]
    else:
        a = 0
        b = 0
        c = 0
    name = "G196_" + str(z)
    qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
    qc.write("\n")
    qc.write(name)
    qc.write(";")
    qc.write("RD0196")
    qc.write(";")
    qc.write(str(z))
    qc.write(";")
    qc.write("0")
    qc.write(";")
    qc.write("0")
    qc.write(";")
    qc.write(str(a))
    qc.write(";")
    qc.write(str(b))
    qc.write(";")
    qc.write(str(c))
    qc.close()
    for d in range (2, 22, 2):
        start=galaxy[:] - [d,d,0]
        end=galaxy[:] + [d,-d,0]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0196")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("0")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                        
        for i in range(18, 90, 18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] - [d,d,-j]
            end = galaxy[:] + [d,-d,-j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                       
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0196")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(i))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1                                
        start=galaxy[:] + [0,-d,d]
        end=galaxy[:] - [0,d,d]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')        
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0196")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("90")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                                    
        for i in range(72, 0, -18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] + [d,-d,j]
            end = galaxy[:] - [d,d,j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')            
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            h = 180 - i
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0196")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(h))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1

#sim 180
            

fn = 'RD0180/RD0180'
ds = yt.load(fn)

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")

galaxy=np.array([10501.472788071803,11192.100157077519,10801.248702515462])

n=401
for k in range (22, 26, 1):
    z = round(k*0.1,1)
    #rayo centro
    start=galaxy[:] - 10
    end=galaxy[:] + [11,10,10] #separamos x_end 1 kpc más línea S II 766 [765.684000 A] error en espectro MemoryError: Unable to allocate 6.10 GiB for an array 
    u = ds.units
    datastart = start*u.kpc
    ray_start=datastart.to('code_length')
    v = ds.units
    dataend = end*v.kpc
    ray_end=dataend.to('code_length')
    ds.derived_field_list
    line_list = 'all'
    ray = trident.make_simple_ray(ds, 
                              start_position=ray_start,
                              end_position=ray_end,  redshift= (z), 
          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_" + str(z)+ "_ray.h5")
    sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
    sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'absorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'spec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'spec.txt')
    sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'HIabsorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'HIspec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'HIspec.txt')
    f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'HIspec.h5')
    flux = f['flux'][:]
    wavelength = f["wavelength"][:]
    f.close()
    HI_parameters = {'name':'HI',
    'f': [.4164],
    'Gamma':[6.265E8],
    'wavelength':[1215.67],
    'numLines':1,
    'maxN': 1E22, 'minN':1E11,
    'maxb': 300, 'minb':1,
    'maxz': 6, 'minz':0,
    'init_b':30,
    'init_N':1E14}
    speciesDicts = {'HI':HI_parameters}
    from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
    orderFits = ['HI']
    fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G180_' + str(z)+ 'HIspecFIT.h5')
    fitted_lines
    if fitted_lines['HI']['N'].size > 0:
        a=fitted_lines['HI']['N'][0]
        b=fitted_lines['HI']['b'][0]
        c=fitted_lines['HI']['z'][0]
    else:
        a = 0
        b = 0
        c = 0
    name = "G180_" + str(z)
    qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
    qc.write("\n")
    qc.write(name)
    qc.write(";")
    qc.write("RD0180")
    qc.write(";")
    qc.write(str(z))
    qc.write(";")
    qc.write("0.4")
    qc.write(";")
    qc.write("0")
    qc.write(";")
    qc.write(str(a))
    qc.write(";")
    qc.write(str(b))
    qc.write(";")
    qc.write(str(c))
    qc.close()
    for d in range (17, 187, 17):
        start=galaxy[:] - [d,d,0]
        end=galaxy[:] + [d,-d,0]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0180")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("0")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                        
        for i in range(18, 90, 18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] - [d,d,-j]
            end = galaxy[:] + [d,-d,-j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                       
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0180")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(i))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1                                
        start=galaxy[:] + [0,-d,d]
        end=galaxy[:] - [0,d,d]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')        
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0180")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("90")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                                    
        for i in range(72, 0, -18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] + [d,-d,j]
            end = galaxy[:] - [d,d,j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')            
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            h = 180 - i
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0180")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(h))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1

#sim 165

fn = 'RD0165/RD0165'
ds = yt.load(fn)

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")

galaxy=np.array([9350.457417869537,9893.696818590259,9571.331842591271])

n=801
for k in range (22, 26, 1):
    z = round(k*0.1,1)
    #rayo centro
    start=galaxy[:] - 10
    end=galaxy[:] + [10,10,10]
    u = ds.units
    datastart = start*u.kpc
    ray_start=datastart.to('code_length')
    v = ds.units
    dataend = end*v.kpc
    ray_end=dataend.to('code_length')
    ds.derived_field_list
    line_list = 'all'
    ray = trident.make_simple_ray(ds, 
                              start_position=ray_start,
                              end_position=ray_end,  redshift= (z), 
          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_" + str(z)+ "_ray.h5")
    sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
    sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'absorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'spec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'spec.txt')
    sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'HIabsorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'HIspec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'HIspec.txt')
    f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'HIspec.h5')
    flux = f['flux'][:]
    wavelength = f["wavelength"][:]
    f.close()
    HI_parameters = {'name':'HI',
    'f': [.4164],
    'Gamma':[6.265E8],
    'wavelength':[1215.67],
    'numLines':1,
    'maxN': 1E22, 'minN':1E11,
    'maxb': 300, 'minb':1,
    'maxz': 6, 'minz':0,
    'init_b':30,
    'init_N':1E14}
    speciesDicts = {'HI':HI_parameters}
    from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
    orderFits = ['HI']
    fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G165_' + str(z)+ 'HIspecFIT.h5')
    fitted_lines
    if fitted_lines['HI']['N'].size > 0:
        a=fitted_lines['HI']['N'][0]
        b=fitted_lines['HI']['b'][0]
        c=fitted_lines['HI']['z'][0]
    else:
        a = 0
        b = 0
        c = 0
    name = "G165_" + str(z)
    qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
    qc.write("\n")
    qc.write(name)
    qc.write(";")
    qc.write("RD0165")
    qc.write(";")
    qc.write(str(z))
    qc.write(";")
    qc.write("0")
    qc.write(";")
    qc.write("0")
    qc.write(";")
    qc.write(str(a))
    qc.write(";")
    qc.write(str(b))
    qc.write(";")
    qc.write(str(c))
    qc.close()
    for d in range (30, 330, 30):
        start=galaxy[:] - [d,d,0]
        end=galaxy[:] + [d,-d,0]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0165")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("0")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                        
        for i in range(18, 90, 18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] - [d,d,-j]
            end = galaxy[:] + [d,-d,-j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                       
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0165")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(i))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1                                
        start=galaxy[:] + [0,-d,d]
        end=galaxy[:] - [0,d,d]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')        
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0165")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("90")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                                    
        for i in range(72, 0, -18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] + [d,-d,j]
            end = galaxy[:] - [d,d,j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')            
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            h = 180 - i
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0165")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(h))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1

#sim 152

fn = 'RD0152/RD0152'
ds = yt.load(fn)

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")

galaxy=np.array([8440.156388155721,8875.605892201605,8606.241907442962])

n=1201
for k in range (22, 26, 1):
    z = round(k*0.1,1)
    #rayo centro
    start=galaxy[:] - 10
    end=galaxy[:] + [10,10,10]
    u = ds.units
    datastart = start*u.kpc
    ray_start=datastart.to('code_length')
    v = ds.units
    dataend = end*v.kpc
    ray_end=dataend.to('code_length')
    ds.derived_field_list
    line_list = 'all'
    ray = trident.make_simple_ray(ds, 
                              start_position=ray_start,
                              end_position=ray_end,  redshift= (z), 
          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_" + str(z)+ "_ray.h5")
    sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
    sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'absorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'spec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'spec.txt')
    sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'HIabsorbers.txt', store_observables=True)
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'HIspec.h5')
    sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'HIspec.txt')
    f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'HIspec.h5')
    flux = f['flux'][:]
    wavelength = f["wavelength"][:]
    f.close()
    HI_parameters = {'name':'HI',
    'f': [.4164],
    'Gamma':[6.265E8],
    'wavelength':[1215.67],
    'numLines':1,
    'maxN': 1E22, 'minN':1E11,
    'maxb': 300, 'minb':1,
    'maxz': 6, 'minz':0,
    'init_b':30,
    'init_N':1E14}
    speciesDicts = {'HI':HI_parameters}
    from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
    orderFits = ['HI']
    fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G152_' + str(z)+ 'HIspecFIT.h5')
    fitted_lines
    if fitted_lines['HI']['N'].size > 0:
        a=fitted_lines['HI']['N'][0]
        b=fitted_lines['HI']['b'][0]
        c=fitted_lines['HI']['z'][0]
    else:
        a = 0
        b = 0
        c = 0
    name = "G152_" + str(z)
    qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
    qc.write("\n")
    qc.write(name)
    qc.write(";")
    qc.write("RD0152")
    qc.write(";")
    qc.write(str(z))
    qc.write(";")
    qc.write("0")
    qc.write(";")
    qc.write("0")
    qc.write(";")
    qc.write(str(a))
    qc.write(";")
    qc.write(str(b))
    qc.write(";")
    qc.write(str(c))
    qc.close()
    for d in range (30, 330, 30):
        start=galaxy[:] - [d,d,0]
        end=galaxy[:] + [d,-d,0]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0152")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("0")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                        
        for i in range(18, 90, 18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] - [d,d,-j]
            end = galaxy[:] + [d,-d,-j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                       
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0152")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(i))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1                                
        start=galaxy[:] + [0,-d,d]
        end=galaxy[:] - [0,d,d]
        u = ds.units
        datastart = start*u.kpc
        ray_start=datastart.to('code_length')
        v = ds.units
        dataend = end*v.kpc
        ray_end=dataend.to('code_length')
        ds.derived_field_list
        line_list = 'all'
        ray = trident.make_simple_ray(ds,
                                  start_position=ray_start,
                                  end_position=ray_end,  redshift= (z),                              
                                  lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")                
        sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
        sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
        sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')        
        f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
        flux = f['flux'][:]
        wavelength = f["wavelength"][:]
        f.close()
        HI_parameters = {'name':'HI',
        'f': [.4164],
        'Gamma':[6.265E8],
        'wavelength':[1215.67],
        'numLines':1,
        'maxN': 1E22, 'minN':1E11,
        'maxb': 300, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14}
        speciesDicts = {'HI':HI_parameters}
        from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
        orderFits = ['HI']
        fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
        fitted_lines
        if fitted_lines['HI']['N'].size > 0:
            a=fitted_lines['HI']['N'][0]
            b=fitted_lines['HI']['b'][0]
            c=fitted_lines['HI']['z'][0]
        else:
            a = 0
            b = 0
            c = 0
        name = "G_" + str(n)
        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
        qc.write("\n")
        qc.write(name)
        qc.write(";")
        qc.write("RD0152")
        qc.write(";")
        qc.write(str(z))
        qc.write(";")
        qc.write(str(d))
        qc.write(";")
        qc.write("90")
        qc.write(";")
        qc.write(str(a))
        qc.write(";")
        qc.write(str(b))
        qc.write(";")
        qc.write(str(c))
        qc.close()
        n += 1                                    
        for i in range(72, 0, -18):
            l=math.radians(i)
            j = d * (math.sin(l) / math.sin(90-l))
            start = galaxy[:] + [d,-d,j]
            end = galaxy[:] - [d,d,j]
            u = ds.units
            datastart = start*u.kpc
            ray_start=datastart.to('code_length')
            v = ds.units
            dataend = end*v.kpc
            ray_end=dataend.to('code_length')
            ds.derived_field_list
            line_list = 'all'
            ray = trident.make_simple_ray(ds,
                                          start_position=ray_start,
                                          end_position=ray_end,  redshift= (z),                              
                                          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_" + str(n)+ "ray.h5")
            sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'absorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'spec.txt')
            sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIabsorbers.txt', store_observables=True)
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.txt')            
            f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspec.h5')
            flux = f['flux'][:]
            wavelength = f["wavelength"][:]
            f.close()
            HI_parameters = {'name':'HI',
            'f': [.4164],
            'Gamma':[6.265E8],
            'wavelength':[1215.67],
            'numLines':1,
            'maxN': 1E22, 'minN':1E11,
            'maxb': 300, 'minb':1,
            'maxz': 6, 'minz':0,
            'init_b':30,
            'init_N':1E14}
            speciesDicts = {'HI':HI_parameters}
            from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
            orderFits = ['HI']
            fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_' + str(n)+ 'HIspecFIT.h5')
            fitted_lines
            if fitted_lines['HI']['N'].size > 0:
                a=fitted_lines['HI']['N'][0]
                b=fitted_lines['HI']['b'][0]
                c=fitted_lines['HI']['z'][0]
            else:
                a = 0
                b = 0
                c = 0
            name = "G_" + str(n)
            h = 180 - i
            qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv","a")
            qc.write("\n")
            qc.write(name)
            qc.write(";")
            qc.write("RD0152")
            qc.write(";")
            qc.write(str(z))
            qc.write(";")
            qc.write(str(d))
            qc.write(";")
            qc.write(str(h))
            qc.write(";")
            qc.write(str(a))
            qc.write(";")
            qc.write(str(b))
            qc.write(";")
            qc.write(str(c))
            qc.close()
            n += 1


            
            
            
            

          



            
            
            
            



            
            
            
            

          


          
            
            

          


