#Code for a 165 DLA at z=2.5
import yt   
from yt.units import kpc
from yt.units import Mpc
from yt import YTArray
import numpy as np
import trident
from yt.mods import *
fn = 'RD0165/RD0165'
ds = yt.load(fn)

Z_Solar= 0.02041
def _metallicity_e(field, data): 
        return data["gas", "metal_density"] / data["gas", "density"]/Z_Solar

ds.add_field(
        ("gas", "metallicity"), 
        function=_metallicity_e, sampling_type="local", force_override=True, display_name="Z/Z$_{\odot}$", take_log=True, units="")

galaxy=np.array([9350.457417869537,9893.696818590259,9571.331842591271])
start=galaxy[:] - 300
end=galaxy[:] + [310,300,300]
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
                              end_position=ray_end,  redshift= 2.5, 
          lines=line_list, fields=['density', 'temperature', 'metallicity'], data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dray.h5")
p = yt.SlicePlot( ds,  "z",  ("gas", "density"),
    center=YTArray([9350.457417869537,9893.696818590259,9571.331842591271],'kpc'),
    width=(30, "kpc"),    buff_size=(10000, 10000))
p.annotate_ray(ray, arrow=True)
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dzray.png')
p = yt.SlicePlot( ds,  "x",  ("gas", "density"),
    center=YTArray([9350.457417869537,9893.696818590259,9571.331842591271],'kpc'),
    width=(30, "kpc"),    buff_size=(10000, 10000))
p.annotate_ray(ray, arrow=True)
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dxray.png')
p = yt.SlicePlot( ds,  "y",  ("gas", "density"),
    center=YTArray([9350.457417869537,9893.696818590259,9571.331842591271],'kpc'),
    width=(30, "kpc"),    buff_size=(10000, 10000))
p.annotate_ray(ray, arrow=True)
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dyray.png')

#PARA PLOT PANEL 3 RAYOS EN OVERLEAF
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

fig = plt.figure()
grid = AxesGrid(fig, (0.065,0.075,0.85,0.85),
 nrows_ncols = (2, 2),
 axes_pad = 1.9,
 label_mode = "all",
 share_all = True,
 cbar_location="right",
 cbar_mode="each",
 cbar_size="3%",
 cbar_pad="0%")
ds.coordinates.x_axis[0] = 2
ds.coordinates.x_axis['x'] = 2
ds.coordinates.y_axis[0] = 1
ds.coordinates.y_axis['x'] = 1
p=yt.SlicePlot(ds,'z',('gas','density'),width=(30, "kpc"),center=YTArray([9350.457417869537,9893.696818590259,9571.331842591271],'kpc'), buff_size=(10000, 10000))
p.annotate_ray(ray, arrow=True)
p.set_font({'family': 'sans-serif', 'size': 15})
plot = p.plots[('gas','density')]
plot.figure = fig
plot.axes = grid[0].axes
plot.cax = grid.cbar_axes[0]
p._setup_plots()
p=yt.SlicePlot(ds,'y',('gas','density'),width=(30, "kpc"),center=YTArray([9350.457417869537,9893.696818590259,9571.331842591271],'kpc'), buff_size=(10000, 10000))
p.annotate_ray(ray, arrow=True)
p.set_font({'family': 'sans-serif', 'size': 15})
plot = p.plots[('gas','density')]
plot.figure = fig
plot.axes = grid[1].axes
plot.cax = grid.cbar_axes[1]
p._setup_plots()
p=yt.SlicePlot(ds,'x',('gas','density'),width=(30, "kpc"),center=YTArray([9350.457417869537,9893.696818590259,9571.331842591271],'kpc'), buff_size=(10000, 10000))
p.annotate_ray(ray, arrow=True)
p.set_font({'family': 'sans-serif', 'size': 15})
plot = p.plots[('gas','density')]
plot.figure = fig
plot.axes = grid[2].axes
plot.cax = grid.cbar_axes[2]
p._setup_plots()
fig.set_size_inches(12.80, 12.80)
fig.delaxes(grid[3].axes)
fig.delaxes(grid.cbar_axes[3])
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/panel165-1-Dray.png')
plt.clf()


sg = trident.SpectrumGenerator(lambda_min=3000, lambda_max=9000, dlambda=0.8)
sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dabsorbers.txt', store_observables=True)
sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dspec.h5')
sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dspec.txt')

## sg.plot_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/165/C165-1-Dspec.png')

##PARA PLOT USAMOS MATPLOT
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 
fig, ax = plt.subplots()
data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dspec.txt", delimiter=' ', unpack=True)
x = data[0] # wavelength
y1 = data[1] 
y2 = data[2] # flux
y3 = data[3] 
plt.figure(figsize=(12.4,4.8))
plt.xticks(fontsize=17) 
plt.yticks(fontsize=17) 
plt.plot(x,y2) 
plt.xlim(3000, 9000) 
plt.ylim(0, ) 
plt.xlabel('Wavelength [$\AA$]', fontsize=20)
plt.ylabel('Relative Flux', fontsize=20)
plt.title('165 DLA z = 2.5', fontsize=20)
plt.show()
plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dspec.png')
plt.close()

column_density_HI = ray.r[('gas', 'H_p0_number_density')] * ray.r[('gas', 'dl')]
print('HI Column Density = %g' % column_density_HI.sum())
## DENSIDAD DE COLUMNA
import yt
import trident
from yt import YTArray
import numpy as np
import h5py
f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Dspec.h5')
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
fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/165/1651DspecFIT.h5')
fitted_lines
