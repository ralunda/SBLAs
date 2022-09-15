#OBTENCIÓN DE UN QSO EN 165 z=3.5
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
ds.derived_field_list
ray = trident.make_compound_ray(fn, simulation_type='Enzo', near_redshift=0, far_redshift=3.5, max_box_fraction= 1.4, lines='all', ftype='gas', fields=['density', 'temperature', 'metallicity'], solution_filename="/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qray.txt",  data_filename="/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qray.h5", seed= 2822863)
p = yt.ProjectionPlot(ds, 'x', 'density')
p.annotate_ray(ray, arrow=True)
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qxray.png')
p = yt.ProjectionPlot(ds, 'y', 'density')
p.annotate_ray(ray, arrow=True)
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qyray.png')
p = yt.ProjectionPlot(ds, 'z', 'density')
p.annotate_ray(ray, arrow=True)
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qzray.png')

#PARA PLOT PANEL 3 RAYOS EN OVERLEAF
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

fig = plt.figure()
grid = AxesGrid(fig, (0.067,0.075,0.85,0.85),
 nrows_ncols = (2, 2),
 axes_pad = 2,
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
p = yt.ProjectionPlot(ds, 'x', 'density')
p.annotate_ray(ray, arrow=True)
p.set_font({'family': 'sans-serif', 'size': 15})
plot = p.plots[('gas','density')]
plot.figure = fig
plot.axes = grid[0].axes
plot.cax = grid.cbar_axes[0]
p._setup_plots()
p = yt.ProjectionPlot(ds, 'y', 'density')
p.annotate_ray(ray, arrow=True)
p.set_font({'family': 'sans-serif', 'size': 15})
plot = p.plots[('gas','density')]
plot.figure = fig
plot.axes = grid[1].axes
plot.cax = grid.cbar_axes[1]
p._setup_plots()
p = yt.ProjectionPlot(ds, 'z', 'density')
p.annotate_ray(ray, arrow=True)
p.set_font({'family': 'sans-serif', 'size': 15})
plot = p.plots[('gas','density')]
plot.figure = fig
plot.axes = grid[2].axes
plot.cax = grid.cbar_axes[2]
p._setup_plots()
fig.set_size_inches(12.90, 12.90)
fig.delaxes(grid[3].axes)
fig.delaxes(grid.cbar_axes[3])
p.save('/home/ralunda/anaconda3/simulaciones/ENZO/165/panel165-1-Qray.png')
plt.clf()


sg = trident.SpectrumGenerator(lambda_min=3000, lambda_max=9000, dlambda=0.8)
sg.make_spectrum(ray, lines='all', output_absorbers_file='/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qabsorbers.txt', store_observables=True)
sg.add_qso_spectrum(emitting_redshift=3.5)
sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qspec.h5')
sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qspec.txt')
##PARA PLOT USAMOS MATPLOT
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 
fig, ax = plt.subplots()
data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qspec.txt", delimiter=' ', unpack=True)
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
plt.title('165 QSO z = 3.5', fontsize=20)
plt.show()
plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qspec.png')

column_density_HI = ray.r[('gas', 'H_p0_number_density')] * ray.r[('gas', 'dl')] 
print('HI Column Density = %g' % column_density_HI.sum())

## 5.37326e+18
## DENSIDAD DE COLUMNA 
import yt
import trident
from yt import YTArray
import numpy as np
import h5py
f = h5py.File('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qspec.h5')
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
'init_N':1E11}
speciesDicts = {'HI':HI_parameters}
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
orderFits = ['HI']
fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file='/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QspecFIT.h5')
fitted_lines
