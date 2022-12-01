import yt
import trident
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
# Ignasi: not sure what we are importing here. Not sure if it's used either
# I'm commenting it for the time being
#from yt.mods import *
import math
import h5py

# some configuration variables used by run_simple_ray
HI_parameters = {
    'name':'HI',
    'f': [.4164],
    'Gamma':[6.265E8],
    'wavelength':[1215.67],
    'numLines':1,
    'maxN': 1E22,
    'minN':1E11,
    'maxb': 300,
    'minb':1,
    'maxz': 6,
    'minz':0,
    'init_b':30,
    'init_N':1E14,
}
u = ds.units
datastart = start*u.kpc
ray_start=datastart.to('code_length')
v = ds.units
dataend = end*v.kpc
ray_end=dataend.to('code_length')
ds.derived_field_list
line_list = 'all'


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
        catalogue_file.write(write(";".join(["Name QSO", "Name Galaxy", "z QSO", "z Galaxy", "d", "Ang", "N [cm^-2]", "Absortion System"])))
    elif quasar:
        catalogue_file.write(";".join(["Name", "Dataset", "z", "seed"]))
    elif galaxy:
        catalogue_file.write(";".join(["Name", "Dataset", "z", "d", "Ang", "N [cm^-2]", "b [km/s]", "zfit"]))
    else:
        catalogue_file.close()
        raise IOError("Unkown catalogue to initialize")

    catalogue_file.write("\n")

    return catalogue_file


def load_snapshot(fn):
    """Load the simulation snapshot

    Arguments
    ---------
    fn: str
    Name of the snapshot

    Return
    ------
    ds: ?
    The loaded snapshot
    """
    ds = yt.load(fn)
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

def run_galaxy_snapshot(fn, galaxy_pos, z_min, z_max, z_step, dist_min, dist_max, dist_step, base_name, catalogue_file, starting_n=1):
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
    for z in range (z_min, z_max, z_step):
        #rayo centro
        base_name_alt = base_name.replace("/G_", "/G" + snapshot_name.replace("RD0", "")) + f"_{z:.1f}"
        if snapshot_name == "RD0180":
            #separamos x_end 1 kpc más línea S II 766 [765.684000 A] error en espectro MemoryError: Unable to allocate 6.10 GiB for an array
            run_simple_ray(z, 0.4, 0, snapshot_name, galaxy_pos, base_name_alt, catalogue_file, start_shift=10, end_shift=[11, 10, 10])
        else:
            run_simple_ray(z, 0, 0, snapshot_name, galaxy_pos, base_name_alt, catalogue_file, start_shift=10, end_shift=[10, 10, 10])

        # loop over distances
        for d in range (dist_min, dist_max, dist_step):
            run_simple_ray(z, d, 0, snapshot_name, galaxy_pos, f"{base_name}{n}", catalogue_file)
            n += 1
            # case: angles up to 90
            for i in range(18, 90, 18):
                run_simple_ray(z, d, i, snapshot_name, galaxy_pos, f"{base_name}{n}", catalogue_file)
                n += 1
            # special case: 90
            run_simple_ray(z, d, 90, snapshot_name, galaxy_pos, f"{base_name}{n}", catalogue_file)
            n += 1
            # case: angles from 90 to 180
            for i in range(72, 0, -18):
                run_simple_ray(z, d, 90, snapshot_name, galaxy_pos, f"{base_name}{n}", catalogue_file, neg_angles=True)
                n += 1

def run_quasar_snapshot(fn, z_min, z_max, z_step, base_name, catalogue_file, starting_n=1):
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
    for i in range(z_min, z_max, z_step):
        for j in range (12061943, 12061953, 1):
            ray = trident.make_compound_ray(
                fn,
                simulation_type='Enzo',
                near_redshift=0,
                far_redshift=(i),
                max_box_fraction= 1.4,
                lines='all', ftype='gas',
                fields=['density', 'temperature', 'metallicity'],
                solution_filename=f"{base_name}{n}ray.txt",
                data_filename=f"{base_name}{n}ray.h5",
                seed= (j))
            sg = trident.SpectrumGenerator(lambda_min=3000, lambda_max=9000, dlambda=0.8)
            sg.make_spectrum(ray, lines='all', 	output_absorbers_file=f"{base_name}{n}absorbers.txt", store_observables=True)
            sg.add_qso_spectrum(emitting_redshift=(i))
            sg.save_spectrum(f"{base_name}{n}spec.h5")
            sg.save_spectrum(f"{base_name}{n}spec.txt")
            name = "QSO_" + str(n)
            catalogue_file.write(name)
            catalogue_file.write(";")
            catalogue_file.write(snapshot_name)
            catalogue_file.write(";")
            catalogue_file.write(str(i))
            catalogue_file.write(";")
            catalogue_file.write(str(j))
            catalogue_file.wirte("\n")

def run_simple_ray(z, d, i, snapshot_name, galaxy_pos, base_name, catalogue_file, neg_angles=False, start_shift=None, end_shift=None):
    """Run simple rays for a specified distance to the center of the galaxy

    Arguments
    ---------
    z: float
    Redshift

    d: float
    Distance to the center of the galaxy

    i: float
    Incidence angle (in degrees)

    snapshot_name: str
    Name of the snapshot used (e.g. "RD0196")

    galaxy: array
    3D position of the galaxy in the snapshot

    base_name: str
    Base name used to name the outputs

    catalogue_name: str
    The name of the catalogue

    neg_angles: bool - Default: False
    Ignasi: Not sure how to call this, but when True what it does is to swap
    start = galaxy_pos[:] - [d,d,-j]
    end = galaxy_pos[:] + [d,-d,-j]
    by
    start = galaxy_pos[:] - [d,-d,j]
    end = galaxy_pos[:] + [d,d,j]
    and save 180-i instead of i

    start_shift: float, array or None - Default: None
    If not None, manually set the start_shift. If a float, all 3 dimensions
    are shifted similarly. If an array, it must have size=3 and each dimension
    will be shifted independently

    end_shift: float, array or None - Default: None
    If not None, manually set the end_shift. If a float, all 3 dimensions
    are shifted similarly. If an array, it must have size=3 and each dimension
    will be shifted independently.
    """
    # Ignasi: I think we are computing the ray trajectory here
    if start_shift is not None and end_shift is not None:
        start = galaxy_pos[:] - start_shift
        end = galaxy_pos[:] + end_shift
    elif start_shift is None and end_shift is None:
        raise Exception("Both or none of 'start_shift' and 'end_shift' must be specified")
    else:
        l = math.radians(i)
        # Ignasi: I think there is a mistake here, the 90-l should be math.pi-l
        j = d * (math.sin(l) / math.sin(math.pi/2.,-l))
        if neg_angles:
            start = galaxy_pos[:] - [d,-d,j]
            end = galaxy_pos[:] + [d,d,j]
        else:
            start = galaxy_pos[:] - [d,d,-j]
            end = galaxy_pos[:] + [d,-d,-j]
    ray = trident.make_simple_ray(
        ds,
        start_position=ray_start,
        end_position=ray_end,  redshift= (z),
        lines=line_list,
        fields=['density', 'temperature', 'metallicity'],
        data_filename=f"{base_name}ray.h5")
    sg = trident.SpectrumGenerator(lambda_min= 3000, lambda_max= 9000, dlambda=0.8)
    sg.make_spectrum(ray, lines='all', output_absorbers_file=f"{base_name}absorbers.txt", store_observables=True)
    sg.save_spectrum(f"{base_name}spec.h5")
    sg.save_spectrum(f"{base_name}spec.txt")
    sg.make_spectrum(ray, lines=['H I 1216'], output_absorbers_file=f"{base_name}HIabsorbers.txt", store_observables=True)
    sg.save_spectrum(f"{base_name}HIspec.h5")
    sg.save_spectrum(f"{base_name}HIspec.txt")
    f = h5py.File(f"{base_name}HIspec.h5")
    flux = f['flux'][:]
    wavelength = f["wavelength"][:]
    f.close()
    speciesDicts = {'HI':HI_parameters}
    orderFits = ['HI']
    fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts, maxLength=9000, output_file=f"{base_name}HIspecFIT.h5")
    fitted_lines
    if fitted_lines['HI']['N'].size > 0:
        a=fitted_lines['HI']['N'][0]
        b=fitted_lines['HI']['b'][0]
        c=fitted_lines['HI']['z'][0]
    else:
        a = 0
        b = 0
        c = 0
    name = base_name.spit("/")[-1] + str(n)
    catalogue_file.write(name)
    catalogue_file.write(";")
    catalogue_file.write(snapshot_name)
    catalogue_file.write(";")
    catalogue_file.write(str(z))
    catalogue_file.write(";")
    catalogue_file.write(str(d))
    catalogue_file.write(";")
    if neg_angles:
        catalogue_file.write(str(180 - i))
    else:
        catalogue_file.write(str(i))
    catalogue_file.write(";")
    catalogue_file.write(str(a))
    catalogue_file.write(";")
    catalogue_file.write(str(b))
    catalogue_file.write(";")
    catalogue_file.write(str(c))
    catalogue_file.write("\n")
