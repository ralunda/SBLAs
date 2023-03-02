import yt
import trident
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
# Ignasi: not sure what we are importing here. Not sure if it's used either
# I'm commenting it for the time being
#from yt.mods import *
import math
import h5py
import numpy as np

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
line_list = 'all'

Z_Solar = 0.02041

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
    ds: ?
    The loaded snapshot
    """
    ds = yt.load(dir+fn)
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
                start_shift = [-dist, -dist, -j]
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
                   dist,
                   start_shift,
                   end_shift,
                   snapshot_name,
                   galaxy_pos,
                   base_name):
    """Run a simple ray from a specified start and end shifts from the centre
    of a galaxy

    Arguments
    ---------
    ds: ?
    The loaded snapshot

    z: float
    The redshift of the ray

    dist: float
    Distance to the centre of the galaxy

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
        lines=line_list,
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
    file = h5py.File(f"{base_name}HIspec.h5")
    flux = file['flux'][:]
    wavelength = file["wavelength"][:]
    file.close()
    speciesDicts = {'HI':HI_parameters}
    orderFits = ['HI']
    fitted_lines, fitted_flux = generate_total_fit(
        wavelength,
        flux,
        orderFits,
        speciesDicts,
        maxLength=9000,
        output_file=f"{base_name}HIspecFIT.h5")
    if fitted_lines['HI']['N'].size > 0:
        nhi = fitted_lines['HI']['N'][0]
        b_kms = fitted_lines['HI']['b'][0]
        zfit = fitted_lines['HI']['z'][0]
    else:
        nhi = 0
        b_kms = 0
        zfit = 0
    name = base_name.split("/")[-1]
    return (
        f"{name}; {snapshot_name}; {z}; {dist}; "
        f"{start_shift[0]}; {start_shift[1]}; {start_shift[2]}; "
        f"{end_shift[0]}; {end_shift[1]}; {end_shift[2]}; "
        f"{nhi}; {b_kms}; {zfit}\n"
    )
