"""
This file contains a script to run the galaxy simulations using a
uniformly distributed random rays through the 3D volume of the galaxy.
"""
import argparse
import os
from itertools import repeat
import multiprocessing
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
import time

from yt.utilities.logger import set_log_level as set_log_level_yt
from yt.config import ytcfg

from utils import (
    fit_lines,
    load_snapshot,
    run_galaxy_snapshot,
    run_simple_ray,
    run_simple_ray_fast
)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def generate_ray(rho, theta_e, theta_r, phi_r, radius):
    """Function to compute the starting and ending points of a ray

    Arguments
    ---------
    rho: array of float
    Distance between the ray and the centre of the galaxy

    theta_e: array of float
    Angle within the cone of possible rays

    theta_r: array of float
    Rotatiom angle theta for the starting and ending points

    phi_r: array of float
    Rotation angle phi for the starting and ending points

    radius: array of float
    Radius of the sphere including the starting and ending point

    Return
    ------
    x_start: array of float
    x coordinate of the starting point

    y_start: array of float
    y coordinate of the starting point

    z_start: array of float
    z coordinate of the starting point

    x_end: array of float
    x coordinate of the ending point

    y_end: array of float
    y coordinate of the ending point

    z_end: array of float
    z coordinate of the ending point
    """
    # ray start
    x_start, y_start, z_start = rotate(-radius, 0, 0, theta_r, phi_r)

    # ray end
    circle_radius = np.sqrt(4*rho**2-4*rho**4/radius**2)
    x_end_norot = radius - 2*rho**2/radius
    y_end_norot = circle_radius*np.cos(theta_e)
    z_end_norot = circle_radius*np.sin(theta_e)
    x_end, y_end, z_end = rotate(x_end_norot, y_end_norot, z_end_norot, theta_r, phi_r)

    return x_start, y_start, z_start, x_end, y_end, z_end

def rotate(x, y, z, theta, phi):
    """Rotate coordinates around the centre.
    Two rotations are performed following:
    - A rotation on the z axis with angle theta
    - A rotation on the y axis with angle phi

    Arguments
    ---------
    x : float
    x coordinate of the point to rotate

    y : float
    y coordinate of the point to rotate

    z : float
    z coordinate of the point to rotate

    theta: float
    Rotation angle along z axis

    phi: float
    Rotation angle along y axis

    Return
    ------
    x2: float
    x coordinate of the rotated point

    y2: float
    y coordinate of the rotated point

    z2: float
    z coordinate of the rotated point
    """
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)

    # rotation on z axis
    x1 = cos_theta*x - sin_theta*y
    y1 = sin_theta*x + cos_theta*y
    z1 = z

    # rotation on y axis
    x2 = cos_phi*x1 + sin_phi*z1
    y2 = y1
    z2 = -sin_phi*x1 + cos_phi*z1

    return x2, y2, z2

def find_rho_max(redshifs, snapshots):
    """Find the maximum rho at a given redshift for the specified snapshots

    Arguments
    ---------
    redshif: array of float
    Chosen redshifts

    snapshots: np.structured_array
    A named array with the snapshots information. Must contain fields "name",
    "rho_max" and "z_max"

    Return
    ------
    rho_max: array of float
    The maximum rho
    """
    rho_max = np.zeros_like(redshifs)
    for index, redshif in enumerate(redshifs):
        pos = np.argwhere((snapshots["z_max"] > redshif))
        rho_max[index] = np.max(snapshots["rho_max"][pos])
    return rho_max

def select_snapshot(redshif, rho, snapshots):
    """Randomly select a snapshot given a choice of z and distance

    Arguments
    ---------
    redshif: float
    Chosen redshift

    rho: float
    Chosen distance

    snapshots: np.structured_array
    A named array with the snapshots information. Must contain fields "name",
    "rho_max" and "z_max"

    Return
    ------
    choice: str
    The position of the chosen
    """
    pos = np.argwhere((snapshots["rho_max"] > rho) & (snapshots["z_max"] > redshif))
    choice = np.random.choice(pos.reshape(pos.shape[0]))
    return choice

def main(args):
    """Run the galaxy simulations with rays with a random distribution

    Arguments
    ---------
    args: argparse.Namespace
    Parsed arguments (parser at the bottom of the file).
    """
    if args.num_processors == 0:
        args.num_processors = (multiprocessing.cpu_count() // 2)
    print("Running with nproc = ", args.num_processors)

    # set yt log level
    if not args.verbose:
        set_log_level_yt("error")
        ytcfg.update({"yt": {"suppress_stream_logging": True}})

    t0_0 = time.time()

    #################################
    # continue previous run:        #
    #  - read catalogue             #
    #  - resume skewers computation #
    #################################
    if args.continue_run:

        print("Continuing with exisiting run")
        print("Loading catalogue")
        t1_0 = time.time()
        # load catalogue
        catalogue = Table.read(f"{args.output_dir}/{args.catalogue_file}")

        # select the entries that were not previously run
        not_run_mask = np.array([
            not (os.path.isfile(f"{args.output_dir}/"+entry["name"]+"_spec_nonoise.fits.gz") and (entry["noise"] < 0.0 or os.path.isfile(f"{args.output_dir}/"+entry["name"]+"_spec.fits.gz")))
            for entry in catalogue
        ])
        not_run_catalogue = catalogue[not_run_mask]
        print(f"{len(not_run_catalogue)} skewer(s) were not previously run")

        # prepare variables to run
        snapshot_names = not_run_catalogue["snapshot_name"]
        redshifts = not_run_catalogue["z"]
        names = not_run_catalogue["name"]
        noise = not_run_catalogue["noise"]
        start_shifts = np.vstack([
            not_run_catalogue["start_shift_x"],
            not_run_catalogue["start_shift_y"],
            not_run_catalogue["start_shift_z"],
        ]).transpose()
        end_shifts = np.vstack([
            not_run_catalogue["end_shift_x"],
            not_run_catalogue["end_shift_y"],
            not_run_catalogue["end_shift_z"],
        ]).transpose()
        galaxy_positions = np.vstack([
            not_run_catalogue["gal_pos_x"],
            not_run_catalogue["gal_pos_y"],
            not_run_catalogue["gal_pos_z"],
        ]).transpose()

        t1_1 = time.time()
        print(f"INFO: Catalogue loaded. Elapsed time: {(t1_1-t1_0)/60.0} minutes")

    ####################################
    # new run:                         #
    #  - compute catalogue and skewers #
    ####################################
    else:
        print("Computing catalogue")
        t1_0 = time.time()
        np.random.seed(args.seed)

        # load snapshots info
        snapshots = np.genfromtxt(args.snapshots, names=True, dtype=None, encoding="UTF-8")
        snapshots_zmax = np.amax(snapshots["z_max"])

        # generate redshift distributions
        ndz = np.genfromtxt(args.z_dist, names=True, encoding="UTF-8")
        z_from_prob = interp1d(ndz["ndz_pdf"], ndz["z"])
        probs = np.random.uniform(0.0, 1.0, size=args.n_points)
        redshifts = z_from_prob(probs)
        pos = np.where(redshifts > snapshots_zmax)
        while pos[0].size > 0:
            print(
                f"WARNING: {pos[0].size} of the selected redshifts are higher "
                "than the largest snaphot redshift. I will now reassign these "
                "redshifs. This means the redshift distribution will be trimmed. "
                "Consider adding snapshots at larger redshifts or changing the "
                "input redshift distribution")
            probs = np.random.uniform(0.0, 1.0, size=pos[0].size)
            redshifts[pos] = z_from_prob(probs)
            pos = np.where(redshifts > snapshots_zmax)

        # generate random list of starting and ending points for the rays
        snapshots_rho_max = find_rho_max(redshifts, snapshots)
        rho = snapshots_rho_max * np.random.uniform(0, 1, size=args.n_points)**(1/3)
        theta_e = np.random.uniform(0, 2*np.pi, size=args.n_points)
        theta_r = np.random.uniform(0, 2*np.pi, size=args.n_points)
        phi_r = np.random.uniform(-np.pi, np.pi, size=args.n_points)
        x_start, y_start, z_start, x_end, y_end, z_end = generate_ray(
            rho, theta_e, theta_r, phi_r, 3*snapshots_rho_max)
        start_shifts = np.vstack([x_start, y_start, z_start]).transpose()
        end_shifts = np.vstack([x_end, y_end, z_end]).transpose()

        # generate noise distributions
        if args.noise_dist is not None:
            # TODO: draw noises
            # This needs to be replaced
            noise = np.zeros_like(redshifts) -1.0
        else:
            noise = np.zeros_like(redshifts) -1.0

        # choose snapshots
        choices = [
            select_snapshot(z_aux, rho_aux, snapshots)
            for z_aux, rho_aux in zip(redshifts, rho)]
        snapshot_names = snapshots["name"][choices]
        galaxy_position_x = snapshots["galaxy_pos_x"][choices]
        galaxy_position_y = snapshots["galaxy_pos_y"][choices]
        galaxy_position_z = snapshots["galaxy_pos_z"][choices]
        galaxy_positions = np.vstack([galaxy_position_x,
                                      galaxy_position_y,
                                      galaxy_position_z]).transpose()

        # get the simulation names
        names = np.array([
            (f"{args.base_name}_{snapshot}_z{z:.4f}_x{xs:.4f}_{xe:.4f}_"
             f"y{ys:.4f}_{ye:.4f}_z{zs:.4f}_{ze:.4f}")
            for snapshot, z, xs, xe, ys, ye, zs, ze in zip(
                snapshot_names,
                redshifts,
                x_start,
                x_end,
                y_start,
                y_end,
                z_start,
                z_end,
            )
        ])

        # save catalogue
        catalogue = Table({
            "name": names,
            "snapshot_name": snapshot_names,
            "z": redshifts,
            "rho": rho,
            "start_shift_x": x_start,
            "start_shift_y": y_start,
            "start_shift_z": z_start,
            "end_shift_x": x_end,
            "end_shift_y": y_end,
            "end_shift_z": z_end,
            "gal_pos_x": galaxy_position_x,
            "gal_pos_y": galaxy_position_y,
            "gal_pos_z": galaxy_position_z,
            "noise": noise,
        })
        catalogue.write(f"{args.output_dir}/{args.catalogue_file}")

        t1_1 = time.time()
        print(f"INFO: Catalogue created. Elapsed time: {(t1_1-t1_0)/60.0} minutes")

    # run the skewers in parallel
    print("Running skewers")
    t2_0 = time.time()
    for snapshot in np.unique(snapshot_names):
        ds = load_snapshot(snapshot, args.snapshots_dir)
        pos = np.where(snapshot_names == snapshot)
        context = multiprocessing.get_context('fork')
        with context.Pool(processes=args.num_processors) as pool:
            arguments = zip(
                repeat(ds),
                redshifts[pos],
                start_shifts[pos],
                end_shifts[pos],
                snapshot_names[pos],
                galaxy_positions[pos],
                names[pos],
                repeat(args.output_dir),
                noise[pos])

            pool.starmap(run_simple_ray_fast, arguments)


    t2_1 = time.time()
    print(f"INFO: Run {len(snapshot_names)} skewers. Elapsed time: {(t2_1-t2_0)/60.0} minutes")


    ###########################
    # fit for NHI, b and z    #
    # (spectra without noise) #
    ###########################
    if args.fit:
        print("Fitting profiles")
        t3_0 = time.time()

        print(catalogue["name"])

        context = multiprocessing.get_context('fork')
        with context.Pool(processes=args.num_processors) as pool:
            arguments = zip(
                catalogue["name"],
                repeat(".fits.gz"),
                repeat(False),
                catalogue["noise"]
            )

            print("Here")
            fit_results_list = pool.starmap(fit_lines, arguments)
            print("There")

            fit_results = np.array(
                fit_results_list,
                dtype=[("N [cm^-2]", float), ("b [km/s]", float), ("zfit", float)]
            )

            print(fit_results_list)
            print(fit_results_list[0])
            print(fit_results_list[:])
            print(fit_results_list[:][0])

            print(fit_results)

            # update catalogue
            catalogue["N [cm^-2]"] = fit_results[:][0]
            catalogue["b [km/s]"] = fit_results[:][1]
            catalogue["zfit"] = fit_results[:][2]
            print(catalogue)
            catalogue.write(args.catalogue_file)

        t3_1 = time.time()
        print(f"INFO: Fits done. Elapsed time: {(t3_1-t3_0)/60.0} minutes")

    t0_1 = time.time()
    print(f"INFO: total elapsed time: {(t0_1-t0_0)/60.0} minutes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=("Generate rays going randomly through a galaxy"))

    parser.add_argument("--base-name",
                        type=str,
                        default=f"galaxy_",
                        help="Base name for output files.")
    parser.add_argument("--catalogue-file",
                        type=str,
                        default=f"{THIS_DIR}/simulations/GCAT/G_catalog.csv",
                        help="Output catalogue filename. Extension should be csv")
    parser.add_argument("--continue-run",
                        action="store_true",
                        help="""Continue a previous run""")
    parser.add_argument("--fit",
                        action="store_true",
                        help="""If passed, fit line parameters on noiseless spectra
                            and add them to the catalogue""")
    parser.add_argument("--n-points",
                        type=int,
                        default=10,
                        help='Number of rays to draw')
    parser.add_argument("--noise-dist",
                        type=str,
                        default=None,
                        help="""File with the noise distribution of objects.
                            Currently this is ignored""")
    parser.add_argument("--num-processors",
                        type=int,
                        default=0,
                        help="""Number of processors to use. If 0 use
                            multiprocessing.cpu_count() // 2)""")
    parser.add_argument("--output-dir",
                        type=str,
                        required=True,
                        help="""Output directory""")
    parser.add_argument("--snapshots",
                        type=str,
                        default=f"{THIS_DIR}/../Data/snapshots.txt",
                        help='File containing the snapshot information')
    parser.add_argument("--snapshots-dir",
                        type=str,
                        default="",
                        help="Directory where the snapshots are placed")
    parser.add_argument("--seed",
                        type=int,
                        default=458467463,
                        help='Seed for the random number generator')
    parser.add_argument("--verbose",
                        action="store_true",
                        help="""If passed, then print all the info from yt and
                            trident""")
    parser.add_argument("--z-dist",
                        type=str,
                        default=f"{THIS_DIR}/../Data/dr16_dla_ndz.txt",
                        help="""File with the redshift distribution of objects.
                            Must have fields z and ndz_pdf""")
    # TODO: fix this so that the value passed here is actually used
    parser.add_argument("--z-sun",
                        type=float,
                        default=0.02041,
                        help="Solar metallicity")

    args = parser.parse_args()

    main(args)
