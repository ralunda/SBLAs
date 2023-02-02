"""
This file contains a script to run the galaxy simulations using a
uniformly distributed random rays through the 3D volume of the galaxy.
"""
import argparse
from itertools import repeat
import multiprocessing
import numpy as np
from scipy.interpolate import interp1d

from utils import (
    initialize_catalogue,
    load_snapshot,
    run_galaxy_snapshot,
    run_simple_ray
)

def compute_rays(ds,
                 z,
                 dist,
                 start_shift,
                 end_shift,
                 snapshot_name,
                 galaxy_pos,
                 base_name):
    """Test function to replace run_simple_ray.

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
    """
    name = base_name.split("/")[-1]
    return (
        f"{name}; {snapshot_name}; {z}; {dist}; "
        f"{start_shift[0]}; {start_shift[1]}; {start_shift[2]}; "
        f"{end_shift[0]}; {end_shift[1]}; {end_shift[2]}; "
        f"{0.0}; {0.0}; {0.0}\n"
    )


def generate_ray(rho, theta_e, theta_r, phi_r, radius):
    """Function to compute the starting and ending points of a ray

    Arguments
    ---------
    rho: float
    Distance between the ray and the centre of the galaxy

    theta_e: float
    Angle within the cone of possible rays

    theta_r: float
    Rotatiom angle theta for the starting and ending points

    phi_r: float
    Rotation angle phi for the starting and ending points

    radius: float
    Radius of the sphere including the starting and ending point

    Return
    ------
    x_start: float
    x coordinate of the starting point

    y_start: float
    y coordinate of the starting point

    z_start: float
    z coordinate of the starting point

    x_end: float
    x coordinate of the ending point

    y_end: float
    y coordinate of the ending point

    z_end: float
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

def select_snapshot(z, rho, snapshots):
    """Randomly select a snapshot given a choice of z and distance

    Arguments
    ---------
    z: float
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
    pos = np.argwhere((snapshots["rho_max"] > rho) & (snapshots["z_max"] > z))
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

    np.random.seed(args.seed)

    # generate random list of starting and ending points for the rays
    rho = args.rho_max * np.random.uniform(0, 1, size=args.n_points)**(1/3)
    theta_e = np.random.uniform(0, 2*np.pi, size=args.n_points)
    theta_r = np.random.uniform(0, 2*np.pi, size=args.n_points)
    phi_r = np.random.uniform(-np.pi, np.pi, size=args.n_points)
    x_start, y_start, z_start, x_end, y_end, z_end = generate_ray(
        rho, theta_e, theta_r, phi_r, 3*args.rho_max)
    start_shifts = np.vstack([x_start, y_start, z_start]).transpose()
    end_shifts = np.vstack([x_end, y_end, z_end]).transpose()

    # generate redshift distributions
    ndz = np.genfromtxt(args.z_dist, names=True, encoding="UTF-8")
    z_from_prob = interp1d(ndz["ndz_pdf"], ndz["z"])
    probs = np.random.uniform(0.0, 1.0, size=args.n_points)
    z = z_from_prob(probs)

    # choose snapshots
    snapshots = np.genfromtxt(args.snapshots, names=True, dtype=None, encoding="UTF-8")
    choices = [select_snapshot(z_aux, rho_aux, snapshots) for z_aux, rho_aux in zip(z, rho)]
    snapshot_names = snapshots["name"][choices]
    galaxy_positions = np.vstack([snapshots["galaxy_pos_x"][choices],
                                  snapshots["galaxy_pos_y"][choices],
                                  snapshots["galaxy_pos_z"][choices]]).transpose()

    # get the simulation names
    names = np.array([f"{args.base_name}{i}" for i in np.arange(args.n_points)])

    # initialize catalogue
    catalogue_file = initialize_catalogue(
        args.catalogue_file,
        quasar=False,
        galaxy=True)

    # run the rays in parallel
    try:
        for snapshot in np.unique(snapshot_names):
            if args.test:
                ds = None
            else:
                ds = load_snapshot(snapshot)
            pos = np.where(snapshot_names == snapshot)
            context = multiprocessing.get_context('fork')
            with context.Pool(processes=args.num_processors) as pool:
                arguments = zip(
                    repeat(ds),
                    z[pos],
                    rho[pos],
                    start_shifts[pos],
                    end_shifts[pos],
                    snapshot_names[pos],
                    galaxy_positions[pos],
                    names[pos])

                if args.test:
                    imap_it = pool.starmap(compute_rays, arguments)
                else:
                    imap_it = pool.starmap(run_simple_ray, arguments)
            for line in imap_it:
                catalogue_file.write(line)

    # close catalogue
    finally:
        catalogue_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=("Generate rays going randomly through a galaxy"))

    parser.add_argument("--base-name",
                        type=str,
                        default="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_",
                        help="Output catalogue filename. Extension should be csv")
    parser.add_argument("--catalogue-file",
                        type=str,
                        default="/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv",
                        help="Output catalogue filename. Extension should be csv")
    parser.add_argument("--rho-max",
                        type=float,
                        default=300,
                        help="Maximum distance to probe (in kpc)")
    parser.add_argument("--n-points",
                        type=int,
                        default=10,
                        help='Number of rays to draw')
    parser.add_argument("--num-processors",
                        type=int,
                        default=0,
                        help="""Number of processors to use. If 0 use
                            multiprocessing.cpu_count() // 2)""")
    parser.add_argument("--z-dist",
                        type=str,
                        default="../Data/dr16_dla_ndz.txt",
                        help="""File with the redshift distribution of objects.
                            Must have fields z and ndz_pdf""")
    parser.add_argument("--seed",
                        type=int,
                        default=458467463,
                        help='Seed for the random number generator')
    parser.add_argument("--snapshots",
                        type=str,
                        default="../Data/snapshots.txt",
                        help='File containing the snapshot information')
    parser.add_argument("--test",
                        action="store_true",
                        help='Use the test function instead of run_simple_ray')
    parser.add_argument("--z-sun",
                        type=float,
                        default=0.02041,
                        help="Solar metallicity")

    args = parser.parse_args()

    main(args)
