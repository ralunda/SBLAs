
import numpy as np
import trident
# Ignasi: not sure what we are importing here. Not sure if it's used either
# I'm commenting it for the time being
#from yt.mods import *
import math
import h5py

from utils import run_snapshot, initialize_catalogue

# output configuration
base_name = "/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_"
catalogue_name = "/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/G_catalog.csv"

# common configuration in all runs
Z_Solar= 0.02041
z_min = 2.2
z_max = 2.6
z_step = 0.1
# TODO: add the angles configuration here
# (Ignasi: I currently don't understand how they're set)

def main():
    """Run the different simulations"""
    catalogue_file = initialize_catalogue(catalogue_name)

    ################
    # run: sim 196 #
    ################
    # specific configuration
    fn = 'RD0196/RD0196'
    galaxy_pos=np.array([11928.448765067616,12822.112568926215,12327.245902774346])
    dist_min = 2.0 # kpc
    dist_max = 22.0 # kpc
    dist_step = 2.0 # kpc
    # do the actual run
    try:
        run_snapshot(fn, galaxy_pos, z_min, z_max, z_step, dist_min, dist_max, dist_step, base_name, catalogue_file, starting_n=1)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    ################
    # run: sim 180 #
    ################
    # specific configuration
    fn = 'RD0180/RD0180'
    galaxy=np.array([10501.472788071803,11192.100157077519,10801.248702515462])
    dist_min = 17.0 # kpc
    dist_max = 187.0 # kpc
    dist_step = 17.0 # kpc
    # do the actual run
    try:
        run_snapshot(fn, galaxy_pos, z_min, z_max, z_step, dist_min, dist_max, dist_step, base_name, catalogue_file, starting_n=401)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    ################
    # run: sim 165 #
    ################
    # specific configuration
    fn = 'RD0165/RD0165'
    galaxy_pos=np.array([9350.457417869537,9893.696818590259,9571.331842591271])
    dist_min = 30.0 # kpc
    dist_max = 330.0 # kpc
    dist_step = 30.0 # kpc
    # do the actual run
    try:
        run_snapshot(fn, galaxy_pos, z_min, z_max, z_step, dist_min, dist_max, dist_step, base_name, catalogue_file, starting_n=801)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    ################
    # run: sim 152 #
    ################
    # specific configuration
    fn = 'RD0152/RD0152'
    galaxy_pos=np.array([8440.156388155721,8875.605892201605,8606.241907442962])
    dist_min = 30.0 # kpc
    dist_max = 330.0 # kpc
    dist_step = 30.0 # kpc
    # do the actual run
    try:
        run_snapshot(fn, galaxy_pos, z_min, z_max, z_step, dist_min, dist_max, dist_step, base_name, catalogue_file, starting_n=1201)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    # end of the run: close the catalogue file
    catalogue_file.close()

if __name__ == "__main__":
    main()
