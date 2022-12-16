from utils import run_quasar_snapshot, initialize_catalogue

# output configuration
base_name = "/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/QSO_"
catalogue_name = "/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/QSO_catalog.csv"

# common configuration in all runs
Z_Solar= 0.02041
z_min = 2.50
z_step = 0.05


def main():
    """Run the quasar simulations"""
    catalogue_file = initialize_catalogue(catalogue_name, quasar=True)

    ################
    # run: sim 196 #
    ################
    # specific configuration
    fn = 'RD0196/RD0196'
    z_max = 2.55
    # do the actual run
    try:
        run_quasar_snapshot(fn, z_min, z_max, z_step, base_name, catalogue_file, starting_n=1)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    ################
    # run: sim 180 #
    ################
    # specific configuration
    fn = 'RD0180/RD0180'
    z_max = 3.05
    # do the actual run
    try:
        run_quasar_snapshot(fn, z_min, z_max, z_step, base_name, catalogue_file, starting_n=11)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    ################
    # run: sim 165 #
    ################
    # specific configuration
    fn = 'RD0165/RD0165'
    z_max = 3.55
    # do the actual run
    try:
        run_quasar_snapshot(fn, z_min, z_max, z_step, base_name, catalogue_file, starting_n=121)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    ################
    # run: sim 152 #
    ################
    # specific configuration
    fn = 'RD0152/RD0152'
    z_max = 4.05
    # do the actual run
    try:
        run_quasar_snapshot(fn, z_min, z_max, z_step, base_name, catalogue_file, starting_n=331)
    # in case of errors, stop the run and close the catalogue file
    except:
        catalogue.close()
        return

    # end of the run: close the catalogue file
    catalogue_file.close()

if __name__ == "__main__":
    main()
