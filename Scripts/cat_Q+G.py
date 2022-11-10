import csv

from utils import initialize_catalogue

# output configuration
catalogue_name = "/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/QSO_catalog.csv"
catQ = "/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/QSO_catalog.csv"
catG = "/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/G_catalog.csv"

def main():
    """Run the quasar+galaxy simulations"""

    catalogue_file = initialize_catalogue(catalogue_name, quasar=True, galaxy=True)

    with open(catQ, "r") as q:
        lectorq = csv.reader(q, delimiter=";")
        # Omitir el encabezado
        next(lectorq, None)
        for fila in lectorq:
            nameq = fila[0]
            datasetq = fila[1]
            zq = float(fila[2])
            seedq = int(fila[3])
            with open(catG, "r") as g:
                lectorg = csv.reader(g, delimiter=";")
                # Omitir el encabezado
                next(lectorg, None)
                for fila in lectorg:
                    nameg = fila[0]
                    datasetg = fila[1]
                    zg = float(fila[2])
                    dg = float(fila[3])
                    ag = int(fila[4])
                    ng = float(fila[5])
                    bg = float(fila[6])
                    zfg = float(fila[7])
                    if  zq>zg:
                            catalogue_file.write("\n")
                            catalogue_file.write(nameq)
                            catalogue_file.write(";")
                            catalogue_file.write(nameg)
                            catalogue_file.write(";")
                            catalogue_file.write(str(zq))
                            catalogue_file.write(";")
                            catalogue_file.write(str(zg))
                            catalogue_file.write(";")
                            catalogue_file.write(str(dg))
                            catalogue_file.write(";")
                            catalogue_file.write(str(ag))
                            catalogue_file.write(";")
                            catalogue_file.write(str(ng))
                            catalogue_file.write(";")
                            if 10**14 <= ng <= 10**16:
                                abs = "SBLA"
                            elif 10**16 < ng < 10**18:
                                abs = "LLS"
                            elif 10**18 <= ng < 2*10**20:
                                abs = "sub-DLA"
                            elif 2*10**20 <= ng:
                                abs = "DLA"
                            else:
                                abs = "insufficient column density"
                            catalogue_file.write(abs)
                            catalogue_file.write("\n")

    # end of the run: close the catalogue file
    catalogue_file.close()

if __name__ == "__main__":
    main()
