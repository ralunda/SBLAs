import csv
qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/Q+G_catalog.csv","a")
qc.write(";".join(["Name QSO", "Name Galaxy", "z QSO", "z Galaxy", "d", "Ang", "N [cm^-2]", "Absortion System"]))
qc.close()
catQ = "/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/QSO_catalog.csv"
catG = "/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/G_catalog.csv"
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
                        qc = open("/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/Q+G_catalog.csv","a")
                        qc.write("\n")
                        qc.write(nameq)
                        qc.write(";")
                        qc.write(nameg)
                        qc.write(";")
                        qc.write(str(zq))
                        qc.write(";")
                        qc.write(str(zg))
                        qc.write(";")
                        qc.write(str(dg))
                        qc.write(";")
                        qc.write(str(ag))
                        qc.write(";")
                        qc.write(str(ng))
                        qc.write(";")
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
                        qc.write(abs)
                        qc.close()


