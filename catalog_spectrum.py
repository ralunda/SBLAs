import csv
import yt   
from yt.units import kpc
from yt.units import Mpc
from yt import YTArray
import trident
from yt.mods import *
import math
import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
cat = "/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/Q+G_catalog.csv"
print ("Welcome to the QSO and Galaxy Spectrum generator")
def catalog_menu(): 
    q = input('Enter QSO name: ')
    with open(cat, "r") as c:
        lector = csv.reader(c, delimiter=";")
        # Omitir el encabezado
        next(lector, None)
        contq = 0
        for fila in lector:
            nameq = fila[0]
            nameg = fila[1]
            zq = float(fila[2])
            zg = float(fila[3])
            dg = float(fila[4])
            ag = int(fila[5])
            ng = float(fila[6])
            abs = fila[7]
            if q == nameq:
                contq += 1
        if contq > 0:
            print("QSO is in the database")
            g = input('Enter Galaxy name: ')
            with open(cat, "r") as c:
                lector = csv.reader(c, delimiter=";")
                next(lector, None)
                contg = 0
                for fila in lector:
                    nameq = fila[0]
                    nameg = fila[1]
                    zq = float(fila[2])
                    zg = float(fila[3])
                    dg = float(fila[4])
                    ag = int(fila[5])
                    ng = float(fila[6])
                    abs = fila[7]
                    if g == nameg:
                        contg += 1
                if contg > 0:
                    print("Galaxy is in the database")
                    with open(cat, "r") as c:
                        lector = csv.reader(c, delimiter=";")
                        next(lector, None)
                        for fila in lector:
                            nameq = fila[0]
                            nameg = fila[1]
                            zq = float(fila[2])
                            zg = float(fila[3])
                            if q in nameq:
                                tq = zq
                                break 
                        for fila in lector:
                            nameq = fila[0]
                            nameg = fila[1]
                            zq = float(fila[2])
                            zg = float(fila[3])
                            if g in nameg:
                                tg = zg
                                break
                        if tq > tg:
                            #qso spectrum
                            rcParams.update({'figure.autolayout': True}) 
                            fig, ax = plt.subplots()
                            data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/" + str(q) +"spec.txt", delimiter=' ', unpack=True)
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
                            plt.title(str(q)+ ' z='+ str(tq), fontsize=20)
                            plt.show()
                            plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/'+ str(q) + 'spec.png')
                            plt.close()
                            #Galaxy spectrum
                            rcParams.update({'figure.autolayout': True}) 
                            fig, ax = plt.subplots()
                            data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/" + str(g) + "spec.txt", delimiter=' ', unpack=True)
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
                            plt.title(str(g) + ' z='+ str(tg), fontsize=20)
                            plt.show()
                            plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/' + str(g) + 'spec.png')
                            plt.close()
                            #spectrum combination
                            array1 = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/QCAT/" + str(q) +"spec.txt") 
                            array2 = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/GCAT/" + str(g) +"spec.txt")
                            array1[:,2] *= array2[:,2] 
                            array1[:,3] *= array2[:,3] 
                            array1[:,1] += array2[:,1] 
                            np.savetxt("/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/" + str(q) + "+" + str(g) + "spec.txt", array1)
                            rcParams.update({'figure.autolayout': True}) 
                            fig, ax = plt.subplots()
                            data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/" + str(q) + "+" + str(g) + "spec.txt", delimiter=' ', unpack=True)
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
                            plt.title(str(q) + '+' + str(g), fontsize=20)
                            plt.show()
                            plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/' + str(q) + '+' + str(g) + 'spec.png')
                            plt.close()
                            sg = trident.load_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/' + str(q) + '+' + str(g) + 'spec.txt')
                            sg.add_gaussian_noise(30)
                            sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/' + str(q) + '+' + str(g) + 'specNOISE.txt')
                            rcParams.update({'figure.autolayout': True}) 
                            fig, ax = plt.subplots()
                            data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/" + str(q) + "+" + str(g) + "specNOISE.txt", delimiter=' ', unpack=True)
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
                            plt.title(str(q) + '+' + str(g) + ' S/N 30', fontsize=20)
                            plt.show()
                            plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/CATQ+G/' + str(q) + '+' + str(g) + 'specNOISE.png')
                            plt.close()
                            print("All spectra have been generated")
                        else:
                            print("Non-valid combination of QSO and Galaxy. Please try again")
                            catalog_menu()
                                    
                else:
                    print("Name is not in the database")
                    catalog_menu()
        else:
            print("Name is not in the database")
            catalog_menu()

            
catalog_menu()
       



