##QSO-SBLA
import numpy as np
array1 = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Qspec.txt") 
array2 = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-Sspec.txt")
array1[:,2] *= array2[:,2] 
array1[:,3] *= array2[:,3] 
array1[:,1] += array2[:,1] 
np.savetxt("/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QSspec.txt", array1)
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 
fig, ax = plt.subplots()
data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QSspec.txt", delimiter=' ', unpack=True)
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
plt.title('165 QSO-SBLA', fontsize=20)
plt.show()
plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QSspec.png')

sg = trident.load_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QSspec.txt')
sg.add_gaussian_noise(30)
sg.save_spectrum('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QSspecNOISE.txt')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 
fig, ax = plt.subplots()
data = np.loadtxt("/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QSspecNOISE.txt", delimiter=' ', unpack=True)
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
plt.title('165 QSO-SBLA SNR 30', fontsize=20)
plt.show()
plt.savefig('/home/ralunda/anaconda3/simulaciones/ENZO/165/165-1-QSspecNOISE.png')
