# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 13:56:45 2022

@author: rober
"""

import numpy as np
import statistics as st

from transformacje import *

el_grs80 = Transformacje(model = "grs80")
plik = 'wsp_inp.txt'

tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)
t = np.genfromtxt(plik, delimiter=',', skip_header = 4)
rows,cols = np.shape(tablica)

blh = np.zeros((rows,cols))
XYZ = np.zeros((rows,cols))
xy2000 = np.zeros((rows,2))
xy92 = np.zeros((rows,2))

neu = np.zeros((rows,cols))
az_elev_dis = np.zeros((rows,2))

odleglosci = np.zeros((rows,2))
a = input("Wykonano transformacje do współrzędnych geodezyjnych. Czy wykonac inne transformacje? (Wprowadz numer: 1 - XYZ, 2 - u2000, 3 - u1992, 4 - NEU")
for i in range(rows):
    blh[i] = el_grs80.hirvonen(tablica[i,0],tablica[i,1],tablica[i,2])    
    if a == '1':
        XYZ[i] = el_grs80.blh2xyz(radians(blh[i,0]),radians(blh[i,1]),radians(blh[i,2]))
    elif a == '2':
        xy2000[i] = el_grs80.u2000(radians(blh[i,0]), radians(blh[i,1]))
    elif a == '3':
        xy92[i] = el_grs80.u1992(radians(blh[i,0]), radians(blh[i,1]))

if a == '4':
    neu = el_grs80.NEU(tablica)

b = input("Czy policzyc azymut i elewacje? (TAK / NIE)")
if b == 'TAK':
    az_elev_dis = el_grs80.Azel(tablica)
else:
    pass
c = input("Czy policzyć odleglosci 2D i 3D (TAK/NIE)")
if c == 'TAK':
    odleglosci = el_grs80.odl2D3D(tablica)
else:
    pass