# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 13:36:10 2022

@author: rober
"""

from math import sqrt, atan, sin, cos, degrees, radians, tan, pi, atan2, asin
import numpy as np
import statistics as st
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2
        
    def hirvonen(self, X, Y, Z):
        """
        

        Parameters
        ----------
        X : wspolrzedna X w ukladzie wspolrzednych kartezjanskich [m]
        Y : wspolrzedna Y w ukladzie wspolrzednych kartezjanskich [m]
        Z : wspolrzedna Z w ukladzie wspolrzednych kartezjanskich [m]

        Returns
        -------
        fi : szerokosc geograficzna [stopnie]
        la : dlugosc geograficzna [stopnie]
        h : wysokosc [m]

        """
        l = atan2(Y, X)
        p = sqrt(X**2 + Y**2)
        fi = atan(Z / (p * (1 - self.ecc2)))
        
        while 1:
            N = self.a / sqrt(1 - self.ecc2 * (sin(fi)**2))
            h = p / cos(fi) - N
            fi_p = fi
            fi = atan(Z / (p * (1 - N * self.ecc2 / (N + h))))
            if abs(fi - fi_p) < (0.000001 / 206265):
                break;
        return degrees(fi), degrees(l), h
    
    def blh2xyz(self, B, L, H):
        """
        

        Parameters
        ----------
        B : szerokosc geograficzna [rad]
        L : dlugosc geograficzna [rad]
        H : wysokosc [m]

        Returns
        -------
        x : wspolrzedna X w ukladzie wspolrzednych kartezjanskich [m]
        y : wspolrzedna Y w ukladzie wspolrzednych kartezjanskich [m]
        z : wspolrzedna Z w ukladzie wspolrzednych kartezjanskich [m]

        """
        N = self.a / sqrt(1 - self.ecc2 * (sin(B)**2))
        x = round((N + H) * cos(B) * cos(L), 3)
        y = round((N + H) * cos(B) * sin(L), 3)
        z = round((N * (1 - self.ecc2) + H) * sin(B), 3)
        return x, y, z
    
    def u2000(self,fi,lam):
        """
        

        Parameters
        ----------
        fi : szerokosc geograficzna [rad]
        lam : dlugosc geograficzna [rad]

        Returns
        -------
        x2000 : wspolrzedna X w ukladzie 2000 [m]
        y2000 : wspolrzedna Y w ukladzie 2000 [m]

        """
        b2 = (self.a**2) * (1 - self.ecc2)
        ep2 = (self.a**2 - b2) / b2
        
        lam = degrees(lam)
        
        if lam > 13.5 and lam <= 16.5:
            L0 = 15
        elif lam > 16.5 and lam <= 19.5:
            L0 = 18
        elif lam > 19.5 and lam <= 22.5:
            L0 = 21
        elif lam > 22.5 and lam <= 25.5:
            L0 = 24
        
        strefa = L0 / 3
        dl = radians(lam - L0)
        
        t = tan(fi)
        n2 = ep2 * (cos(fi)**2)
        N = self.a / sqrt(1 - self.ecc2 * (sin(fi)**2))
        
        A0 = 1 - (self.ecc2 / 4) - (3 / 64) * (self.ecc2**2) - (5 / 256) * (self.ecc2**3)
        A2 = (3 / 8) * (self.ecc2 + (1 / 4) * (self.ecc2**2) + (15 / 128) * (self.ecc2**3))
        A4 = (15 / 256) * (self.ecc2**2 + (3 / 4) * (self.ecc2**3))
        A6 = (35 / 3072) * (self.ecc2**3)
        
        si = self.a * (A0 * fi - A2 * sin(2*fi) + A4 * sin(4*fi) - A6 * sin(6*fi))

        xgk = si + (1 / 2) * (dl**2) * N * sin(fi) * cos(fi) * (1 + (1 / 12) * (dl**2) * (cos(fi)**2) * (5 - (t**2) + 9 * n2 + 4 * (n2**2) + (1 / 360) * (dl**4) * (cos(fi)**4) * (61 - 58 * (t**2) + (t**4) + 270 * n2 - 330 * n2 * (t**2))))
        ygk = dl * N * cos(fi) * (1 + (1 / 6) * (dl**2) * (cos(fi)**2) * (1 - (t**2) + n2) + (1 / 120) * (dl**4) * (cos(fi)**4) * (5 - 18 * (t**2) + (t**4) + 14 * n2 - 58 * n2 * (t**2)))
        
        x2000 = round(xgk * 0.999923, 3) 
        y2000 = round(ygk * 0.999923 + strefa * 1000000 + 500000, 3)
        
        return x2000, y2000
    
    def u1992(self,fi,lam):
        """
        

        Parameters
        ----------
        fi : szerokosc geograficzna [rad]
        lam : dlugosc geograficzna [rad]

        Returns
        -------
        x1992 : wspolrzedna X w ukladzie 1992 [m]
        y1992 : wspolrzedna Y w ukladzie 1992 [m]

        """
        b2 = (self.a**2) * (1 - self.ecc2)
        ep2 = (self.a**2 - b2) / b2
        
        L0 = 19 * pi / 180
        dl = lam - L0
        
        t = tan(fi)
        n2 = ep2 * (cos(fi)**2)
        N = self.a / sqrt(1 - self.ecc2 * (sin(fi)**2))
        
        A0 = 1 - (self.ecc2 / 4) - (3 / 64) * (self.ecc2**2) - (5 / 256) * (self.ecc2**3)
        A2 = (3 / 8) * (self.ecc2 + (1 / 4) * (self.ecc2**2) + (15 / 128) * (self.ecc2**3))
        A4 = (15 / 256) * (self.ecc2**2 + (3 / 4) * (self.ecc2**3))
        A6 = (35 / 3072) * (self.ecc2**3)
        
        si = self.a * (A0 * fi - A2 * sin(2*fi) + A4 * sin(4*fi) - A6 * sin(6*fi))

        xgk = si + (1 / 2) * (dl**2) * N * sin(fi) * cos(fi) * (1 + (1 / 12) * (dl**2) * (cos(fi)**2) * (5 - (t**2) + 9 * n2 + 4 * (n2**2) + (1 / 360) * (dl**4) * (cos(fi)**4) * (61 - 58 * (t**2) + (t**4) + 270 * n2 - 330 * n2 * (t**2))))
        ygk = dl * N * cos(fi) * (1 + (1 / 6) * (dl**2) * (cos(fi)**2) * (1 - (t**2) + n2) + (1 / 120) * (dl**4) * (cos(fi)**4) * (5 - 18 * (t**2) + (t**4) + 14 * n2 - 58 * n2 * (t**2)))

        x1992 = round(xgk * 0.9993 - 5300000, 3) 
        y1992 = round(ygk * 0.9993 + 500000, 3)
        
        return x1992, y1992
    
    def Azel(self,tablica):
        """
        
        Parameters
        ----------
        tablica : tablica ze współrzędnymi X Y Z [m]

        Returns
        -------
        Azel: tablica z wartosciami azymutu i elewacji [stopnie]

        """
        neu = self.NEU(tablica)
        
        n = np.zeros((len(neu),1))
        e = np.zeros((len(neu),1))
        u = np.zeros((len(neu),1))
        
        i = 0
        while i < len(neu):
            n[i,0] = neu[i,0]
            e[i,0] = neu[i,1]
            u[i,0] = neu[i,2]
            i += 1
        
        sredni_X = st.mean(tablica[:,0])                                           
        sredni_Y = st.mean(tablica[:,1])
        sredni_Z = st.mean(tablica[:,2])
        
        odl = sqrt(sredni_X ** 2 + sredni_Y ** 2 + sredni_Z ** 2)
        Azel = np.zeros((len(neu),2))
        i = 0
        for i in range(0,len(tablica)):
            Azel[i,0] = degrees(atan2(e[i], n[i]))
            Azel[i,1] = degrees(asin(u[i] / sqrt(n[i]**2 + e[i]**2 + u[i]**2)))
            i +=1
        
        return Azel
    
    def NEU(self,tablica):
        """
        
        Parameters
        ----------
        tablica : tablica ze współrzędnymi X Y Z [m]

        Returns
        -------
        neu : tablica ze współrzędnymi w układzie NEU

        """

        sredni_X = st.mean(tablica[:,0])                                           
        sredni_Y = st.mean(tablica[:,1])
        sredni_Z = st.mean(tablica[:,2])
    
    
        f_sr, l_sr, h_sr = self.hirvonen (sredni_X, sredni_Y, sredni_Z)   
    
        R = np.array([[-np.sin(f_sr)*np.cos(l_sr), -np.sin(l_sr), np.cos(f_sr)*np.cos(l_sr)],
                      [-np.sin(f_sr)*np.sin(l_sr), np.cos(l_sr),np.cos(f_sr)*np.sin(l_sr)],
                      [np.cos(f_sr), 0,np.sin(f_sr)]])
    
    
        R_T = np.transpose(R)                                                      
    
        i = 0
        roznice_xyz = np.zeros((len(tablica), 3))                                  
        while i < len(tablica):
            roznice_xyz[i,0] = tablica[i,0] - sredni_X
            roznice_xyz[i,1] = tablica[i,1] - sredni_Y
            roznice_xyz[i,2] = tablica[i,2] - sredni_Z
            i += 1                                                            
            
        E = []
        N = []
        U = []
        
        i = 0
        for j in roznice_xyz:
            n = np.dot(R_T[0,:], j)                                                
            e = np.dot(R_T[1,:], j)
            u = np.dot(R_T[2,:], j)
            i += 1
            N.append(n)                                                            
            E.append(e)
            U.append(u)
    
        neu = np.zeros((len(N),3))                                                 
        i = 0
        while i < len(N):
            neu[i,0] = N[i]
            neu[i,1] = E[i]
            neu[i,2] = U[i]
            i += 1                         
                                        
        return neu
    
    def odl2D3D(self,tablica):
        """
        
        Parameters
        ----------
        tablica : tablica ze współrzędnymi X Y Z [m]

        Returns
        -------
        odleglosci : tablica z odległosciami 2D (pierwsza kolumna) i 3D (druga kolumna) [m]

        """
        sredni_X = st.mean(tablica[:,0])                                           
        sredni_Y = st.mean(tablica[:,1])
        sredni_Z = st.mean(tablica[:,2])
    
        odleglosci = np.zeros((len(tablica), 2))
        
        for i in range(0, len(tablica)):
            odl2d = sqrt((sredni_X - tablica[i][0])**2 + (sredni_Y - tablica[i][1])**2 )
            odleglosci[i,0] = round(odl2d, 3)  
            odl3d = sqrt((sredni_X - tablica[i][0])**2 + (sredni_Y - tablica[i][1])**2 + (sredni_Z - tablica[i][2])**2 )
            odleglosci[i,1] = round(odl3d, 3)    
            
        return odleglosci