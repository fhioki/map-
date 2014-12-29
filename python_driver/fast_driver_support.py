from mapsys import *
import matplotlib.pyplot as plt
from optparse import OptionParser
import numpy as np
import csv


class Vessel:
    time = []
    x = []
    y = []
    z = []
    phi = []
    the = []
    psi = []
    out = []


def get_vessel_column_index(name,out_chanel):
    index = [0,0,0,0,0,0,0]
    fp = open(name) 
    for i,line in enumerate(fp):
        if i==6:
            words = line = line.strip().split("\t")
            for j in range(0,len(words)):
                if words[j].strip()=='PtfmSurge':
                    index[0] = j;
                if words[j].strip()=='PtfmSway':
                    index[1] = j;
                if words[j].strip()=='PtfmHeave':
                    index[2] = j;
                if words[j].strip()=='PtfmRoll':
                    index[3] = j;
                if words[j].strip()=='PtfmPitch':
                    index[4] = j;
                if words[j].strip()=='PtfmYaw':
                    index[5] = j;
                if words[j].strip()==out_chanel:
                    index[6] = j;
    fp.close()
    return index


def set_vessel_prescribed_motion(table,index):
    vessel = Vessel()
    N = len(table)
    vessel.time = [float(table[i][0]) for i in range(8,N)]
    vessel.x = [float(table[i][index[0]]) for i in range(8,N)]
    vessel.y = [float(table[i][index[1]]) for i in range(8,N)]
    vessel.z = [float(table[i][index[2]]) for i in range(8,N)]
    vessel.phi = [float(table[i][index[3]]) for i in range(8,N)]
    vessel.the = [float(table[i][index[4]]) for i in range(8,N)]
    vessel.psi = [float(table[i][index[5]]) for i in range(8,N)]    
    vessel.out = [float(table[i][index[6]]) for i in range(8,N)]    
    return vessel


def get_line_tension(mooring, vessel, i):
    mooring.displace_vessel(vessel.x[i],vessel.y[i],vessel.z[i],vessel.phi[i],vessel.the[i],vessel.psi[i])
    mooring.update_states(vessel.time[i],0)
    fx,fy,fz = mooring.get_fairlead_force_3d(4)
    return np.sqrt(fx**2 + fy**2 + fz**2)



