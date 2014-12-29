#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Copyright (C) 2014 mdm                                     
  marco[dot]masciola[at]gmail                                
                                                             
Licensed to the Apache Software Foundation (ASF) under one   
or more contributor license agreements.  See the NOTICE file 
distributed with this work for additional information        
regarding copyright ownership.  The ASF licenses this file   
to you under the Apache License, Version 2.0 (the            
"License"); you may not use this file except in compliance   
with the License.  You may obtain a copy of the License at   
                                                             
  http://www.apache.org/licenses/LICENSE-2.0                 
                                                             
Unless required by applicable law or agreed to in writing,   
software distributed under the License is distributed on an  
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       
KIND, either express or implied.  See the License for the    
specific language governing permissions and limitations            
under the License.                                             
'''  


# use like this:
# 
# $ ./fast_driver.py -f ../test/Test22.out -n T[1]



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



if __name__ == '__main__':      
    # command line argument processing
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('-f', '--file', dest='file_name', help='read data from FILENAME')
    parser.add_option('-n', '--output stream', dest='chanel', help='read output stream from file')
    (options, args) = parser.parse_args()

    # index is the surge, ... , yaw column number in the FAST output file
    index = get_vessel_column_index(options.file_name,options.chanel) 

    mooring = Map()
    mooring.map_set_sea_depth(150) # 150 for barge
    mooring.map_set_gravity(9.81) 
    mooring.map_set_sea_density(1020.0)
    mooring.read_file('../test/NREL_Barge.map') # 100 m depth
    mooring.summary_file('barge.sum.txt')
    mooring.init()

    # read column data
    with open(options.file_name, 'rb') as csv_file:        
        reader = csv.reader(csv_file,delimiter='\t')
        table = list(reader)

    # time marching, vessel is displaced here
    vessel = set_vessel_prescribed_motion(table,index)
    tension = [get_line_tension(mooring,vessel,i) for i in range(0,len(vessel.time))]

    mooring.end()

    plt.figure(1)
    plt.title('Tension T[5]')
    plt.xlabel('Time [sec]')
    plt.ylabel('Force [kN]')
    plt.plot(vessel.time,vessel.out,'b',lw=3,alpha=0.5,label='FAST Output')
    plt.plot(vessel.time,tension,'k',label='MAP++ Python Driver')
    plt.legend(loc=2)
    plt.show()    
