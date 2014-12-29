#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Copyright (C) 2014 mdm                                     
  http://www.apache.org/licenses/LICENSE-2.0                 

  use like this:
  $ ./fast_driver.py -f ../test/Test22.out -n T[1]
'''  

from fast_driver_support import *

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
