'''
  Copyright (C) 2014 mdm and others
  Use at your own risk                                                             
  http://www.apache.org/licenses/LICENSE-2.0                 
'''  


if __name__ == '__main__':      
    from mapsys import *
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    np.set_printoptions(precision=2)
    np.set_printoptions(suppress=True)

    mooring_1 = Map( )
    
    mooring_1.map_set_sea_depth(100) # ficticious depth
    mooring_1.map_set_gravity(9.81)
    mooring_1.map_set_sea_density(1020.0)
    
    mooring_1.read_file("../test/OMAERB_Mooring.dat") # 100 m depth

    # mooring_1.summary_file('name_me.txt')
    mooring_1.init( )
    

    epsilon = 1e-5
    K = mooring_1.linear(epsilon)    
    print "\nHere is the linearized stiffness matrix with zero vessel displacement:"
    print np.array(K)

    # H,V = mooring_1.get_fairlead_force_2d(0)    
    # print H, "  ", V
      
    # fx,fy,fz = mooring_1.get_fairlead_force_3d(0)    
    # print fx, "  ", fy, "  ", fz

    fig = plt.figure()
    ax = Axes3D(fig)
    for i in range(0,mooring_1.size_lines()):
        x = mooring_1.plot_x( i, 10 )
        y = mooring_1.plot_y( i, 10 )
        z = mooring_1.plot_z( i, 10 )        
        ax.plot(x,y,z,'b-')
     
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')        
    ax.set_xlim([-3.0,3])        
    ax.set_ylim([-3.0,3])        
    ax.set_zlim([-3.0,0])        
     
    plt.show()
    
    mooring_1.end( )
    
