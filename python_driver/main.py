if __name__ == '__main__':      
    from mapsys import *
    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    # import numpy as np
    # np.set_printoptions(precision=2)
    # np.set_printoptions(suppress=True)

    mooring_1 = Map( )
    
    mooring_1.map_set_sea_depth(100)
    mooring_1.map_set_gravity(9.81)
    mooring_1.map_set_sea_density(1020.0)

    mooring_1.read_file("../test/baseline_1.map") # 100 m depth    
    # mooring_1.read_file("../test/baseline_2.map") # 100 m depth
    # mooring_1.read_file("../test/NRELOffshrBsline5MW_Platform_OC3Hywind_segmented.map") # 320 m depth
    
    mooring_1.init( )

    
    # epsilon = 1e-3
    # K = mooring_1.linear(epsilon)    
    # print "\nHere is the linearized stiffness matrix with zero vessel displacement:"
    # print np.array(K)
    # 
    # H,V = mooring_1.get_fairlead_force_2d(0)    
    # print H, "  ", V
    # 
    # fx,fy,fz = mooring_1.get_fairlead_force_3d(0)    
    # print fx, "  ", fy, "  ", fz
    # 
    # ''' 
    # function residual at (hopefully) the solution
    # '''
    # # print mooring_1.funch(0) 
    # # print mooring_1.funcl(0)
    # # print mooring_1.funch(1)
    # # print mooring_1.funcl(1)
    # 
    # '''
    # derivatives at solution
    # '''
    # # print mooring_1.dxdh(0)
    # # print mooring_1.dxdv(0)    
    # # print mooring_1.dzdh(0)
    # # print mooring_1.dzdv(0)
    # # 
    # # print mooring_1.dxdh(1)
    # # print mooring_1.dxdv(1)    
    # # print mooring_1.dzdh(1)
    # # print mooring_1.dzdv(1)
    # 
    # x0 = mooring_1.plot_x( 0, 50 )
    # y0 = mooring_1.plot_y( 0, 50 )
    # z0 = mooring_1.plot_z( 0, 50 )    
    # 
    # x1 = mooring_1.plot_x( 1, 50 )
    # y1 = mooring_1.plot_y( 1, 50 )
    # z1 = mooring_1.plot_z( 1, 50 )
    # 
    # x2 = mooring_1.plot_x( 2, 50 )
    # y2 = mooring_1.plot_y( 2, 50 )
    # z2 = mooring_1.plot_z( 2, 50 )
    # 
    # x3 = mooring_1.plot_x( 3, 50 )
    # y3 = mooring_1.plot_y( 3, 50 )
    # z3 = mooring_1.plot_z( 3, 50 )    
    # 
    # x4 = mooring_1.plot_x( 4, 50 )
    # y4 = mooring_1.plot_y( 4, 50 )
    # z4 = mooring_1.plot_z( 4, 50 )
    # 
    # x5 = mooring_1.plot_x( 5, 50 )
    # y5 = mooring_1.plot_y( 5, 50 )
    # z5 = mooring_1.plot_z( 5, 50 )
    #     
    # fig = plt.figure()
    # ax = Axes3D(fig)
    # 
    # ax.plot(x0,y0,z0,'b-')
    # ax.plot(x1,y1,z1,'r-')
    # 
    # ax.plot(x2,y2,z2,'b-')
    # ax.plot(x3,y3,z3,'r-')
    # 
    # ax.plot(x4,y4,z4,'b-')
    # ax.plot(x5,y5,z5,'r-')
    # 
    # ax.set_xlabel('X [m]')
    # ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Z [m]')        
    # #ax.set_xlim([-50.0,50])        
    # #ax.set_ylim([-50.0,50])        
    # #ax.set_zlim([-30.0,-12])        
    # 
    # plt.show()
    #
    # mooring_1.end( )
    
